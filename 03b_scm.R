#!/usr/bin/env Rscript
# ==============================================================================
# Stage 3b: Abadie-style Synthetic Control Method (SCM)
#
# Implements a proper state-level SCM for California's TANF lifetime limit
# extension (effective 2021).  Collapses the HH-year panel to state-year means,
# constructs a synthetic California from donor states, and runs:
#   - Main SCM (EPAY and unconditional TPAYWK, separately)
#   - Permutation / placebo inference (rank-based p-value)
#   - Leave-one-out donor sensitivity
#   - In-time placebo (fake treatment at 2017)
#
# All outputs land in BASE_DIR.  2020 is excluded from optimization but its
# gap is still reported and plotted.
# ==============================================================================

rm(list = ls())

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# ---- Section 0: Config -------------------------------------------------------

read_env_config <- function(path) {
  if (!file.exists(path)) stop("Missing shared config: ", path)
  x <- trimws(readLines(path, warn = FALSE))
  x <- x[nzchar(x) & substr(x, 1, 1) != "#"]
  out <- list()
  for (ln in x) {
    p <- strsplit(ln, "=", fixed = TRUE)[[1]]
    if (length(p) >= 2) out[[trimws(p[1])]] <- trimws(paste(p[-1], collapse = "="))
  }
  out
}

script_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) normalizePath("scripts/03b_scm.R", winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base       <- cfg$BASE_DIR
panel_file <- file.path(base, cfg$ANNUAL_PANEL_FILE)
out_dir    <- base

cat("======================================================================\n")
cat("[config] panel:", panel_file, "\n")
if (!file.exists(panel_file)) stop("Panel file missing: ", panel_file)

# SCM parameters
treat_state  <- as.integer(cfg$SCM_TREAT_STATE     %||% "6")
pre_start    <- as.integer(cfg$SCM_PRE_START        %||% "2014")
pre_end      <- as.integer(cfg$SCM_PRE_END          %||% "2019")
post_start   <- as.integer(cfg$SCM_POST_START       %||% "2021")
min_hh       <- as.integer(cfg$SCM_MIN_HH_PER_STATE_YEAR %||% "5")
excl_donors  <- as.integer(strsplit(cfg$SCM_EXCLUDE_DONORS %||% "6,9,10,44", ",", fixed = TRUE)[[1]])
exclude_analysis_years <- as.integer(strsplit(cfg$EXCLUDE_YEARS %||% "2020", ",", fixed = TRUE)[[1]])

# Predictor variables (non-outcome)
predictor_vars <- trimws(strsplit(cfg$CONTROLS %||%
  "ERN,WRTH,N_children,TMWKHRS_SUM_PARENTS,HHMAXEEDUC_HS", ",", fixed = TRUE)[[1]])

# Outcomes to run SCM on
scm_outcomes <- c("EPAY", "TPAYWK")   # TPAYWK is unconditional (0 for non-payers)

cat("[config] treat_state:", treat_state, "\n")
cat("[config] pre-period: ", pre_start, "-", pre_end, "\n", sep = "")
cat("[config] post-period:", post_start, "+\n")
cat("[config] excl_donors:", paste(excl_donors, collapse = ","), "\n")
cat("[config] excl_analysis_years:", paste(exclude_analysis_years, collapse = ","), "\n")
cat("[config] outcomes:  ", paste(scm_outcomes, collapse = ", "), "\n")
cat("[config] predictors:", paste(predictor_vars, collapse = ", "), "\n")
cat("======================================================================\n")

# ---- Section 1: Load and prepare panel ---------------------------------------

D <- read.csv(panel_file, check.names = FALSE, stringsAsFactors = FALSE)
num_vars <- unique(c("TEHC_ST", "CALYR", "TREATED_STATE", "N_children",
                     scm_outcomes, predictor_vars, "SSUID", "TPAYWK", "EPAY"))
for (v in intersect(num_vars, names(D))) D[[v]] <- suppressWarnings(as.numeric(D[[v]]))
D <- D[!is.na(D$N_children) & D$N_children > 0, , drop = FALSE]
if (!nrow(D)) stop("No rows after child restriction.")

# Ensure TPAYWK_ASINH exists but is not used as SCM outcome (unconditional TPAYWK is)
if ("TPAYWK" %in% names(D)) {
  tp <- D$TPAYWK
  tp[!is.na(tp) & tp < 0] <- NA_real_
  D$TPAYWK <- tp   # cleaned version
}

cat(sprintf("[data] Loaded %d HH-year rows across %d states\n",
            nrow(D), length(unique(D$TEHC_ST[!is.na(D$TEHC_ST)]))))

# ---- Section 2: Collapse to state-year means ---------------------------------

collapse_to_state_year <- function(D, outcome_vars, pred_vars) {
  vars_use  <- intersect(unique(c(outcome_vars, pred_vars)), names(D))
  states    <- sort(unique(D$TEHC_ST[!is.na(D$TEHC_ST)]))
  all_years <- sort(unique(D$CALYR[!is.na(D$CALYR)]))
  rows <- list()
  for (st in states) {
    for (yr in all_years) {
      m <- !is.na(D$TEHC_ST) & D$TEHC_ST == st & !is.na(D$CALYR) & D$CALYR == yr
      if (!any(m)) next
      row <- data.frame(
        TEHC_ST = st,
        CALYR   = yr,
        N_hh    = length(unique(as.character(D$SSUID[m]))),
        stringsAsFactors = FALSE
      )
      for (v in vars_use) row[[v]] <- mean(D[[v]][m], na.rm = TRUE)
      rows[[length(rows) + 1]] <- row
    }
  }
  if (!length(rows)) return(data.frame())
  SY <- do.call(rbind, rows)
  SY
}

cat("[collapse] Building state-year means...\n")
SY <- collapse_to_state_year(D, scm_outcomes, predictor_vars)
cat(sprintf("[collapse] %d state-year cells across %d states\n",
            nrow(SY), length(unique(SY$TEHC_ST))))

# ---- Section 3: Donor pool eligibility ---------------------------------------

get_eligible_donors <- function(SY, treat_st, excl, pre_yrs, min_hh_per_yr) {
  candidates <- setdiff(unique(SY$TEHC_ST), excl)
  keep <- candidates[sapply(candidates, function(st) {
    ss <- SY[SY$TEHC_ST == st, , drop = FALSE]
    years_present <- ss$CALYR
    # Must have data in every pre-period year
    if (!all(pre_yrs %in% years_present)) return(FALSE)
    # Must meet minimum HH threshold in every pre-period year
    hh_counts <- ss$N_hh[match(pre_yrs, ss$CALYR)]
    all(!is.na(hh_counts) & hh_counts >= min_hh_per_yr)
  })]
  keep
}

pre_years <- seq.int(pre_start, pre_end)
donors    <- get_eligible_donors(SY, treat_state, excl_donors, pre_years, min_hh)
cat(sprintf("[donors] Eligible donor states (%d): %s\n",
            length(donors), paste(donors, collapse = ", ")))

# Check treated unit has full pre-period coverage
treat_pre <- SY[SY$TEHC_ST == treat_state & SY$CALYR %in% pre_years, ]
if (nrow(treat_pre) < length(pre_years)) {
  cat(sprintf("[WARN] Treated state %d missing %d pre-period years. Proceeding with available years.\n",
              treat_state, length(pre_years) - nrow(treat_pre)))
  pre_years <- sort(intersect(pre_years, treat_pre$CALYR))
  cat(sprintf("[donors] Revised pre-years: %s\n", paste(pre_years, collapse = ",")))
}

if (length(donors) < 2) stop("Fewer than 2 donor states qualify. Cannot estimate SCM.")

# ---- Section 4: SCM weight optimizer -----------------------------------------
# Minimise ||x_treat - X_donors %*% w||^2 subject to w >= 0, sum(w) = 1
# Parameterised via softmax: w = exp(theta) / sum(exp(theta))

make_predictor_vec <- function(SY, state_id, pre_yrs, outcome_vars, pred_vars) {
  ss <- SY[SY$TEHC_ST == state_id, , drop = FALSE]
  ss <- ss[match(pre_yrs, ss$CALYR), , drop = FALSE]
  # Outcome lags: each pre-year value separately (most important predictors)
  vec <- c()
  for (ov in outcome_vars) {
    if (ov %in% names(ss)) vec <- c(vec, as.numeric(ss[[ov]]))
  }
  # Predictor means over pre-period
  for (pv in pred_vars) {
    if (pv %in% names(ss)) vec <- c(vec, mean(as.numeric(ss[[pv]]), na.rm = TRUE))
  }
  # Replace non-finite values with column mean (shouldn't happen often)
  bad <- !is.finite(vec)
  if (any(bad) && any(!bad)) vec[bad] <- mean(vec[!bad])
  vec
}

fit_scm_weights <- function(SY, treat_st, donor_sts, pre_yrs, outcome_vars, pred_vars,
                             n_restarts = 5) {
  x_t <- make_predictor_vec(SY, treat_st,   pre_yrs, outcome_vars, pred_vars)
  X_d <- sapply(donor_sts, function(st)
    make_predictor_vec(SY, st, pre_yrs, outcome_vars, pred_vars))
  if (is.null(dim(X_d))) X_d <- matrix(X_d, ncol = 1)

  K <- length(donor_sts)
  obj <- function(theta) {
    exp_t <- exp(theta - max(theta))
    w     <- exp_t / sum(exp_t)
    fit   <- as.numeric(X_d %*% w)
    sum((x_t - fit)^2, na.rm = TRUE)
  }
  grad <- function(theta) {
    exp_t <- exp(theta - max(theta))
    w     <- exp_t / sum(exp_t)
    fit   <- as.numeric(X_d %*% w)
    resid <- x_t - fit
    # Jacobian of w wrt theta_k: dw_k/dtheta_j = w_k*(I(j==k) - w_j)
    dLdw  <- -2 * as.numeric(t(X_d) %*% resid)
    # Chain rule: dL/dtheta_j = sum_k dL/dw_k * dw_k/dtheta_j
    #                          = sum_k dL/dw_k * w_k * (I(k==j) - w_j)
    #                          = dL/dw_j * w_j - w_j * sum_k(dL/dw_k * w_k)
    dLdtheta <- dLdw * w - w * sum(dLdw * w)
    dLdtheta
  }

  best_opt <- NULL
  for (i in seq_len(n_restarts)) {
    theta0 <- if (i == 1) rep(0, K) else rnorm(K, 0, 1)
    opt <- tryCatch(
      optim(theta0, obj, grad, method = "BFGS",
            control = list(maxit = 2000, reltol = 1e-10)),
      error = function(e) NULL
    )
    if (!is.null(opt) && (is.null(best_opt) || opt$value < best_opt$value))
      best_opt <- opt
  }
  if (is.null(best_opt)) return(NULL)

  exp_t <- exp(best_opt$par - max(best_opt$par))
  w     <- exp_t / sum(exp_t)
  names(w) <- as.character(donor_sts)

  # RMSPE over pre-period (fit quality)
  fit_pre <- as.numeric(X_d %*% w)
  # only outcome-lag portion
  n_out_lags <- length(pre_yrs) * length(outcome_vars)
  rmspe_pre <- sqrt(mean((x_t[seq_len(n_out_lags)] - fit_pre[seq_len(n_out_lags)])^2,
                         na.rm = TRUE))
  list(weights = w, rmspe_pre = rmspe_pre, converged = (best_opt$convergence == 0))
}

# ---- Section 5: Gap series computation ---------------------------------------

compute_gap_series <- function(SY, treat_st, donor_sts, w, outcome_var, all_years) {
  treat_row <- match(all_years, SY$CALYR[SY$TEHC_ST == treat_st])
  actual <- sapply(all_years, function(yr) {
    v <- SY[[outcome_var]][SY$TEHC_ST == treat_st & SY$CALYR == yr]
    if (!length(v)) NA_real_ else v[1]
  })
  synthetic <- sapply(all_years, function(yr) {
    vals <- sapply(donor_sts, function(st) {
      v <- SY[[outcome_var]][SY$TEHC_ST == st & SY$CALYR == yr]
      if (!length(v)) NA_real_ else v[1]
    })
    sum(w * vals, na.rm = FALSE)  # NA if any donor missing
  })
  data.frame(CALYR = all_years, actual = actual, synthetic = synthetic,
             gap = actual - synthetic, stringsAsFactors = FALSE)
}

# ---- Section 6: RMSPE ratio (post/pre) ---------------------------------------

rmspe_ratio <- function(gaps_df, pre_yrs, post_yrs) {
  pre_gaps  <- gaps_df$gap[gaps_df$CALYR %in% pre_yrs]
  post_gaps <- gaps_df$gap[gaps_df$CALYR %in% post_yrs]
  rmspe_pre  <- sqrt(mean(pre_gaps^2,  na.rm = TRUE))
  rmspe_post <- sqrt(mean(post_gaps^2, na.rm = TRUE))
  if (!is.finite(rmspe_pre) || rmspe_pre == 0) return(NA_real_)
  rmspe_post / rmspe_pre
}

# ---- Section 7: Main estimation loop (one per outcome) -----------------------

all_years_data    <- sort(unique(SY$CALYR))
# Years used for analysis (optimization): exclude contaminated years
analysis_years    <- setdiff(all_years_data, exclude_analysis_years)
post_years        <- analysis_years[analysis_years >= post_start]

cat("\n======================================================================\n")
cat("  MAIN SCM ESTIMATION\n")
cat("======================================================================\n")

results <- list()

for (outcome_var in scm_outcomes) {
  cat(sprintf("\n--- Outcome: %s ---\n", outcome_var))

  # Use only the current outcome as the "lag" predictor
  cur_outcomes <- outcome_var
  cur_preds    <- predictor_vars

  # Filter donors: re-check with this specific outcome (NaN check)
  donors_ok <- donors[sapply(donors, function(st) {
    ss <- SY[SY$TEHC_ST == st & SY$CALYR %in% pre_years, , drop = FALSE]
    !any(is.nan(ss[[outcome_var]]))
  })]
  cat(sprintf("  Donors after NaN check: %d\n", length(donors_ok)))

  scm_fit <- fit_scm_weights(SY, treat_state, donors_ok, pre_years,
                              cur_outcomes, cur_preds)
  if (is.null(scm_fit)) {
    cat("  FAILED: weight optimization returned NULL\n")
    next
  }

  cat(sprintf("  Converged: %s | Pre-RMSPE: %.6f\n",
              scm_fit$converged, scm_fit$rmspe_pre))
  cat("  Weights (>0.01):\n")
  w <- scm_fit$weights
  for (i in order(w, decreasing = TRUE)) {
    if (w[i] > 0.005) cat(sprintf("    FIPS %3d: %.4f\n", donors_ok[i], w[i]))
  }

  # Gap series over ALL years (including excluded 2020)
  gaps <- compute_gap_series(SY, treat_state, donors_ok, w, outcome_var, all_years_data)
  gaps$excluded_year <- as.integer(gaps$CALYR %in% exclude_analysis_years)
  gaps$pre_period    <- as.integer(gaps$CALYR %in% pre_years)
  gaps$post_period   <- as.integer(gaps$CALYR %in% post_years)

  ratio <- rmspe_ratio(gaps[!(gaps$CALYR %in% exclude_analysis_years), ],
                       pre_years, post_years)
  cat(sprintf("  RMSPE ratio (post/pre): %.4f\n", ratio))

  # ---- Permutation inference -------------------------------------------------
  cat("  Running permutation tests...\n")
  placebo_gaps <- list()
  for (pl_st in donors_ok) {
    pl_donors <- setdiff(donors_ok, pl_st)
    if (length(pl_donors) < 2) next
    pl_fit <- tryCatch(
      fit_scm_weights(SY, pl_st, pl_donors, pre_years, cur_outcomes, cur_preds,
                      n_restarts = 3),
      error = function(e) NULL
    )
    if (is.null(pl_fit)) next
    pl_gaps <- compute_gap_series(SY, pl_st, pl_donors, pl_fit$weights,
                                  outcome_var, all_years_data)
    pl_rmspe_pre  <- rmspe_ratio(pl_gaps, pre_years, post_years)
    placebo_gaps[[as.character(pl_st)]] <- data.frame(
      TEHC_ST = pl_st,
      pl_gaps[, c("CALYR", "gap")],
      rmspe_pre  = sqrt(mean(pl_gaps$gap[pl_gaps$CALYR %in% pre_years]^2, na.rm = TRUE)),
      rmspe_post = sqrt(mean(pl_gaps$gap[pl_gaps$CALYR %in% post_years]^2, na.rm = TRUE)),
      rmspe_ratio = pl_rmspe_pre,
      stringsAsFactors = FALSE
    )
  }
  pl_df <- if (length(placebo_gaps)) do.call(rbind, placebo_gaps) else data.frame()

  # Permutation p-value: share of units (CA included) with RMSPE ratio >= CA
  ca_ratio <- ratio
  all_ratios <- c(ca_ratio,
                  if (nrow(pl_df)) unique(pl_df$rmspe_ratio[is.finite(pl_df$rmspe_ratio)]) else numeric(0))
  pval <- if (length(all_ratios) > 0 && is.finite(ca_ratio)) {
    mean(all_ratios >= ca_ratio, na.rm = TRUE)
  } else NA_real_
  cat(sprintf("  Permutation p-value: %.4f (%d placebo runs)\n",
              pval, length(placebo_gaps)))

  # ---- Leave-one-out sensitivity ---------------------------------------------
  cat("  Running leave-one-out...\n")
  loo_gaps <- list()
  for (drop_st in donors_ok) {
    loo_donors <- setdiff(donors_ok, drop_st)
    if (length(loo_donors) < 2) next
    loo_fit <- tryCatch(
      fit_scm_weights(SY, treat_state, loo_donors, pre_years, cur_outcomes, cur_preds,
                      n_restarts = 3),
      error = function(e) NULL
    )
    if (is.null(loo_fit)) next
    loo_g <- compute_gap_series(SY, treat_state, loo_donors, loo_fit$weights,
                                outcome_var, all_years_data)
    loo_gaps[[as.character(drop_st)]] <- data.frame(
      dropped_state = drop_st, loo_g[, c("CALYR", "gap")],
      stringsAsFactors = FALSE
    )
  }
  loo_df <- if (length(loo_gaps)) do.call(rbind, loo_gaps) else data.frame()

  # ---- In-time placebo (fake treatment 2017) ---------------------------------
  cat("  Running in-time placebo (fake treat = 2017)...\n")
  fake_pre_end   <- 2016L
  fake_pre_years <- seq.int(pre_start, fake_pre_end)
  itp_fit <- tryCatch(
    fit_scm_weights(SY, treat_state, donors_ok, fake_pre_years,
                    cur_outcomes, cur_preds, n_restarts = 3),
    error = function(e) NULL
  )
  itp_gaps <- if (!is.null(itp_fit)) {
    ig <- compute_gap_series(SY, treat_state, donors_ok, itp_fit$weights,
                             outcome_var, analysis_years)
    ig$fake_post <- as.integer(ig$CALYR >= 2017 & ig$CALYR <= 2019)
    ig
  } else NULL
  if (!is.null(itp_gaps)) {
    fake_rmspe_pre  <- sqrt(mean(itp_gaps$gap[itp_gaps$CALYR %in% fake_pre_years]^2, na.rm = TRUE))
    fake_rmspe_post <- sqrt(mean(itp_gaps$gap[itp_gaps$fake_post == 1]^2, na.rm = TRUE))
    cat(sprintf("  In-time placebo RMSPE pre=%.4f post(2017-19)=%.4f ratio=%.4f\n",
                fake_rmspe_pre, fake_rmspe_post,
                ifelse(fake_rmspe_pre > 0, fake_rmspe_post / fake_rmspe_pre, NA_real_)))
  }

  # ---- 2020 gap report -------------------------------------------------------
  gap_2020 <- gaps$gap[gaps$CALYR == 2020]
  if (length(gap_2020) && is.finite(gap_2020)) {
    cat(sprintf("  2020 gap (COVID year, excluded from optim): %.4f\n", gap_2020))
  }

  # ---- Collect results -------------------------------------------------------
  results[[outcome_var]] <- list(
    weights      = w,
    donors       = donors_ok,
    gaps         = gaps,
    placebo_gaps = pl_df,
    loo_gaps     = loo_df,
    itp_gaps     = itp_gaps,
    rmspe_ratio  = ratio,
    pval         = pval
  )
}

# ---- Section 8: Write CSVs ---------------------------------------------------

cat("\n======================================================================\n")
cat("  WRITING OUTPUTS\n")
cat("======================================================================\n")

for (outcome_var in names(results)) {
  res <- results[[outcome_var]]
  tag <- tolower(outcome_var)

  # Weights
  w_df <- data.frame(
    TEHC_ST = as.integer(names(res$weights)),
    weight  = as.numeric(res$weights),
    stringsAsFactors = FALSE
  )
  w_df <- w_df[order(w_df$weight, decreasing = TRUE), ]
  fname <- file.path(out_dir, sprintf("scm_%s_weights.csv", tag))
  write.csv(w_df, fname, row.names = FALSE)
  cat("Saved:", basename(fname), "\n")

  # Gap series (including 2020)
  fname <- file.path(out_dir, sprintf("scm_%s_gaps.csv", tag))
  write.csv(res$gaps, fname, row.names = FALSE)
  cat("Saved:", basename(fname), "\n")

  # Placebo gaps
  if (nrow(res$placebo_gaps) > 0) {
    fname <- file.path(out_dir, sprintf("scm_%s_placebo_gaps.csv", tag))
    write.csv(res$placebo_gaps, fname, row.names = FALSE)
    cat("Saved:", basename(fname), "\n")
  }

  # Leave-one-out gaps
  if (nrow(res$loo_gaps) > 0) {
    fname <- file.path(out_dir, sprintf("scm_%s_loo_gaps.csv", tag))
    write.csv(res$loo_gaps, fname, row.names = FALSE)
    cat("Saved:", basename(fname), "\n")
  }

  # In-time placebo gaps
  if (!is.null(res$itp_gaps)) {
    fname <- file.path(out_dir, sprintf("scm_%s_intime_placebo.csv", tag))
    write.csv(res$itp_gaps, fname, row.names = FALSE)
    cat("Saved:", basename(fname), "\n")
  }
}

# Summary table
summary_rows <- lapply(names(results), function(ov) {
  res <- results[[ov]]
  data.frame(
    outcome         = ov,
    n_donors        = length(res$donors),
    rmspe_ratio     = res$rmspe_ratio,
    pval_permutation = res$pval,
    top_donor_fips  = if (length(res$weights)) as.integer(names(which.max(res$weights))) else NA_integer_,
    top_donor_weight = if (length(res$weights)) max(res$weights) else NA_real_,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(out_dir, "scm_summary.csv"), row.names = FALSE)
cat("Saved: scm_summary.csv\n")

# ---- Section 9: Plots --------------------------------------------------------

cat("\nGenerating plots...\n")

for (outcome_var in names(results)) {
  res <- results[[outcome_var]]
  tag <- tolower(outcome_var)
  gaps <- res$gaps

  # ── (a) Actual vs synthetic path ─────────────────────────────────────────
  fname_path <- file.path(out_dir, sprintf("scm_%s_path.png", tag))
  tryCatch({
    png(fname_path, width = 900, height = 500)
    yrng <- range(c(gaps$actual, gaps$synthetic), na.rm = TRUE)
    if (!all(is.finite(yrng))) yrng <- c(0, 1)
    pad <- diff(yrng) * 0.12
    yrng <- yrng + c(-pad, pad)

    plot(gaps$CALYR, gaps$actual, type = "b", lwd = 2.5, pch = 16,
         col = "#1f77b4", ylim = yrng,
         xlab = "Year", ylab = sprintf("Mean %s", outcome_var),
         main = sprintf("SCM: Actual vs Synthetic California (%s)", outcome_var))
    lines(gaps$CALYR, gaps$synthetic, type = "b", lwd = 2.5, pch = 17,
          col = "#d62728", lty = 2)

    # Shade excluded years
    excl_yrs <- gaps$CALYR[gaps$excluded_year == 1]
    for (ey in excl_yrs) {
      rect(ey - 0.45, yrng[1], ey + 0.45, yrng[2],
           col = rgb(0.9, 0.9, 0.1, 0.25), border = NA)
    }
    abline(v = post_start - 0.5, lty = 3, lwd = 1.5, col = "gray30")
    abline(v = pre_start  - 0.5, lty = 3, lwd = 1,   col = "gray60")
    legend("topleft",
           legend = c("Actual CA", "Synthetic CA",
                      "Treatment onset",
                      if (length(excl_yrs)) sprintf("Excl. %s", paste(excl_yrs, collapse = ",")) else NULL),
           col = c("#1f77b4", "#d62728", "gray30",
                   if (length(excl_yrs)) rgb(0.9, 0.9, 0.1, 0.6) else NULL),
           lty = c(1, 2, 3, NA),
           pch = c(16, 17, NA, NA),
           fill = c(NA, NA, NA, if (length(excl_yrs)) rgb(0.9, 0.9, 0.1, 0.25) else NULL),
           border = c(NA, NA, NA, if (length(excl_yrs)) "gray60" else NULL),
           bty = "n", cex = 0.85)
    dev.off()
    cat("Saved:", basename(fname_path), "\n")
  }, error = function(e) cat("  Path plot failed:", e$message, "\n"))

  # ── (b) Gap plot with placebo gaps overlaid ───────────────────────────────
  fname_gap <- file.path(out_dir, sprintf("scm_%s_gap_placebo.png", tag))
  tryCatch({
    png(fname_gap, width = 900, height = 500)
    pl_df <- res$placebo_gaps

    # Gather all gaps for y-range
    all_gaps <- gaps$gap
    if (nrow(pl_df) > 0) all_gaps <- c(all_gaps, pl_df$gap)
    grng <- range(all_gaps, na.rm = TRUE)
    if (!all(is.finite(grng))) grng <- c(-0.5, 0.5)
    pad  <- diff(grng) * 0.12
    grng <- grng + c(-pad, pad)

    plot(gaps$CALYR, gaps$gap, type = "n", ylim = grng,
         xlab = "Year", ylab = sprintf("Gap: Actual - Synthetic (%s)", outcome_var),
         main = sprintf("SCM Gap Plot with Placebo Tests (%s)\np=%.3f",
                        outcome_var, res$pval))

    # Placebo lines (thin, gray)
    if (nrow(pl_df) > 0) {
      for (pl_st in unique(pl_df$TEHC_ST)) {
        sub <- pl_df[pl_df$TEHC_ST == pl_st, ]
        lines(sub$CALYR, sub$gap, col = rgb(0.6, 0.6, 0.6, 0.5), lwd = 0.8)
      }
    }

    # Shade excluded years
    excl_yrs <- gaps$CALYR[gaps$excluded_year == 1]
    for (ey in excl_yrs) {
      rect(ey - 0.45, grng[1], ey + 0.45, grng[2],
           col = rgb(0.9, 0.9, 0.1, 0.25), border = NA)
    }

    # CA gap line (bold, blue)
    lines(gaps$CALYR, gaps$gap, type = "b", lwd = 2.5, pch = 16, col = "#1f77b4")
    abline(h = 0, lty = 2, col = "black", lwd = 1)
    abline(v = post_start - 0.5, lty = 3, lwd = 1.5, col = "gray30")
    legend("topleft",
           legend = c("CA gap", "Placebo gaps"),
           col = c("#1f77b4", rgb(0.6, 0.6, 0.6, 0.8)),
           lwd = c(2.5, 0.8), pch = c(16, NA), bty = "n", cex = 0.85)
    dev.off()
    cat("Saved:", basename(fname_gap), "\n")
  }, error = function(e) cat("  Gap/placebo plot failed:", e$message, "\n"))

  # ── (c) Leave-one-out sensitivity ────────────────────────────────────────
  fname_loo <- file.path(out_dir, sprintf("scm_%s_loo.png", tag))
  loo_df <- res$loo_gaps
  if (nrow(loo_df) > 0) {
    tryCatch({
      png(fname_loo, width = 900, height = 500)
      # Show synthetic paths for each LOO run vs actual
      loo_synth <- lapply(unique(loo_df$dropped_state), function(ds) {
        sub <- loo_df[loo_df$dropped_state == ds, ]
        # synthetic = actual - gap
        data.frame(CALYR = sub$CALYR,
                   synthetic_loo = gaps$actual[match(sub$CALYR, gaps$CALYR)] - sub$gap,
                   dropped = ds, stringsAsFactors = FALSE)
      })
      loo_synth_df <- do.call(rbind, loo_synth)

      yrng2 <- range(c(gaps$actual, loo_synth_df$synthetic_loo), na.rm = TRUE)
      if (!all(is.finite(yrng2))) yrng2 <- c(0, 1)
      pad2 <- diff(yrng2) * 0.12
      yrng2 <- yrng2 + c(-pad2, pad2)

      plot(gaps$CALYR, gaps$actual, type = "b", lwd = 2.5, pch = 16,
           col = "#1f77b4", ylim = yrng2,
           xlab = "Year", ylab = sprintf("Mean %s", outcome_var),
           main = sprintf("SCM Leave-One-Out Sensitivity (%s)", outcome_var))

      for (ds in unique(loo_synth_df$dropped)) {
        sub <- loo_synth_df[loo_synth_df$dropped == ds, ]
        lines(sub$CALYR, sub$synthetic_loo, col = rgb(0.8, 0.3, 0.3, 0.4), lwd = 0.9)
      }
      # Main synthetic
      lines(gaps$CALYR, gaps$synthetic, type = "b", lwd = 2, pch = 17,
            col = "#d62728", lty = 2)
      abline(v = post_start - 0.5, lty = 3, lwd = 1.5, col = "gray30")
      legend("topleft",
             legend = c("Actual CA", "Main synthetic", "LOO synthetics"),
             col = c("#1f77b4", "#d62728", rgb(0.8, 0.3, 0.3, 0.6)),
             lwd = c(2.5, 2, 0.9), pch = c(16, 17, NA), bty = "n", cex = 0.85)
      dev.off()
      cat("Saved:", basename(fname_loo), "\n")
    }, error = function(e) cat("  LOO plot failed:", e$message, "\n"))
  }
}

# ---- Final summary -----------------------------------------------------------

cat("\n======================================================================\n")
cat("  SCM FINAL SUMMARY\n")
cat("======================================================================\n")
print(summary_df)
cat("\nDone stage 3b (SCM).\n")
