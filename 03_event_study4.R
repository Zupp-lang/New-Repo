#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════
# Stage 3: DiD Event Study with Extensive-Margin Primary Outcome
#
# Key changes from original:
#   - EPAY (any childcare spending) as primary outcome
#   - Outcomes configurable from pipeline_config.env
#   - NO lag variables in regression controls
#   - Relaxed sample construction (per-outcome, impute missing controls)
#   - Defensive Wald pre-trend tests
#   - Comprehensive pre-regression diagnostics
# ══════════════════════════════════════════════════════════════════════

rm(list = ls())

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

# ── Section 0: Load config ──────────────────────────────────────────

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
  error = function(e) normalizePath("scripts/03_event_study.R",
                                    winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base       <- cfg$BASE_DIR
file_panel <- file.path(base, cfg$ANNUAL_PANEL_FILE)
out_folder <- base

cat("══════════════════════════════════════════════════════════════\n")
cat("[config] panel:", file_panel, "\n")
if (!file.exists(file_panel)) stop("Panel file missing: ", file_panel)

if (!requireNamespace("fixest", quietly = TRUE)) stop("Install package: fixest")
library(fixest)

# ── Configurable parameters ─────────────────────────────────────────

# Outcomes: read from config, easy to swap
primary_outcome    <- cfg$PRIMARY_OUTCOME
secondary_outcomes <- strsplit(cfg$SECONDARY_OUTCOMES, ",", fixed = TRUE)[[1]]
all_outcomes       <- unique(c(primary_outcome, secondary_outcomes))

# Controls: NO LAG VARIABLES.
controls_base <- strsplit(cfg$CONTROLS, ",", fixed = TRUE)[[1]]

# Balancing / robustness methods
balance_methods <- strsplit(cfg$BALANCE_METHODS %||% "none", ",", fixed = TRUE)[[1]]
balance_covars  <- strsplit(cfg$BALANCE_COVARS %||% cfg$CONTROLS, ",", fixed = TRUE)[[1]]
smd_warn        <- as.numeric(cfg$SMD_WARN_THRESHOLD %||% "0.10")
scm_treat_state <- as.integer(cfg$SCM_TREAT_STATE %||% "6")
scm_pre_start   <- as.integer(cfg$SCM_PRE_START %||% "2014")
scm_pre_end     <- as.integer(cfg$SCM_PRE_END %||% "2019")
scm_min_hh      <- as.integer(cfg$SCM_MIN_HH_PER_STATE_YEAR %||% "5")

# Event study parameters
max_event_abs <- as.integer(cfg$MAX_EVENT_ABS %||% "3")
ref_event     <- as.integer(cfg$REF_EVENT %||% "-1")

# EXCLUDE_YEARS (COVID contamination)
exclude_years <- integer(0)
if (!is.null(cfg$EXCLUDE_YEARS) && nzchar(cfg$EXCLUDE_YEARS)) {
  exclude_years <- as.integer(strsplit(cfg$EXCLUDE_YEARS, ",", fixed = TRUE)[[1]])
  exclude_years <- exclude_years[!is.na(exclude_years)]
}

# Guardrails
MIN_CELL_SIZE <- 5
MIN_CLUSTERS  <- as.integer(cfg$MIN_CLUSTERS %||% "5")
MIN_HTE_N     <- 50

# Identifiers
cluster_var <- "TEHC_ST"
id_var      <- "SSUID"
time_var    <- "CALYR"

RESTRICT_CONTROL_YEARS <- TRUE

cat("[config] primary outcome: ", primary_outcome, "\n")
cat("[config] secondary:       ", paste(secondary_outcomes, collapse = ", "), "\n")
cat("[config] controls:        ", paste(controls_base, collapse = ", "), "\n")
cat("[config] balance methods: ", paste(balance_methods, collapse = ", "), "\n")
cat("[config] max_event_abs:   ", max_event_abs, "\n")
cat("[config] ref_event:       ", ref_event, "\n")
cat("[config] scm treat/pre:   state ", scm_treat_state, ", ", scm_pre_start, "-", scm_pre_end, "\n", sep = "")
if (length(exclude_years)) {
  cat("[config] EXCLUDE_YEARS:   ", paste(exclude_years, collapse = ","), "\n")
} else {
  cat("[config] EXCLUDE_YEARS:    (none)\n")
}
cat("══════════════════════════════════════════════════════════════\n")

# ── Helper functions ────────────────────────────────────────────────

build_formula <- function(lhs, rhs_terms, fe_terms) {
  rhs <- if (length(rhs_terms)) paste(rhs_terms, collapse = " + ") else "1"
  as.formula(sprintf("%s ~ %s | %s", lhs, rhs, paste(fe_terms, collapse = " + ")))
}

# ── RELAXED sample construction ─────────────────────────────────────
# Only require non-missing on: outcome + treatment indicators + cluster.
# Impute missing controls at within-year median.

build_analysis_sample <- function(D, outcome_var, ctrls, cluster_var, id_var, time_var) {
  # Step 1: require non-missing on outcome + structural variables
  required <- c(outcome_var, "TREATED_STATE", "TREAT", cluster_var, id_var, time_var)
  required <- intersect(required, names(D))
  mask <- complete.cases(D[, required, drop = FALSE]) & is.finite(as.numeric(D[[outcome_var]]))
  D_sub <- D[mask, , drop = FALSE]
  
  if (!nrow(D_sub)) return(D_sub)
  
  # Step 2: impute missing controls at within-year median
  for (ctrl in ctrls) {
    if (!(ctrl %in% names(D_sub))) next
    n_miss <- sum(is.na(D_sub[[ctrl]]))
    if (n_miss > 0) {
      # Impute at year-level median
      for (yr in unique(D_sub[[time_var]])) {
        yr_mask <- D_sub[[time_var]] == yr
        med_val <- median(D_sub[[ctrl]][yr_mask], na.rm = TRUE)
        if (is.na(med_val)) med_val <- median(D_sub[[ctrl]], na.rm = TRUE)
        fix_mask <- yr_mask & is.na(D_sub[[ctrl]])
        if (any(fix_mask)) D_sub[[ctrl]][fix_mask] <- med_val
      }
      cat(sprintf("    Imputed %d missing values for %s\n", n_miss, ctrl))
    }
  }
  
  # Step 3: report
  n_cl <- length(unique(D_sub[[cluster_var]]))
  n_tr <- sum(D_sub$TREATED_STATE == 1, na.rm = TRUE)
  n_ct <- sum(D_sub$TREATED_STATE == 0, na.rm = TRUE)
  cat(sprintf("    Sample: N=%d, clusters=%d, treated-state=%d, control=%d\n",
              nrow(D_sub), n_cl, n_tr, n_ct))
  
  return(D_sub)
}

# ── Defensive Wald pre-trend test ───────────────────────────────────

pretrend_ftest_safe <- function(fit, pre_coef_names) {
  pre_coef_names <- intersect(pre_coef_names, names(coef(fit)))
  q <- length(pre_coef_names)
  
  if (q == 0) {
    return(list(F_stat = NA_real_, p_value = NA_real_, df1 = 0L,
                note = "No pre-period coefficients found"))
  }
  if (q == 1) {
    # Single coefficient: just report its t-test
    ct <- coeftable(fit)
    if (pre_coef_names[1] %in% rownames(ct)) {
      pv <- ct[pre_coef_names[1], "Pr(>|t|)"]
      return(list(F_stat = ct[pre_coef_names[1], "t value"]^2, p_value = pv,
                  df1 = 1L, note = "Single pre-period coefficient (t-test squared)"))
    }
  }
  
  # Check residual DoF
  n_obs    <- nobs(fit)
  n_params <- length(coef(fit))
  dof_res  <- n_obs - n_params
  if (dof_res < q) {
    return(list(F_stat = NA_real_, p_value = NA_real_, df1 = q,
                note = sprintf("Insufficient residual DoF (%d) for %d pre-trend coefs", dof_res, q)))
  }
  
  tryCatch({
    wt <- wald(fit, keep = pre_coef_names)
    if (!is.finite(wt$stat) || wt$stat > 1e6) {
      return(list(F_stat = NA_real_, p_value = NA_real_, df1 = q,
                  note = "Numerically unstable"))
    }
    list(F_stat = wt$stat, p_value = wt$p, df1 = q, note = "OK")
  }, error = function(e) {
    list(F_stat = NA_real_, p_value = NA_real_, df1 = q,
         note = paste("wald() failed:", e$message))
  })
}

# ── SMD helper ──────────────────────────────────────────────────────

compute_smd <- function(xt, xc, wt = NULL, wc = NULL) {
  if (is.null(wt)) wt <- rep(1, length(xt))
  if (is.null(wc)) wc <- rep(1, length(xc))
  ok_t <- !is.na(xt) & is.finite(wt); ok_c <- !is.na(xc) & is.finite(wc)
  xt <- xt[ok_t]; wt <- wt[ok_t]; xc <- xc[ok_c]; wc <- wc[ok_c]
  if (length(xt) < 2 || length(xc) < 2) return(NA_real_)
  wt <- wt / sum(wt); wc <- wc / sum(wc)
  mt <- sum(wt * xt); mc <- sum(wc * xc)
  vt <- sum(wt * (xt - mt)^2); vc <- sum(wc * (xc - mc)^2)
  denom <- sqrt(0.5 * (vt + vc))
  if (!is.finite(denom) || denom == 0) NA_real_ else (mt - mc) / denom
}

# ── Synthetic control weights (state-year, Abadie-style) ─────────────

compute_scm_state_weights <- function(D, outcome_vars, predictor_vars,
                                      treat_state = 6,
                                      pre_start = 2014, pre_end = 2019,
                                      min_hh_per_state_year = 5) {
  req <- intersect(c("TEHC_ST", "CALYR", "SSUID", outcome_vars, predictor_vars), names(D))
  X <- D[, req, drop = FALSE]
  
  # Collapse to state-year means + HH counts
  states <- sort(unique(X$TEHC_ST[!is.na(X$TEHC_ST)]))
  years  <- sort(unique(X$CALYR[!is.na(X$CALYR)]))
  sy_rows <- list()
  idx <- 1L
  for (st in states) {
    for (yr in years) {
      m <- X$TEHC_ST == st & X$CALYR == yr
      m[is.na(m)] <- FALSE
      if (!any(m)) next
      row <- data.frame(TEHC_ST = st, CALYR = yr,
                        N_hh = length(unique(as.character(X$SSUID[m]))),
                        stringsAsFactors = FALSE)
      for (v in unique(c(outcome_vars, predictor_vars))) {
        if (v %in% names(X)) row[[v]] <- mean(as.numeric(X[[v]][m]), na.rm = TRUE)
      }
      sy_rows[[idx]] <- row
      idx <- idx + 1L
    }
  }
  if (!length(sy_rows)) return(NULL)
  SY <- do.call(rbind, sy_rows)
  
  pre_years <- seq.int(pre_start, pre_end)
  if (length(pre_years) < 3)
    cat(sprintf("    NOTE: Only %d SCM pre-period years — interpret cautiously.\n", length(pre_years)))
  SY_pre <- SY[SY$CALYR %in% pre_years, , drop = FALSE]
  if (!nrow(SY_pre)) return(NULL)
  
  donor_states <- setdiff(unique(SY_pre$TEHC_ST), treat_state)
  donor_keep <- donor_states[sapply(donor_states, function(st) {
    ss <- SY_pre[SY_pre$TEHC_ST == st, , drop = FALSE]
    all(pre_years %in% ss$CALYR) && all(ss$N_hh[match(pre_years, ss$CALYR)] >= min_hh_per_state_year)
  })]
  
  treat_pre <- SY_pre[SY_pre$TEHC_ST == treat_state, , drop = FALSE]
  if (!nrow(treat_pre) || !all(pre_years %in% treat_pre$CALYR) || !length(donor_keep)) return(NULL)
  
  # Build predictors: outcome lags by year + average predictors over pre-period
  outcome_use <- intersect(outcome_vars, names(SY_pre))
  pred_use <- intersect(predictor_vars, names(SY_pre))
  
  make_vec <- function(state_id) {
    ss <- SY_pre[SY_pre$TEHC_ST == state_id, , drop = FALSE]
    ss <- ss[match(pre_years, ss$CALYR), , drop = FALSE]
    vec <- c()
    for (ov in outcome_use) vec <- c(vec, as.numeric(ss[[ov]]))
    for (pv in pred_use) vec <- c(vec, mean(as.numeric(ss[[pv]]), na.rm = TRUE))
    vec[!is.finite(vec)] <- mean(vec[is.finite(vec)], na.rm = TRUE)
    vec
  }
  
  x_t <- make_vec(treat_state)
  X_d <- sapply(donor_keep, make_vec)
  if (is.null(dim(X_d))) X_d <- matrix(X_d, ncol = 1)
  
  # Minimize ||x_t - X_d w||^2 with simplex weights via softmax parameterization
  obj <- function(theta) {
    exp_t <- exp(theta - max(theta))
    w <- exp_t / sum(exp_t)
    fit <- as.numeric(X_d %*% w)
    sum((x_t - fit)^2)
  }
  opt <- tryCatch(optim(rep(0, ncol(X_d)), obj, method = "BFGS"), error = function(e) NULL)
  if (is.null(opt)) return(NULL)
  exp_t <- exp(opt$par - max(opt$par))
  w <- exp_t / sum(exp_t)
  names(w) <- as.character(donor_keep)
  w
}

# ── Overlap trimming via propensity score ───────────────────────────

compute_overlap_weights <- function(D, covars, trim_alpha = 0.05) {
  # Estimate propensity score, trim outside [alpha, 1-alpha]
  # Returns: logical mask of rows to keep + propensity scores
  covars_use <- intersect(covars, names(D))
  covars_use <- covars_use[sapply(covars_use, function(v) {
    x <- as.numeric(D[[v]]); sum(!is.na(x)) > 10 && sd(x, na.rm = TRUE) > 0
  })]
  if (length(covars_use) == 0) {
    return(list(keep = rep(TRUE, nrow(D)), pscore = rep(0.5, nrow(D))))
  }
  
  cc <- complete.cases(D[, c("TREATED_STATE", covars_use), drop = FALSE])
  if (sum(cc) < 20) {
    return(list(keep = rep(TRUE, nrow(D)), pscore = rep(0.5, nrow(D))))
  }
  
  D_fit <- D[cc, , drop = FALSE]
  fml <- as.formula(paste("TREATED_STATE ~", paste(covars_use, collapse = " + ")))
  fit <- tryCatch(
    glm(fml, data = D_fit, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(keep = rep(TRUE, nrow(D)), pscore = rep(0.5, nrow(D))))
  }
  
  ps <- rep(0.5, nrow(D))
  ps_fit <- tryCatch(predict(fit, newdata = D_fit, type = "response"), error = function(e) NULL)
  if (!is.null(ps_fit)) ps[cc] <- as.numeric(ps_fit)
  ps[is.na(ps)] <- 0.5
  
  keep <- ps >= trim_alpha & ps <= (1 - trim_alpha)
  cat(sprintf("    Overlap trim: dropped %d of %d obs outside [%.2f, %.2f]\n",
              sum(!keep), length(keep), trim_alpha, 1 - trim_alpha))
  list(keep = keep, pscore = ps)
}

# ── Balance diagnostics table builder ───────────────────────────────

build_balance_table <- function(D, covars, w = NULL, method_label = "none") {
  covars_use <- intersect(covars, names(D))
  treat_mask <- D$TREATED_STATE == 1
  rows <- list()
  for (v in covars_use) {
    xt <- as.numeric(D[[v]][treat_mask])
    xc <- as.numeric(D[[v]][!treat_mask])
    wt <- if (!is.null(w)) w[treat_mask] else NULL
    wc <- if (!is.null(w)) w[!treat_mask] else NULL
    smd_val <- compute_smd(xt, xc, wt, wc)
    rows[[length(rows) + 1]] <- data.frame(
      method = method_label, variable = v,
      mean_treat = mean(xt, na.rm = TRUE),
      mean_ctrl_wtd = if (is.null(wc)) mean(xc, na.rm = TRUE) else {
        wc2 <- wc[!is.na(xc)]; xc2 <- xc[!is.na(xc)]
        if (sum(wc2) > 0) sum(wc2 * xc2) / sum(wc2) else NA_real_
      },
      SMD = smd_val,
      flag = ifelse(is.na(smd_val), "NA", ifelse(abs(smd_val) > smd_warn, "WARN", "OK")),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

# ══════════════════════════════════════════════════════════════════════
# Load and prepare data
# ══════════════════════════════════════════════════════════════════════

M <- read.csv(file_panel, check.names = FALSE, stringsAsFactors = FALSE)

# Coerce types
numeric_vars <- unique(c(cluster_var, time_var, "TREATED_STATE", "EVENTTIME",
                         "POST", "TREAT", "N_children", controls_base,
                         balance_covars, "TREAT_ONSET", all_outcomes, "TPAYWK"))
for (v in numeric_vars) {
  if (v %in% names(M)) M[[v]] <- as.numeric(M[[v]])
}
M[[id_var]] <- as.character(M[[id_var]])

# Keep only HH with children
M <- M[!is.na(M$N_children) & M$N_children > 0, , drop = FALSE]
if (!nrow(M)) stop("No rows after child restriction")

# ── Apply EXCLUDE_YEARS ──────────────────────────────────────────────
if (length(exclude_years)) {
  n_before <- nrow(M)
  M <- M[!(M[[time_var]] %in% exclude_years), , drop = FALSE]
  cat(sprintf("[EXCLUDE_YEARS] Dropped %d rows for years: %s\n",
              n_before - nrow(M), paste(exclude_years, collapse = ",")))
  cat(sprintf("[EXCLUDE_YEARS] Remaining: %s\n",
              paste(sort(unique(M[[time_var]])), collapse = ",")))
  if (!nrow(M)) stop("No rows remain after year exclusion.")
}

# Ensure TPAYWK_ASINH exists
if (!("TPAYWK_ASINH" %in% names(M)) && "TPAYWK" %in% names(M)) {
  tp <- M$TPAYWK; tp[tp < 0] <- NA; M$TPAYWK_ASINH <- asinh(tp)
}

# Restrict control-group years to overlap with treated group
if (RESTRICT_CONTROL_YEARS && any(M$TREATED_STATE == 1, na.rm = TRUE)) {
  tr_yrs <- unique(M[[time_var]][M$TREATED_STATE == 1])
  M <- M[M[[time_var]] %in% tr_yrs | M$TREATED_STATE == 1, , drop = FALSE]
}

# ── Event-time bins ─────────────────────────────────────────────────

et <- as.numeric(M$EVENTTIME)
et[!is.na(et) & et < -max_event_abs] <- -max_event_abs
et[!is.na(et) & et >  max_event_abs] <-  max_event_abs
M$EVENT_BIN     <- et
M$EVENT_BIN_USE <- ifelse(M$TREATED_STATE == 1, M$EVENT_BIN, 0)
M$EVENT_BIN_USE[is.na(M$EVENT_BIN_USE)] <- 0

# Compressed event component
M$EVENT_COMP <- ifelse(M$EVENT_BIN_USE <= -2, "pre",
                       ifelse(M$EVENT_BIN_USE == -1, "ref",
                              ifelse(M$EVENT_BIN_USE == 0, "impact", "post")))
M$EVENT_COMP[M$TREATED_STATE == 0] <- "control"
M$EVENT_COMP <- factor(M$EVENT_COMP, levels = c("control", "pre", "ref", "impact", "post"))

# Factor variables
M$SSUID_f  <- factor(M[[id_var]])
M$CALYR_f  <- factor(M[[time_var]])
M$STATE_f  <- factor(M[[cluster_var]])
M$NKID_CAT <- ifelse(M$N_children == 1, "1",
                     ifelse(M$N_children == 2, "2", "3plus"))

# ══════════════════════════════════════════════════════════════════════
# PRE-REGRESSION DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  PRE-REGRESSION DIAGNOSTICS\n")
cat("══════════════════════════════════════════════════════════════\n")

cat("\nTotal rows (children > 0):", nrow(M), "\n")
cat("Unique HH:", length(unique(M[[id_var]])), "\n")
cat("Unique states:", length(unique(M[[cluster_var]])), "\n")
cat("Treated-state rows:", sum(M$TREATED_STATE == 1, na.rm = TRUE), "\n")
cat("Control-state rows:", sum(M$TREATED_STATE == 0, na.rm = TRUE), "\n")

cat("\nOutcome variable summaries:\n")
for (ov in all_outcomes) {
  if (!(ov %in% names(M))) { cat(sprintf("  %-15s: NOT IN DATA\n", ov)); next }
  vals <- as.numeric(M[[ov]])
  n_nm <- sum(!is.na(vals))
  n_pos <- sum(vals > 0, na.rm = TRUE)
  cat(sprintf("  %-15s: %d non-NA (%.1f%%), >0: %d, mean: %.3f\n",
              ov, n_nm, 100 * n_nm / nrow(M), n_pos,
              ifelse(n_nm > 0, mean(vals, na.rm = TRUE), NA_real_)))
}

cat("\nControl variable coverage:\n")
for (cv in controls_base) {
  if (!(cv %in% names(M))) { cat(sprintf("  %-25s: NOT IN DATA\n", cv)); next }
  n_nm <- sum(!is.na(M[[cv]]))
  cat(sprintf("  %-25s: %d non-NA (%.1f%%)\n", cv, n_nm, 100 * n_nm / nrow(M)))
}

cat("\nEvent-bin cell sizes (treated only):\n")
cell_tab <- table(M$EVENT_BIN_USE[M$TREATED_STATE == 1])
print(cell_tab)
write.csv(data.frame(EventBin = as.numeric(names(cell_tab)),
                     N_treated = as.integer(cell_tab)),
          file.path(out_folder, "eventstudy_cell_sizes.csv"), row.names = FALSE)

cat("\nYear × Treatment status:\n")
print(table(Year = M[[time_var]], Treated = M$TREATED_STATE))

# ══════════════════════════════════════════════════════════════════════
# MULTI-METHOD ESTIMATION LOOP
# For each method in BALANCE_METHODS, compute weights, then run
# Simple DiD + Compressed ES + Full ES for every outcome.
# Collect all results into a robustness grid.
# ══════════════════════════════════════════════════════════════════════

robustness_grid    <- list()   # one row per (outcome, method, spec)
all_balance_diag   <- list()   # balance tables across methods
all_coef_tables    <- list()   # event-study coef tables
all_pretrend       <- list()   # pretrend test results
all_did            <- list()   # simple DiD results

for (method in balance_methods) {
  if (method == "entropy_balance") {
    cat("  note: method 'entropy_balance' is deprecated; using synthetic_control instead.\n")
    method <- "synthetic_control"
  }
  
  cat("\n\n")
  cat("################################################################\n")
  cat(sprintf("  METHOD: %s\n", toupper(method)))
  cat("################################################################\n")
  
  # ── Compute weights for this method ─────────────────────────────────
  
  M$w_bal <- 1.0  # default: unit weights
  
  if (method == "synthetic_control") {
    cat("\n  Computing synthetic-control weights on state-year pre-period...\n")
    out_for_scm <- intersect(c(primary_outcome, "TPAYWK"), names(M))
    pred_for_scm <- intersect(c(balance_covars, controls_base), names(M))
# Auto-narrow SCM pre-period to years actually in data
    actual_years <- sort(unique(M_method$CALYR[!is.na(M_method$CALYR)]))
    scm_pre_avail <- intersect(seq.int(scm_pre_start, scm_pre_end), actual_years)
    if (length(scm_pre_avail) < length(seq.int(scm_pre_start, scm_pre_end))) {
      cat(sprintf("    WARNING: SCM pre-period narrowed from %d-%d to {%s}\n",
                  scm_pre_start, scm_pre_end, paste(scm_pre_avail, collapse=",")))
    }
    scm_ps <- if (length(scm_pre_avail) >= 2) min(scm_pre_avail) else scm_pre_start
    scm_pe <- if (length(scm_pre_avail) >= 2) max(scm_pre_avail) else scm_pre_end
    
    scm_w <- compute_scm_state_weights(
      M,
      outcome_vars = out_for_scm,
      predictor_vars = pred_for_scm,
      treat_state = scm_treat_state,
      pre_start = scm_ps,
      pre_end = scm_pe,
      min_hh_per_state_year = scm_min_hh
    )
    
    if (!is.null(scm_w) && length(scm_w)) {
      donor_states <- as.numeric(names(scm_w))
      M$w_bal <- ifelse(M$TEHC_ST == scm_treat_state, 1.0,
                        ifelse(M$TEHC_ST %in% donor_states,
                               as.numeric(scm_w[as.character(M$TEHC_ST)]),
                               NA_real_))
      cat(sprintf("    SCM donor states: %d, min=%.4f, median=%.4f, max=%.4f\n",
                  length(scm_w), min(scm_w), median(scm_w), max(scm_w)))
    } else {
      cat("    SCM weight construction failed; leaving unit weights.\n")
      M$w_bal <- 1.0
    }
    
  } else if (method == "overlap_trim") {
    cat("\n  Computing overlap-trimmed sample...\n")
    ov <- compute_overlap_weights(M, balance_covars)
    M$w_bal <- ifelse(ov$keep, 1.0, NA_real_)
  }
  
  # ── Balance diagnostics for this method ──────────────────────────────
  bal_tbl <- build_balance_table(M[!is.na(M$w_bal), ], balance_covars,
                                 w = M$w_bal[!is.na(M$w_bal)],
                                 method_label = method)
  if (nrow(bal_tbl)) {
    all_balance_diag[[length(all_balance_diag) + 1]] <- bal_tbl
    max_smd <- max(abs(bal_tbl$SMD), na.rm = TRUE)
    cat(sprintf("  Balance: max |SMD| = %.4f %s\n", max_smd,
                ifelse(max_smd > smd_warn, "[WARN > threshold]", "[OK]")))
  }
  
  # ── Working sample for this method ───────────────────────────────────
  M_method <- M[!is.na(M$w_bal) & M$w_bal > 0, , drop = FALSE]
  if (!nrow(M_method)) {
    cat("  No observations after weighting/trimming. Skipping method.\n")
    next
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # Stage 1: Simple DiD (per method)
  # ══════════════════════════════════════════════════════════════════════
  
  cat(sprintf("\n  [Stage 1] Simple DiD (%s)\n", method))
  
  for (yvar in all_outcomes) {
    if (!(yvar %in% names(M_method))) next
    Y <- as.numeric(M_method[[yvar]]); Y[Y < 0] <- NA_real_; M_method$Y_did <- Y
    ctrls <- intersect(controls_base, names(M_method))
    M_sub <- build_analysis_sample(M_method, "Y_did", ctrls, cluster_var, id_var, time_var)
    n_cl <- length(unique(M_sub[[cluster_var]])); n_t <- sum(M_sub$TREAT == 1, na.rm = TRUE)
    if (nrow(M_sub) == 0 || n_cl < 2 || n_t == 0) next
    
    rhs <- c("TREAT", ctrls)
    use_wts <- if (method %in% c("synthetic_control")) TRUE else FALSE
    fit <- tryCatch({
      if (use_wts) {
        feols(build_formula("Y_did", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, weights = ~w_bal, cluster = ~STATE_f, notes = FALSE)
      } else {
        feols(build_formula("Y_did", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, cluster = ~STATE_f, notes = FALSE)
      }
    }, error = function(e) { cat("    feols error:", e$message, "\n"); NULL })
    if (is.null(fit)) next
    
    ct <- coeftable(fit)
    if ("TREAT" %in% rownames(ct)) {
      row <- data.frame(
        outcome = yvar, method = method, spec = "did",
        Coef = ct["TREAT", "Estimate"], SE = ct["TREAT", "Std. Error"],
        tstat = ct["TREAT", "t value"], pval = ct["TREAT", "Pr(>|t|)"],
        N = nobs(fit), N_clust = n_cl, pretrend_p = NA_real_,
        stringsAsFactors = FALSE
      )
      all_did[[paste(yvar, method, sep = "|")]] <- row
      robustness_grid[[length(robustness_grid) + 1]] <- row
      cat(sprintf("    %s DiD: %.4f (SE %.4f, p=%.3f)\n",
                  yvar, row$Coef, row$SE, row$pval))
    }
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # Stage 2: Compressed ES (per method)
  # ══════════════════════════════════════════════════════════════════════
  
  cat(sprintf("\n  [Stage 2] Compressed ES (%s)\n", method))
  
  for (yvar in all_outcomes) {
    if (!(yvar %in% names(M_method))) next
    Y <- as.numeric(M_method[[yvar]]); Y[Y < 0] <- NA; M_method$Y_comp <- Y
    ctrls <- intersect(controls_base, names(M_method))
    M_sub <- build_analysis_sample(M_method, "Y_comp", ctrls, cluster_var, id_var, time_var)
    if (nrow(M_sub) == 0) next
    
    rhs <- c("i(EVENT_COMP, TREATED_STATE, ref = 'ref')", ctrls)
    use_wts <- method == "synthetic_control"
    fit <- tryCatch({
      if (use_wts) {
        feols(build_formula("Y_comp", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, weights = ~w_bal, cluster = ~STATE_f, notes = FALSE)
      } else {
        feols(build_formula("Y_comp", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, cluster = ~STATE_f, notes = FALSE)
      }
    }, error = function(e) NULL)
    if (is.null(fit)) next
    
    ct <- coeftable(fit)
    ev_terms <- grepl("EVENT_COMP", rownames(ct))
    if (any(ev_terms)) {
      out <- data.frame(outcome = yvar, method = method,
                        term = rownames(ct)[ev_terms],
                        Coef = ct[ev_terms, "Estimate"],
                        SE = ct[ev_terms, "Std. Error"],
                        tstat = ct[ev_terms, "t value"],
                        pval = ct[ev_terms, "Pr(>|t|)"],
                        stringsAsFactors = FALSE)
      # Save compressed results tagged by method
      write.csv(out, file.path(out_folder, sprintf("compressed_%s_%s.csv", yvar, method)),
                row.names = FALSE)
    }
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # Stage 3: Full Event Study (per method)
  # ══════════════════════════════════════════════════════════════════════
  
  cat(sprintf("\n  [Stage 3] Full ES (%s)\n", method))
  
  for (yvar in all_outcomes) {
    if (!(yvar %in% names(M_method))) next
    Y <- as.numeric(M_method[[yvar]]); Y[Y < 0] <- NA; M_method$Y_es <- Y
    is_primary <- (yvar == primary_outcome)
    cat(sprintf("\n    --- %s%s ---\n", yvar, ifelse(is_primary, " [PRIMARY]", "")))
    
    ctrls <- intersect(controls_base, names(M_method))
    M_sub <- build_analysis_sample(M_method, "Y_es", ctrls, cluster_var, id_var, time_var)
    n_cl <- length(unique(M_sub[[cluster_var]]))
    n_t <- sum(M_sub$TREATED_STATE == 1, na.rm = TRUE)
    n_c <- sum(M_sub$TREATED_STATE == 0, na.rm = TRUE)
    if (nrow(M_sub) == 0 || n_t == 0 || n_c == 0 || n_cl < 2) {
      cat(sprintf("      SKIPPED: N=%d, treated=%d, control=%d\n", nrow(M_sub), n_t, n_c))
      next
    }
    
    ref_use <- ref_event
    treated_bins <- sort(unique(M_sub$EVENT_BIN_USE[M_sub$TREATED_STATE == 1]))
    treated_bins <- treated_bins[is.finite(treated_bins)]
    if (!(ref_use %in% treated_bins)) {
      cand <- treated_bins[treated_bins < 0]
      ref_use <- if (length(cand)) max(cand) else 0
      cat(sprintf("      note: ref_event=%d absent; using ref=%d for this fit\n", ref_event, ref_use))
    }
    
    rhs <- c(sprintf("i(EVENT_BIN_USE, TREATED_STATE, ref = %d)", ref_use), ctrls)
    use_wts <- method == "synthetic_control"
    fit <- tryCatch({
      if (use_wts) {
        feols(build_formula("Y_es", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, weights = ~w_bal, cluster = ~STATE_f, notes = FALSE)
      } else {
        feols(build_formula("Y_es", rhs, c("SSUID_f", "CALYR_f")),
              data = M_sub, cluster = ~STATE_f, notes = FALSE)
      }
    }, error = function(e) { cat("      feols error:", e$message, "\n"); NULL })
    if (is.null(fit)) next
    
    ct <- coeftable(fit)
    ev_rows <- grepl("EVENT_BIN_USE::", rownames(ct)) & grepl("TREATED_STATE", rownames(ct))
    ev <- ct[ev_rows, , drop = FALSE]
    if (!nrow(ev)) { cat("      No event-time coefficients\n"); next }
    
    bins <- as.numeric(gsub(".*EVENT_BIN_USE::(-?[0-9]+).*", "\\1", rownames(ev)))
    coef_table <- data.frame(
      outcome = yvar, method = method, EventBin = bins,
      Coef = ev[, "Estimate"], SE = ev[, "Std. Error"],
      CI_lo = ev[, "Estimate"] - 1.96 * ev[, "Std. Error"],
      CI_hi = ev[, "Estimate"] + 1.96 * ev[, "Std. Error"],
      tstat = ev[, "t value"], pval = ev[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    coef_table <- coef_table[order(coef_table$EventBin), ]
    key <- paste(yvar, method, sep = "|")
    all_coef_tables[[key]] <- coef_table
    
    write.csv(coef_table,
              file.path(out_folder, sprintf("eventstudy_%s_%s_coefs.csv", yvar, method)),
              row.names = FALSE)
    
    cat("      Coefficients:\n")
    for (r in seq_len(nrow(coef_table))) {
      cat(sprintf("        t=%+d: %.4f (SE %.4f, p=%.3f)\n",
                  coef_table$EventBin[r], coef_table$Coef[r],
                  coef_table$SE[r], coef_table$pval[r]))
    }
    
    # Pre-trend test
    pre_names <- rownames(ev)[bins < ref_event]
    ft <- pretrend_ftest_safe(fit, pre_names)
    all_pretrend[[key]] <- data.frame(
      outcome = yvar, method = method, F_stat = ft$F_stat,
      df1 = ft$df1, p_value = ft$p_value, note = ft$note,
      stringsAsFactors = FALSE
    )
    cat(sprintf("      Pre-trend: F=%.3f, p=%.4f (%s)\n",
                ifelse(is.na(ft$F_stat), NaN, ft$F_stat),
                ifelse(is.na(ft$p_value), NaN, ft$p_value), ft$note))
    
    # Add to robustness grid (use avg post coef as summary)
    post_coefs <- coef_table$Coef[coef_table$EventBin >= 0]
    avg_post <- if (length(post_coefs)) mean(post_coefs, na.rm = TRUE) else NA_real_
    robustness_grid[[length(robustness_grid) + 1]] <- data.frame(
      outcome = yvar, method = method, spec = "event_study",
      Coef = avg_post,
      SE = mean(coef_table$SE[coef_table$EventBin >= 0], na.rm = TRUE),
      tstat = NA_real_, pval = NA_real_,
      N = nobs(fit), N_clust = n_cl, pretrend_p = ft$p_value,
      stringsAsFactors = FALSE
    )
    
    # Coefficient plot (method-tagged)
    tryCatch({
      ylim_range <- range(c(coef_table$CI_lo, coef_table$CI_hi), na.rm = TRUE)
      if (!all(is.finite(ylim_range))) ylim_range <- c(-1, 1)
      pad <- diff(ylim_range) * 0.1
      ylim_range <- ylim_range + c(-pad, pad)
      
      png(file.path(out_folder, sprintf("eventstudy_%s_%s_plot.png", yvar, method)),
          width = 900, height = 500)
      plot(coef_table$EventBin, coef_table$Coef, type = "o", pch = 16,
           ylim = ylim_range, xlab = "Event time (years relative to policy)",
           ylab = sprintf("Coefficient (%s)", yvar),
           main = sprintf("Event Study: %s [%s]", yvar, method))
      arrows(coef_table$EventBin, coef_table$CI_lo,
             coef_table$EventBin, coef_table$CI_hi,
             angle = 90, code = 3, length = 0.05, col = "gray40")
      abline(h = 0, lty = 2, col = "red")
      abline(v = -0.5, lty = 3, col = "blue")
      grid()
      legend("topleft", legend = c("Point estimate", "95% CI", "Zero", "Policy onset"),
             col = c("black", "gray40", "red", "blue"),
             lty = c(1, 1, 2, 3), pch = c(16, NA, NA, NA), bty = "n", cex = 0.8)
      dev.off()
    }, error = function(e) cat(sprintf("      Plot failed: %s\n", e$message)))
  }
  
}  # end method loop

# ══════════════════════════════════════════════════════════════════════
# Export balance diagnostics
# ══════════════════════════════════════════════════════════════════════

if (length(all_balance_diag)) {
  bal_all <- do.call(rbind, all_balance_diag)
  write.csv(bal_all, file.path(out_folder, "eventstudy_balance_diagnostics.csv"),
            row.names = FALSE)
  cat("\nBalance diagnostics saved: eventstudy_balance_diagnostics.csv\n")
}

# ══════════════════════════════════════════════════════════════════════
# Export robustness grid
# ══════════════════════════════════════════════════════════════════════

if (length(robustness_grid)) {
  grid_df <- do.call(rbind, robustness_grid)
  write.csv(grid_df, file.path(out_folder, "eventstudy_robustness_grid.csv"),
            row.names = FALSE)
  cat("Robustness grid saved: eventstudy_robustness_grid.csv\n")
}

# ══════════════════════════════════════════════════════════════════════
# Export combined event-study coefs and pretrend tests
# ══════════════════════════════════════════════════════════════════════

if (length(all_coef_tables)) {
  write.csv(do.call(rbind, all_coef_tables),
            file.path(out_folder, "eventstudy_all_coefs.csv"), row.names = FALSE)
}
if (length(all_pretrend)) {
  write.csv(do.call(rbind, all_pretrend),
            file.path(out_folder, "pretrend_tests.csv"), row.names = FALSE)
}

# ══════════════════════════════════════════════════════════════════════
# HTE: Family size heterogeneity (unweighted baseline only)
# ══════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  [HTE] Family Size Heterogeneity\n")
cat("══════════════════════════════════════════════════════════════\n")

hte_outcome <- primary_outcome
if (hte_outcome %in% names(M)) {
  Y_hte <- as.numeric(M[[hte_outcome]]); Y_hte[Y_hte < 0] <- NA
  M$Y_hte <- Y_hte
  ctrls <- intersect(controls_base, names(M))
  M_sub <- build_analysis_sample(M, "Y_hte", ctrls, cluster_var, id_var, time_var)
  
  if (nrow(M_sub) >= MIN_HTE_N) {
    rhs <- c("TREAT * NKID_CAT", ctrls)
    fit <- tryCatch(
      feols(build_formula("Y_hte", rhs, c("SSUID_f", "CALYR_f")),
            data = M_sub, cluster = ~STATE_f, notes = FALSE),
      error = function(e) { cat("  HTE feols error:", e$message, "\n"); NULL }
    )
    if (!is.null(fit)) {
      ct <- coeftable(fit); tr_rows <- grepl("TREAT", rownames(ct))
      hte_df <- data.frame(variable = rownames(ct)[tr_rows],
                           Coef = ct[tr_rows, "Estimate"],
                           SE = ct[tr_rows, "Std. Error"],
                           tstat = ct[tr_rows, "t value"],
                           pval = ct[tr_rows, "Pr(>|t|)"])
      write.csv(hte_df, file.path(out_folder, "eventstudy_hte_familysize.csv"),
                row.names = FALSE)
      cat("  HTE results:\n"); print(hte_df)
    }
  } else {
    cat(sprintf("  SKIPPED: N=%d < MIN_HTE_N=%d\n", nrow(M_sub), MIN_HTE_N))
  }
}

# ══════════════════════════════════════════════════════════════════════
# Save all results and final summary
# ══════════════════════════════════════════════════════════════════════

saveRDS(list(
  did_results    = all_did,
  coef_tables    = all_coef_tables,
  pretrend       = all_pretrend,
  balance_diag   = all_balance_diag,
  robustness     = robustness_grid
), file.path(out_folder, "eventstudy_all_results.rds"))

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  FINAL SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n")

if (length(robustness_grid)) {
  grid_df <- do.call(rbind, robustness_grid)
  cat("\nRobustness grid (DiD rows):\n")
  did_rows <- grid_df[grid_df$spec == "did", ]
  if (nrow(did_rows)) {
    for (i in seq_len(nrow(did_rows))) {
      cat(sprintf("  %-15s [%-17s]: coef=%7.4f  SE=%6.4f  p=%5.3f  N=%d\n",
                  did_rows$outcome[i], did_rows$method[i],
                  did_rows$Coef[i], did_rows$SE[i],
                  ifelse(is.na(did_rows$pval[i]), NaN, did_rows$pval[i]),
                  did_rows$N[i]))
    }
  }
}

if (length(all_pretrend)) {
  pt_df <- do.call(rbind, all_pretrend)
  cat("\nPre-trend tests:\n")
  for (i in seq_len(nrow(pt_df))) {
    cat(sprintf("  %-15s [%-17s]: F=%6.3f  p=%5.3f  [%s]\n",
                pt_df$outcome[i], pt_df$method[i],
                ifelse(is.na(pt_df$F_stat[i]), NaN, pt_df$F_stat[i]),
                ifelse(is.na(pt_df$p_value[i]), NaN, pt_df$p_value[i]),
                pt_df$note[i]))
  }
}

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Done stage 3.\n")
