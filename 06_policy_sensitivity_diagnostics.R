#!/usr/bin/env Rscript
# Stage 6: Pre-estimation diagnostics for potential design fixes
# - Year exclusion effects (e.g., drop 2020)
# - Pscore/prognostic-score feasibility
# - Entropy balancing diagnostics (via optional WeightIt fallback)
# - Overlap trimming impact
# - Weighting aptitude + sensitivity-grid template

rm(list = ls())

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

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
  error = function(e) normalizePath("scripts/06_policy_sensitivity_diagnostics.R", winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base <- cfg$BASE_DIR
panel_file <- file.path(base, cfg$ANNUAL_PANEL_FILE)
out_dir <- base
if (!file.exists(panel_file)) stop("Missing annual panel: ", panel_file)

primary_outcome <- cfg$PRIMARY_OUTCOME %||% "EPAY"
secondary_outcomes <- trimws(strsplit(cfg$SECONDARY_OUTCOMES %||% "TPAYWK,TPAYWK_ASINH,EDAYCARE", ",", fixed = TRUE)[[1]])
outcomes <- unique(c(primary_outcome, secondary_outcomes))
controls <- trimws(strsplit(cfg$CONTROLS %||% "ERN,WRTH,N_children,TMWKHRS_SUM_PARENTS,HHMAXEEDUC_HS", ",", fixed = TRUE)[[1]])

D <- read.csv(panel_file, check.names = FALSE, stringsAsFactors = FALSE)

num_vars <- unique(c("TREATED_STATE", "EVENTTIME", "POST", "TREAT", "CALYR", "TEHC_ST", "SSUID", "TREAT_ONSET", controls, outcomes, "TPAYWK"))
for (v in intersect(num_vars, names(D))) D[[v]] <- suppressWarnings(as.numeric(D[[v]]))
D <- D[!is.na(D$N_children) & D$N_children > 0, , drop = FALSE]

if (!"TPAYWK_ASINH" %in% names(D) && "TPAYWK" %in% names(D)) {
  tp <- D$TPAYWK
  tp[!is.na(tp) & tp < 0] <- NA_real_
  D$TPAYWK_ASINH <- asinh(tp)
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# -------------------------------------------------------------------
# 1) Year-shock diagnostics (2019 vs 2021 and role of 2020)
# -------------------------------------------------------------------
year_rows <- list()
for (yvar in outcomes) {
  if (!(yvar %in% names(D))) next
  y <- as.numeric(D[[yvar]])
  y[!is.na(y) & y < 0] <- NA_real_
  
  for (g in c(0, 1)) {
    m2019 <- mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == g])
    m2020 <- mean_or_na(y[D$CALYR == 2020 & D$TREATED_STATE == g])
    m2021 <- mean_or_na(y[D$CALYR == 2021 & D$TREATED_STATE == g])
    year_rows[[length(year_rows) + 1]] <- data.frame(
      outcome = yvar,
      treated_state = g,
      mean_2019 = m2019,
      mean_2020 = m2020,
      mean_2021 = m2021,
      delta_2021_vs_2019 = m2021 - m2019,
      delta_2020_vs_2019 = m2020 - m2019,
      delta_2021_vs_2020 = m2021 - m2020,
      stringsAsFactors = FALSE
    )
  }
  
  did_19_21 <- (mean_or_na(y[D$CALYR == 2021 & D$TREATED_STATE == 1]) - mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 1])) -
    (mean_or_na(y[D$CALYR == 2021 & D$TREATED_STATE == 0]) - mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 0]))
  did_19_20 <- (mean_or_na(y[D$CALYR == 2020 & D$TREATED_STATE == 1]) - mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 1])) -
    (mean_or_na(y[D$CALYR == 2020 & D$TREATED_STATE == 0]) - mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 0]))
  
  year_rows[[length(year_rows) + 1]] <- data.frame(
    outcome = yvar,
    treated_state = 9,
    mean_2019 = NA_real_,
    mean_2020 = NA_real_,
    mean_2021 = NA_real_,
    delta_2021_vs_2019 = did_19_21,
    delta_2020_vs_2019 = did_19_20,
    delta_2021_vs_2020 = did_19_21 - did_19_20,
    stringsAsFactors = FALSE
  )
}
if (length(year_rows)) {
  year_diag <- do.call(rbind, year_rows)
  year_diag$group <- ifelse(year_diag$treated_state == 0, "control", ifelse(year_diag$treated_state == 1, "treated", "did_contrast"))
  write.csv(year_diag, file.path(out_dir, "sensitivity_year_shock_summary.csv"), row.names = FALSE)
}

# -------------------------------------------------------------------
# 2) Propensity/prognostic score feasibility
# -------------------------------------------------------------------
pre <- D[D$CALYR <= 2019, , drop = FALSE]
ps_covars <- intersect(controls, names(pre))
ps_covars <- ps_covars[sapply(pre[, ps_covars, drop = FALSE], function(z) sum(!is.na(z)) > 10)]

ps_out <- data.frame()
prog_out <- data.frame()
if (length(ps_covars) > 0 && "TREATED_STATE" %in% names(pre)) {
  ps_mask <- complete.cases(pre[, c("TREATED_STATE", ps_covars), drop = FALSE])
  ps_data <- pre[ps_mask, c("TREATED_STATE", ps_covars), drop = FALSE]
  
  if (nrow(ps_data) > 100 && length(unique(ps_data$TREATED_STATE)) == 2) {
    f_ps <- as.formula(paste("TREATED_STATE ~", paste(ps_covars, collapse = " + ")))
    ps_fit <- tryCatch(glm(f_ps, data = ps_data, family = binomial()), error = function(e) NULL)
    
    if (!is.null(ps_fit)) {
      ps <- pmin(pmax(as.numeric(predict(ps_fit, type = "response")), 1e-4), 1 - 1e-4)
      ps_data$ps <- ps
      
      tr <- ps_data$TREATED_STATE == 1
      rng_t <- range(ps[tr], na.rm = TRUE)
      rng_c <- range(ps[!tr], na.rm = TRUE)
      lo <- max(rng_t[1], rng_c[1])
      hi <- min(rng_t[2], rng_c[2])
      in_overlap <- ps >= lo & ps <= hi
      
      lps <- qlogis(ps)
      cal <- 0.2 * sd(lps, na.rm = TRUE)
      nn_dist <- rep(NA_real_, sum(tr))
      ctrl_lps <- lps[!tr]
      if (length(ctrl_lps) > 0) {
        tr_idx <- which(tr)
        for (i in seq_along(tr_idx)) nn_dist[i] <- min(abs(lps[tr_idx[i]] - ctrl_lps), na.rm = TRUE)
      }
      
      ps_out <- data.frame(
        n_total = nrow(ps_data),
        n_treated = sum(tr),
        n_control = sum(!tr),
        overlap_low = lo,
        overlap_high = hi,
        treated_in_overlap = sum(tr & in_overlap),
        control_in_overlap = sum((!tr) & in_overlap),
        treated_matchable_caliper = sum(nn_dist <= cal, na.rm = TRUE),
        treated_unmatchable_caliper = sum(nn_dist > cal, na.rm = TRUE),
        caliper_logit = cal,
        stringsAsFactors = FALSE
      )
      
      # Prognostic score (control-only fit of primary outcome)
      if (primary_outcome %in% names(pre)) {
        y <- as.numeric(pre[[primary_outcome]])
        y[!is.na(y) & y < 0] <- NA_real_
        prog_mask <- complete.cases(pre[, ps_covars, drop = FALSE]) & !is.na(y)
        prog_data <- pre[prog_mask, , drop = FALSE]
        y_prog <- y[prog_mask]
        
        c_idx <- prog_data$TREATED_STATE == 0
        if (sum(c_idx) > 50 && sum(!c_idx) > 20) {
          f_prog <- as.formula(paste("y_prog ~", paste(ps_covars, collapse = " + ")))
          prog_fit <- tryCatch(lm(f_prog, data = prog_data[c_idx, , drop = FALSE]), error = function(e) NULL)
          if (!is.null(prog_fit)) {
            prog_score <- as.numeric(predict(prog_fit, newdata = prog_data))
            tr2 <- prog_data$TREATED_STATE == 1
            ctrl_prog <- prog_score[!tr2]
            tr_prog <- prog_score[tr2]
            nn_prog <- sapply(tr_prog, function(v) min(abs(v - ctrl_prog), na.rm = TRUE))
            cal_prog <- 0.2 * sd(prog_score, na.rm = TRUE)
            prog_out <- data.frame(
              outcome = primary_outcome,
              n_total = nrow(prog_data),
              n_treated = sum(tr2),
              n_control = sum(!tr2),
              sd_progscore = sd(prog_score, na.rm = TRUE),
              caliper_prog = cal_prog,
              treated_matchable_prog = sum(nn_prog <= cal_prog, na.rm = TRUE),
              treated_unmatchable_prog = sum(nn_prog > cal_prog, na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}
write.csv(ps_out, file.path(out_dir, "sensitivity_pscore_feasibility.csv"), row.names = FALSE)
write.csv(prog_out, file.path(out_dir, "sensitivity_progscore_feasibility.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 3) Entropy-balancing diagnostics (optional package-assisted)
# -------------------------------------------------------------------
entropy_diag <- data.frame(method = "none", status = "not_run", note = "Install WeightIt package to compute entropy balancing weights", stringsAsFactors = FALSE)
if (length(ps_covars) > 0 && requireNamespace("WeightIt", quietly = TRUE)) {
  ew_mask <- complete.cases(pre[, c("TREATED_STATE", ps_covars), drop = FALSE])
  ew <- pre[ew_mask, c("TREATED_STATE", ps_covars), drop = FALSE]
  if (nrow(ew) > 100 && length(unique(ew$TREATED_STATE)) == 2) {
    wobj <- tryCatch(WeightIt::weightit(as.formula(paste("TREATED_STATE ~", paste(ps_covars, collapse = " + "))),
                                        data = ew,
                                        method = "ebal",
                                        estimand = "ATT"), error = function(e) NULL)
    if (!is.null(wobj)) {
      w <- as.numeric(wobj$weights)
      tr <- ew$TREATED_STATE == 1
      smd_before <- sapply(ps_covars, function(v) {
        x <- ew[[v]]
        d <- sqrt(0.5 * (var(x[tr], na.rm = TRUE) + var(x[!tr], na.rm = TRUE)))
        if (!is.finite(d) || d == 0) return(NA_real_)
        (mean(x[tr], na.rm = TRUE) - mean(x[!tr], na.rm = TRUE)) / d
      })
      smd_after <- sapply(ps_covars, function(v) {
        x <- ew[[v]]
        m1 <- weighted.mean(x[tr], w[tr], na.rm = TRUE)
        m0 <- weighted.mean(x[!tr], w[!tr], na.rm = TRUE)
        d <- sqrt(0.5 * (var(x[tr], na.rm = TRUE) + var(x[!tr], na.rm = TRUE)))
        if (!is.finite(d) || d == 0) return(NA_real_)
        (m1 - m0) / d
      })
      
      entropy_diag <- data.frame(
        method = "ebal",
        status = "ok",
        covariate = ps_covars,
        abs_smd_before = abs(smd_before),
        abs_smd_after = abs(smd_after),
        stringsAsFactors = FALSE
      )
      
      write.csv(data.frame(weight = w, treated = ew$TREATED_STATE), file.path(out_dir, "sensitivity_entropy_weights_raw.csv"), row.names = FALSE)
    }
  }
}
write.csv(entropy_diag, file.path(out_dir, "sensitivity_entropy_balance_diagnostics.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4) Overlap trimming impact + extreme region detection
# -------------------------------------------------------------------
trim_diag <- data.frame()
if (nrow(ps_out) > 0) {
  ps_mask <- complete.cases(pre[, c("TREATED_STATE", ps_covars), drop = FALSE])
  ps_data <- pre[ps_mask, c("TREATED_STATE", ps_covars), drop = FALSE]
  f_ps <- as.formula(paste("TREATED_STATE ~", paste(ps_covars, collapse = " + ")))
  ps_fit <- tryCatch(glm(f_ps, data = ps_data, family = binomial()), error = function(e) NULL)
  if (!is.null(ps_fit)) {
    ps <- pmin(pmax(as.numeric(predict(ps_fit, type = "response")), 1e-4), 1 - 1e-4)
    tr <- ps_data$TREATED_STATE == 1
    qs <- c(0.01, 0.05, 0.10)
    trim_diag <- do.call(rbind, lapply(qs, function(q) {
      lo <- quantile(ps, q, na.rm = TRUE)
      hi <- quantile(ps, 1 - q, na.rm = TRUE)
      keep <- ps >= lo & ps <= hi
      data.frame(
        trim_quantile = q,
        lower = as.numeric(lo),
        upper = as.numeric(hi),
        kept_total = sum(keep),
        kept_treated = sum(keep & tr),
        kept_control = sum(keep & !tr),
        dropped_total = sum(!keep),
        dropped_treated = sum(!keep & tr),
        dropped_control = sum(!keep & !tr),
        stringsAsFactors = FALSE
      )
    }))
    
    extreme_flag <- data.frame(
      metric = c("ps_lt_0.05", "ps_gt_0.95", "ps_lt_0.10", "ps_gt_0.90"),
      share_total = c(mean(ps < 0.05), mean(ps > 0.95), mean(ps < 0.10), mean(ps > 0.90)),
      share_treated = c(mean(ps[tr] < 0.05), mean(ps[tr] > 0.95), mean(ps[tr] < 0.10), mean(ps[tr] > 0.90)),
      share_control = c(mean(ps[!tr] < 0.05), mean(ps[!tr] > 0.95), mean(ps[!tr] < 0.10), mean(ps[!tr] > 0.90)),
      stringsAsFactors = FALSE
    )
    write.csv(extreme_flag, file.path(out_dir, "sensitivity_overlap_extremes.csv"), row.names = FALSE)
  }
}
write.csv(trim_diag, file.path(out_dir, "sensitivity_overlap_trimming.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 5) Weighting aptitude (IPW stability proxy)
# -------------------------------------------------------------------
weighting_diag <- data.frame()
if (nrow(ps_out) > 0) {
  ps_mask <- complete.cases(pre[, c("TREATED_STATE", ps_covars), drop = FALSE])
  ps_data <- pre[ps_mask, c("TREATED_STATE", ps_covars), drop = FALSE]
  f_ps <- as.formula(paste("TREATED_STATE ~", paste(ps_covars, collapse = " + ")))
  ps_fit <- tryCatch(glm(f_ps, data = ps_data, family = binomial()), error = function(e) NULL)
  if (!is.null(ps_fit)) {
    ps <- pmin(pmax(as.numeric(predict(ps_fit, type = "response")), 1e-4), 1 - 1e-4)
    tr <- ps_data$TREATED_STATE == 1
    p_t <- mean(tr)
    sw <- ifelse(tr, p_t / ps, (1 - p_t) / (1 - ps))
    ess <- (sum(sw)^2) / sum(sw^2)
    
    weighting_diag <- data.frame(
      n = length(sw),
      ess = ess,
      ess_ratio = ess / length(sw),
      w_mean = mean(sw),
      w_sd = sd(sw),
      w_cv = sd(sw) / mean(sw),
      w_p99 = as.numeric(quantile(sw, 0.99, na.rm = TRUE)),
      w_max = max(sw, na.rm = TRUE),
      apt_for_weighting = ifelse((ess / length(sw)) >= 0.5 && max(sw, na.rm = TRUE) < 25, "yes", "caution"),
      stringsAsFactors = FALSE
    )
  }
}
write.csv(weighting_diag, file.path(out_dir, "sensitivity_weighting_aptitude.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 6) Sensitivity-grid design scaffold (ready to merge with stage-3 outputs)
# -------------------------------------------------------------------
methods <- c("baseline", "drop2020", "ps_match", "prog_match", "ps_prog_hybrid", "entropy_balance", "overlap_trim", "ipw", "entropy_plus_trim")
grid <- expand.grid(outcome = outcomes[outcomes %in% names(D)], method = methods, stringsAsFactors = FALSE)
grid$requires_same_units <- grid$method %in% c("ps_match", "prog_match", "ps_prog_hybrid")
grid$requires_overlap <- grid$method %in% c("ps_match", "ps_prog_hybrid", "ipw", "overlap_trim", "entropy_plus_trim")
grid$requires_weights_stable <- grid$method %in% c("ipw", "entropy_balance", "entropy_plus_trim")

if (nrow(ps_out) > 0) {
  grid$ps_matchable_treated <- ps_out$treated_matchable_caliper[1]
  grid$ps_unmatchable_treated <- ps_out$treated_unmatchable_caliper[1]
  grid$ps_overlap_treated <- ps_out$treated_in_overlap[1]
}
if (nrow(prog_out) > 0) {
  grid$prog_matchable_treated <- prog_out$treated_matchable_prog[1]
  grid$prog_unmatchable_treated <- prog_out$treated_unmatchable_prog[1]
}
if (nrow(weighting_diag) > 0) {
  grid$weighting_aptitude <- weighting_diag$apt_for_weighting[1]
  grid$weighting_ess_ratio <- weighting_diag$ess_ratio[1]
}

grid$feasibility_flag <- "CHECK"
if ("ps_matchable_treated" %in% names(grid)) {
  grid$feasibility_flag[grid$method %in% c("ps_match", "ps_prog_hybrid") & grid$ps_matchable_treated < 25] <- "LOW_MATCH_SUPPORT"
}
if ("prog_matchable_treated" %in% names(grid)) {
  grid$feasibility_flag[grid$method %in% c("prog_match", "ps_prog_hybrid") & grid$prog_matchable_treated < 25] <- "LOW_MATCH_SUPPORT"
}
if ("weighting_aptitude" %in% names(grid)) {
  grid$feasibility_flag[grid$method %in% c("ipw", "entropy_balance", "entropy_plus_trim") & grid$weighting_aptitude == "caution"] <- "WEIGHT_INSTABILITY_RISK"
}

write.csv(grid, file.path(out_dir, "sensitivity_grid_design.csv"), row.names = FALSE)

cat("Saved policy sensitivity diagnostics to:", out_dir, "\n")
cat("- sensitivity_year_shock_summary.csv\n")
cat("- sensitivity_pscore_feasibility.csv\n")
cat("- sensitivity_progscore_feasibility.csv\n")
cat("- sensitivity_entropy_balance_diagnostics.csv\n")
cat("- sensitivity_overlap_trimming.csv\n")
cat("- sensitivity_overlap_extremes.csv\n")
cat("- sensitivity_weighting_aptitude.csv\n")
cat("- sensitivity_grid_design.csv\n")
