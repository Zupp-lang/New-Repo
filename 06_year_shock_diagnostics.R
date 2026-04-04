#!/usr/bin/env Rscript
# Stage 6: Year-shock diagnostics (2019 vs 2020 vs 2021)
#
# CHANGED: Gutted from original 06_policy_sensitivity_diagnostics.R.
# Retained only Section 1 (year-shock summary) -- the rest (propensity-score
# feasibility, entropy balancing, overlap trimming impact, IPW stability,
# sensitivity grid) referenced matching/weighting designs that are no longer
# used. The causal strategy is now Abadie-style SCM (see 03b_scm.R).

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
  error = function(e) normalizePath("scripts/06_year_shock_diagnostics.R", winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base <- cfg$BASE_DIR
panel_file <- file.path(base, cfg$ANNUAL_PANEL_FILE)
out_dir <- base
if (!file.exists(panel_file)) stop("Missing annual panel: ", panel_file)

primary_outcome    <- cfg$PRIMARY_OUTCOME %||% "EPAY"
secondary_outcomes <- trimws(strsplit(cfg$SECONDARY_OUTCOMES %||% "TPAYWK,TPAYWK_ASINH,EDAYCARE", ",", fixed = TRUE)[[1]])
outcomes           <- unique(c(primary_outcome, secondary_outcomes))
controls           <- trimws(strsplit(cfg$CONTROLS %||% "ERN,WRTH,N_children,TMWKHRS_SUM_PARENTS,HHMAXEEDUC_HS", ",", fixed = TRUE)[[1]])

D <- read.csv(panel_file, check.names = FALSE, stringsAsFactors = FALSE)

num_vars <- unique(c("TREATED_STATE", "EVENTTIME", "POST", "TREAT", "CALYR", "TEHC_ST",
                     "SSUID", "TREAT_ONSET", "N_children", controls, outcomes, "TPAYWK"))
for (v in intersect(num_vars, names(D))) D[[v]] <- suppressWarnings(as.numeric(D[[v]]))
D <- D[!is.na(D$N_children) & D$N_children > 0, , drop = FALSE]

if (!"TPAYWK_ASINH" %in% names(D) && "TPAYWK" %in% names(D)) {
  tp <- D$TPAYWK
  tp[!is.na(tp) & tp < 0] <- NA_real_
  D$TPAYWK_ASINH <- asinh(tp)
}

mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# -------------------------------------------------------------------
# Year-shock diagnostics (2019 vs 2020 vs 2021)
# Useful for understanding COVID contamination and whether 2020 should
# be excluded from the analysis window.
# -------------------------------------------------------------------
cat("== Year-shock diagnostics: 2019 / 2020 / 2021 ==\n\n")

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

  did_19_21 <- (mean_or_na(y[D$CALYR == 2021 & D$TREATED_STATE == 1]) -
                mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 1])) -
               (mean_or_na(y[D$CALYR == 2021 & D$TREATED_STATE == 0]) -
                mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 0]))

  did_19_20 <- (mean_or_na(y[D$CALYR == 2020 & D$TREATED_STATE == 1]) -
                mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 1])) -
               (mean_or_na(y[D$CALYR == 2020 & D$TREATED_STATE == 0]) -
                mean_or_na(y[D$CALYR == 2019 & D$TREATED_STATE == 0]))

  year_rows[[length(year_rows) + 1]] <- data.frame(
    outcome = yvar,
    treated_state = 9,   # sentinel for DiD contrast row
    mean_2019 = NA_real_,
    mean_2020 = NA_real_,
    mean_2021 = NA_real_,
    delta_2021_vs_2019 = did_19_21,
    delta_2020_vs_2019 = did_19_20,
    delta_2021_vs_2020 = did_19_21 - did_19_20,
    stringsAsFactors = FALSE
  )

  cat(sprintf("%-20s  DiD(2021v2019)=%+.4f  DiD(2020v2019)=%+.4f\n",
              yvar, did_19_21, did_19_20))
}

if (length(year_rows)) {
  year_diag <- do.call(rbind, year_rows)
  year_diag$group <- ifelse(year_diag$treated_state == 0, "control",
                            ifelse(year_diag$treated_state == 1, "treated", "did_contrast"))
  write.csv(year_diag, file.path(out_dir, "sensitivity_year_shock_summary.csv"), row.names = FALSE)
  cat("\nSaved: sensitivity_year_shock_summary.csv\n")
}

cat("\nDone stage 6.\n")
