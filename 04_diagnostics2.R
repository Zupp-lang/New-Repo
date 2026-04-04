#!/usr/bin/env Rscript
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
  error = function(e) normalizePath("scripts/04_diagnostics.R", winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base <- cfg$BASE_DIR
file_panel <- file.path(base, cfg$ANNUAL_PANEL_FILE)
out_folder <- base
cat("[config] panel:", file_panel, "\n")
if (!file.exists(file_panel)) stop("Missing annual panel: ", file_panel)

# New-model aligned config
exclude_years <- integer(0)
if (!is.null(cfg$EXCLUDE_YEARS) && nzchar(cfg$EXCLUDE_YEARS)) {
  exclude_years <- as.integer(strsplit(cfg$EXCLUDE_YEARS, ",", fixed = TRUE)[[1]])
  exclude_years <- exclude_years[!is.na(exclude_years)]
}
scm_treat_state <- as.integer(cfg$SCM_TREAT_STATE %||% "6")
scm_pre_start   <- as.integer(cfg$SCM_PRE_START %||% "2014")
scm_pre_end     <- as.integer(cfg$SCM_PRE_END %||% "2019")
scm_min_hh      <- as.integer(cfg$SCM_MIN_HH_PER_STATE_YEAR %||% "5")

covars_balance <- c(
  "PARENT_AGE_MEAN", "KDAG", "KD2", "KD3", "N_children", "LK3", "UN3", "DIS3",
  "TMWKHRS_SUM_PARENTS", "ERN", "WRTH", "HHMAXEEDUC_HS", "HHMAXEEDUC_POSTSEC", "HHMAXEEDUC_MAUP"
)
outcomes_plot <- unique(c(cfg$PRIMARY_OUTCOME %||% "EPAY", "TPAYWK", "EDAYCARE"))
binary_outcomes <- c("EPAY", "EDAYCARE", "EFAM", "EGRAN", "ENREL", "ELIST", "EWORKMORE", "ETIMELOST")
controls_base <- trimws(strsplit(cfg$CONTROLS %||% "ERN,WRTH,N_children,TMWKHRS_SUM_PARENTS,HHMAXEEDUC_HS", ",", fixed = TRUE)[[1]])

smd <- function(xt, xc) {
  mt <- mean(xt, na.rm = TRUE); mc <- mean(xc, na.rm = TRUE)
  vt <- var(xt, na.rm = TRUE); vc <- var(xc, na.rm = TRUE)
  d <- sqrt(0.5 * (vt + vc))
  if (!is.finite(d) || d == 0) NA_real_ else (mt - mc) / d
}
recode_binary_12 <- function(x) {
  x <- as.numeric(x)
  ifelse(x == 1, 1, ifelse(x == 2, 0, NA_real_))
}

model_ready_mask <- function(D, yvar) {
  y <- as.numeric(D[[yvar]])
  y[!is.na(y) & y < 0] <- NA_real_
  required <- unique(c("SSUID", "CALYR", "TEHC_ST", "TREATED_STATE", controls_base, yvar))
  required <- intersect(required, names(D))
  base_mask <- complete.cases(D[, required, drop = FALSE])
  base_mask & !is.na(y)
}

D <- read.csv(file_panel, check.names = FALSE, stringsAsFactors = FALSE)
for (v in unique(c("TREATED_STATE", "EVENTTIME", "CALYR", "N_children", "TEHC_ST", "TREAT_ONSET",
                   covars_balance, outcomes_plot, binary_outcomes, controls_base, "TPAYWK"))) {
  if (v %in% names(D)) D[[v]] <- as.numeric(D[[v]])
}
D <- D[!is.na(D$N_children) & D$N_children > 0, , drop = FALSE]
if (length(exclude_years)) D <- D[!(D$CALYR %in% exclude_years), , drop = FALSE]

for (v in intersect(binary_outcomes, names(D))) D[[v]] <- recode_binary_12(D[[v]])

if (!"TREAT_ONSET" %in% names(D) || all(is.na(D$TREAT_ONSET[D$TREATED_STATE == 1]))) {
  earliest_onset <- min(D$CALYR[D$TREATED_STATE == 1], na.rm = TRUE)
} else {
  earliest_onset <- min(D$TREAT_ONSET[D$TREATED_STATE == 1], na.rm = TRUE)
}
D$PRE_COMMON <- as.numeric(D$CALYR < earliest_onset)
D_pre <- D[D$PRE_COMMON == 1, , drop = FALSE]

# 1) Balance table (single consolidated CSV)
covars_present <- intersect(covars_balance, names(D_pre))
if (!length(covars_present)) {
  balance_tbl <- data.frame()
} else {
  balance_tbl <- do.call(rbind, lapply(covars_present, function(v) {
    xt <- D_pre[[v]][D_pre$TREATED_STATE == 1]
    xc <- D_pre[[v]][D_pre$TREATED_STATE == 0]
    data.frame(
      Variable = v,
      Mean_Treat = round(mean(xt, na.rm = TRUE), 4),
      Mean_Ctrl = round(mean(xc, na.rm = TRUE), 4),
      SMD = round(smd(xt, xc), 4),
      N_Treat = sum(!is.na(xt)),
      N_Ctrl = sum(!is.na(xc)),
      stringsAsFactors = FALSE
    )
  }))
}

# 2) Trends and non-missing (all outcomes in one file)
trend_rows <- list(); k <- 1L
for (yvar in outcomes_plot) {
  if (!(yvar %in% names(D))) next
  y <- as.numeric(D[[yvar]])
  if (yvar %in% binary_outcomes) y <- recode_binary_12(y) else y[!is.na(y) & y < 0] <- NA_real_
  
  for (yr in sort(unique(D$CALYR))) {
    mt <- mean(y[D$CALYR == yr & D$TREATED_STATE == 1], na.rm = TRUE)
    mc <- mean(y[D$CALYR == yr & D$TREATED_STATE == 0], na.rm = TRUE)
    trend_rows[[k]] <- data.frame(
      outcome = yvar, CALYR = yr,
      mean_treat = ifelse(is.nan(mt), NA_real_, mt),
      mean_ctrl = ifelse(is.nan(mc), NA_real_, mc),
      nonmiss_treat = sum(!is.na(y[D$CALYR == yr & D$TREATED_STATE == 1])),
      nonmiss_ctrl = sum(!is.na(y[D$CALYR == yr & D$TREATED_STATE == 0])),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
trends_tbl <- if (length(trend_rows)) do.call(rbind, trend_rows) else data.frame()

# 3) Sample counts by year
yr_list <- sort(unique(D$CALYR))
counts_tbl <- do.call(rbind, lapply(yr_list, function(yr) {
  data.frame(
    CALYR = yr,
    N_treat = sum(D$CALYR == yr & D$TREATED_STATE == 1, na.rm = TRUE),
    N_ctrl = sum(D$CALYR == yr & D$TREATED_STATE == 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

# 4) TPAYWK extensive-margin summary
ext_tbl <- data.frame()
if ("TPAYWK" %in% names(D)) {
  tp <- as.numeric(D$TPAYWK); tp[tp < 0] <- NA_real_
  ext_tbl <- do.call(rbind, lapply(yr_list, function(yr) {
    data.frame(
      CALYR = yr,
      share_pos_treat = mean(tp[D$CALYR == yr & D$TREATED_STATE == 1] > 0, na.rm = TRUE),
      share_pos_ctrl = mean(tp[D$CALYR == yr & D$TREATED_STATE == 0] > 0, na.rm = TRUE),
      nonmiss_tpaywk_treat = sum(!is.na(tp[D$CALYR == yr & D$TREATED_STATE == 1])),
      nonmiss_tpaywk_ctrl = sum(!is.na(tp[D$CALYR == yr & D$TREATED_STATE == 0])),
      stringsAsFactors = FALSE
    )
  }))
}

# 5) Model-ready counts + by-year + by-state + SCM donor feasibility in one file
tracked_outcomes <- unique(c(cfg$PRIMARY_OUTCOME %||% "EPAY", trimws(strsplit(cfg$SECONDARY_OUTCOMES %||% "", ",", fixed = TRUE)[[1]])))
tracked_outcomes <- tracked_outcomes[nzchar(tracked_outcomes) & tracked_outcomes %in% names(D)]

ready_rows <- list(); r <- 1L
for (yvar in tracked_outcomes) {
  m <- model_ready_mask(D, yvar)
  ready_rows[[r]] <- data.frame(section = "overall", outcome = yvar, key1 = "all", key2 = "all",
                                value1 = sum(m), value2 = sum(m & D$TREATED_STATE == 1, na.rm = TRUE), value3 = sum(m & D$TREATED_STATE == 0, na.rm = TRUE),
                                stringsAsFactors = FALSE); r <- r + 1L
  for (yr in sort(unique(D$CALYR))) {
    ready_rows[[r]] <- data.frame(section = "by_year", outcome = yvar, key1 = as.character(yr), key2 = "all",
                                  value1 = sum(m & D$CALYR == yr), value2 = sum(m & D$CALYR == yr & D$TREATED_STATE == 1, na.rm = TRUE), value3 = sum(m & D$CALYR == yr & D$TREATED_STATE == 0, na.rm = TRUE),
                                  stringsAsFactors = FALSE); r <- r + 1L
  }
}

if ("TEHC_ST" %in% names(D)) {
  y0 <- if (length(tracked_outcomes)) tracked_outcomes[1] else NA_character_
  if (!is.na(y0)) {
    m0 <- model_ready_mask(D, y0)
    states <- sort(unique(D$TEHC_ST[!is.na(D$TEHC_ST)]))
    for (st in states) {
      ready_rows[[r]] <- data.frame(section = "by_state", outcome = y0, key1 = as.character(st), key2 = as.character(as.integer(any(D$TREATED_STATE[D$TEHC_ST == st] == 1, na.rm = TRUE))),
                                    value1 = sum(m0 & D$TEHC_ST == st, na.rm = TRUE), value2 = NA_real_, value3 = NA_real_,
                                    stringsAsFactors = FALSE); r <- r + 1L
    }
  }
}

# SCM donor feasibility diagnostics
pre_years <- seq.int(scm_pre_start, scm_pre_end)
if ("TEHC_ST" %in% names(D)) {
  states <- sort(unique(D$TEHC_ST[!is.na(D$TEHC_ST)]))
  for (st in states) {
    by_yr <- sapply(pre_years, function(yr) length(unique(as.character(D$SSUID[D$TEHC_ST == st & D$CALYR == yr]))))
    all_years_present <- all(pre_years %in% D$CALYR[D$TEHC_ST == st])
    donor_ok <- as.integer(st != scm_treat_state && all_years_present && all(by_yr >= scm_min_hh))
    ready_rows[[r]] <- data.frame(section = "scm_donor_check", outcome = "SCM", key1 = as.character(st), key2 = paste0("treat=", as.integer(st == scm_treat_state)),
                                  value1 = donor_ok, value2 = sum(by_yr, na.rm = TRUE), value3 = min(by_yr, na.rm = TRUE),
                                  stringsAsFactors = FALSE); r <- r + 1L
  }
}
model_ready_tbl <- if (length(ready_rows)) do.call(rbind, ready_rows) else data.frame()

# 6) Single multi-page trend plot file
plot_file <- file.path(out_folder, "diag_parallel_trends.pdf")
pdf(plot_file, width = 10, height = 6)
if (nrow(trends_tbl)) {
  for (yvar in unique(trends_tbl$outcome)) {
    DD <- trends_tbl[trends_tbl$outcome == yvar, , drop = FALSE]
    yrng <- range(c(DD$mean_treat, DD$mean_ctrl), na.rm = TRUE)
    if (!all(is.finite(yrng))) yrng <- c(0, 1)
    plot(DD$CALYR, DD$mean_treat, type = "b", lwd = 2, pch = 16, col = "#1f77b4",
         ylim = yrng, xlab = "Reference year", ylab = paste("Mean", yvar),
         main = paste("Trend diagnostics:", yvar))
    lines(DD$CALYR, DD$mean_ctrl, type = "b", lwd = 2, pch = 17, col = "#ff7f0e")
    legend("topleft", legend = c("Treated", "Control"), col = c("#1f77b4", "#ff7f0e"),
           lwd = 2, pch = c(16, 17), bty = "n")
    abline(v = earliest_onset - 0.5, lty = 3, col = "gray40")
    grid()
  }
} else {
  plot.new(); title("No trend data available")
}
dev.off()

# Write consolidated outputs (<= 15 files total; here 6 files)
write.csv(balance_tbl, file.path(out_folder, "diag_balance_pre.csv"), row.names = FALSE)
write.csv(trends_tbl, file.path(out_folder, "diag_trends_nonmissing.csv"), row.names = FALSE)
write.csv(counts_tbl, file.path(out_folder, "diag_sample_counts_by_year.csv"), row.names = FALSE)
write.csv(ext_tbl, file.path(out_folder, "diag_tpaywk_extensive_margin.csv"), row.names = FALSE)
write.csv(model_ready_tbl, file.path(out_folder, "diag_model_ready_and_scm.csv"), row.names = FALSE)

manifest <- data.frame(
  file = c(
    "diag_balance_pre.csv",
    "diag_trends_nonmissing.csv",
    "diag_sample_counts_by_year.csv",
    "diag_tpaywk_extensive_margin.csv",
    "diag_model_ready_and_scm.csv",
    "diag_parallel_trends.pdf"
  ),
  purpose = c(
    "Pre-period covariate balance (SMD)",
    "Outcome means + non-missing counts by year and treatment",
    "Raw sample size by year and treatment",
    "TPAYWK extensive margin and non-missing counts",
    "Model-ready counts (overall/by-year/by-state) + SCM donor feasibility",
    "All trend plots in one multi-page file"
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(out_folder, "diag_manifest.csv"), row.names = FALSE)

cat("Diagnostics stage complete. Outputs written (7 files total, consolidated).\n")
