#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════
# Stage 2: Build HH-year panel from per-year SIPP extracts
#
# Key changes from original:
#   - Dual eligibility mode (ever_tanf / lifetime_limit) via config
#   - Comprehensive diagnostics printed at end
#   - Lag variables computed but NOT used in regressions
# ══════════════════════════════════════════════════════════════════════

rm(list = ls())

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Section 0: Load config ──────────────────────────────────────────

read_env_config <- function(path) {
  if (!file.exists(path)) stop("Missing shared config: ", path)
  x <- readLines(path, warn = FALSE)
  x <- trimws(x)
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
  error = function(e) normalizePath("scripts/02_transform_annual.R",
                                    winslash = "/", mustWork = FALSE)
)
cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))

base_folder <- cfg$BASE_DIR
years        <- as.integer(strsplit(cfg$YEARS, ",", fixed = TRUE)[[1]])
out_csv      <- file.path(base_folder, cfg$ANNUAL_PANEL_FILE)
clean_tmpl   <- cfg$CLEAN_EXTRACT_TEMPLATE
sample_mode  <- cfg$SAMPLE_MODE %||% "ever_tanf"

# Parse treatment state onsets: "6:2021,9:2021,10:2021,44:2022"
treat_raw <- strsplit(cfg$TREAT_STATE_ONSET, ",", fixed = TRUE)[[1]]
TREAT_FIPS <- as.integer(sub(":.*", "", treat_raw))
TREAT_ONSET <- as.integer(sub(".*:", "", treat_raw))
names(TREAT_ONSET) <- as.character(TREAT_FIPS)

CHILD_AGE_MAX  <- as.integer(cfg$CHILD_AGE_MAX %||% "13")
PARENT_AGE_MIN <- as.integer(cfg$PARENT_AGE_MIN %||% "15")

# Parse EXCLUDE_YEARS (e.g. "2020" or "2020,2021")
exclude_years <- integer(0)
if (!is.null(cfg$EXCLUDE_YEARS) && nzchar(cfg$EXCLUDE_YEARS)) {
  exclude_years <- as.integer(strsplit(cfg$EXCLUDE_YEARS, ",", fixed = TRUE)[[1]])
  exclude_years <- exclude_years[!is.na(exclude_years)]
}

cat("══════════════════════════════════════════════════════════════\n")
cat("[config] base:        ", base_folder, "\n")
cat("[config] years:       ", paste(years, collapse = ","), "\n")
cat("[config] panel:       ", out_csv, "\n")
cat("[config] sample_mode: ", sample_mode, "\n")
cat("[config] treat states:", paste(TREAT_FIPS, TREAT_ONSET, sep = "→", collapse = ", "), "\n")
if (length(exclude_years)) {
  cat("[config] EXCLUDE_YEARS:", paste(exclude_years, collapse = ","), "\n")
} else {
  cat("[config] EXCLUDE_YEARS: (none)\n")
}
cat("══════════════════════════════════════════════════════════════\n")

# ── Constants ───────────────────────────────────────────────────────

BINARY_CHILDCARE <- c(
  "EPAY", "EDAYCARE", "EFAM", "EGRAN", "ENREL", "ELIST", "EWORKMORE",
  "ETIMELOST", "EPAYHELP", "EDAYHS", "EHEADST", "EKIDSELF", "ENUR",
  "ENURHS", "EOTHR", "EPAR", "EPROG", "ESELF", "ESIB15",
  "EWHOPAID1", "EWHOPAID2", "EWHOPAID3", "EWHOPAID4", "EWHOPAID5"
)
CHILDCARE_OUTCOMES <- unique(c("TPAYWK", BINARY_CHILDCARE))

REL_PARENT_PRIMARY <- c(1, 2, 5, 6, 7, 18)
REL_GRANDPARENT    <- c(8)
REL_EXTENDED       <- c(16, 17)

# ── Helper functions ────────────────────────────────────────────────

max_or_na   <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else max(x) }
sum_omit_na <- function(x) if (all(is.na(x))) 0 else sum(x, na.rm = TRUE)
sum_or_na   <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
mean_or_na  <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

combine_tristate <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) -1 else if (any(x == 1)) 1 else if (any(x == 0)) 0 else -1
}

recode_3         <- function(x) ifelse(x == 1, 1, ifelse(x == 2, 0, -1))
recode_binary_12 <- function(x) { x <- as.numeric(x); ifelse(x == 1, 1, ifelse(x == 2, 0, NA_real_)) }
safe_get         <- function(df, nm) if (nm %in% names(df)) as.numeric(df[[nm]]) else rep(NA_real_, nrow(df))

# ══════════════════════════════════════════════════════════════════════
# Pass 1: Identify eligible households
# ══════════════════════════════════════════════════════════════════════

cat("\n[Pass 1] Identifying eligible HH (mode:", sample_mode, ")\n")

hh_tanf  <- new.env(parent = emptyenv())
hh_child <- new.env(parent = emptyenv())

# Also track the broader TANF-receipt flag separately
hh_ever_tanf <- new.env(parent = emptyenv())

for (yr in years) {
  fin <- file.path(base_folder, gsub("\\{year\\}", as.character(yr), clean_tmpl))
  if (!file.exists(fin)) { cat("[pass1] missing", fin, "\n"); next }
  
  T <- read.csv(fin, check.names = FALSE, stringsAsFactors = FALSE,
                colClasses = c(SSUID = "character"))
  if (!nrow(T)) next
  
  ss <- as.character(T$SSUID)
  
  # Lifetime-limit flag (original strict filter)
  tanf_limit <- (as.numeric(T$TTANF_ERSN) == 3) | (as.numeric(T$RTANF_LCYR) == 3)
  tanf_limit[is.na(tanf_limit)] <- FALSE
  
  # Ever-received-TANF flag (broad filter; coded consistently across years)
  tanf_receipt <- as.numeric(T$TTANF_AMT) > 0
  tanf_receipt[is.na(tanf_receipt)] <- FALSE
  
  # Child flag
  child_hit <- !is.na(as.numeric(T$TAGE)) & as.numeric(T$TAGE) <= CHILD_AGE_MAX
  
  for (i in seq_len(nrow(T))) {
    if (isTRUE(tanf_limit[i]))   hh_tanf[[ss[i]]]      <- TRUE
    if (isTRUE(tanf_receipt[i])) hh_ever_tanf[[ss[i]]]  <- TRUE
    if (isTRUE(child_hit[i]))    hh_child[[ss[i]]]      <- TRUE
  }
}

# Report both pools before applying filter
n_limit     <- length(ls(hh_tanf))
n_ever_tanf <- length(ls(hh_ever_tanf))
n_child     <- length(ls(hh_child))

cat("\n[Pass 1 diagnostics]\n")
cat("  HH hitting lifetime limit (TTANF_ERSN=3 | RTANF_LCYR=3):", n_limit, "\n")
cat("  HH with any TANF amount   (TTANF_AMT>0):                ", n_ever_tanf, "\n")
cat("  HH with children <= 13:                                  ", n_child, "\n")

# Apply selected eligibility mode
if (sample_mode == "ever_tanf") {
  tanf_pool <- ls(hh_ever_tanf)
} else if (sample_mode == "lifetime_limit") {
  tanf_pool <- ls(hh_tanf)
} else {
  stop("Unknown SAMPLE_MODE: ", sample_mode)
}

eligible <- intersect(tanf_pool, ls(hh_child))
keep_env <- new.env(parent = emptyenv())
for (id in eligible) keep_env[[id]] <- TRUE

cat("  Eligible HH (", sample_mode, " ∩ children):", length(eligible), "\n")

if (sample_mode == "ever_tanf") {
  # Also report how many of those are lifetime-limit HH (for comparison)
  n_also_limit <- length(intersect(eligible, ls(hh_tanf)))
  cat("  Of which, also hit lifetime limit:", n_also_limit, "\n")
}

if (!length(eligible)) stop("No eligible HH found. Check data paths and variable coding.")

# ══════════════════════════════════════════════════════════════════════
# Pass 2: Build HH-year panel
# ══════════════════════════════════════════════════════════════════════

cat("\n[Pass 2] Building HH-year panel\n")
HH_all <- data.frame()

for (yr in years) {
  fin <- file.path(base_folder, gsub("\\{year\\}", as.character(yr), clean_tmpl))
  if (!file.exists(fin)) next
  
  T <- read.csv(fin, check.names = FALSE, stringsAsFactors = FALSE,
                colClasses = c(SSUID = "character"))
  if (!nrow(T)) next
  
  # Keep only eligible HH
  T <- T[vapply(as.character(T$SSUID),
                function(id) exists(id, envir = keep_env, inherits = FALSE),
                logical(1)), , drop = FALSE]
  if (!nrow(T)) next
  
  T$CALYR_USE <- if ("RHCALYR" %in% names(T)) as.numeric(T$RHCALYR) else yr
  T <- T[!is.na(T$CALYR_USE), , drop = FALSE]
  if (!nrow(T)) next
  
  # ── Person-year intermediate ─────────────────────────────────────
  person_id <- if ("SPUID" %in% names(T)) as.character(T$SPUID) else rep(NA_character_, nrow(T))
  bad_pid   <- is.na(person_id) | !nzchar(person_id)
  if (any(bad_pid) && "PNUM" %in% names(T))
    person_id[bad_pid] <- paste0("PNUM_", as.character(T$PNUM[bad_pid]))
  person_id[is.na(person_id) | !nzchar(person_id)] <-
    paste0("ROW_", which(is.na(person_id) | !nzchar(person_id)))
  T$PERSON_ID <- person_id
  
  T$PYKEY  <- paste(T$SSUID, T$CALYR_USE, T$PERSON_ID, sep = "\x01")
  py_keys  <- unique(T$PYKEY)
  py_id    <- match(T$PYKEY, py_keys)
  
  PY <- data.frame(PYKEY = py_keys, stringsAsFactors = FALSE)
  pparts     <- strsplit(py_keys, "\x01", fixed = TRUE)
  PY$SSUID     <- vapply(pparts, `[`, character(1), 1)
  PY$CALYR     <- as.numeric(vapply(pparts, `[`, character(1), 2))
  PY$PERSON_ID <- vapply(pparts, function(z) paste(z[-c(1, 2)], collapse = "_"), character(1))
  
  rep_vars <- c("TAGE", "ERELRPE", "TPEARN", "TNETWORTH", "TMWKHRS", "EEDUC",
                "ENJ_LKWRK", "ENJ_NOWRK2", "EJOBCANT",
                "EHEARING", "ESEEING", "ECOGNIT", "EAMBULAT", "ESELFCARE",
                "EERRANDS", "TEHC_ST")
  for (v in rep_vars) {
    vals     <- safe_get(T, v)
    PY[[v]]  <- as.numeric(tapply(vals, py_id, max_or_na))
  }
  PY$EEDUC[!is.na(PY$EEDUC) & PY$EEDUC < 0] <- NA_real_
  
  # ── HH-year from person-year ─────────────────────────────────────
  PY$HHKEY <- paste(PY$SSUID, PY$CALYR, sep = "\x01")
  hkeys    <- unique(PY$HHKEY)
  hid      <- match(PY$HHKEY, hkeys)
  
  HH <- data.frame(HHKEY = hkeys, stringsAsFactors = FALSE)
  hparts    <- strsplit(hkeys, "\x01", fixed = TRUE)
  HH$SSUID  <- vapply(hparts, `[`, character(1), 1)
  HH$CALYR  <- as.numeric(vapply(hparts, `[`, character(1), 2))
  HH$TEHC_ST <- as.numeric(tapply(PY$TEHC_ST, hid, function(x) {
    x <- x[!is.na(x)]; if (!length(x)) NA_real_ else as.numeric(names(which.max(table(x))))
  }))
  
  # Children
  age      <- PY$TAGE
  is_child <- !is.na(age) & age <= CHILD_AGE_MAX
  is_adult <- !is.na(age) & age >= PARENT_AGE_MIN
  rel      <- PY$ERELRPE
  is_parent <- is_adult & rel %in% c(REL_PARENT_PRIMARY, REL_GRANDPARENT, REL_EXTENDED)
  
  HH$N_children <- as.numeric(tapply(as.numeric(is_child), hid, sum_omit_na))
  HH$KD2  <- as.numeric(HH$N_children %in% c(1, 2))
  HH$KD3  <- as.numeric(HH$N_children > 2)
  
  child_age <- ifelse(is_child, age, NA_real_)
  youngest  <- as.numeric(tapply(child_age, hid,
                                 function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)))
  HH$KDAG <- ifelse(HH$N_children > 0 & !is.na(youngest), youngest, 0)
  
  # Earnings / wealth
  HH$SUM_TPEARN    <- as.numeric(tapply(PY$TPEARN, hid, sum_omit_na))
  HH$SUM_TNETWORTH <- as.numeric(tapply(PY$TNETWORTH, hid, sum_omit_na))
  HH$ERN  <- asinh(HH$SUM_TPEARN)
  HH$WRTH <- asinh(HH$SUM_TNETWORTH)
  
  # Parent demographics
  HH$N_PARENTS <- as.numeric(tapply(as.numeric(is_parent), hid, sum_omit_na))
  age_parent   <- ifelse(is_parent, age, NA_real_)
  age_adult    <- ifelse(is_adult, age, NA_real_)
  p_age <- as.numeric(tapply(age_parent, hid, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)))
  a_age <- as.numeric(tapply(age_adult, hid, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)))
  HH$PARENT_AGE_MEAN <- ifelse(HH$N_PARENTS > 0, p_age, a_age)
  
  # Hours worked
  hrs_parent <- ifelse(is_parent, PY$TMWKHRS, NA_real_)
  hrs_adult  <- ifelse(is_adult, PY$TMWKHRS, NA_real_)
  HH$TMWKHRS_SUM_PARENTS <- ifelse(
    HH$N_PARENTS > 0,
    as.numeric(tapply(hrs_parent, hid, sum_or_na)),
    as.numeric(tapply(hrs_adult, hid, sum_or_na))
  )
  
  # Job search / disability tristate indicators
  LK3 <- recode_3(PY$ENJ_LKWRK)
  N2  <- recode_3(PY$ENJ_NOWRK2)
  JC  <- recode_3(PY$EJOBCANT)
  UN3 <- ifelse(N2 == 1 | JC == 1, 1, ifelse(N2 == 0 | JC == 0, 0, -1))
  
  DIS3 <- rep(-1, nrow(PY))
  dmat <- cbind(PY$EHEARING, PY$ESEEING, PY$ECOGNIT, PY$EAMBULAT, PY$ESELFCARE, PY$EERRANDS)
  for (i in seq_len(nrow(PY))) {
    v <- dmat[i, ]
    if (any(v == 1, na.rm = TRUE)) DIS3[i] <- 1
    else if (any(v %in% c(1, 2), na.rm = TRUE)) DIS3[i] <- 0
  }
  
  lk_par <- ifelse(is_parent, LK3, NA_real_)
  un_par <- ifelse(is_parent, UN3, NA_real_)
  ds_par <- ifelse(is_parent, DIS3, NA_real_)
  
  HH$LK3  <- ifelse(HH$N_PARENTS > 0,
                    as.numeric(tapply(lk_par, hid, combine_tristate)),
                    as.numeric(tapply(LK3, hid, combine_tristate)))
  HH$UN3  <- ifelse(HH$N_PARENTS > 0,
                    as.numeric(tapply(un_par, hid, combine_tristate)),
                    as.numeric(tapply(UN3, hid, combine_tristate)))
  HH$DIS3 <- ifelse(HH$N_PARENTS > 0,
                    as.numeric(tapply(ds_par, hid, combine_tristate)),
                    as.numeric(tapply(DIS3, hid, combine_tristate)))
  
  # Education
  HH$HHMAXEEDUC <- as.numeric(tapply(PY$EEDUC, hid,
                                     function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
  HH$HHMAXEEDUC_HS      <- as.numeric(!is.na(HH$HHMAXEEDUC) & HH$HHMAXEEDUC >= 39)
  HH$HHMAXEEDUC_POSTSEC <- as.numeric(!is.na(HH$HHMAXEEDUC) & HH$HHMAXEEDUC >= 40 & HH$HHMAXEEDUC <= 43)
  HH$HHMAXEEDUC_MAUP    <- as.numeric(!is.na(HH$HHMAXEEDUC) & HH$HHMAXEEDUC >= 44)
  
  # ── Childcare outcomes from raw T ────────────────────────────────
  HHKEY_raw  <- paste(T$SSUID, T$CALYR_USE, sep = "\x01")
  hh_raw_id  <- match(HHKEY_raw, HH$HHKEY)
  
  for (v in CHILDCARE_OUTCOMES) {
    if (!(v %in% names(T))) { HH[[v]] <- NA_real_; next }
    
    if (v %in% BINARY_CHILDCARE) {
      vals    <- recode_binary_12(as.numeric(T[[v]]))
      HH[[v]] <- as.numeric(tapply(vals, hh_raw_id, max_or_na))
      next
    }
    
    if (v == "TPAYWK") {
      vals <- as.numeric(T[[v]])
      vals[!is.na(vals) & vals < 0] <- NA_real_
      
      # For 2018 and earlier, TPAYWK references a typical week this month;
      # use fall months only for comparability with later years.
      if ("MONTHCODE" %in% names(T)) {
        monthcode <- as.numeric(T$MONTHCODE)
        is_fall_pre2019 <- !is.na(T$CALYR_USE) & T$CALYR_USE <= 2018 & monthcode %in% c(9, 10, 11, 12)
        vals[!is_fall_pre2019 & !is.na(T$CALYR_USE) & T$CALYR_USE <= 2018] <- NA_real_
      }
      
      # Avoid implicit person-weighting: first collapse to one TPAYWK value per
      # HH-year-month, then average across months within HH-year.
      if ("MONTHCODE" %in% names(T)) {
        monthcode <- as.numeric(T$MONTHCODE)
        hh_month_key <- paste(HHKEY_raw, monthcode, sep = "")
        month_keys <- unique(hh_month_key)
        month_id <- match(hh_month_key, month_keys)
        tp_month <- as.numeric(tapply(vals, month_id, max_or_na))
        month_hh <- sub(".*$", "", month_keys)
        month_hh_id <- match(month_hh, HH$HHKEY)
        HH[[v]] <- as.numeric(tapply(tp_month, month_hh_id, mean_or_na))
      } else {
        HH[[v]] <- as.numeric(tapply(vals, hh_raw_id, mean_or_na))
      }
      
      # Make TPAYWK unconditional (extensive + intensive margin).
      if ("EPAY" %in% names(HH)) {
        HH[[v]][!is.na(HH$EPAY) & HH$EPAY == 0] <- 0
      }
      HH[[v]][is.na(HH[[v]])] <- 0
      next
    }
    
    vals <- as.numeric(T[[v]])
    vals[!is.na(vals) & vals < 0] <- NA_real_
    HH[[v]] <- as.numeric(tapply(vals, hh_raw_id, max_or_na))
  }
  
  cat(sprintf("  %d: %d HH-year rows\n", yr, nrow(HH)))
  HH_all <- rbind(HH_all, HH)
}

# ══════════════════════════════════════════════════════════════════════
# Deduplicate and assign treatment
# ══════════════════════════════════════════════════════════════════════

HH_all <- HH_all[!duplicated(HH_all$HHKEY), , drop = FALSE]
if (!nrow(HH_all)) {
  write.csv(HH_all, out_csv, row.names = FALSE, quote = TRUE)
  stop("No HH-year rows after filtering.")
}

# ── Apply EXCLUDE_YEARS filter ─────────────────────────────────────
if (length(exclude_years)) {
  n_before <- nrow(HH_all)
  HH_all <- HH_all[!(HH_all$CALYR %in% exclude_years), , drop = FALSE]
  n_dropped <- n_before - nrow(HH_all)
  cat(sprintf("\n[EXCLUDE_YEARS] Dropped %d HH-year rows for years: %s\n",
              n_dropped, paste(exclude_years, collapse = ",")))
  cat(sprintf("[EXCLUDE_YEARS] Remaining years: %s\n",
              paste(sort(unique(HH_all$CALYR)), collapse = ",")))
  if (!nrow(HH_all)) stop("No rows remain after year exclusion.")
}

# Treatment assignment using config-driven onset mapping
onset <- rep(NA_real_, nrow(HH_all))
for (i in seq_along(TREAT_FIPS)) {
  mask <- HH_all$TEHC_ST == TREAT_FIPS[i]
  mask[is.na(mask)] <- FALSE
  onset[mask] <- TREAT_ONSET[i]
}

HH_all$TREATED_STATE <- as.numeric(!is.na(onset))
HH_all$TREAT_ONSET   <- onset
HH_all$EVENTTIME     <- HH_all$CALYR - onset
HH_all$POST          <- as.numeric(!is.na(onset) & HH_all$CALYR >= onset)
HH_all$TREAT         <- as.numeric(HH_all$TREATED_STATE == 1 & HH_all$POST == 1)

# Sort by HH and year
HH_all <- HH_all[order(HH_all$SSUID, HH_all$CALYR), ]

# ── Lag variables (DESCRIPTIVE ONLY — do NOT use in regressions) ───
# Including lags in regression controls requires consecutive-year
# non-missing values. TPAYWK is reported for only ~8-30% of HH, so
# requiring TPAYWK_LAG1 destroys the sample. HH fixed effects already
# absorb within-unit level differences.

lag_outcomes <- intersect(
  c("TPAYWK", "EPAY", "EDAYCARE", "EFAM", "EGRAN", "ENREL", "ELIST", "EWORKMORE"),
  names(HH_all)
)
for (v in lag_outcomes) {
  lag  <- rep(NA_real_, nrow(HH_all))
  vals <- as.numeric(HH_all[[v]])
  for (i in 2:nrow(HH_all)) {
    if (HH_all$SSUID[i] == HH_all$SSUID[i - 1] &&
        HH_all$CALYR[i] == HH_all$CALYR[i - 1] + 1) {
      lag[i] <- vals[i - 1]
    }
  }
  HH_all[[paste0(v, "_LAG1")]] <- lag
}

# Asinh transform of TPAYWK
if ("TPAYWK" %in% names(HH_all)) {
  tp <- as.numeric(HH_all$TPAYWK)
  tp[!is.na(tp) & tp < 0] <- NA_real_
  HH_all$TPAYWK_ASINH <- asinh(tp)
}

# Drop internal key before writing
HH_all$HHKEY <- NULL

# ══════════════════════════════════════════════════════════════════════
# Write output
# ══════════════════════════════════════════════════════════════════════

write.csv(HH_all, out_csv, row.names = FALSE, quote = TRUE)
cat("\nSaved annual panel:", out_csv, "\n")

# ══════════════════════════════════════════════════════════════════════
# COMPREHENSIVE DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════════════════════════\n")
cat("  DIAGNOSTIC SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n")

cat("\n--- PANEL SIZE ---\n")
cat("  Total HH-year rows:  ", nrow(HH_all), "\n")
cat("  Unique HH (SSUID):   ", length(unique(HH_all$SSUID)), "\n")
cat("  Year range:          ", paste(range(HH_all$CALYR), collapse = " - "), "\n")
cat("  Sample mode:         ", sample_mode, "\n")

cat("\n--- TREATED vs CONTROL ---\n")
cat("  Treated-state HH-years:  ", sum(HH_all$TREATED_STATE == 1, na.rm = TRUE), "\n")
cat("  Treated-state unique HH: ", length(unique(HH_all$SSUID[HH_all$TREATED_STATE == 1])), "\n")
cat("  Control-state HH-years:  ", sum(HH_all$TREATED_STATE == 0, na.rm = TRUE), "\n")
cat("  Control-state unique HH: ", length(unique(HH_all$SSUID[HH_all$TREATED_STATE == 0])), "\n")
cat("  TREAT==1 (post × treated):", sum(HH_all$TREAT == 1, na.rm = TRUE), "\n")

cat("\n--- OUTCOME VARIABLE COVERAGE ---\n")
for (ov in c("EPAY", "TPAYWK", "TPAYWK_ASINH", "EDAYCARE", "EFAM", "EGRAN", "ENREL")) {
  if (ov %in% names(HH_all)) {
    vals  <- as.numeric(HH_all[[ov]])
    n_nm  <- sum(!is.na(vals))
    pct   <- round(100 * n_nm / nrow(HH_all), 1)
    n_pos <- sum(vals > 0, na.rm = TRUE)
    cat(sprintf("  %-15s: %d non-NA (%5.1f%%), >0: %d, mean(non-NA): %.2f\n",
                ov, n_nm, pct, n_pos, mean(vals, na.rm = TRUE)))
  }
}

# Lag coverage (showing why they shouldn't be in regressions)
cat("\n--- LAG VARIABLE COVERAGE (descriptive only) ---\n")
for (lv in paste0(c("TPAYWK", "EPAY"), "_LAG1")) {
  if (lv %in% names(HH_all)) {
    n_nm <- sum(!is.na(HH_all[[lv]]))
    pct  <- round(100 * n_nm / nrow(HH_all), 1)
    cat(sprintf("  %-15s: %d non-NA (%5.1f%%) — would drop %d rows if required\n",
                lv, n_nm, pct, nrow(HH_all) - n_nm))
  }
}

cat("\n--- YEAR × TREATMENT STATUS ---\n")
print(table(Year = HH_all$CALYR, Treated = HH_all$TREATED_STATE))

cat("\n--- STATE DISTRIBUTION (top 20) ---\n")
st_tab <- sort(table(HH_all$TEHC_ST), decreasing = TRUE)
print(head(st_tab, 20))

cat("\n--- TREATED STATE DETAIL ---\n")
for (i in seq_along(TREAT_FIPS)) {
  fips <- TREAT_FIPS[i]
  n_hh <- sum(HH_all$TEHC_ST == fips, na.rm = TRUE)
  n_id <- length(unique(HH_all$SSUID[HH_all$TEHC_ST == fips & !is.na(HH_all$TEHC_ST)]))
  cat(sprintf("  FIPS %2d (onset %d): %d HH-years, %d unique HH\n",
              fips, TREAT_ONSET[i], n_hh, n_id))
}

cat("\n--- CHILDREN DISTRIBUTION ---\n")
print(table(N_children = HH_all$N_children))

# Guardrail warnings
min_hh <- as.integer(cfg$MIN_HH_TOTAL %||% "50")
min_tr <- as.integer(cfg$MIN_HH_TREATED %||% "10")
n_unique_hh <- length(unique(HH_all$SSUID))
n_unique_tr <- length(unique(HH_all$SSUID[HH_all$TREATED_STATE == 1]))

cat("\n--- GUARDRAIL CHECKS ---\n")
if (n_unique_hh < min_hh) {
  cat(sprintf("  ⚠ WARNING: Only %d unique HH (threshold: %d)\n", n_unique_hh, min_hh))
} else {
  cat(sprintf("  OK: %d unique HH >= threshold %d\n", n_unique_hh, min_hh))
}
if (n_unique_tr < min_tr) {
  cat(sprintf("  ⚠ WARNING: Only %d unique treated HH (threshold: %d)\n", n_unique_tr, min_tr))
} else {
  cat(sprintf("  OK: %d unique treated HH >= threshold %d\n", n_unique_tr, min_tr))
}

cat("\n══════════════════════════════════════════════════════════════\n")
cat("Done stage 2.\n")
