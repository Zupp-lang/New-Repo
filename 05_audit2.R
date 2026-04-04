#!/usr/bin/env Rscript
local({
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
  
  read_env_config <- function(path) {
    if (!file.exists(path)) stop("Missing shared config: ", path)
    x <- trimws(readLines(path, warn = FALSE)); x <- x[nzchar(x) & substr(x,1,1)!="#"]
    out <- list(); for (ln in x) {p <- strsplit(ln, "=", fixed=TRUE)[[1]]; if(length(p)>=2) out[[trimws(p[1])]] <- trimws(paste(p[-1], collapse="="))}; out
  }
  script_path <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE), error = function(e) normalizePath("scripts/05_audit.R", winslash = "/", mustWork = FALSE))
  cfg <- read_env_config(file.path(dirname(script_path), "pipeline_config.env"))
  base <- cfg$BASE_DIR
  f_panel <- file.path(base, cfg$ANNUAL_PANEL_FILE)
  out_summary <- file.path(base, "pipeline_audit_summary.csv")
  out_hh <- file.path(base, "pipeline_audit_hh_support.csv")
  out_event <- file.path(base, "pipeline_audit_event_support.csv")
  out_treated <- file.path(base, "audit_treated_state_support.csv")
  cat("[config] panel:", f_panel, "\n")
  if (!file.exists(f_panel)) stop("Missing annual panel: ", f_panel)
  
  D <- read.csv(f_panel, check.names = FALSE, stringsAsFactors = FALSE)
  for (v in c("TPAYWK","EPAY","EDAYCARE","EVENTTIME","TREATED_STATE","CALYR","TEHC_ST","N_children","TREAT",
              "ERN","WRTH","TMWKHRS_SUM_PARENTS","HHMAXEEDUC_HS","SSUID","TREAT_ONSET")) if (v %in% names(D)) D[[v]] <- suppressWarnings(as.numeric(D[[v]]))
  
  metric <- character(); value <- character(); add <- function(k,v){metric <<- c(metric,k); value <<- c(value, as.character(v))}
  
  add("total_hh_year_rows", nrow(D)); add("unique_ssuid", length(unique(D$SSUID))); add("duplicate_hh_year_keys", sum(duplicated(paste(D$SSUID,D$CALYR,sep="_"))))
  
  controls_base <- trimws(strsplit(cfg$CONTROLS %||% "ERN,WRTH,N_children,TMWKHRS_SUM_PARENTS,HHMAXEEDUC_HS", ",", fixed = TRUE)[[1]])
  fe_vars <- c("SSUID","CALYR","TEHC_ST","TREATED_STATE")
  outcomes <- c("TPAYWK","TPAYWK_ASINH","EPAY","EDAYCARE","EFAM","EGRAN","ENREL","ELIST","EWORKMORE")
  
  audit_rows <- list()
  for (dv in outcomes) {
    if (!(dv %in% names(D))) next
    y <- as.numeric(D[[dv]]); y[!is.na(y) & y < 0] <- NA
    lag_name <- paste0(gsub("_ASINH$","",dv), "_LAG1")
    rhs <- intersect(c(controls_base, lag_name, "EVENTTIME"), names(D))
    cols <- unique(c(rhs, fe_vars))
    mask <- complete.cases(D[, cols, drop = FALSE]) & !is.na(y)
    n_ready <- sum(mask)
    n_t <- sum(mask & D$TREATED_STATE == 1, na.rm = TRUE)
    n_c <- sum(mask & D$TREATED_STATE == 0, na.rm = TRUE)
    add(paste0("ready_rows_", dv), n_ready)
    add(paste0("ready_treated_", dv), n_t)
    add(paste0("ready_control_", dv), n_c)
    
    sample_flag <- if (n_ready == 0 || n_t == 0 || n_c == 0) {
      "FAIL_treated_or_control_missing"
    } else {
      tr_share <- n_t / n_ready
      if (tr_share > 0.95 || tr_share < 0.05) "WARN_near_single_group" else "PASS"
    }
    add(paste0("ready_sample_flag_", dv), sample_flag)
    
    et <- as.numeric(D$EVENTTIME)
    et_ok <- mask & D$TREATED_STATE == 1 & is.finite(et)
    if (any(et_ok)) {
      et_tab <- table(et[et_ok])
      min_cell <- min(as.integer(et_tab))
      add(paste0("min_event_bin_", dv), min_cell)
      add(paste0("event_bins_below5_", dv), sum(as.integer(et_tab) < 5))
    } else {
      add(paste0("min_event_bin_", dv), NA)
      add(paste0("event_bins_below5_", dv), NA)
    }
    
    tr_state <- sort(unique(D$TEHC_ST[D$TREATED_STATE == 1 & mask]))
    if (length(tr_state) > 0) {
      state_counts <- sapply(tr_state, function(st) sum(mask & D$TREATED_STATE == 1 & D$TEHC_ST == st, na.rm = TRUE))
      total_tr <- sum(state_counts)
      shares <- if (total_tr > 0) state_counts / total_tr else rep(NA_real_, length(state_counts))
      for (i in seq_along(tr_state)) {
        audit_rows[[length(audit_rows) + 1]] <- data.frame(
          outcome = dv,
          TEHC_ST = tr_state[i],
          treated_rows = as.integer(state_counts[i]),
          treated_share = as.numeric(shares[i]),
          support_flag = ifelse(state_counts[i] == 0, "FAIL_zero_model_ready", ifelse(shares[i] > 0.8, "WARN_concentrated", "PASS")),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(audit_rows)) {
    treated_support <- do.call(rbind, audit_rows)
    write.csv(treated_support, out_treated, row.names = FALSE)
    fail_zero <- any(treated_support$support_flag == "FAIL_zero_model_ready")
    warn_conc <- any(treated_support$support_flag == "WARN_concentrated")
    add("treated_state_support_fail", fail_zero)
    add("treated_state_support_warn", warn_conc)
  }
  
  if ("TPAYWK" %in% names(D)) {
    tp <- as.numeric(D$TPAYWK); tp[!is.na(tp) & tp < 0] <- NA
    hh <- as.character(D$SSUID)
    yh <- split(tp, hh)
    sd_within <- sapply(yh, function(z){z <- z[!is.na(z)]; if(length(z)<2) NA_real_ else sd(z)})
    nn <- sapply(yh, function(z) sum(!is.na(z)))
    uq <- sapply(yh, function(z) length(unique(z[!is.na(z)])))
    write.csv(data.frame(SSUID=names(nn), outcome_nonmissing=as.integer(nn), outcome_unique_values=as.integer(uq), within_sd=as.numeric(sd_within), estimable_hh=as.integer(nn>=2 & uq>=2)), out_hh, row.names = FALSE)
    
    et <- as.numeric(D$EVENTTIME); tr <- D$TREATED_STATE==1; et_ok <- is.finite(et) & tr
    ev_out <- if (any(et_ok)) {
      tab <- table(et[et_ok]); nn2 <- tapply(!is.na(tp[et_ok]), et[et_ok], sum)
      data.frame(event_bin=as.numeric(names(tab)), hh_year_rows=as.integer(tab), tpaywk_nonmissing=as.integer(nn2[names(tab)]))
    } else data.frame(event_bin=numeric(0), hh_year_rows=integer(0), tpaywk_nonmissing=integer(0))
    write.csv(ev_out[order(ev_out$event_bin),], out_event, row.names = FALSE)
  }
  
  out <- data.frame(metric=metric, value=value, stringsAsFactors = FALSE)
  write.csv(out, out_summary, row.names = FALSE)
  cat("Saved audit outputs:\n", out_summary, "\n", out_hh, "\n", out_event, "\n", out_treated, "\n")
})
