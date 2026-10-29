# ==============================================================================
# R/steps/fantasia.R -- Result 4: young vs older adults (Fantasia), stability of the
# representation, and respiration / blood-pressure context. Called by analysis/datasets/06_fantasia.R.
#
#   fantasia_summary_table()        one row per recording: latent summaries + conventional HRV
#   fantasia_group_contrasts()      older - young, latent and HRV families (BH within family)
#   fantasia_incremental_info()     leave-one-out AUC / R^2 of latent summaries vs HRV
#   fantasia_stability()            test-retest of two non-overlapping windows
#   fantasia_physiology_table()     respiratory coupling of the states, baroreflex alpha, BP
#   fantasia_physiology_analysis()  group contrasts, paired tests and associations on that table
# ==============================================================================

# Latent summaries compared between groups (column names of the summary table)
FANTASIA_LATENT_VARIABLES <- c("Nu_0", "Tau_P", "Tau_S", "lambda_P", "Sigma_sys2",
                               "rr_jitter_ms", "sd_x_p", "sd_x_s", "sd_ratio", "log_kappa_ratio", "share_p")
FANTASIA_HRV_VARIABLES <- c("mean_hr", "SDNN", "RMSSD", "LF_HF")


# ---- Summary table ------------------------------------------------------------------------------
#' Latent summaries of one recording from its fit.
fantasia_latent_row <- function(id, fits) {
  summary <- fits$summary[SubjectID == id]
  params <- fits$params[[id]]
  states <- fits$states[SubjectID == id]
  var_p <- var(states$X_p)
  var_s <- var(states$X_s)
  data.table(
    SubjectID = id, Nu_0 = summary$Nu_0, Tau_P = summary$Tau_P, Tau_S = summary$Tau_S,
    lambda_P = params$lp, Sigma_sys2 = summary$Sigma_sys2,
    rr_jitter_ms = params$rr_jitter_ms, sd_x_p = sqrt(var_p), sd_x_s = sqrt(var_s),
    sd_ratio = sqrt(var_p / var_s), log_kappa_ratio = log(params$kp / params$ks),
    share_p = var_p / (var_p + var_s)
  )
}

#' One row per recording: group, latent summaries, conventional HRV (and log HRV for the models).
fantasia_summary_table <- function(fits, rr_by_recording, cfg) {
  ids <- fits$index$SubjectID
  latent <- merge(fits$index[, .(SubjectID, Group)],
                  rbindlist(lapply(ids, fantasia_latent_row, fits = fits)), by = "SubjectID")
  hrv <- rbindlist(lapply(ids, function(id) cbind(SubjectID = id, hrv_metrics(rr_by_recording[[id]], cfg$hrv))))
  table <- merge(latent, hrv, by = "SubjectID")
  for (name in c("SDNN", "RMSSD", "LF_HF")) table[, paste0("log_", name) := log(get(name))]
  table
}


# ---- Group contrasts and information over HRV ------------------------------------------------------------
#' Older - young contrasts for the latent and the HRV families, BH-adjusted within each family.
fantasia_group_contrasts <- function(table, cfg) {
  contrast_family <- function(variables, family) {
    contrasts <- rbindlist(lapply(variables, function(v) group_contrast(table, v, n_boot = cfg$stats$n_boot)),
                           fill = TRUE)
    contrasts[, `:=`(family = family, p_BH = p.adjust(p, "BH"))]
  }
  rbind(contrast_family(FANTASIA_LATENT_VARIABLES, "latent"),
        contrast_family(FANTASIA_HRV_VARIABLES, "HRV"), fill = TRUE)
}

#' Leave-one-out AUC (young vs older) with HRV, latent or both; leave-one-out R^2 with which
#' HRV predicts each latent summary. Exploratory with small groups.
fantasia_incremental_info <- function(table) {
  is_older <- as.integer(table$Group == "older")
  hrv_predictors <- table[, .(log_SDNN, log_RMSSD, log_LF_HF)]
  latent_predictors <- table[, .(log_tau_p = log(Tau_P), log_tau_s = log(Tau_S),
                                 log_lambda_P = log(lambda_P), log_sigma2 = log(Sigma_sys2))]
  list(
    auc_hrv = loo_auc(is_older, hrv_predictors),
    auc_latent = loo_auc(is_older, latent_predictors),
    auc_both = loo_auc(is_older, cbind(hrv_predictors, latent_predictors)),
    r2_latent_from_hrv = sapply(as.list(latent_predictors), function(v) loo_r2(v, hrv_predictors))
  )
}


# ---- Stability across windows ------------------------------------------------------------------------------
#' Refit one window of a recording; returns its latent summaries and HRV.
fantasia_window_fit <- function(id, window, dataset, cfg) {
  rec <- window_rec(dataset$rec_full[[id]], start_s = window[1], end_s = window[2])
  if (length(rec$rr) < cfg$io$min_beats) stop("window too short: ", id)
  fit <- fit_recording(rec, id, "fantasia", spark_priors(), spark_control())
  cbind(
    data.table(SubjectID = id, lambda_P = fit$lp, sd_x_p = sd(fit$states$X_p), sd_x_s = sd(fit$states$X_s)),
    fit$summary[, .(Nu_0, Tau_P, Tau_S, Sigma_sys2)],
    hrv_metrics(rec$rr, cfg$hrv)[, .(SDNN, RMSSD, LF_HF)]
  )
}

#' Test-retest of window A vs window B (Spearman and ICC of the log values); NULL if fewer than
#' 4 recordings have both windows.
fantasia_stability <- function(ids, dataset, cfg) {
  windows <- cfg$windows$fantasia
  log_msg("Stability: refitting windows A=%s s and B=%s s",
          paste(windows$stability_A, collapse = "-"), paste(windows$stability_B, collapse = "-"))
  fits_a <- map_recordings(ids, fantasia_window_fit, cfg$cores, window = windows$stability_A, dataset = dataset, cfg = cfg)
  fits_b <- map_recordings(ids, fantasia_window_fit, cfg$cores, window = windows$stability_B, dataset = dataset, cfg = cfg)
  both_ok <- !(vapply(fits_a, is_err, logical(1)) | vapply(fits_b, is_err, logical(1)))
  if (sum(both_ok) < 4) {
    log_msg("Stability analysis skipped: fewer than 4 recordings with both windows.")
    return(NULL)
  }
  a <- rbindlist(fits_a[both_ok])
  b <- rbindlist(fits_b[both_ok])

  agreement <- function(variable) {
    log_a <- log(a[[variable]]); log_b <- log(b[[variable]])
    finite <- is.finite(log_a) & is.finite(log_b)
    data.table(variable = variable, n = sum(finite),
               spearman = suppressWarnings(cor(log_a[finite], log_b[finite], method = "spearman")),
               icc_log = icc_consistency(log_a[finite], log_b[finite]),
               rel_change_median = median(abs(b[[variable]][finite] - a[[variable]][finite]) / abs(a[[variable]][finite])))
  }
  stability <- rbindlist(lapply(setdiff(names(a), "SubjectID"), agreement))
  stability[, family := ifelse(variable %in% c("SDNN", "RMSSD", "LF_HF"), "HRV", "latent")]
  stability[]
}


# ---- Respiration and blood pressure ---------------------------------------------------------------------------
#' Respiratory frequency and the respiratory-band content of x_p, x_s and RR (recordings with respiration).
respiratory_coupling_row <- function(states, cfg) {
  physio <- cfg$physio
  f_resp <- resp_frequency(states$Time, states$resp, physio$resp_band)
  coupling <- function(series) {
    resp_coupling(states$Time, series, states$resp, f_resp[["f_resp"]], physio$coh_halfwidth_hz)
  }
  x_p <- coupling(states$X_p); x_s <- coupling(states$X_s); rr <- coupling(states$RR_Interval)
  row <- data.table(
    f_resp = f_resp[["f_resp"]], resp_peak_frac = f_resp[["peak_frac"]],
    resp_in_HF = f_resp[["f_resp"]] >= cfg$hrv$hf[1] & f_resp[["f_resp"]] <= cfg$hrv$hf[2],
    coh_x_p = x_p[["coh"]], coh_x_s = x_s[["coh"]], coh_rr = rr[["coh"]],
    resp_sd_x_p = x_p[["band_sd"]], resp_sd_x_s = x_s[["band_sd"]],
    resp_share_x_p = x_p[["share_total"]], resp_share_x_s = x_s[["share_total"]]
  )
  row[, `:=`(coh_diff = coh_x_p - coh_x_s, resp_sd_ratio = resp_sd_x_p / resp_sd_x_s)]
  row
}

#' Blood-pressure summary and baroreflex alpha index (recordings with blood pressure).
blood_pressure_row <- function(states, cfg) {
  map <- if (column_available(states, "map", cfg$io$min_cov_coverage)) states$map else NULL
  summary <- bp_summary(states$systolic, states$diastolic, map)
  alpha <- brs_alpha(states$Time, states$systolic, states$RR_Interval, cfg$physio$brs_bands, cfg$physio$brs_min_coh)
  cbind(summary, as.data.table(as.list(alpha)))
}

fantasia_physiology_row <- function(id, fit_states, cfg) {
  states <- fit_states[SubjectID == id]
  coverage <- cfg$io$min_cov_coverage
  row <- data.table(SubjectID = id)
  if (column_available(states, "resp", coverage)) row <- cbind(row, respiratory_coupling_row(states, cfg))
  if (column_available(states, "systolic", coverage) && column_available(states, "diastolic", coverage)) {
    row <- cbind(row, blood_pressure_row(states, cfg))
  }
  row
}

#' One row per recording with its respiratory / blood-pressure measures (missing where the
#' signal is absent) next to the latent and HRV summaries used in the associations.
fantasia_physiology_table <- function(fits, summary_table, cfg) {
  rows <- rbindlist(lapply(fits$index$SubjectID, fantasia_physiology_row, fit_states = fits$states, cfg = cfg),
                    fill = TRUE)
  table <- merge(fits$index[, .(SubjectID, Group)], rows, by = "SubjectID")
  latent_columns <- summary_table[, .(SubjectID, tau_p = Tau_P, tau_s = Tau_S, lambda_P, sd_x_p, sd_x_s,
                                      share_p, RMSSD, SDNN)]
  merge(table, latent_columns, by = "SubjectID")
}

#' Group contrasts, paired tests (vagal-like vs sympathetic-like state) and associations
#' on the physiology table. Variables absent from the table are skipped.
fantasia_physiology_analysis <- function(table, cfg) {
  present <- function(variables) intersect(variables, names(table))
  contrast_variables <- present(c("f_resp", "resp_sd_x_p", "resp_sd_x_s", "resp_share_x_p", "resp_share_x_s",
                                  "coh_x_p", "coh_x_s", "sbp", "map", "pulse_pressure", "alpha_LF", "alpha_HF"))
  contrasts <- rbindlist(lapply(contrast_variables, function(v) group_contrast(table, v, n_boot = cfg$stats$n_boot)),
                         fill = TRUE)
  if (nrow(contrasts) && "p" %in% names(contrasts)) contrasts[, p_BH := p.adjust(p, "BH")]

  has <- function(column) column %in% names(table)
  tests <- list(
    coh_p_vs_s = if (has("coh_x_p")) paired_comparison(table, "coh_x_p", "coh_x_s"),
    resp_sd_p_vs_s = if (has("resp_sd_x_p")) paired_comparison(table, "resp_sd_x_p", "resp_sd_x_s"),
    resp_share_p_vs_s = if (has("resp_share_x_p")) paired_comparison(table, "resp_share_x_p", "resp_share_x_s"),
    resp_in_HF = if (has("resp_in_HF")) list(
      n = sum(is.finite(table$f_resp)), frac = mean(table$resp_in_HF[is.finite(table$f_resp)]),
      f_resp = median_range(table$f_resp))
  )
  associations <- spearman_table(table, present(c("alpha_LF", "alpha_HF", "resp_sd_x_p", "f_resp")),
                                 c("tau_p", "tau_s", "lambda_P", "sd_x_p", "sd_x_s", "share_p", "RMSSD"))
  list(contrasts = contrasts, tests = tests, associations = associations)
}
