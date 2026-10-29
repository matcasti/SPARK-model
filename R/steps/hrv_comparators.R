# ==============================================================================
# R/steps/hrv_comparators.R -- conventional HRV (mean HR, SDNN, RMSSD, LF/HF) as comparators of the latent
# summaries. Called by analysis/validation/10_hrv_comparators.R. Descriptive only: no claim of superiority.
#
#   hrv_recording_table()       latent summaries + HRV of every recording of a dataset
#   hrv_phase_table()           HRV in the rest/baseline and exercise/tilt windows
#   hrv_latent_associations()   Spearman correlations latent vs HRV, per dataset
#   hrv_loo_r2_table()          how well HRV (linear) predicts each latent summary, per dataset
#   hrv_response_associations() phase-to-phase change of HRV vs latent response
# ==============================================================================

HRV_METRICS <- c("mean_hr", "SDNN", "RMSSD", "LF_HF")
LATENT_FOR_COMPARISON <- c("tau_p", "tau_s", "lambda_P", "sigma2", "sd_x_p", "sd_x_s", "share_p")

# Labels of the two phases compared in each perturbation dataset
PHASE_LABELS <- list(prcp = c("baseline", "tilt"), exercise = c("rest", "exercise"))


# ---- Per recording -------------------------------------------------------------------------------------
latent_summary_row <- function(id, fits) {
  params <- fits$params[[id]]
  states <- fits$states[SubjectID == id]
  var_p <- var(states$X_p)
  var_s <- var(states$X_s)
  data.table(SubjectID = id, tau_p = 1 / params$kp, tau_s = 1 / params$ks, lambda_P = params$lp,
             sigma2 = params$sig2, sd_x_p = sqrt(var_p), sd_x_s = sqrt(var_s), share_p = var_p / (var_p + var_s))
}

#' One row per recording of `dataset`: latent summaries merged with conventional HRV.
hrv_recording_table <- function(dataset, fits, cfg) {
  ids <- fits$summary$SubjectID
  latent <- rbindlist(lapply(ids, latent_summary_row, fits = fits))
  hrv <- rbindlist(lapply(ids, function(id) {
    cbind(SubjectID = id, hrv_metrics(fits$states[SubjectID == id, RR_Interval], cfg$hrv))
  }))
  cbind(Dataset = dataset, merge(latent, hrv, by = "SubjectID"))
}


# ---- Per phase -----------------------------------------------------------------------------------------------
#' Events with a defined onset: SubjectID, t_onset, t_end (the end of the tilt or of exercise; NA if none) and,
#' for tilts, t_base_from (start of the event's own baseline). One row per recording for exercise, one per
#' analysed tilt for PRCP. NULL for datasets without perturbation events.
phase_events <- function(dataset, cfg) {
  if (dataset == "prcp") {
    tilts <- read_result("prcp", cfg, FALSE)$tilts
    if (is.null(tilts)) return(NULL)
    events <- as.data.table(tilts)[analysed == TRUE, .(SubjectID, t_onset, t_end = t_tilt_end, t_base_from)]
  } else if (dataset == "exercise") {
    info <- read_result("exercise", cfg, FALSE)$info
    if (is.null(info)) return(NULL)
    events <- as.data.table(info)[usable == TRUE, .(SubjectID, t_onset, t_end)]
  } else {
    return(NULL)
  }
  events[is.finite(t_onset)]
}

#' HRV in the reference window (before the onset) and the perturbation window of each event; a recording
#' with several tilts contributes the median over its tilts, so that it counts once.
hrv_phase_table <- function(dataset, states, events, cfg) {
  labels <- PHASE_LABELS[[dataset]]
  rows <- list()
  for (i in seq_len(nrow(events))) {
    id <- events$SubjectID[i]
    recording <- states[SubjectID == id]
    onset <- events$t_onset[i]
    end <- if (is.finite(events$t_end[i])) events$t_end[i] else Inf
    base_from <- if ("t_base_from" %in% names(events)) events$t_base_from[i] else cfg$exercise$rest_skip_s
    windows <- list(
      recording[Time >= base_from & Time < onset - 5, RR_Interval],
      recording[Time >= onset + 10 & Time < end, RR_Interval]
    )
    for (k in 1:2) {
      if (length(windows[[k]]) < 30) next
      rows[[length(rows) + 1L]] <- cbind(data.table(Dataset = dataset, SubjectID = id, phase = labels[k]),
                                         hrv_metrics(windows[[k]], cfg$hrv))
    }
  }
  if (!length(rows)) return(list())
  list(rbindlist(rows)[, lapply(.SD, median, na.rm = TRUE), by = .(Dataset, SubjectID, phase)])
}


# ---- Across recordings ------------------------------------------------------------------------------------------------
#' Spearman correlation of each latent summary (log) with each HRV metric (log for LF/HF), per dataset.
hrv_latent_associations <- function(per_recording) {
  per_dataset <- function(dataset) {
    x <- per_recording[Dataset == dataset]
    if (nrow(x) < 5) return(NULL)
    pairs <- CJ(latent = LATENT_FOR_COMPARISON, hrv = HRV_METRICS, sorted = FALSE)
    rows <- lapply(seq_len(nrow(pairs)), function(i) {
      a <- log(x[[pairs$latent[i]]])
      b <- if (pairs$hrv[i] == "LF_HF") log(x[[pairs$hrv[i]]]) else x[[pairs$hrv[i]]]
      finite <- is.finite(a) & is.finite(b)
      if (sum(finite) < 5) return(NULL)
      data.table(Dataset = dataset, latent = pairs$latent[i], hrv = pairs$hrv[i], n = sum(finite),
                 spearman = suppressWarnings(cor(a[finite], b[finite], method = "spearman")))
    })
    rbindlist(rows)
  }
  rbindlist(lapply(unique(per_recording$Dataset), per_dataset))
}

#' Leave-one-out R^2 of (log SDNN, log RMSSD, mean HR) predicting each latent summary, per dataset
#' (datasets with fewer than 8 recordings are skipped).
hrv_loo_r2_table <- function(per_recording) {
  per_dataset <- function(dataset) {
    x <- per_recording[Dataset == dataset]
    if (nrow(x) < 8) return(NULL)
    predictors <- data.frame(SDNN = log(x$SDNN), RMSSD = log(x$RMSSD), mean_hr = x$mean_hr)
    rbindlist(lapply(LATENT_FOR_COMPARISON, function(latent) {
      data.table(Dataset = dataset, latent = latent, loo_r2 = loo_r2(log(x[[latent]]), predictors))
    }))
  }
  rbindlist(lapply(unique(per_recording$Dataset), per_dataset))
}

#' Spearman correlation between the HRV change (perturbation - reference window) and the latent
#' response, per dataset. NULL when no dataset has enough recordings with both.
hrv_response_associations <- function(per_phase, cfg) {
  if (is.null(per_phase)) return(NULL)
  wide <- dcast(per_phase, Dataset + SubjectID ~ phase, value.var = HRV_METRICS)
  per_dataset <- function(dataset) {
    labels <- PHASE_LABELS[[dataset]]
    results <- read_result(dataset, cfg, FALSE)
    latent_response <- if (is.null(results)) NULL else if (dataset == "prcp") results$response else results$rest_response
    if (is.null(latent_response)) return(NULL)                    # no events were defined for this dataset
    latent_response <- as.data.table(latent_response)[, .(SubjectID, z_d_x_p, z_d_x_s)]

    hrv <- wide[Dataset == dataset]
    if (!nrow(hrv)) return(NULL)
    change <- function(metric) hrv[[paste0(metric, "_", labels[2])]] - hrv[[paste0(metric, "_", labels[1])]]
    hrv_change <- data.table(SubjectID = hrv$SubjectID, d_mean_hr = change("mean_hr"),
                             d_SDNN = change("SDNN"), d_RMSSD = change("RMSSD"))
    merged <- merge(latent_response, hrv_change, by = "SubjectID")
    if (nrow(merged) < 5) return(NULL)

    pairs <- CJ(latent = c("z_d_x_p", "z_d_x_s"), hrv = c("d_mean_hr", "d_SDNN", "d_RMSSD"), sorted = FALSE)
    rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
      data.table(Dataset = dataset, latent = pairs$latent[i], hrv = pairs$hrv[i], n = nrow(merged),
                 spearman = suppressWarnings(cor(merged[[pairs$latent[i]]], merged[[pairs$hrv[i]]],
                                                 method = "spearman", use = "complete.obs")))
    }))
  }
  rbindlist(lapply(c("prcp", "exercise"), per_dataset))
}
