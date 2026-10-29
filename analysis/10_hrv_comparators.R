# ==============================================================================
# 10_hrv_comparators.R -- conventional HRV (SDNN, RMSSD, LF/HF, mean HR) as comparators.
#
# Reads   results/fit_<dataset>.rds    (from 01_fit_all.R)
#         results/prcp.rds, results/exercise.rds   (event times, latent responses; optional)
# Writes  results/hrv.rds   per_recording, per_phase, association, loo_r2, response_association
# Method  Per recording and per detected phase; association with the latent summaries; and
#         how well HRV alone predicts each latent summary (leave-one-out R^2: low values mean the
#         summary is not directly encoded by these metrics). Descriptive; no claim of superiority.
#         The steps are in R/steps/hrv_comparators.R.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. HRV and latent summaries per recording, and HRV per phase --------------------------------------------
per_recording_rows <- list()
phase_rows <- list()
for (dataset in cfg$datasets) {
  fits <- read_result(paste0("fit_", dataset), cfg)
  per_recording_rows[[dataset]] <- hrv_recording_table(dataset, fits, cfg)
  events <- phase_events(dataset, cfg)
  if (!is.null(events)) phase_rows <- c(phase_rows, hrv_phase_table(dataset, fits$states, events, cfg))
}
per_recording <- rbindlist(per_recording_rows, fill = TRUE)
per_phase <- rbindlist(phase_rows, fill = TRUE)
if (!nrow(per_phase)) per_phase <- NULL

# ---- 2. Associations with the latent summaries --------------------------------------------------------------------
association <- hrv_latent_associations(per_recording)
loo_r2_table <- hrv_loo_r2_table(per_recording)
response_association <- hrv_response_associations(per_phase, cfg)
if (!nrow(association)) association <- NULL
if (!nrow(loo_r2_table)) loo_r2_table <- NULL

# ---- 3. Save -------------------------------------------------------------------------------------------------------------------
write_result(list(per_recording = per_recording, per_phase = per_phase, association = association,
                  loo_r2 = loo_r2_table, response_association = response_association), "hrv", cfg)

log_msg("HRV comparators saved (%d recordings, %d phase windows).", nrow(per_recording),
        if (is.null(per_phase)) 0L else nrow(per_phase))
if (!is.null(loo_r2_table)) {
  print(loo_r2_table[, .(median_loo_r2 = round(median(loo_r2, na.rm = TRUE), 2)), by = Dataset])
}
if (!is.null(association)) print(association[latent == "sd_x_p" & hrv == "RMSSD"])
if (!is.null(response_association)) print(response_association)
