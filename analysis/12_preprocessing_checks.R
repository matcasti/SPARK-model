# ==============================================================================
# 12_preprocessing_checks.R -- checks of the data-preparation choices on simulated recordings with
# known truth, written through emulations of the earlier ("legacy") and the current tools/prep-data.R.
#
# Reads   nothing (simulated recordings, seeded)
# Writes  results/prep_checks.rds   per_recording, summary, cleaning, cleaning_summary
# Checks  (a) join artifacts of the earlier script: repeated rri values, lost alignment of rri with time
#         (b) fits from the legacy rri column, from the RR rebuilt from its `time`, and from files
#             written by the current script (Fantasia-like and PRCP-like)
#         (c) effect of the cleaning gate on beats, SDNN and RMSSD (current 3.5 MAD x 2 vs earlier 2 MAD x 5)
#         The checks are in R/steps/checks.R.
# ==============================================================================
source("analysis/00_setup.R")
n_simulations <- cfg$simulation$n_prep
directory <- tempfile()
dir.create(directory)

# ---- 1. Simulate, write both formats, read back and fit ---------------------------------------------------------------
log_msg("Preprocessing checks on %d simulated recordings per kind", n_simulations)
results <- drop_failed(map_recordings(seq_len(n_simulations), prep_check_one_simulation, cfg$cores,
                                      cfg = cfg, directory = directory), "simulation")
per_recording <- rbindlist(lapply(results, `[[`, "src"))
cleaning <- rbindlist(lapply(results, `[[`, "clean"))

# ---- 2. Summarise -----------------------------------------------------------------------------------------------------------
summary <- per_recording[, .(
  n_recordings = .N, consistency = median(consistency), dup_frac = median(dup_frac),
  true_dup_frac = median(true_dup_frac), sd_rr_ms = median(sd_rr_ms), true_sd_rr_ms = median(true_sd_rr_ms),
  jitter_ms = median(jitter_ms), tau_p = median(tau_p), true_tau_p = median(true_tau_p),
  abs_ln_err_tau_p = median(abs(log(tau_p / true_tau_p)))), by = .(kind, script, src)]
cleaning_summary <- cleaning[, .(frac_changed = mean(frac_changed), rmssd_ratio = mean(rmssd_ratio),
                                 sdnn_ratio = mean(sdnn_ratio)), by = .(group, n_mad, n_iter)]

# ---- 3. Save --------------------------------------------------------------------------------------------------------------------
write_result(list(per_recording = per_recording, summary = summary, cleaning = cleaning,
                  cleaning_summary = cleaning_summary), "prep_checks", cfg)
print(summary)
print(cleaning_summary)
