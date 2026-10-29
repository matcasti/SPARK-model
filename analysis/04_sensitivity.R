# ==============================================================================
# 04_sensitivity.R -- sensitivity of kinetics, variance partition and reconstructed trajectories
# to the physiological hyperparameters; leave-one-prior-out; prior-vs-posterior contrast.
#
# Reads   the CSV recordings (a subset of each dataset: cfg$sensitivity$max_recordings)
# Writes  results/sensitivity.rds   sensitivity (one row per recording x setting), summary,
#                                   prior_posterior
# Method  Reference model M3; each setting is refitted from the engine start values.
#         The steps are in R/steps/sensitivity.R.
# ==============================================================================
source("analysis/00_setup.R")
n_settings <- sum(lengths(SENSITIVITY_SETTINGS)) + 6            # hyperparameters + leave-one-out + M0 + reference

# ---- 1. Analyse a subset of the recordings of every dataset -------------------------------------------------
per_dataset <- list()
for (dataset in cfg$datasets) {
  data <- load_dataset(dataset, cfg)
  ids <- pick_subset(data$index, cfg$sensitivity$max_recordings, cfg$seed)
  log_msg("Sensitivity | %s: %d recordings x %d settings", dataset, length(ids), n_settings)
  jobs <- lapply(ids, function(id) list(ds = dataset, id = id, rr = data$rr[[id]]))
  per_dataset[[dataset]] <- drop_failed(map_recordings(jobs, sensitivity_analyse_recording, cfg$cores))
}
results <- unlist(per_dataset, recursive = FALSE)

# ---- 2. Combine and summarise ---------------------------------------------------------------------------------------
sens <- rbindlist(lapply(results, `[[`, "sens"), fill = TRUE)
prior_posterior <- rbindlist(lapply(results, `[[`, "prior_post"), fill = TRUE)
sens_summary <- sensitivity_summarise(sens)

# ---- 3. Save -------------------------------------------------------------------------------------------------------------
write_result(list(sensitivity = sens, summary = sens_summary, prior_posterior = prior_posterior),
             "sensitivity", cfg)

print(sens_summary[order(kind, setting)], nrows = 60)
print(prior_posterior[, .(med_z_tau_p = median(z_tau_p), med_z_tau_s = median(z_tau_s),
                          frac_in_p = mean(in_prior_tau_p), frac_in_s = mean(in_prior_tau_s),
                          med_contr_p = median(contraction_tau_p), med_contr_s = median(contraction_tau_s)),
                      by = Dataset])
