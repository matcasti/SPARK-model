# ==============================================================================
# 03_identifiability.R -- Results 1-2: is the latent decomposition identified by RR alone, and do
# physiological constraints reduce the ambiguity?
#
# Reads   the CSV recordings (a subset of each dataset: cfg$identifiability$max_recordings)
# Writes  results/identifiability.rds   metrics (per recording x specification), starts (every
#                                       fit), hessian, surfaces (profile over tau_p, tau_s),
#                                       trajectories (first recordings), summary
# Method  Dispersed multi-start fits of M0-M3 (gate off): same starting points for all
#         specifications; ambiguity metrics; Hessian of M0 and M3; profile surfaces.
#         The per-recording steps are in R/model/identifiability.R.
# ==============================================================================
source("analysis/00_setup.R")
settings <- cfg$identifiability

# ---- 1. Analyse a subset of the recordings of every dataset ---------------------------------------------------
per_dataset <- list()
for (dataset in cfg$datasets) {
  data <- load_dataset(dataset, cfg)
  ids <- pick_subset(data$index, settings$max_recordings, cfg$seed)
  log_msg("Identifiability | %s: %d of %d recordings, %d starts x %d specs", dataset, length(ids),
          nrow(data$index), settings$n_starts, length(IDENTIFIABILITY_SPECS))
  jobs <- lapply(seq_along(ids), function(i) {
    first_recordings <- i <= settings$n_surface_recordings            # surfaces and trajectories: first ones only
    list(ds = dataset, id = ids[i], rr = data$rr[[ids[i]]], surface = first_recordings, keep_traj = first_recordings)
  })
  per_dataset[[dataset]] <- drop_failed(map_recordings(jobs, identifiability_analyse_recording, cfg$cores, cfg = cfg))
}
results <- unlist(per_dataset, recursive = FALSE)

# ---- 2. Combine --------------------------------------------------------------------------------------------------------
stack <- function(field) {
  rbindlist(unlist(lapply(results, function(r) r[[field]]), recursive = FALSE), fill = TRUE)
}
metrics <- stack("metrics")
starts <- stack("starts")
hessian <- stack("hess")
surfaces <- stack("surf")
trajectories <- unlist(lapply(results, function(r) if (length(r$traj)) list(r$traj) else NULL), recursive = FALSE)
names(trajectories) <- NULL
summary <- identifiability_summarise(metrics)

# ---- 3. Save -------------------------------------------------------------------------------------------------------------
write_result(list(metrics = metrics, starts = starts, hessian = hessian, surfaces = surfaces,
                  trajectories = trajectories, summary = summary), "identifiability", cfg)

print(summary[, .(Dataset, spec, n_ok, n_solutions, n_equiv, sd_ln_tau_p, sd_ln_tau_s, collapse_frac,
                  frac_in_phys_range, traj_cor_p_equiv, best_rmse_ms)])
if (nrow(hessian)) {
  print(hessian[, .(median_cond = median(cond), median_sd_ln_tau_s = median(sd_ln_tau_s)), by = .(Dataset, spec)])
}
