# ==============================================================================
# 01_fit_all.R -- fit every recording with M3 (time-scale and separation priors,
# spectral anchor and timing-jitter prior).
#
# Reads   the CSV recordings (cfg$paths$prcp / fantasia / exercise)
# Writes  results/config_used.rds          the configuration of this run
#         results/fit_<dataset>.rds        summary, states, pacf_matrix, params, spectral, index
#         results/df_population_summary.rds, population_summary.csv   one row per recording
# Feeds   fit/innovation diagnostics and kinetics (Results), and every later script.
# ==============================================================================
source("analysis/00_setup.R")

# The configuration used is stored with the results (the manuscript quotes windows and settings from it)
write_result(cfg[setdiff(names(cfg), "is_placeholder")], "config_used", cfg)

# ---- 1. Fit each dataset ----------------------------------------------------------------------------------------
summaries <- list()
for (dataset in cfg$datasets) {
  data <- load_dataset(dataset, cfg)
  fits <- fit_dataset(data$rec, dataset, cfg)
  fits$index <- data$index[SubjectID %in% fits$summary$SubjectID]
  write_result(fits, paste0("fit_", dataset), cfg)
  summaries[[dataset]] <- fits$summary
  log_msg("%s: %d/%d converged | RMSE %s ms | MAPE %s %%", dataset, sum(fits$summary$Converged),
          nrow(fits$summary), format_median_range(fits$summary$RMSE_ms, 2), format_median_range(fits$summary$MAPE_pct, 2))
}

# ---- 2. Population table (one row per recording) -------------------------------------------------------------------
population <- rbindlist(summaries, fill = TRUE)
if ("fantasia" %in% cfg$datasets) {
  groups <- read_result("fit_fantasia", cfg)$index[, .(SubjectID, Group)]
  population <- merge(population, groups, by = "SubjectID", all.x = TRUE, sort = FALSE)
}
write_result(as.data.frame(population), "df_population_summary", cfg)
fwrite(population, file.path(cfg$paths$results, "population_summary.csv"))
log_msg("Saved results/df_population_summary.rds (%d recordings).", nrow(population))
