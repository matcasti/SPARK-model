# ==============================================================================
# 13_identifiability_boundary.R -- where can the vagal / sympathetic decomposition be recovered? Simulated recordings with known drives over the
# time-scale ratio tau_s / tau_p (tau_p fixed), the variance ratio lambda_P and the recording length; every specification (M0-M3) fitted to each.
#
# Reads   nothing (simulation)
# Writes  results/ident_boundary.rds   rows (one per simulation x specification), summary by cell, collapse by regime
# Method  R/steps/extended_validation.R (boundary_jobs, boundary_job); grid in cfg$ident_boundary.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. Simulate and fit ---------------------------------------------------------------------------------------------
jobs <- boundary_jobs(cfg$ident_boundary, cfg$seed)
log_msg("Identifiability boundary: %d simulated recordings, 4 specifications each", length(jobs))
rows <- rbindlist(drop_failed(map_recordings(jobs, boundary_job, cfg$cores, settings = cfg$ident_boundary, h = cfg$simulation$h), "simulation"), fill = TRUE)
rows <- rows[ok == TRUE]
rows[, collapse := est_share < cfg$identifiability$collapse_share[1] | est_share > cfg$identifiability$collapse_share[2]]

# ---- 2. Summaries ------------------------------------------------------------------------------------------------------------
by_cell <- rows[, .(n = .N, cor_x_p = median(cor_x_p), cor_x_s = median(cor_x_s), collapse = mean(collapse), abs_ln_err_tau_p = median(abs(ln_err_tau_p)),
                    abs_ln_err_tau_s = median(abs(ln_err_tau_s))), by = .(axis, ratio, lambda, duration, spec)]
rows[, regime := fifelse(axis == "time-scale ratio" & ratio < 10, "ratio < 10", fifelse(axis == "time-scale ratio" & ratio >= 10 & duration >= 600, "ratio >= 10, >= 10 min", NA_character_))]
by_regime <- rows[!is.na(regime), .(n = .N, collapse = mean(collapse), cor_x_s = median(cor_x_s), cor_x_p = median(cor_x_p)), by = .(regime, spec)]

# ---- 3. Save -------------------------------------------------------------------------------------------------------------------
write_result(list(rows = rows, by_cell = by_cell, by_regime = by_regime, settings = cfg$ident_boundary), "ident_boundary", cfg)
print(by_regime[order(regime, spec)])
