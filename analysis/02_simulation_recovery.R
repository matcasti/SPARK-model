# ==============================================================================
# 02_simulation_recovery.R -- parameter and state recovery under known truth, for M0-M3.
#
# Reads   nothing (simulated recordings, seeded)
# Writes  results/sim_recovery.rds   table (one row per simulation x specification), summary
# Scenarios (R/steps/simulation_recovery.R)
#   in_prior   true time scales inside the prior ranges
#   off_prior  true time scales outside them (what do the priors do?)
#   vagal_dom / symp_dom   contrasting true variance partition (stress on the spectral anchor)
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. One job per simulated recording ------------------------------------------------------------------------
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(scenario) {
  lapply(seq_len(cfg$simulation$n_per_scenario), function(i) {
    list(scenario = scenario, sim = i, seed = id_seed(paste(scenario, i), cfg$seed))
  })
}), recursive = FALSE)
log_msg("Simulation recovery: %d recordings x 4 specs", length(jobs))

# ---- 2. Simulate and fit ----------------------------------------------------------------------------------------------
results <- drop_failed(map_recordings(jobs, recovery_one_simulation, cfg$cores, cfg = cfg), "simulation")
table <- rbindlist(results)
summary <- recovery_summarise(table)

# ---- 3. Save -------------------------------------------------------------------------------------------------------------
write_result(list(table = table, summary = summary), "sim_recovery", cfg)
print(summary[order(scenario, spec)])
