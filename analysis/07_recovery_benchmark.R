# ==============================================================================
# 07_recovery_benchmark.R -- empirical choice of the recovery-time definition and of the
# exercise-end rule, on simulated recordings with KNOWN recovery constants.
#
# Reads   nothing (simulated recordings, seeded)
# Writes  results/recovery_choice.rds   choice (definition, end_method, performance, aggregate,
#                                       end_rule, n_profiles), table (every metric)
# Criteria, per state (R/steps/benchmarks.R)
#   success      share of recordings with a defined, uncensored estimate
#   truth rank   Spearman(estimate, true recovery constant)
#   retest rank  Spearman between two noise realizations of the same true profile
#   score = success * (truth + retest) / 2, averaged over x_p and x_s.
#   The exercise-end rule with the smallest median absolute error against the true end wins.
# ==============================================================================
source("analysis/00_setup.R")
n_profiles <- if (!is.null(cfg$recovery_benchmark$n)) cfg$recovery_benchmark$n else 40L

# ---- 1. Simulate, fit and measure every profile ----------------------------------------------------------------
log_msg("Recovery benchmark: %d simulated profiles x 2 realizations", n_profiles)
seeds <- id_seed(seq_len(n_profiles), cfg$seed)
results <- drop_failed(map_recordings(seeds, benchmark_one_profile, cfg$cores, cfg = cfg), "profile")
table <- rbindlist(results, fill = TRUE)

# ---- 2. Choose ------------------------------------------------------------------------------------------------------------
choice <- choose_recovery_definition(table, n_profiles)

# ---- 3. Save --------------------------------------------------------------------------------------------------------------
write_result(list(choice = choice, table = table), "recovery_choice", cfg)
log_msg("Chosen recovery definition: %s | exercise-end rule: %s", choice$definition, choice$end_method)
print(choice$aggregate)
print(choice$end_rule)
