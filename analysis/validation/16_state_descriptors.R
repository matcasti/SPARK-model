# ==============================================================================
# 16_state_descriptors.R -- perturbation-response descriptors (maximum displacement, recovery time, probability of return, regime transition) computed from the
# reconstructed states with their uncertainty: first on simulated recordings with known states (can they be recovered?), then on the exercise recordings and
# the PRCP tilts and stand-ups.
#
# Reads   results/exercise.rds, results/prcp.rds (states, events)
# Writes  results/state_descriptors.rds   sim (one row per simulated recording), exercise and tilt (one row per perturbation)
# Method  R/signals/state_descriptors.R, R/signals/regimes.R, R/steps/extended_validation.R; settings in cfg$state_descriptors and cfg$descriptor_sims.
# ==============================================================================
source("analysis/00_setup.R")
settings <- cfg$state_descriptors

# ---- 1. Simulations with known states ----------------------------------------------------------------------------------------------
jobs <- unlist(lapply(c("exercise", "tilt", "shift"), function(sc) lapply(seq_len(cfg$descriptor_sims$n_per_scenario), function(i)
  list(scenario = sc, sim = i, seed = id_seed(paste("stdesc", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Perturbation descriptors: %d simulated recordings", length(jobs))
sim <- rbindlist(drop_failed(map_recordings(jobs, descriptor_sim_job, cfg$cores, settings = settings), "simulation"), fill = TRUE)

# ---- 2. Real recordings ---------------------------------------------------------------------------------------------------------------------
exercise <- read_result("exercise", cfg); prcp <- read_result("prcp", cfg)
ex_jobs <- lapply(names(exercise$series_all), function(id) { s <- exercise$series_all[[id]]
  list(id = id, subject = id, kind = "exercise", s = s, rest = c(cfg$exercise$rest_skip_s, s$t_on - 5), pert = c(s$t_on, s$t_end), H = settings$horizon_exercise,
       t_max = NA_real_, whole = TRUE, seed = id_seed(paste("stex", id), cfg$seed)) })
tilts <- prcp$tilts[usable == TRUE]
tl_jobs <- lapply(seq_len(nrow(tilts)), function(i) { e <- tilts[i]
  list(id = paste(e$SubjectID, e$event, sep = "_"), subject = e$SubjectID, kind = e$kind, s = prcp$series_all[[e$SubjectID]], rest = c(e$t_base_from, e$t_onset - 5),
       pert = c(e$t_onset, e$t_tilt_end), H = settings$horizon_tilt, t_max = e$t_rec_to, whole = FALSE, seed = id_seed(paste("sttl", e$SubjectID, e$event), cfg$seed)) })
log_msg("Perturbation descriptors: %d exercise recordings, %d PRCP maneuvers", length(ex_jobs), length(tl_jobs))
ex <- rbindlist(drop_failed(map_recordings(ex_jobs, descriptor_real_job, cfg$cores, settings = settings), "exercise"), fill = TRUE)
tl <- rbindlist(drop_failed(map_recordings(tl_jobs, descriptor_real_job, cfg$cores, settings = settings), "maneuver"), fill = TRUE)

# Recovery times of the main analysis (chosen definition) for comparison
rec <- exercise$recovery[definition == exercise$summary$recovery_definition & status == "recovered", .(id = SubjectID, state, time)]
recovery_main <- dcast(rec, id ~ state, value.var = "time")

# ---- 3. Save --------------------------------------------------------------------------------------------------------------------------------------
write_result(list(sim = sim, exercise = ex, tilt = tl, recovery_main = recovery_main, settings = settings), "state_descriptors", cfg)
print(sim[, .(n = .N, Trec_true = median(Trec_true, na.rm = TRUE), Trec = median(Trec, na.rm = TRUE), transition = mean(transition, na.rm = TRUE)), by = scenario])
