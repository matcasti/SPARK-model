# ==============================================================================
# 09_exercise.R -- Result 5: latent dynamics during exercise and recovery (rest -> exercise
# -> recovery) in older adults performing the two-minute step test.
#
# Reads   results/fit_exercise.rds                 (from 01_fit_all.R)
#         results/recovery_choice.rds              (optional; from 07_recovery_benchmark.R)
# Writes  results/exercise.rds   info, rest_response, recovery, geometry, km, summary,
#                                example and series_all (for figure 9)
# Method  Phases are detected from heart rate. No subject characteristics or clinical
#         outcomes exist for this cohort: the analysis is descriptive.
#         The steps are in R/steps/exercise.R.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. Load and choose the recovery definition ------------------------------------------------
fits <- read_result("fit_exercise", cfg)
settings <- exercise_settings(cfg, read_result("recovery_choice", cfg, required = FALSE))
log_msg("Exercise analysis | recovery definition: %s | exercise-end rule: %s",
        settings$recovery_definition, settings$end_method)

# ---- 2. Analyse every recording --------------------------------------------------------------------------
results <- map_recordings(fits$index$SubjectID, exercise_analyse_recording, cfg$cores,
                          states = fits$states, params = fits$params, settings = settings, cfg = cfg)
results <- drop_failed(results)

# ---- 3. Combine and summarise -------------------------------------------------------------------------------
info     <- bind_results(results, "info")
base     <- bind_results(results, "base")
recovery <- bind_results(results, "recovery")
geometry <- bind_results(results, "geometry")

km        <- exercise_km_table(recovery, settings$recovery_definition)
agreement <- exercise_definition_agreement(recovery)
summary   <- exercise_summarise(info, base, geometry, km, agreement, settings)

# Example recording for figure 9: the largest standardized response
example_id <- base[which.max(abs(z_d_x_p) + abs(z_d_x_s)), SubjectID]
with_series <- Filter(function(r) !is.null(r$series), results)
series_all <- setNames(lapply(with_series, `[[`, "series"), vapply(with_series, function(r) r$info$SubjectID, character(1)))
example <- series_all[[example_id]]
landscape_modes <- pooled_phase_landscapes(series_all, "exercise", cfg$figures$landscape)$modes

# ---- 4. Save ------------------------------------------------------------------------------------------------------
write_result(list(info = info, rest_response = base, recovery = recovery, geometry = geometry, km = km,
                  summary = summary, example = example, example_id = example_id, series_all = series_all,
                  landscape_modes = landscape_modes),
             "exercise", cfg)

log_msg("Exercise: usable %d/%d", summary$n_usable, summary$n)
print(base[, .(SubjectID, rest_hr = round(rest_hr), d_hr = round(d_hr),
               z_d_x_p = round(z_d_x_p, 1), z_d_x_s = round(z_d_x_s, 1))])
print(recovery[definition == settings$recovery_definition,
               .(SubjectID, state, time = round(time, 1), status)][order(state, SubjectID)])
print(km)
print(geometry[, .(SubjectID, path_length = round(path_length), loop_area = round(loop_area), U_max = round(U_max))])
