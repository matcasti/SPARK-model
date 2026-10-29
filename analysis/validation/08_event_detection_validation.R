# ==============================================================================
# 08_event_detection_validation.R -- validation of heart-rate-based event detection on simulated
# recordings with known timing.
#
# Reads   nothing (simulated recordings, seeded)
# Writes  results/event_detection.rds   tilt, repeated (per tilt), exercise, rest (per simulation),
#                                       false_detection (by minimum-shift rule), summary
# Reports detection rate and onset error for tilt-like recordings (single and repeated tilts), onset / end errors for
#         exercise-like recordings (three end rules), and the false-detection rate on stationary
#         recordings as a function of the minimum-shift rule (R/steps/benchmarks.R).
# ==============================================================================
source("analysis/00_setup.R")
n_simulations <- if (!is.null(cfg$simulation$n_detect)) cfg$simulation$n_detect else 20L
seeds <- id_seed(seq_len(n_simulations), cfg$seed)

# ---- 1. Simulate and detect -------------------------------------------------------------------------------------------
stack <- function(rows) rbindlist(Filter(Negate(is_err), rows))
tilt     <- stack(map_recordings(seeds, detection_tilt_row, cfg$cores, cfg = cfg))
repeated <- stack(map_recordings(seeds, detection_repeated_row, cfg$cores, cfg = cfg))
exercise <- stack(map_recordings(seeds, detection_exercise_row, cfg$cores, cfg = cfg))
rest     <- stack(map_recordings(seeds, detection_rest_rows, cfg$cores, cfg = cfg))

# ---- 2. Summarise ----------------------------------------------------------------------------------------------------------
false_detection <- rest[, .(false_detection_rate = mean(detected), n = .N), by = min_shift_hz]
detection_rate <- rbindlist(lapply(MIN_SHIFT_CANDIDATES_HZ, detection_rate_for_shift, seeds = seeds, cfg = cfg))
end_error <- function(column) median(abs(exercise[[column]]), na.rm = TRUE)
summary <- list(
  tilt_detected = mean(tilt$detected), tilt_onset_err = median_range(tilt$onset_err),
  repeated_recall = mean(repeated$matched), repeated_onset_err = median_range(repeated$onset_err),
  repeated_extra = mean(repeated[, .(extra = n_extra[1]), by = seed]$extra),
  ex_onset_err_cp = median_range(exercise$onset_err_cp), ex_onset_err_profile = median_range(exercise$onset_err_profile),
  ex_end_mae = c(peak = end_error("end_err_peak"), cp = end_error("end_err_cp"), profile = end_error("end_err_profile"))
)

# ---- 3. Save --------------------------------------------------------------------------------------------------------------------
write_result(list(tilt = tilt, repeated = repeated, exercise = exercise, rest = rest,
                  false_detection = merge(false_detection, detection_rate, by = "min_shift_hz"), summary = summary),
             "event_detection", cfg)
log_msg(paste("Event detection: tilt detected %.0f%%, onset error median %.1f s | repeated tilts: %.0f%% found, %.1f extra events per recording |",
              "exercise end MAE (peak/cp/profile): %s s | false detection at %.2f Hz: %.0f%%"),
        100 * mean(tilt$detected), median(tilt$onset_err, na.rm = TRUE),
        100 * summary$repeated_recall, summary$repeated_extra,
        paste(round(summary$ex_end_mae, 1), collapse = "/"), cfg$segmentation$min_shift_hz,
        100 * false_detection[min_shift_hz == cfg$segmentation$min_shift_hz, false_detection_rate])
print(merge(false_detection, detection_rate, by = "min_shift_hz"))
