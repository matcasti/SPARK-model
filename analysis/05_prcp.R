# ==============================================================================
# 05_prcp.R -- Result 3: do the latent states respond coherently to orthostatic perturbation?
#
# Reads   results/fit_prcp.rds                       (from 01_fit_all.R)
# Writes  results/prcp.rds   recordings (one row each), tilts (one row per tilt: times, windows, analysed),
#                            response / bp / lead (one row per recording, median over its tilts) and
#                            response_events / bp_events / lead_events (one row per tilt), lagged correlation,
#                            dose-response, validation of the RR-based events (per angle-derived tilt),
#                            summary numbers, landscape_modes and the series of every recording (series_all)
#                            from which the figures take the showcased example
# Method  A recording holds several maneuvers. Tilts come from the tilt angle when present (labelled slow / rapid from the protocol
#         annotations, data/prcp_events.csv), else from the annotations, else from heart rate; the stand-ups of the protocol, which
#         the angle cannot show, come from the annotations. Each maneuver is analysed against its own baseline and spaces the next one
#         (R/signals/events.R, tilt_windows()). The pooled analyses use the tilts; responses by kind (slow tilt, rapid tilt, stand-up) are in
#         `by_kind` and `kind_contrasts`. Blood pressure
#         (when available) gives an independent response to compare with the latent states. Between-recording
#         statistics use one row per recording (median over its tilts). Steps: R/steps/prcp.R.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. Load -------------------------------------------------------------------------------
fits <- read_result("fit_prcp", cfg)
prcp_events <- read_prcp_events(cfg)
if (is.null(prcp_events)) {
  log_msg("No protocol annotations (%s): events from the tilt angle and heart rate only, no stand-ups, tilts not labelled slow / rapid.", if (is.null(cfg$paths$prcp_events)) "not configured" else cfg$paths$prcp_events)
} else if (!isTRUE(cfg$tilt$use_annotations)) {
  log_msg("Protocol annotations found but cfg$tilt$use_annotations is FALSE: not used.")
} else {
  log_msg("Protocol annotations: %d maneuvers in %d records (%s)", nrow(prcp_events), uniqueN(prcp_events$id), paste(names(table(prcp_events$kind)), table(prcp_events$kind), collapse = ", "))
}

# ---- 2. Analyse every recording -------------------------------------------------------------------
results <- map_recordings(fits$index$SubjectID, prcp_analyse_recording, cfg$cores,
                          states = fits$states, cfg = cfg, maneuvers = prcp_events)
results <- drop_failed(results)

recordings <- bind_results(results, "recording")
tilts      <- bind_results(results, "events")                  # every maneuver (tilts and stand-ups)
response_events_all <- bind_results(results, "response")
is_tilt <- function(x) if (is.null(x) || !"kind" %in% names(x)) x else x[kind != "stand-up"]
response_events <- is_tilt(response_events_all)                 # pooled analyses: the tilts (stand-ups are reported by kind below)
if (is.null(response_events) || !nrow(response_events)) {      # nothing to analyze: keep the events tables and let the pipeline continue
  write_result(list(events = recordings, tilts = tilts), "prcp", cfg)
  log_msg(paste("No perturbation event could be analysed in any PRCP recording (no usable tilt angle and",
                "no heart-rate change detected, or no event with enough baseline); only the events tables were saved.",
                "Result 3 cannot be populated."))
  quit(save = "no", status = 0)
}

# ---- 3. Combine ---------------------------------------------------------------------------------------
phase_events_all <- bind_results(results, "phase"); lead_events_all <- bind_results(results, "lead"); bp_events_all <- bind_results(results, "bp")
phase_events <- is_tilt(phase_events_all); lead_events <- is_tilt(lead_events_all); bp_events <- is_tilt(bp_events_all)
xcorr        <- bind_results(results, "xcorr")
dose         <- bind_results(results, "dose")

response <- prcp_by_recording(response_events)       # one row per recording: median over its tilts
lead     <- prcp_by_recording(lead_events)
bp       <- prcp_by_recording(bp_events)
if (!"n_analysed" %in% names(recordings)) recordings[, n_analysed := 0L]
recordings[is.na(n_analysed), n_analysed := 0L]

by_kind        <- prcp_kind_summary(response_events_all, bp_events_all, lead_events_all)
kind_contrasts <- prcp_kind_contrasts(response_events_all)
validation  <- prcp_event_validation(tilts)
dose_cor    <- prcp_dose_correlations(dose)
dose_levels <- prcp_dose_levels(dose)
summary <- prcp_summarise(recordings, response, bp, xcorr, dose_cor, lead, validation, tilts, response_events)

# Series of every recording (the figures choose which one to showcase: cfg$figures$example$prcp)
with_series <- Filter(function(r) !is.null(r$series), results)
series_all <- setNames(lapply(with_series, `[[`, "series"), vapply(with_series, function(r) r$recording$SubjectID, character(1)))
landscape_modes <- pooled_phase_landscapes(series_all, "tilt", cfg$figures$landscape)$modes

# Default example (largest standardized response), kept for backward compatibility
example_id <- response[which.max(abs(z_d_x_p) + abs(z_d_x_s)), SubjectID]
example <- series_all[[example_id]]

# ---- 4. Save ----------------------------------------------------------------------------------------------
write_result(list(events = recordings, tilts = tilts, response = response, bp = bp, lead = lead,
                  response_events = response_events, bp_events = bp_events, lead_events = lead_events,
                  xcorr = xcorr, dose = dose, dose_cor = dose_cor, dose_levels = dose_levels, validation = validation,
                  phase = phase_events, summary = summary, example = example, example_id = example_id,
                  by_kind = by_kind, kind_contrasts = kind_contrasts, response_events_all = response_events_all, bp_events_all = bp_events_all,
                  annotations_used = !is.null(prcp_events) && isTRUE(cfg$tilt$use_annotations),
                  series_all = series_all, landscape_modes = landscape_modes), "prcp", cfg)

log_msg("PRCP: %d recordings, %d tilts and %d stand-ups analysed | tilts from angle %d, annotations %d, RR %d, none %d | RR events matched %s, extra %d | BP in %d | dose-response in %d",
        summary$n, summary$n_events, summary$n_stand, summary$n_angle, summary$n_annotation, summary$n_rr, summary$n - summary$n_detected,
        if (is.finite(summary$val_recall)) sprintf("%.0f%%", 100 * summary$val_recall) else "n/a", summary$val_n_extra,
        summary$n_bp, summary$n_dose)
if (!is.null(by_kind)) print(by_kind[, .(kind, n_events, n_recordings, both = round(frac_both_expected, 2), z_d_x_p = round(z_d_x_p_med, 1), z_d_x_s = round(z_d_x_s_med, 1), lat_p = round(lat_p_med), lat_s = round(lat_s_med))])
print(tilts[, .(SubjectID, event, kind, source, t_onset = round(t_onset), t_return = round(t_return), onset_rr = round(onset_rr),
                baseline_s = round(baseline_s), analysed)])
print(response_events[, .(SubjectID, event, z_d_x_p = round(z_d_x_p, 1), z_d_x_s = round(z_d_x_s, 1),
                          lat_p = round(lat_p), lat_s = round(lat_s),
                          tau_p = round(tau_p_resp, 1), tau_s = round(tau_s_resp, 1))])
if (!is.null(bp_events)) print(bp_events[, .(SubjectID, event, d_map = round(d_map, 1), d_pp = round(d_pp, 1),
                                             tau_map = round(tau_map, 1), lat_map = round(lat_map))])
if (!is.null(xcorr)) print(xcorr[, .(SubjectID, state, peak_r = round(peak_r, 2), lag_s, p = round(p, 2))])
if (!is.null(dose_cor)) print(dose_cor)
