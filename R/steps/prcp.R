# ==============================================================================
# R/steps/prcp.R -- Result 3: orthostatic perturbation (PRCP tilt-table recordings).
# Called by analysis/datasets/05_prcp.R. A recording holds several tilts; each tilt is one event, analysed against its
# own baseline (tilt_windows()). prcp_analyse_recording() chains the steps below; the per-event tables are
# condensed to one row per recording (median over its events) before any between-recording statistic.
#
#   prcp_define_events()       every tilt: from the tilt angle, else from heart rate; their windows
#   prcp_state_response()      change of x_p, x_s between phases + step-response fits (one event)
#   prcp_placebo_response()    negative control: the same contrast between the two halves of the resting baseline (one event)
#   prcp_pressure_response()   same for blood pressure (when available)
#   prcp_lagged_correlation()  states vs MAP after removing slow trends (whole recording)
#   prcp_dose_response()       mean states / HR / MAP over tilt-angle plateaus (whole recording)
#   prcp_cwt_timing()          first departure from baseline: states vs CWT band power (one event)
#   prcp_by_recording()        per-event table -> one row per recording
# ==============================================================================


# ---- Step 1: events -----------------------------------------------------------------------------
#' Every maneuver of a recording. Tilts come from the tilt angle when it is present (labelled slow / rapid from the protocol
#' annotations when there are any), else from the annotations, else from heart rate; the stand-ups of the protocol, which the table angle
#' cannot show, come from the annotations. The heart-rate events are always detected as well and matched to the angle / annotation events
#' (validation of the RR-only detection). Every maneuver is an event with its own baseline, and every maneuver spaces the next baseline.
#' @param maneuvers annotated maneuvers of this recording (annotated_maneuvers()), or NULL
#' @return list(recording = one-row table, events = one row per maneuver (kind, times, matched heart-rate times, windows),
#'              windows = tilt_windows() of the events used, with their kind)
prcp_define_events <- function(states, id, cfg, maneuvers = NULL) {
  t_end <- tail(states$Time, 1)
  if (!isTRUE(cfg$tilt$use_annotations)) maneuvers <- NULL
  from_rr <- detect_events(states$RR_Interval, "tilt", cfg$segmentation)
  from_angle <- if (column_available(states, "angle", cfg$io$min_cov_coverage)) {
    tilt_events(states$Time, states$angle, cfg$tilt)
  } else {
    list(detected = FALSE, n_events = 0L, base_angle = NA_real_, events = NULL)
  }
  tilts_annotated <- if (is.null(maneuvers)) NULL else maneuvers[maneuvers$kind != "stand-up", ]
  stands <- if (is.null(maneuvers)) NULL else maneuvers[maneuvers$kind == "stand-up", ]
  has_annotated_tilts <- !is.null(tilts_annotated) && nrow(tilts_annotated) > 0
  has_stands <- !is.null(stands) && nrow(stands) > 0

  source <- if (from_angle$detected) "angle" else if (has_annotated_tilts) "annotation" else if (from_rr$detected) "rr" else if (has_stands) "annotation" else "none"
  empty <- data.table(event = integer(), kind = character(), source = character(), t_onset = numeric(), t_return = numeric(), peak_angle = numeric())
  tilt_part <- switch(source,
    angle = {
      tilts <- from_angle$events[, .(event, t_onset, t_return, peak_angle)]
      cbind(tilts[, .(event)], kind = label_tilt_kinds(tilts$t_onset, maneuvers, cfg$tilt$annotation_slack_s), source = "angle", tilts[, .(t_onset, t_return, peak_angle)])
    },
    annotation = if (has_annotated_tilts) {
      annotated <- data.table(event = seq_len(nrow(tilts_annotated)), kind = tilts_annotated$kind, source = "annotation", t_onset = tilts_annotated$onset_s,
                              t_return = ifelse(is.finite(tilts_annotated$return_end_s), tilts_annotated$return_end_s, tilts_annotated$return_s), peak_angle = NA_real_)
      annotated[is.finite(t_onset) & is.finite(t_return)]                  # a tilt whose return was not annotated cannot define a window
    } else empty[0],
    rr = {
      rr_events <- from_rr$events[, .(event, t_onset, t_return = t_recovery)]
      if (has_stands) rr_events <- rr_events[vapply(t_onset, function(t) all(abs(t - stands$onset_s) > cfg$tilt$match_tol_s), logical(1))]   # those are the annotated stand-ups
      cbind(rr_events[, .(event)], kind = "tilt", source = "rr", rr_events[, .(t_onset, t_return)], peak_angle = NA_real_)
    },
    empty[0])
  stand_part <- if (has_stands) data.table(event = seq_len(nrow(stands)), kind = "stand-up", source = "annotation", t_onset = stands$onset_s,
                                           t_return = stands$return_s, peak_angle = NA_real_) else empty[0]
  used <- rbind(tilt_part, stand_part)
  setorder(used, t_onset)
  used[, event := seq_len(.N)]

  if (!nrow(used)) {
    return(list(recording = data.table(SubjectID = id, source = "none", detected = FALSE, n_events = 0L, n_usable = 0L,
                                       n_angle = from_angle$n_events, n_rr = from_rr$n_events, n_rr_matched = NA_integer_,
                                       n_rr_extra = NA_integer_, K_rr = from_rr$K, base_angle = from_angle$base_angle, n_stand = 0L),
                events = NULL, windows = tilt_windows(numeric(), numeric(), 0, t_end, cfg$tilt)))
  }
  reference <- used$source != "rr"                              # events whose timing does not come from heart rate
  matched <- if (any(reference)) match_events(used$t_onset[reference], from_rr$events$t_onset, cfg$tilt$match_tol_s) else NULL
  windows <- tilt_windows(used$t_onset, used$t_return, 0, t_end, cfg$tilt)
  windows[, kind := used$kind]
  events <- cbind(SubjectID = rep(id, nrow(used)), used[, .(source, kind, event, t_onset, t_return, peak_angle)])
  events[, `:=`(onset_rr = NA_real_, return_rr = NA_real_)]
  if (!is.null(matched)) {
    rows <- which(reference)[!is.na(matched$detected)]
    detected <- matched$detected[!is.na(matched$detected)]
    events[rows, `:=`(onset_rr = from_rr$events$t_onset[detected], return_rr = from_rr$events$t_recovery[detected])]
  }
  events[source == "rr", `:=`(onset_rr = t_onset, return_rr = t_return)]
  events <- cbind(events, windows[, .(t_base_from, t_tilt_end, t_rec_to, baseline_s, usable)])

  recording <- data.table(
    SubjectID = id, source = source, detected = TRUE, n_events = nrow(used), n_usable = sum(windows$usable),
    n_angle = from_angle$n_events, n_rr = from_rr$n_events,
    n_rr_matched = if (is.null(matched)) NA_integer_ else sum(!is.na(matched$detected)),
    n_rr_extra = if (is.null(matched)) NA_integer_ else from_rr$n_events - sum(!is.na(matched$detected)),
    K_rr = from_rr$K, base_angle = from_angle$base_angle, n_stand = sum(used$kind == "stand-up"))
  list(recording = recording, events = events, windows = windows)
}

# ---- Step 2: response of the latent states ---------------------------------------------------------------
#' Phase contrast (tilt - baseline, in baseline SDs) and step-response fits of x_p and x_s for one event `w`
#' (one row of tilt_windows()). Only the beats of the event's own windows are used.
#' @return list(response = one-row table, phase_summary = per-phase means), or NULL if a phase is too short
prcp_state_response <- function(states, id, w) {
  d <- states[Time >= w$t_base_from & Time < w$t_rec_to]
  d[, phase := factor(fcase(Time < w$t_onset, "baseline", Time < w$t_tilt_end, "tilt", default = "recovery"),
                      levels = c("baseline", "tilt", "recovery"))]
  phase_means <- phase_summary(data.table(t = d$Time, x_p = d$X_p, x_s = d$X_s, phase = d$phase))
  baseline <- phase_means[phase == "baseline"]
  tilt <- phase_means[phase == "tilt"]
  if (!nrow(baseline) || !nrow(tilt) || baseline$n < 20 || tilt$n < 20) return(NULL)

  baseline_sd <- c(p = sd(d[phase == "baseline", X_p]), s = sd(d[phase == "baseline", X_s]))
  change_p <- tilt$mean_x_p - baseline$mean_x_p
  change_s <- tilt$mean_x_s - baseline$mean_x_s
  until_return <- d[Time < w$t_tilt_end]
  step_fit <- lapply(c(p = "X_p", s = "X_s"), function(column) fit_step_response(until_return$Time, until_return[[column]], w$t_onset))

  response <- data.table(
    SubjectID = id, event = w$event, d_x_p = change_p, d_x_s = change_s,
    z_d_x_p = change_p / baseline_sd[["p"]], z_d_x_s = change_s / baseline_sd[["s"]],
    expected_p = change_p < 0, expected_s = change_s > 0,
    hr_baseline = 60 / mean(d[phase == "baseline", RR_Interval]),
    hr_tilt = 60 / mean(d[phase == "tilt", RR_Interval]),
    t_on_p = field_or_na(step_fit$p, "t_on"), t_on_s = field_or_na(step_fit$s, "t_on"),
    tau_p_resp = field_or_na(step_fit$p, "tau"), tau_s_resp = field_or_na(step_fit$s, "tau"),
    amp_p = field_or_na(step_fit$p, "D"), amp_s = field_or_na(step_fit$s, "D"),
    r2_p = field_or_na(step_fit$p, "r2"), r2_s = field_or_na(step_fit$s, "r2"),
    bound_p = isTRUE(step_fit$p$at_bound), bound_s = isTRUE(step_fit$s$at_bound)
  )
  response[, `:=`(lag_s_minus_p = t_on_s - t_on_p, tau_ratio_s_over_p = tau_s_resp / tau_p_resp,
                  lat_p = t_on_p - w$t_onset, lat_s = t_on_s - w$t_onset)]
  list(response = response, phase_summary = cbind(SubjectID = id, event = w$event, phase_means))
}


#' Negative control for the phase contrast: the resting baseline of an event is split in two halves and the second half is contrasted with
#' the first exactly as the tilt is contrasted with the baseline (difference of the means in SDs of the first part; expected signs: x_p
#' down, x_s up). No maneuver separates the two halves, so any "expected" pattern here arises from drift, noise or from the way
#' the contrast is built, and not from a perturbation.
#' @return one-row table with the same columns as the main contrast, or NULL when a half has fewer than `min_beats` beats
prcp_placebo_response <- function(states, id, w, min_beats = 20L) {
  middle <- w$t_base_from + (w$t_onset - w$t_base_from) / 2
  first <- states[Time >= w$t_base_from & Time < middle]
  second <- states[Time >= middle & Time < w$t_onset]
  if (nrow(first) < min_beats || nrow(second) < min_beats) return(NULL)
  change_p <- mean(second$X_p) - mean(first$X_p)
  change_s <- mean(second$X_s) - mean(first$X_s)
  data.table(SubjectID = id, event = w$event, z_d_x_p = change_p / sd(first$X_p), z_d_x_s = change_s / sd(first$X_s),
             expected_p = change_p < 0, expected_s = change_s > 0)
}


# ---- Step 3: blood pressure -------------------------------------------------------------------------------
#' Pressure change between the event's baseline and its tilt (first 30 s after the onset excluded) and step
#' response of MAP.
#' @return one-row table or NULL when the recording has no blood pressure
prcp_pressure_response <- function(states, id, w, response, min_coverage) {
  has_pressure <- all(vapply(c("map", "systolic", "diastolic"), column_available, logical(1),
                             states = states, min_coverage = min_coverage))
  if (!has_pressure) return(NULL)

  d <- states[Time >= w$t_base_from & Time < w$t_tilt_end]
  d[, pulse_pressure := systolic - diastolic]
  in_tilt <- d$Time >= w$t_onset + 30
  in_baseline <- d$Time < w$t_onset - 5
  change <- function(column) mean(d[[column]][in_tilt], na.rm = TRUE) - mean(d[[column]][in_baseline], na.rm = TRUE)
  map_fit <- fit_step_response(d$Time, d$map, w$t_onset)

  pressure <- data.table(
    SubjectID = id, event = w$event, d_sbp = change("systolic"), d_dbp = change("diastolic"),
    d_map = change("map"), d_pp = change("pulse_pressure"),
    z_d_map = change("map") / sd(d$map[in_baseline], na.rm = TRUE),
    t_on_map = field_or_na(map_fit, "t_on"), tau_map = field_or_na(map_fit, "tau"),
    r2_map = field_or_na(map_fit, "r2")
  )
  pressure[, `:=`(lat_map = t_on_map - w$t_onset,
                  lag_map_minus_s = t_on_map - response$t_on_s,
                  lag_map_minus_p = t_on_map - response$t_on_p)]
  pressure
}

#' Peak lagged correlation of each state with MAP after removing slow trends
#' (positive lag: MAP lags the state), with a circular-shift surrogate p-value. Whole recording.
prcp_lagged_correlation <- function(states, id, cfg) {
  one_state <- function(column) {
    xcorr_surrogate(states$Time, states[[column]], states$map, cfg$physio$max_lag_s, 60,
                    cfg$physio$n_surrogates, id_seed(paste(id, column), cfg$seed))
  }
  correlations <- lapply(c(X_p = "X_p", X_s = "X_s"), one_state)
  rows <- lapply(names(correlations), function(state) {
    r <- correlations[[state]]
    if (is.null(r)) return(NULL)
    data.table(SubjectID = id, state = state, peak_r = r$peak_r, lag_s = r$lag_s, p = r$p, r0 = r$r0)
  })
  rbindlist(rows)
}


# ---- Step 4: dose-response over tilt-angle plateaus --------------------------------------------------------
#' Mean x_p, x_s, HR (and MAP) over the second half of each angle plateau; NULL without plateaus.
prcp_dose_response <- function(states, id, cfg) {
  if (!column_available(states, "angle", cfg$io$min_cov_coverage)) return(NULL)
  variables <- list(x_p = states$X_p, x_s = states$X_s, hr = 60 / states$RR_Interval)
  if (column_available(states, "map", cfg$io$min_cov_coverage)) variables$map <- states$map
  plateaus <- plateau_means(states$Time, states$angle, variables, cfg$tilt)
  if (is.null(plateaus) || !nrow(plateaus)) NULL else cbind(SubjectID = id, plateaus)
}


# ---- Step 5: temporal resolution vs CWT -----------------------------------------------------------------------
#' First time each series leaves its baseline band, relative to the onset: filtered states (causal) versus
#' the CWT moving LF / HF band power of the same recording, for one event `w`. The search stops at the return.
#' @param cwt cwt_band_series() of the whole recording (same clock as the states)
#' @return one-row table
prcp_cwt_timing <- function(states, cwt, id, w) {
  baseline_window <- c(max(w$t_base_from, w$t_onset - 150), w$t_onset - 15)
  departure <- function(time, series) first_departure(time, series, baseline_window, t_max = w$t_tilt_end)

  timing <- data.table(
    SubjectID = id, event = w$event, onset = w$t_onset,
    dep_x_p = departure(states$Time, states$X_p), dep_x_s = departure(states$Time, states$X_s),
    dep_LF = departure(cwt$t, cwt$LF), dep_HF = departure(cwt$t, cwt$HF)
  )
  timing[, `:=`(lead_x_p = dep_x_p - onset, lead_x_s = dep_x_s - onset,
                lead_LF = dep_LF - onset, lead_HF = dep_HF - onset)]
}


# ---- One recording: chain the steps -----------------------------------------------------------------------------
#' Analyse every usable tilt of one recording.
#' @param id recording id; @param states the fitted states of the whole dataset (all recordings)
#' @return list of per-recording tables (`events` one row per tilt, `recording`, and per-event `response`, `bp`,
#'         `phase`, `lead`; per-recording `xcorr`, `dose`) and the `series` used for the example figures;
#'         only `recording` and `events` when no event could be analysed.
prcp_analyse_recording <- function(id, states, cfg, maneuvers = NULL) {
  recording <- states[SubjectID == id]
  min_coverage <- cfg$io$min_cov_coverage

  ev <- prcp_define_events(recording, id, cfg, annotated_maneuvers(maneuvers, id))
  if (!ev$recording$detected) return(list(recording = ev$recording, events = ev$events))

  cwt <- cwt_band_series(recording$RR_Interval, cfg$hrv, tt = recording$Time)      # same clock as the states
  one_event <- function(i) {
    w <- ev$windows[i]
    state_response <- prcp_state_response(recording, id, w)
    if (is.null(state_response)) return(NULL)
    list(response = state_response$response, phase = state_response$phase_summary,
         bp = prcp_pressure_response(recording, id, w, state_response$response, min_coverage),
         lead = prcp_cwt_timing(recording, cwt, id, w), placebo = prcp_placebo_response(recording, id, w))
  }
  per_event <- setNames(lapply(which(ev$windows$usable), one_event), which(ev$windows$usable))
  per_event <- Filter(Negate(is.null), per_event)
  ev$events[, analysed := event %in% as.integer(names(per_event))]
  analysed_kinds <- ev$events[analysed == TRUE, kind]
  ev$recording[, `:=`(n_analysed = sum(analysed_kinds != "stand-up"), n_stand_analysed = sum(analysed_kinds == "stand-up"), detected = length(per_event) > 0)]
  if (!length(per_event)) return(list(recording = ev$recording, events = ev$events))
  tag_kind <- function(x) { if (!is.null(x) && nrow(x) && "event" %in% names(x)) x[, kind := ev$events$kind[match(event, ev$events$event)]]; x }

  lagged <- if (column_available(recording, "map", min_coverage)) prcp_lagged_correlation(recording, id, cfg) else NULL
  tilt_windows_used <- ev$windows[event %in% as.integer(names(per_event)) & kind != "stand-up"]   # baselines for the standardization of the figures: tilts only
  phases <- assign_tilt_phases(recording$Time, tilt_windows_used)
  first_tilt <- tilt_windows_used[1]
  list(
    recording = ev$recording, events = ev$events,
    response = tag_kind(rbindlist(lapply(per_event, `[[`, "response"))), bp = tag_kind(rbindlist(lapply(per_event, `[[`, "bp"))),
    phase = tag_kind(rbindlist(lapply(per_event, `[[`, "phase"))), lead = tag_kind(rbindlist(lapply(per_event, `[[`, "lead"))),
    placebo = tag_kind(rbindlist(lapply(per_event, `[[`, "placebo"))),
    xcorr = lagged, dose = prcp_dose_response(recording, id, cfg),
    series = list(
      t = recording$Time, x_p = recording$X_p, x_s = recording$X_s,
      vp = recording$X_p_var, vs = recording$X_s_var, rr = recording$RR_Interval,
      map = if (column_available(recording, "map", min_coverage)) recording$map else NULL,
      angle = if (column_available(recording, "angle", min_coverage)) recording$angle else NULL,
      ev = list(windows = ev$windows[event %in% as.integer(names(per_event))],
                t_onset = if (nrow(first_tilt)) first_tilt$t_onset else NA_real_,
                t_recovery = if (nrow(first_tilt)) first_tilt$t_return else NA_real_),
      cwt = cwt,
      sd_p = sd(recording$X_p[which(phases$phase == "baseline")]), sd_s = sd(recording$X_s[which(phases$phase == "baseline")])
    )
  )
}

#' One row per recording from a per-event table: the median of every numeric column over its events
#' (`event` excluded), the share of events for logical columns and the number of events. Between-recording
#' statistics use these rows, so that a recording with more tilts does not weigh more.
prcp_by_recording <- function(per_event) {
  if (is.null(per_event) || !nrow(per_event)) return(NULL)
  numeric_columns <- setdiff(names(per_event)[vapply(per_event, is.numeric, logical(1))], c("event"))
  logical_columns <- names(per_event)[vapply(per_event, is.logical, logical(1))]
  safe_median <- function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
  rows <- per_event[, c(lapply(.SD, safe_median), list(n_events = .N)), by = SubjectID, .SDcols = numeric_columns]
  if (length(logical_columns)) {
    shares <- per_event[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), by = SubjectID, .SDcols = logical_columns]
    rows <- merge(rows, shares, by = "SubjectID", sort = FALSE)
    rows[, (logical_columns) := lapply(.SD, function(x) x >= 0.5), .SDcols = logical_columns]
  }
  if ("d_x_p" %in% names(rows)) rows[, `:=`(expected_p = d_x_p < 0, expected_s = d_x_s > 0)]
  rows[]
}


# ---- Across recordings -----------------------------------------------------------------------------------------------
#' Accuracy of the heart-rate-based events relative to the angle- or annotation-derived ones (tilts and stand-ups), one row per such maneuver.
#' `matched` is FALSE when no heart-rate event lies within `match_tol_s` of the onset.
prcp_event_validation <- function(tilts) {
  if (is.null(tilts)) return(data.table(SubjectID = character(), event = integer(), matched = logical(),
                                        err_onset = numeric(), err_return = numeric()))
  tilts[source %in% c("angle", "annotation"), .(SubjectID, event, kind, matched = is.finite(onset_rr),
                             err_onset = onset_rr - t_onset, err_return = return_rr - t_return)]
}

#' Within-recording rank correlation between tilt angle and each variable (recordings with
#' at least 4 levels). NULL when no recording qualifies (e.g. a single-step protocol).
prcp_dose_correlations <- function(dose) {
  if (is.null(dose)) return(NULL)
  variables <- intersect(c("x_p", "x_s", "hr", "map"), names(dose))
  per_recording <- function(d) {
    if (uniqueN(d$level) < 4) return(NULL)
    rbindlist(lapply(variables, function(v) {
      data.table(SubjectID = d$SubjectID[1], variable = v, n_levels = uniqueN(d$level),
                 rho = suppressWarnings(cor(d$level, d[[v]], method = "spearman")))
    }))
  }
  correlations <- rbindlist(lapply(split(dose, dose$SubjectID), per_recording), fill = TRUE)
  if (nrow(correlations)) correlations else NULL
}

#' Mean of each variable by tilt-angle level (across recordings).
prcp_dose_levels <- function(dose) {
  if (is.null(dose)) return(NULL)
  variables <- intersect(c("x_p", "x_s", "hr", "map"), names(dose))
  dose[, lapply(.SD, mean, na.rm = TRUE), by = level, .SDcols = variables][order(level)]
}

#' Negative control against the real contrast, for the tilts: share of events with the expected pattern in the placebo contrast (halves
#' of the baseline) and in the tilt contrast, with exact 95% intervals (descriptive: events of one recording are not independent), the
#' median absolute standardized change of each state in both, and the paired comparison over recordings (median over the recording's
#' events; exact Wilcoxon signed-rank). Only events that have both contrasts are used.
#' @param response_events,placebo_events per-event tables of the tilts (columns SubjectID, event, z_d_x_p, z_d_x_s, expected_p, expected_s)
prcp_placebo_summary <- function(response_events, placebo_events) {
  if (is.null(placebo_events) || !nrow(placebo_events)) return(NULL)
  keys <- c("SubjectID", "event")
  both <- merge(response_events[, c(keys, "z_d_x_p", "z_d_x_s", "expected_p", "expected_s"), with = FALSE],
                placebo_events[, c(keys, "z_d_x_p", "z_d_x_s", "expected_p", "expected_s"), with = FALSE], by = keys, suffixes = c("", "_placebo"))
  if (!nrow(both)) return(NULL)
  n <- nrow(both)
  k_real <- sum(both$expected_p & both$expected_s); k_placebo <- sum(both$expected_p_placebo & both$expected_s_placebo)
  per_recording <- both[, .(abs_p = median(abs(z_d_x_p)), abs_s = median(abs(z_d_x_s)),
                            abs_p_placebo = median(abs(z_d_x_p_placebo)), abs_s_placebo = median(abs(z_d_x_s_placebo))), by = SubjectID]
  signed_rank <- function(a, b) if (nrow(per_recording) >= 4) suppressWarnings(stats::wilcox.test(a, b, paired = TRUE, alternative = "greater", exact = TRUE)$p.value) else NA_real_
  list(n_events = n, n_recordings = nrow(per_recording),
       k_real = k_real, k_placebo = k_placebo, frac_real = k_real / n, frac_placebo = k_placebo / n,
       ci_real = exact_ci(k_real, n), ci_placebo = exact_ci(k_placebo, n),
       abs_z_p = median(abs(both$z_d_x_p)), abs_z_s = median(abs(both$z_d_x_s)),
       abs_z_p_placebo = median(abs(both$z_d_x_p_placebo)), abs_z_s_placebo = median(abs(both$z_d_x_s_placebo)),
       p_paired_p = signed_rank(per_recording$abs_p, per_recording$abs_p_placebo), p_paired_s = signed_rank(per_recording$abs_s, per_recording$abs_s_placebo))
}

#' Numbers quoted in the manuscript (median [2.5th, 97.5th percentile] across recordings, counts, shares).
#' `recordings`, `response`, `bp`, `lead`: one row per recording (prcp_by_recording()); `tilts` and
#' `response_events`: one row per tilt, for counts and the within-recording consistency.
prcp_summarise <- function(recordings, response, bp, xcorr, dose_cor, lead, validation, tilts, response_events) {
  n_or_zero <- function(x, column) if (is.null(x)) 0L else uniqueN(x[[column]])
  matched <- validation[matched == TRUE]
  list(
    n = nrow(recordings), n_detected = sum(recordings$detected),
    n_angle = sum(recordings$source == "angle"), n_annotation = sum(recordings$source == "annotation"), n_rr = sum(recordings$source == "rr"),
    n_events = if (is.null(tilts)) 0L else sum(tilts$analysed & tilts$kind != "stand-up", na.rm = TRUE),
    n_stand = if (is.null(tilts)) 0L else sum(tilts$analysed & tilts$kind == "stand-up", na.rm = TRUE),
    events_per_recording = median_range(recordings$n_analysed),
    frac_expected_p = mean(response$expected_p), frac_expected_s = mean(response$expected_s),
    frac_both_expected = mean(response$expected_p & response$expected_s),
    ci_both_expected = exact_ci(sum(response$expected_p & response$expected_s), nrow(response)),
    frac_events_both_expected = mean(response_events$expected_p & response_events$expected_s),
    z_d_x_p = median_range(response$z_d_x_p), z_d_x_s = median_range(response$z_d_x_s),
    lag_s_minus_p = median_range(response$lag_s_minus_p),
    tau_p_resp = median_range(response$tau_p_resp), tau_s_resp = median_range(response$tau_s_resp),
    lat_p = median_range(response$lat_p), lat_s = median_range(response$lat_s),
    frac_fit_at_bound_p = mean(response$bound_p), frac_fit_at_bound_s = mean(response$bound_s),
    lead_x_s = median_range(lead$lead_x_s), lead_LF = median_range(lead$lead_LF),
    lead_x_p = median_range(lead$lead_x_p), lead_HF = median_range(lead$lead_HF),
    frac_LF_leads_x_s = mean(lead$lead_LF < lead$lead_x_s, na.rm = TRUE),
    n_bp = if (is.null(bp)) 0L else nrow(bp),
    d_map = median_range(bp$d_map), d_sbp = median_range(bp$d_sbp),
    d_dbp = median_range(bp$d_dbp), d_pp = median_range(bp$d_pp),
    tau_map = median_range(bp$tau_map), lat_map = median_range(bp$lat_map),
    lag_map_minus_s = median_range(bp$lag_map_minus_s), lag_map_minus_p = median_range(bp$lag_map_minus_p),
    n_xcorr = n_or_zero(xcorr, "SubjectID"),
    xcorr = if (is.null(xcorr)) NULL else
      xcorr[, .(peak_r = median(peak_r), lag_s = median(lag_s), frac_p05 = mean(p < 0.05), n = .N), by = state],
    n_val = nrow(validation), val_recall = if (nrow(validation)) mean(validation$matched) else NA_real_,
    val_n_extra = sum(recordings$n_rr_extra, na.rm = TRUE),
    val_onset = median_range(matched$err_onset), val_return = median_range(matched$err_return),
    val_frac_within_15s = mean(abs(matched$err_onset) <= 15),
    n_dose = n_or_zero(dose_cor, "SubjectID"),
    dose_rho = if (is.null(dose_cor)) NULL else
      dose_cor[, .(rho_median = median(rho, na.rm = TRUE), n = .N), by = variable]
  )
}


# ---- By type of maneuver ---------------------------------------------------------------------------------------------------
#' Responses by kind of maneuver (slow tilt, rapid tilt, stand-up; "tilt" when the protocol annotations are missing), one row per
#' kind. Counts are events (a recording contributes up to two events per kind); statistics are the median and quartiles over events
#' and are descriptive. `response_events`, `bp_events`: per-event tables with a `kind` column.
prcp_kind_summary <- function(response_events, bp_events = NULL, lead_events = NULL) {
  if (is.null(response_events) || !nrow(response_events) || !"kind" %in% names(response_events)) return(NULL)
  q3 <- function(x) { x <- x[is.finite(x)]; if (!length(x)) rep(NA_real_, 3) else stats::quantile(x, c(0.5, 0.25, 0.75), names = FALSE) }
  add <- function(row, name, x) { q <- q3(x); row[[paste0(name, "_med")]] <- q[1]; row[[paste0(name, "_q25")]] <- q[2]; row[[paste0(name, "_q75")]] <- q[3]; row }
  kinds <- intersect(c("slow tilt", "rapid tilt", "tilt", "stand-up"), unique(response_events$kind))
  rows <- lapply(kinds, function(k) {
    r <- response_events[response_events$kind == k]
    b <- if (is.null(bp_events) || !"kind" %in% names(bp_events)) NULL else bp_events[bp_events$kind == k]
    row <- list(kind = k, n_events = nrow(r), n_recordings = uniqueN(r$SubjectID),
                frac_expected_p = mean(r$expected_p), frac_expected_s = mean(r$expected_s), frac_both_expected = mean(r$expected_p & r$expected_s))
    for (v in c("z_d_x_p", "z_d_x_s", "lat_p", "lat_s", "tau_p_resp", "tau_s_resp")) row <- add(row, v, r[[v]])
    row <- add(row, "d_map", if (is.null(b)) NA_real_ else b$d_map)
    row$n_map <- if (is.null(b)) 0L else sum(is.finite(b$d_map))
    as.data.table(row)
  })
  rbindlist(rows, fill = TRUE)
}

#' Paired within-recording contrasts between kinds (exploratory; recording medians over its events of each kind):
#' rapid - slow tilt, and stand-up - tilt (slow and rapid together). Median of the paired differences and Wilcoxon signed-rank
#' p-value (exact when there are no ties), Benjamini-Hochberg across the table.
prcp_kind_contrasts <- function(response_events, variables = c("z_d_x_p", "z_d_x_s", "lat_p", "lat_s")) {
  if (is.null(response_events) || !nrow(response_events) || !"kind" %in% names(response_events)) return(NULL)
  per_recording <- function(kinds) response_events[response_events$kind %in% kinds, lapply(.SD, function(x) if (any(is.finite(x))) median(x[is.finite(x)]) else NA_real_),
                                                   by = SubjectID, .SDcols = variables]
  contrast <- function(label, a, b) {
    A <- per_recording(a); B <- per_recording(b)
    if (!nrow(A) || !nrow(B)) return(NULL)
    both <- merge(A, B, by = "SubjectID", suffixes = c(".a", ".b"))
    rbindlist(lapply(variables, function(v) {
      d <- both[[paste0(v, ".a")]] - both[[paste0(v, ".b")]]
      d <- d[is.finite(d)]
      data.table(contrast = label, variable = v, n_pairs = length(d), median_diff = if (length(d)) median(d) else NA_real_,
                 p = if (length(d) >= 4 && any(d != 0)) suppressWarnings(stats::wilcox.test(d, exact = length(unique(abs(d))) == length(d))$p.value) else NA_real_)
    }))
  }
  out <- rbindlist(list(contrast("rapid - slow tilt", "rapid tilt", "slow tilt"), contrast("stand-up - tilts", "stand-up", c("slow tilt", "rapid tilt", "tilt"))))
  if (!is.null(out) && nrow(out)) out[, p_BH := p.adjust(p, "BH")]
  out
}
