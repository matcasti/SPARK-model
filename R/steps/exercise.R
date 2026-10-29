# ==============================================================================
# R/steps/exercise.R -- Result 5: exercise and recovery (two-minute step test, TMST).
# Called by analysis/09_exercise.R. Only RR is available, so the rest / exercise / recovery
# phases are detected from heart rate (R/signals/events.R); no subject characteristics or clinical
# outcomes exist for this cohort, so no association analyses are performed.
#
#   exercise_settings()           recovery definition and exercise-end rule (auto = benchmark choice)
#   exercise_define_events()      onset and end of exercise from heart rate
#   exercise_rest_response()      rest level, response at the end of exercise, peak excursion
#   exercise_recovery()           recovery times of x_p and x_s (all candidate definitions)
#   exercise_analyse_recording()  chains the steps for one recording
#   exercise_km_table()           Kaplan-Meier median recovery time, with censoring
#   exercise_summarise()          numbers quoted in the manuscript
# ==============================================================================


# ---- Settings -------------------------------------------------------------------------------------
#' Recovery definition and exercise-end rule. "auto" reads the choice made by
#' analysis/07_recovery_benchmark.R (falls back to 'expfit' / 'profile' when it has not run).
exercise_settings <- function(cfg, recovery_choice) {
  setting <- function(configured, benchmark_field, fallback) {
    if (!identical(configured, "auto")) return(configured)
    if (is.null(recovery_choice)) return(fallback)
    recovery_choice$choice[[benchmark_field]]
  }
  if (identical(cfg$exercise$recovery_definition, "auto") && is.null(recovery_choice)) {
    log_msg("recovery_choice.rds not found: run analysis/07_recovery_benchmark.R first; using 'expfit'.")
  }
  list(recovery_definition = setting(cfg$exercise$recovery_definition, "definition", "expfit"),
       end_method = setting(cfg$exercise$end_method, "end_method", "profile"))
}


# ---- Step 1: events --------------------------------------------------------------------------------------
#' Onset and end of exercise from heart rate; `info` carries the model parameters of the fit.
exercise_define_events <- function(states, id, params, end_method, cfg) {
  events <- detect_events(states$RR_Interval, "exercise", cfg$segmentation)
  t_on <- if (is.finite(events$t_onset_profile)) events$t_onset_profile else events$t_onset
  t_end <- switch(end_method, peak = events$t_end_peak, cp = events$t_end_cp, profile = events$t_end_profile)
  info <- data.table(
    SubjectID = id, detected = events$detected, K = events$K, t_onset = t_on, t_end = t_end,
    profile_r2 = if (!is.null(events$profile)) events$profile$r2 else NA_real_,
    tau_p = 1 / params$kp, tau_s = 1 / params$ks, lambda_P = params$lp, nu0 = params$nu0,
    sigma2 = params$sig2
  )
  list(info = info, t_on = t_on, t_end = t_end)
}


# ---- Step 2: rest level and response --------------------------------------------------------------------------
exercise_rest_response <- function(states, id, t_on, t_end, cfg) {
  ex <- cfg$exercise
  rest <- states[Time >= ex$rest_skip_s & Time < t_on - 5]
  end_of_exercise <- states[Time >= t_end - ex$end_window_s & Time <= t_end]
  smoothed <- function(series) to_grid(states$Time, series, smooth_s = ex$smooth_s)
  hr_smoothed <- smoothed(60 / states$RR_Interval)

  response <- data.table(
    SubjectID = id, rest_x_p = mean(rest$X_p), rest_x_s = mean(rest$X_s),
    rest_sd_x_p = sd(rest$X_p), rest_sd_x_s = sd(rest$X_s),
    rest_hr = 60 / mean(rest$RR_Interval), peak_hr = max(hr_smoothed$x),
    d_x_p = mean(end_of_exercise$X_p) - mean(rest$X_p),
    d_x_s = mean(end_of_exercise$X_s) - mean(rest$X_s)
  )
  response[, `:=`(z_d_x_p = d_x_p / rest_sd_x_p, z_d_x_s = d_x_s / rest_sd_x_s, d_hr = peak_hr - rest_hr)]

  during_exercise <- function(grid) grid[t >= t_on & t <= t_end, x]
  response[, `:=`(peak_d_x_p = min(during_exercise(smoothed(states$X_p))) - rest_x_p,
                  peak_d_x_s = max(during_exercise(smoothed(states$X_s))) - rest_x_s)]
  response
}


# ---- Step 3: recovery ------------------------------------------------------------------------------------------------
#' Recovery metrics of x_p and x_s under every candidate definition (see R/signals/state_summaries.R).
exercise_recovery <- function(states, id, t_on, t_end, cfg) {
  rest_window <- c(cfg$exercise$rest_skip_s, t_on - 5)
  per_state <- lapply(c("p", "s"), function(state) {
    metrics <- recovery_metrics(states$Time, states[[paste0("X_", state)]], rest_window, t_end, cfg$exercise)
    if (is.null(metrics)) NULL else cbind(SubjectID = id, state = state, metrics)
  })
  rbindlist(per_state)
}


# ---- One recording -------------------------------------------------------------------------------------------------------
#' @return list(info, base, recovery, geometry, series); only `info` (usable = FALSE) when the
#'         events could not be defined or the rest period is too short.
exercise_analyse_recording <- function(id, states, params, settings, cfg) {
  recording <- states[SubjectID == id]
  ev <- exercise_define_events(recording, id, params[[id]], settings$end_method, cfg)
  rest_too_short <- ev$t_on - cfg$exercise$rest_skip_s < 30
  if (!is.finite(ev$t_on) || !is.finite(ev$t_end) || rest_too_short) {
    return(list(info = ev$info[, usable := FALSE]))
  }

  response <- exercise_rest_response(recording, id, ev$t_on, ev$t_end, cfg)
  recovery <- exercise_recovery(recording, id, ev$t_on, ev$t_end, cfg)
  geometry <- path_geometry(recording$Time, recording$X_p, recording$X_s,
                            response$rest_sd_x_p, response$rest_sd_x_s)
  list(
    info = ev$info[, usable := TRUE], base = response, recovery = recovery,
    geometry = if (!is.null(geometry)) cbind(SubjectID = id, geometry),
    series = list(t = recording$Time, x_p = recording$X_p, x_s = recording$X_s,
                  vp = recording$X_p_var, vs = recording$X_s_var, rr = recording$RR_Interval,
                  t_on = ev$t_on, t_end = ev$t_end,
                  sd_p = response$rest_sd_x_p, sd_s = response$rest_sd_x_s)
  )
}


# ---- Across recordings -------------------------------------------------------------------------------------------------------
#' Kaplan-Meier median recovery time per state (recordings that did not recover are right-censored
#' at the end of the recording; recordings without an excursion are excluded).
exercise_km_table <- function(recovery, recovery_definition) {
  chosen <- recovery[definition == recovery_definition]
  per_state <- lapply(c("p", "s"), function(which_state) {
    rows <- chosen[state == which_state & status != "no_excursion"]
    if (nrow(rows) < 3) return(NULL)
    time <- ifelse(rows$censored, rows$cens_time, rows$time)
    km <- km_median(time, as.integer(!rows$censored))
    data.table(state = which_state, n = nrow(rows), frac_censored = mean(rows$censored),
               median = km[["median"]], lo = km[["lo"]], hi = km[["hi"]])
  })
  rbindlist(per_state)
}

#' Spearman agreement between the candidate recovery definitions (sympathetic-like state; descriptive).
exercise_definition_agreement <- function(recovery) {
  wide <- dcast(recovery[state == "s"], SubjectID ~ definition, value.var = "time")
  if (!all(RECOVERY_DEFINITIONS %in% names(wide)) || nrow(wide) < 4) return(NULL)
  suppressWarnings(cor(wide[, ..RECOVERY_DEFINITIONS], method = "spearman", use = "pairwise.complete.obs"))
}

exercise_summarise <- function(info, base, geometry, km, agreement, settings) {
  list(
    n = nrow(info), n_usable = sum(info$usable),
    recovery_definition = settings$recovery_definition, end_method = settings$end_method,
    z_d_x_p = median_range(base$z_d_x_p), z_d_x_s = median_range(base$z_d_x_s), d_hr = median_range(base$d_hr),
    path_length = median_range(geometry$path_length), loop_area = median_range(geometry$loop_area),
    U_max = median_range(geometry$U_max), km = km, definition_agreement = agreement
  )
}
