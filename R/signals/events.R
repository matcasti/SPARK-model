# ==============================================================================
# R/signals/events.R -- event detection from heart rate alone (rest / perturbation / recovery).
# Piecewise-constant change points on binned heart rate, chosen by BIC with a minimum-shift rule;
# for exercise, also a parametric rise-and-decay profile. Detected times are estimates, not truth.
# Used where no event information exists (TMST, PRCP recordings without a tilt angle) and validated
# against the tilt angle where both exist (analysis/05_prcp.R) and on simulations (analysis/08_*.R).
#
#   detect_events()          onset / recovery / end-of-exercise candidates for one recording
#                            (tilt mode: every perturbation of the recording)
#   pair_perturbations()     pairs the upward and downward change points into perturbation events
#   match_events()           one-to-one matching of detected to reference event times
#   tilt_windows()           baseline / tilt / recovery window of every tilt event
#   assign_phase()           phase label of each beat time (one event);  assign_tilt_phases(): all events
#   segment_heart_rate()     the change-point segmentation;   fit_piecewise_constant()  its exact solver
#   fit_exercise_profile()   parametric exercise profile
# ==============================================================================


# ---- Binning and change points ---------------------------------------------------------------------------
#' Heart rate (Hz) interpolated at 1 Hz and averaged in non-overlapping bins of `bin_s` seconds.
bin_heart_rate <- function(rr, bin_s = 5) {
  beat_time <- cumsum(rr)
  grid <- seq(ceiling(beat_time[1]), floor(beat_time[length(beat_time)]), by = 1)
  heart_rate <- approx(beat_time, 1 / rr, xout = grid, rule = 2)$y
  bin <- floor((grid - grid[1]) / bin_s)
  bins <- data.table(b = bin, t = grid, hr = heart_rate)[, .(t = mean(t), hr = mean(hr), n = .N), by = b]
  bins[n >= max(1, bin_s %/% 2)][, n := NULL][]
}

#' Optimal K-segment piecewise-constant fit by dynamic programming (exact).
#' @return list(sse, starts) with `starts` the first index of each segment; NULL if the series is too short
fit_piecewise_constant <- function(y, K, min_len) {
  n <- length(y)
  if (n < K * min_len) return(NULL)
  cum_y <- c(0, cumsum(y))
  cum_y2 <- c(0, cumsum(y^2))
  segment_sse <- function(i, j) {                       # sum of squares about the mean of y[i..j]
    m <- j - i + 1
    total <- cum_y[j + 1] - cum_y[i]
    (cum_y2[j + 1] - cum_y2[i]) - total^2 / m
  }
  best <- matrix(Inf, K, n)
  best_start <- matrix(NA_integer_, K, n)
  last_ends <- min_len:n
  best[1, last_ends] <- segment_sse(1, last_ends)
  if (K > 1) {
    for (k in 2:K) for (j in (k * min_len):n) {
      starts <- ((k - 1) * min_len + 1):(j - min_len + 1)
      if (!length(starts)) next
      cost <- best[k - 1, starts - 1] + segment_sse(starts, j)
      winner <- which.min(cost)
      best[k, j] <- cost[winner]
      best_start[k, j] <- starts[winner]
    }
  }
  starts <- integer(K)                                  # backtrack
  j <- n
  for (k in K:1) {
    if (k == 1) { starts[1] <- 1L; break }
    starts[k] <- best_start[k, j]
    j <- starts[k] - 1L
  }
  list(sse = best[K, n], starts = starts)
}

#' Segment heart rate with K chosen by BIC (K = 1 .. max_changepoints + 1), subject to the
#' minimum-shift rule: adjacent segment means must differ by at least `min_shift_hz`.
segment_heart_rate <- function(rr, seg_cfg, max_changepoints = seg_cfg$max_changepoints) {
  bins <- bin_heart_rate(rr, seg_cfg$bin_s)
  y <- bins$hr
  n <- length(y)
  min_len <- max(2L, ceiling(seg_cfg$min_seg_s / seg_cfg$bin_s))
  fits <- lapply(seq_len(max_changepoints + 1L), function(K) fit_piecewise_constant(y, K, min_len))
  feasible <- which(!vapply(fits, is.null, logical(1)))
  bic <- vapply(feasible, function(K) n * log(max(fits[[K]]$sse, 1e-12) / n) + (2 * K - 1) * log(n), numeric(1))

  segment_means <- function(K) {
    starts <- fits[[K]]$starts
    vapply(seq_len(K), function(k) mean(y[starts[k]:(if (k < K) starts[k + 1] - 1L else n)]), numeric(1))
  }
  chosen <- 1L
  for (K in feasible[order(bic)]) {                     # best BIC first, subject to the shift rule
    shifts <- abs(diff(segment_means(K)))
    if (K == 1L || all(shifts >= seg_cfg$min_shift_hz)) { chosen <- K; break }
  }
  starts <- fits[[chosen]]$starts
  change_time <- if (chosen > 1) (bins$t[starts[-1] - 1L] + bins$t[starts[-1]]) / 2 else numeric()
  list(K = chosen, cp_t = change_time, means = segment_means(chosen), bins = bins, bic = setNames(bic, feasible))
}


# ---- Parametric exercise profile -----------------------------------------------------------------------------
#' Shape g(t): 0 before t_on, 1 - exp(-(t - t_on)/tau_on) until t_off, then exponential decay with tau_off.
exercise_profile_shape <- function(t, t_on, t_off, tau_on, tau_off) {
  level_at_end <- 1 - exp(-(t_off - t_on) / tau_on)
  ifelse(t < t_on, 0,
         ifelse(t < t_off, 1 - exp(-(t - t_on) / tau_on), level_at_end * exp(-(t - t_off) / tau_off)))
}

#' Fit HR(t) = h0 + A g(t) to binned heart rate: h0 and A by least squares (variable projection),
#' (t_on, t_off, tau_on, tau_off) by a coarse grid search refined with Nelder-Mead.
fit_exercise_profile <- function(bins, min_ex_s = 40, max_ex_s = 360) {
  t <- bins$t
  y <- bins$hr
  span <- t[length(t)] - t[1]

  # Residual sum of squares (and h0, A) for p = (t_on, t_off, ln tau_on, ln tau_off)
  profile_fit <- function(p) {
    t_on <- p[1]; t_off <- p[2]
    duration <- t_off - t_on
    if (duration < min_ex_s || duration > max_ex_s || t_on < t[1] || t_off > t[length(t)] + 1e-9) {
      return(list(sse = Inf))
    }
    shape <- exercise_profile_shape(t, t_on, t_off, exp(p[3]), exp(p[4]))
    design <- cbind(1, shape)
    coefficients <- tryCatch(solve(crossprod(design), crossprod(design, y)), error = function(e) NULL)
    if (is.null(coefficients) || coefficients[2] <= 0) return(list(sse = Inf))
    list(sse = sum((y - design %*% coefficients)^2), h0 = coefficients[1], A = coefficients[2])
  }

  step <- max(5, round(span / 80))                      # coarse grid: about 80 points over the recording
  grid <- expand.grid(t_on = seq(t[1] + 10, t[length(t)] - min_ex_s, by = step),
                      dt = seq(min_ex_s, min(max_ex_s, span), by = step),
                      ton = log(c(5, 12, 25, 50)), toff = log(c(10, 25, 60, 150)))
  grid$t_off <- grid$t_on + grid$dt
  grid <- grid[grid$t_off <= t[length(t)], ]
  if (!nrow(grid)) return(NULL)

  grid_sse <- vapply(seq_len(nrow(grid)), function(i) {
    profile_fit(c(grid$t_on[i], grid$t_off[i], grid$ton[i], grid$toff[i]))$sse
  }, numeric(1))
  if (!any(is.finite(grid_sse))) return(NULL)
  best <- grid[which.min(grid_sse), ]
  start <- c(best$t_on, best$t_off, best$ton, best$toff)

  refined <- optim(start, function(p) { v <- profile_fit(p)$sse; if (is.finite(v)) v else 1e12 },
                   method = "Nelder-Mead", control = list(maxit = 300))
  fit <- profile_fit(refined$par)
  if (!is.finite(fit$sse)) { fit <- profile_fit(start); refined$par <- start }
  list(t_on = refined$par[1], t_off = refined$par[2], tau_on = exp(refined$par[3]), tau_off = exp(refined$par[4]),
       h0 = fit$h0, A = fit$A, r2 = 1 - fit$sse / sum((y - mean(y))^2))
}


# ---- Events of one recording --------------------------------------------------------------------------------------
#' Pair the change points of a segmentation into perturbation events. Walking through time, an upward shift
#' starts an event; a downward shift ends it only if heart rate falls back at least halfway to the pre-event
#' level (a dip inside a perturbation does not split it). A further upward shift inside an event (stairs) and a
#' downward shift outside one are ignored. A last event that never ends has `t_recovery = NA`.
#' @return data.table(event, t_onset, t_recovery, shift_hz) (zero rows if no upward shift)
pair_perturbations <- function(change_t, means) {
  rows <- list()
  open <- FALSE
  for (j in seq_along(change_t)) {
    if (!open && means[j + 1] > means[j]) {
      open <- TRUE
      onset <- change_t[j]
      pre_level <- means[j]
      top_level <- means[j + 1]
    } else if (open && means[j + 1] > top_level) {
      top_level <- means[j + 1]
    } else if (open && means[j + 1] < means[j] && means[j + 1] <= pre_level + 0.5 * (top_level - pre_level)) {
      rows[[length(rows) + 1L]] <- data.table(t_onset = onset, t_recovery = change_t[j], shift_hz = top_level - pre_level)
      open <- FALSE
    }
  }
  if (open) rows[[length(rows) + 1L]] <- data.table(t_onset = onset, t_recovery = NA_real_, shift_hz = top_level - pre_level)
  if (!length(rows)) return(data.table(event = integer(), t_onset = numeric(), t_recovery = numeric(), shift_hz = numeric()))
  cbind(event = seq_along(rows), rbindlist(rows))
}

#' Detect rest / perturbation / recovery events from the RR series.
#'  mode "tilt":     every perturbation (`$events`): up to `max_tilt_events` rises of heart rate, each ended by a
#'                   later fall. `t_onset` / `t_recovery` repeat the first one.
#'  mode "exercise": one bout: onset as above; the end of exercise is returned under three rules so that the rule
#'                   can be benchmarked (analysis/07_recovery_benchmark.R): the smoothed heart-rate peak
#'                   after the onset (t_end_peak), the next change point (t_end_cp), the fitted
#'                   profile's offset (t_end_profile).
#' Heart-rate rises alone cannot tell a tilt from another manoeuvre: the events are candidates.
detect_events <- function(rr, mode = c("tilt", "exercise"), seg_cfg) {
  mode <- match.arg(mode)
  max_changepoints <- if (mode == "tilt") 2L * seg_cfg$max_tilt_events else seg_cfg$max_changepoints
  segmentation <- segment_heart_rate(rr, seg_cfg, max_changepoints)
  none <- data.table(event = integer(), t_onset = numeric(), t_recovery = numeric(), shift_hz = numeric())
  events <- list(mode = mode, K = segmentation$K, t_onset = NA_real_, t_recovery = NA_real_, t_end_cp = NA_real_,
                 t_end_peak = NA_real_, t_onset_profile = NA_real_, t_end_profile = NA_real_, profile = NULL,
                 detected = FALSE, n_events = 0L, events = none, seg = segmentation)
  if (segmentation$K < 2) return(events)
  paired <- pair_perturbations(segmentation$cp_t, segmentation$means)
  if (mode == "exercise") paired <- paired[seq_len(min(.N, 1L))]
  if (!nrow(paired)) return(events)

  events$events <- paired
  events$n_events <- nrow(paired)
  events$t_onset <- paired$t_onset[1]
  events$t_recovery <- paired$t_recovery[1]
  events$t_end_cp <- paired$t_recovery[1]
  events$detected <- TRUE

  bins <- segmentation$bins
  window <- seg_cfg$smooth_bins
  smoothed <- as.numeric(stats::filter(bins$hr, rep(1 / window, window), sides = 2))
  smoothed[is.na(smoothed)] <- bins$hr[is.na(smoothed)]
  after_onset <- which(bins$t >= events$t_onset)
  if (length(after_onset)) events$t_end_peak <- bins$t[after_onset[which.max(smoothed[after_onset])]]

  if (mode == "exercise") {
    profile <- tryCatch(fit_exercise_profile(bins, seg_cfg$min_exercise_s, seg_cfg$max_exercise_s),
                        error = function(e) NULL)
    if (!is.null(profile)) {
      events$profile <- profile
      events$t_onset_profile <- profile$t_on
      events$t_end_profile <- profile$t_off
    }
  }
  events
}

#' One-to-one matching of detected to reference event times: closest pairs first, each time used once,
#' pairs further apart than `tol_s` never match.
#' @return data.table(reference, detected, error): index of the matched detected time (NA if none) and
#'         detected - reference (s)
match_events <- function(reference, detected, tol_s) {
  result <- data.table(reference = seq_along(reference), detected = NA_integer_, error = NA_real_)
  if (!length(reference) || !length(detected)) return(result)
  pairs <- CJ(i_reference = seq_along(reference), i_detected = seq_along(detected))
  pairs[, error := detected[i_detected] - reference[i_reference]]
  pairs <- pairs[is.finite(error) & abs(error) <= tol_s][order(abs(error))]
  taken_reference <- integer()
  taken_detected <- integer()
  for (i in seq_len(nrow(pairs))) {
    if (pairs$i_reference[i] %in% taken_reference || pairs$i_detected[i] %in% taken_detected) next
    set(result, pairs$i_reference[i], c("detected", "error"), list(pairs$i_detected[i], pairs$error[i]))
    taken_reference <- c(taken_reference, pairs$i_reference[i])
    taken_detected <- c(taken_detected, pairs$i_detected[i])
  }
  result[]
}


# ---- Windows of repeated tilts ------------------------------------------------------------------------------------
#' Baseline / tilt / recovery window of every tilt event of a recording (times in the recording's clock).
#'  baseline: from `washout_s` after the previous return (the recording start for the first event) to the onset,
#'            at most `baseline_max_s` long;   tilt: onset to return (to the end if the tilt never returns);
#'  recovery: return to `washout_s` later, or the next onset if earlier.
#' An event is `usable` when its baseline lasts at least `baseline_min_s`.
#' @param onset,return_ event times (return_ may be NA for a last tilt that never returns)
#' @return data.table(event, t_base_from, t_onset, t_return, t_tilt_end, t_rec_to, baseline_s, usable)
tilt_windows <- function(onset, return_, t_start, t_end, tilt_cfg) {
  n <- length(onset)
  if (!n) return(data.table(event = integer(), t_base_from = numeric(), t_onset = numeric(), t_return = numeric(),
                            t_tilt_end = numeric(), t_rec_to = numeric(), baseline_s = numeric(), usable = logical()))
  tilt_end <- ifelse(is.finite(return_), return_, t_end)
  previous_end <- c(-Inf, tilt_end[-n])
  base_from <- pmax(t_start + tilt_cfg$edge_skip_s, previous_end + tilt_cfg$washout_s, onset - tilt_cfg$baseline_max_s)
  next_onset <- c(onset[-1], Inf)
  rec_to <- pmin(tilt_end + tilt_cfg$washout_s, next_onset, t_end)
  rec_to[!is.finite(return_)] <- tilt_end[!is.finite(return_)]
  data.table(event = seq_len(n), t_base_from = base_from, t_onset = onset, t_return = return_, t_tilt_end = tilt_end,
             t_rec_to = rec_to, baseline_s = onset - base_from, usable = (onset - base_from) >= tilt_cfg$baseline_min_s)
}

#' Phase of every beat time with respect to all usable tilt events: baseline / tilt / recovery inside the
#' event windows, NA elsewhere (washout and manoeuvres that are not analysed).
#' @param windows tilt_windows() result
#' @return data.table(phase, event)
assign_tilt_phases <- function(t, windows) {
  phase <- rep(NA_character_, length(t))
  event <- rep(NA_integer_, length(t))
  for (i in which(windows$usable)) {
    w <- windows[i]
    label <- function(from, to, name) {
      inside <- t >= from & t < to
      phase[inside] <<- name
      event[inside] <<- w$event
    }
    label(w$t_base_from, w$t_onset, "baseline")
    label(w$t_onset, w$t_tilt_end, "tilt")
    label(w$t_tilt_end, w$t_rec_to, "recovery")
  }
  data.table(phase = factor(phase, levels = c("baseline", "tilt", "recovery")), event = event)
}

#' Phase label of each beat time. Phases: baseline / tilt / recovery, or rest / exercise / recovery
#' (the end of exercise follows `end_method`: "peak", "cp" or "profile").
assign_phase <- function(t, events, end_method = "peak") {
  labels <- if (events$mode == "tilt") c("baseline", "tilt", "recovery") else c("rest", "exercise", "recovery")
  t_end <- if (events$mode == "tilt") events$t_recovery else {
    switch(end_method, cp = events$t_end_cp, profile = events$t_end_profile, events$t_end_peak)
  }
  if (!isTRUE(events$detected)) return(factor(rep(NA_character_, length(t)), levels = labels))
  phase <- rep(labels[1], length(t))
  phase[t >= events$t_onset] <- labels[2]
  if (is.finite(t_end)) phase[t >= t_end] <- labels[3]
  factor(phase, levels = labels)
}


# ---- Protocol annotations (PRCP) ---------------------------------------------------------------------------
#' Maneuvers of one recording from the events table (data/prcp_events.csv); NULL when there are none.
annotated_maneuvers <- function(events, subject) {
  if (is.null(events) || !nrow(events)) return(NULL)
  rows <- which(events$id == subject)                                # (a column is named `id`: do not use a variable of that name inside [ ])
  if (length(rows)) events[rows] else NULL
}

#' Kind of every angle-derived tilt: the annotated tilt whose ascent ([onset_s, onset_end_s], +- slack) contains the onset found
#' from the angle ("tilt" when no annotation matches or there are no annotations). Stand-ups are never matched: the angle does not change.
label_tilt_kinds <- function(onset, maneuvers, slack_s = 5) {
  kinds <- rep("tilt", length(onset))
  if (is.null(maneuvers)) return(kinds)
  tilts <- maneuvers[maneuvers$kind != "stand-up", ]
  if (!nrow(tilts)) return(kinds)
  ascent_end <- ifelse(is.finite(tilts$onset_end_s), tilts$onset_end_s, tilts$onset_s + 60)
  for (i in seq_along(onset)) {
    hit <- which(onset[i] >= tilts$onset_s - slack_s & onset[i] <= ascent_end + slack_s)
    if (length(hit)) kinds[i] <- tilts$kind[hit[1]]
  }
  kinds
}
