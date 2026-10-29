# ==============================================================================
# R/model/simulation.R -- generative simulator of the SPARK model (OU drives -> integrate-and-fire
# phase) with deterministic drive profiles that mimic perturbations. Used for recovery studies,
# calibration checks, method benchmarks and the test suite. (The order of the random draws is
# part of the seeded results: do not reorder them.)
#
#   simulate_rr()   core simulator: latent drives -> beat times, RR, true drives at the beats
#   sim_rest()      resting recording, young- or older-like parameters
#   sim_tilt()      tilt-like step response;  sim_tilt_angle()  drives follow a tilt-angle protocol
#   sim_exercise()  rest -> exercise -> recovery;  draw_exercise_pars()  its true profile
#   prof_step(), prof_exercise(), drive_from_angle(), angle_protocol(), rsa_drive()   drive profiles
# ==============================================================================


# ---- Core simulator -----------------------------------------------------------------------------------
#' Simulate RR intervals (s) with known latent drives.
#'
#' The net pacing rate is nu(t) = nu0 - (x_p + m_p(t)) + (x_s + m_s(t)); the phase accumulates it
#' plus threshold noise of variance sig2 * lR per second and resets at each beat.
#' @param m_p,m_s  optional functions of time (s) giving deterministic drive offsets (Hz)
#' @return data.table(beat, t, rr, x_p, x_s, x_p_ou, x_s_ou); attribute "truth" holds the parameters
simulate_rr <- function(duration_s, nu0 = 1.2, kp = 1 / 4, ks = 1 / 150, sig2 = 9e-4, lp = 2.5,
                        lR = 0.02, m_p = NULL, m_s = NULL, h = 0.01, seed = NULL, nu_min = 0.3) {
  if (!is.null(seed)) set.seed(seed)
  n <- ceiling(duration_s / h)
  time <- (0:n) * h

  # Exact OU transition on the grid: x[t + h] = a x[t] + sqrt(V (1 - a^2)) eps
  simulate_ou <- function(k, variance) {
    a <- exp(-k * h)
    noise <- rnorm(n + 1, 0, sqrt(variance * (1 - a^2)))
    as.numeric(stats::filter(noise, a, method = "recursive", init = rnorm(1, 0, sqrt(variance))))
  }
  x_p <- simulate_ou(kp, lp * sig2)
  x_s <- simulate_ou(ks, sig2)
  x_p_total <- x_p + if (is.null(m_p)) 0 else m_p(time)
  x_s_total <- x_s + if (is.null(m_s)) 0 else m_s(time)

  rate <- pmax(nu0 - x_p_total + x_s_total, nu_min)
  increment <- rate * h + sqrt(sig2 * lR * h) * rnorm(n + 1)
  phase <- cumsum(increment)

  # Beat k: first grid time at which the phase accumulated since the last beat reaches 1
  beats <- numeric()
  index <- 1L
  phase_at_last_beat <- 0
  look_ahead <- ceiling(6 / h)
  repeat {
    window_end <- min(index + look_ahead, n + 1L)
    if (index >= window_end) break
    window <- (index + 1L):window_end
    first_crossing <- which(phase[window] - phase_at_last_beat >= 1)[1]
    if (is.na(first_crossing)) break
    i <- window[first_crossing]
    fraction <- min(max((1 - (phase[i - 1L] - phase_at_last_beat)) / increment[i], 0), 1)
    beats <- c(beats, time[i - 1L] + fraction * h)
    phase_at_last_beat <- phase[i - 1L] + fraction * increment[i]
    index <- i
  }

  at_beats <- function(x) approx(time, x, xout = beats, rule = 2)$y
  out <- data.table(beat = seq_along(beats), t = beats, rr = diff(c(0, beats)),
                    x_p = at_beats(x_p_total), x_s = at_beats(x_s_total),
                    x_p_ou = at_beats(x_p), x_s_ou = at_beats(x_s))
  out <- out[-1]
  attr(out, "truth") <- list(nu0 = nu0, kp = kp, ks = ks, sig2 = sig2, lp = lp, lR = lR)
  out
}


# ---- Drive profiles ----------------------------------------------------------------------------------------
#' Exponential approach to `amp` from `t0` with time constant `tau`.
prof_step <- function(t0, amp, tau) {
  function(t) ifelse(t < t0, 0, amp * (1 - exp(-(t - t0) / tau)))
}

#' Rise towards `amp` between t_on and t_off, exponential decay afterwards.
prof_exercise <- function(t_on, t_off, amp, tau_on, tau_off) {
  amplitude_at_end <- amp * (1 - exp(-(t_off - t_on) / tau_on))
  function(t) {
    ifelse(t < t_on, 0,
           ifelse(t < t_off, amp * (1 - exp(-(t - t_on) / tau_on)),
                  amplitude_at_end * exp(-(t - t_off) / tau_off)))
  }
}

#' First-order low-pass drive m' = (target - m) / tau with target = amp * sin(angle(t)).
#' Returns a function of an equally spaced time vector starting at 0.
drive_from_angle <- function(angle_fun, amp, tau) {
  function(t) {
    h <- t[2] - t[1]
    a <- exp(-h / tau)
    target <- amp * sin(angle_fun(t) * pi / 180)
    as.numeric(stats::filter((1 - a) * target, a, method = "recursive"))
  }
}

#' Tilt-table angle protocols (degrees): "step" 0 -> 70 at t0 (held); "stairs" 0, 30, 60, 80, back to 0;
#' "repeated" `n_tilts` tilts to 70 degrees of `tilt_s` seconds, the first at t0, separated by `rest_s` seconds flat.
angle_protocol <- function(protocol, t0, duration_s, n_tilts = 3L, tilt_s = 180, rest_s = 240) {
  switch(protocol,
    step = function(t) ifelse(t < t0, 0, 70),
    stairs = function(t) {
      ifelse(t < t0, 0, ifelse(t < t0 + 60, 30, ifelse(t < t0 + 120, 60, ifelse(t < t0 + 180, 80, 0))))
    },
    repeated = function(t) {
      onsets <- t0 + (seq_len(n_tilts) - 1L) * (tilt_s + rest_s)
      in_tilt <- Reduce(`|`, lapply(onsets, function(o) t >= o & t < o + tilt_s))
      ifelse(in_tilt, 70, 0)
    },
    stop("unknown protocol")
  )
}

#' Respiratory sinus arrhythmia as a deterministic vagal offset at the respiratory frequency.
rsa_drive <- function(f_resp, amp, phase = 0) function(t) -amp * sin(2 * pi * f_resp * t + phase)


# ---- Recordings ---------------------------------------------------------------------------------------------
#' Resting recording with young- or older-like parameters (Fantasia-like groups).
sim_rest <- function(seed, group = c("young", "older"), duration_s = 600, ...) {
  group <- match.arg(group)
  set.seed(seed)
  parameters <- if (group == "young") {
    list(nu0 = runif(1, 1.05, 1.25), kp = 1 / runif(1, 3, 6), ks = 1 / runif(1, 100, 300),
         sig2 = runif(1, 8e-4, 1.4e-3), lp = runif(1, 2, 4))
  } else {
    list(nu0 = runif(1, 1.0, 1.2), kp = 1 / runif(1, 3.5, 7), ks = 1 / runif(1, 100, 300),
         sig2 = runif(1, 2.5e-4, 5e-4), lp = runif(1, 1, 2.2))
  }
  recording <- do.call(simulate_rr, c(list(duration_s = duration_s, lR = 0.02), parameters, list(...)))
  attr(recording, "group") <- group
  recording
}

#' Tilt-like recording: vagal withdrawal (fast) followed by sympathetic rise (slow), as steps at t0.
sim_tilt <- function(seed, duration_s = 420, kp = 1 / 4, ks = 1 / 150, ...) {
  set.seed(seed)
  t0 <- runif(1, 0.28, 0.45) * duration_s
  amp_p <- runif(1, 0.08, 0.14)
  amp_s <- runif(1, 0.06, 0.12)
  tau_p <- runif(1, 4, 9)
  tau_s <- runif(1, 18, 40)
  recording <- simulate_rr(duration_s, kp = kp, ks = ks, m_p = prof_step(t0, -amp_p, tau_p),
                           m_s = prof_step(t0, amp_s, tau_s), ...)
  attr(recording, "events") <- list(t_onset = t0, amp_p = -amp_p, amp_s = amp_s, tau_p_on = tau_p, tau_s_on = tau_s)
  recording
}

#' Tilt recording whose drives follow a tilt-angle protocol (vagal withdrawal fast, sympathetic slow).
#' attr "angle" is the angle function; attr "events" holds the true onset and return of the first tilt and
#' `tilts`, the true onset / return of every tilt (protocol "repeated": `n_tilts` tilts, the first at 150-180 s;
#' choose `duration_s` of at least 150 + n_tilts * (tilt_s + rest_s)).
sim_tilt_angle <- function(seed, duration_s = 600, protocol = c("step", "stairs", "repeated", "mixed"), kp = 1 / 4, ks = 1 / 150,
                           n_tilts = 3L, tilt_s = 180, rest_s = 240, stand_s = 120, ...) {
  protocol <- match.arg(protocol)
  set.seed(seed)
  t0 <- if (protocol %in% c("repeated", "mixed")) runif(1, 150, 180) else runif(1, 0.25, 0.35) * duration_s
  amp_p <- runif(1, 0.10, 0.18)
  amp_s <- runif(1, 0.08, 0.15)
  tau_p <- runif(1, 4, 9)
  tau_s <- runif(1, 18, 40)
  # "mixed": slow tilt, stand-up (the drives respond, the table angle stays at 0), rapid tilt; `maneuvers` are the protocol annotations
  starts <- t0 + (0:2) * (tilt_s + rest_s)
  maneuvers <- if (protocol == "mixed") data.table(kind = c("slow tilt", "stand-up", "rapid tilt"), onset_s = starts,
                                                   onset_end_s = c(starts[1] + 50, NA, starts[3] + 2),
                                                   return_s = c(starts[1] + tilt_s, starts[2] + stand_s, starts[3] + tilt_s),
                                                   return_end_s = c(starts[1] + tilt_s + 50, NA, starts[3] + tilt_s + 2)) else NULL
  if (protocol == "mixed") {
    in_interval <- function(t, from, to) t >= from & t < to
    angle <- function(t) ifelse(in_interval(t, starts[1], starts[1] + tilt_s) | in_interval(t, starts[3], starts[3] + tilt_s), 70, 0)    # table angle
    drive_angle <- function(t) ifelse(angle(t) > 0 | in_interval(t, starts[2], starts[2] + stand_s), 70, 0)                          # what the drives respond to
  } else {
    angle <- angle_protocol(protocol, t0, duration_s, n_tilts, tilt_s, rest_s)
    drive_angle <- angle
  }
  recording <- simulate_rr(duration_s, kp = kp, ks = ks, m_p = drive_from_angle(drive_angle, -amp_p, tau_p),
                           m_s = drive_from_angle(drive_angle, amp_s, tau_s), ...)
  tilts <- switch(protocol,
    step = data.table(t_onset = t0, t_return = NA_real_),
    stairs = data.table(t_onset = t0, t_return = t0 + 180),
    repeated = data.table(t_onset = t0 + (seq_len(n_tilts) - 1L) * (tilt_s + rest_s),
                          t_return = t0 + (seq_len(n_tilts) - 1L) * (tilt_s + rest_s) + tilt_s),
    mixed = data.table(t_onset = starts[c(1, 3)], t_return = starts[c(1, 3)] + tilt_s))
  attr(recording, "angle") <- angle
  attr(recording, "events") <- list(t_onset = t0, t_return = tilts$t_return[1],
                                    amp_p = -amp_p, amp_s = amp_s, tau_p_on = tau_p, tau_s_on = tau_s,
                                    protocol = protocol, tilts = tilts, maneuvers = maneuvers)
  recording
}

#' Draw the true drive profile of a rest - exercise - recovery recording.
draw_exercise_pars <- function(seed, duration_s = 420) {
  set.seed(seed)
  t_on <- runif(1, 0.18, 0.26) * duration_s
  list(t_on = t_on, t_off = t_on + runif(1, 100, 130),
       amp_p = runif(1, 0.25, 0.40), amp_s = runif(1, 0.35, 0.60),
       ton_p = runif(1, 6, 12), ton_s = runif(1, 15, 30),
       toff_p = exp(runif(1, log(12), log(60))), toff_s = exp(runif(1, log(25), log(100))))
}

#' Rest - exercise - recovery recording with known onset, end and recovery constants. `pars` (from
#' draw_exercise_pars) can be reused with a different `noise_seed` to realize the same profile twice
#' (test-retest).
sim_exercise <- function(seed, duration_s = 420, kp = 1 / 4, ks = 1 / 150, pars = NULL, noise_seed = seed, ...) {
  if (is.null(pars)) pars <- draw_exercise_pars(seed, duration_s)
  set.seed(noise_seed)
  recording <- simulate_rr(
    duration_s, kp = kp, ks = ks,
    m_p = prof_exercise(pars$t_on, pars$t_off, -pars$amp_p, pars$ton_p, pars$toff_p),
    m_s = prof_exercise(pars$t_on, pars$t_off, pars$amp_s, pars$ton_s, pars$toff_s), ...)
  attr(recording, "events") <- list(t_onset = pars$t_on, t_end = pars$t_off,
                                    tau_rec_p = pars$toff_p, tau_rec_s = pars$toff_s,
                                    amp_p = -pars$amp_p, amp_s = pars$amp_s)
  recording
}
