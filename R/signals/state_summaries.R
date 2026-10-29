# ==============================================================================
# R/signals/state_summaries.R -- summaries of the reconstructed latent states.
#
#   to_grid()            resample a beat-wise series to a regular grid (optionally smoothed)
#   phase_summary()      mean / SD of the states by detected phase
#   fit_step_response()  x(t) = x0 + D (1 - exp(-(t - t_on)/tau)): onset, time constant, amplitude
#   recovery_metrics()   recovery time under several candidate definitions (RECOVERY_DEFINITIONS)
#   path_geometry()      length, loop area and energy of the (x_p, x_s) path
# ==============================================================================


# ---- Resampling and phase means --------------------------------------------------------------------------------
#' Resample a series onto a regular grid (default 1 Hz) and smooth it with a centered moving average.
to_grid <- function(t, x, by = 1, smooth_s = 0) {
  grid <- seq(ceiling(t[1]), floor(t[length(t)]), by = by)
  y <- approx(t, x, xout = grid, rule = 2)$y
  window <- max(1L, round(smooth_s / by))
  if (window > 1) {
    smoothed <- as.numeric(stats::filter(y, rep(1 / window, window), sides = 2))
    y <- ifelse(is.na(smoothed), y, smoothed)
  }
  data.table(t = grid, x = y)
}

#' Mean and SD of the states by detected phase; `d` has columns t, x_p, x_s, phase.
phase_summary <- function(d) {
  d[!is.na(phase), .(n = .N, mean_x_p = mean(x_p), mean_x_s = mean(x_s), sd_x_p = sd(x_p), sd_x_s = sd(x_s)),
    by = phase]
}


# ---- Step response ----------------------------------------------------------------------------------------------
#' Step-response fit x(t) = x0 + D (1 - exp(-(t - t_on)/tau)) for t >= t_on. The onset t_on is
#' searched within +/- `search_s` of `t_ref` and tau on a log grid; x0 and D are least squares.
#' `at_bound` flags fits at the edge of the search ranges.
#' @return list(t_on, tau, x0, D, r2, at_bound), or NULL if no candidate could be fitted
fit_step_response <- function(t, x, t_ref, search_s = 40, tau_range = c(2, 250), n_tau = 25,
                              step_s = 2, smooth_s = 5) {
  grid <- to_grid(t, x, smooth_s = smooth_s)
  onsets <- seq(t_ref - search_s, t_ref + search_s, by = step_s)
  taus <- exp(seq(log(tau_range[1]), log(tau_range[2]), length.out = n_tau))

  best <- list(sse = Inf)
  for (onset in onsets) for (tau in taus) {
    shape <- ifelse(grid$t < onset, 0, 1 - exp(-(grid$t - onset) / tau))
    if (sd(shape) < 1e-6) next
    design <- cbind(1, shape)
    coefficients <- solve(crossprod(design), crossprod(design, grid$x))
    sse <- sum((grid$x - design %*% coefficients)^2)
    if (sse < best$sse) best <- list(sse = sse, t_on = onset, tau = tau, x0 = coefficients[1], D = coefficients[2])
  }
  if (!is.finite(best$sse)) return(NULL)

  best$r2 <- 1 - best$sse / sum((grid$x - mean(grid$x))^2)
  best$at_bound <- best$tau <= tau_range[1] * 1.001 || best$tau >= tau_range[2] * 0.999 ||
                   abs(best$t_on - t_ref) >= search_s - 1e-9
  best$sse <- NULL
  best
}


# ---- Recovery definitions ---------------------------------------------------------------------------------------------
RECOVERY_DEFINITIONS <- c("sustained", "first", "frac90", "expfit")

#' Common setting of the definitions: the smoothed state `grid` (t, x), its rest level and SD, the end
#' of exercise. The excursion `dz` is signed to be positive (direction of the end-of-exercise change).
recovery_setup <- function(grid, rest_mean, rest_sd, t_end, settings) {
  end_window <- grid[t >= t_end - settings$end_window_s & t <= t_end]
  if (!nrow(end_window)) return(NULL)
  direction <- sign(median(end_window$x) - rest_mean)
  if (direction == 0) direction <- 1
  list(dz = direction * (grid$x - rest_mean), t = grid$t, A_end = direction * (median(end_window$x) - rest_mean),
       sd0 = rest_sd, t_end = t_end, dirn = direction)
}

#' Time (s after the end of exercise) at which a recovery criterion is met; otherwise the series is
#' right-censored at the end of the recording. Definitions:
#'   sustained  within tol_sd rest SDs for hold_s seconds     first   first time within tol_sd rest SDs
#'   frac90     excursion below 10% of its end-of-exercise value
#'   expfit     time constant of an exponential decay fitted to the excursion
recovery_time <- function(setup, definition, settings) {
  after_end <- setup$t >= setup$t_end
  time <- setup$t[after_end]
  dz <- setup$dz[after_end]
  censor_time <- tail(time, 1) - setup$t_end
  censored <- list(time = NA_real_, censored = TRUE, status = "censored", cens_time = censor_time)
  if (setup$A_end <= settings$tol_sd * setup$sd0) {
    return(list(time = NA_real_, censored = FALSE, status = "no_excursion", cens_time = censor_time))
  }
  first_true <- function(condition) {
    i <- which(condition)[1]
    if (is.na(i)) censored else list(time = time[i] - setup$t_end, censored = FALSE, status = "recovered",
                                     cens_time = censor_time)
  }
  tolerance <- settings$tol_sd * setup$sd0

  switch(definition,
    sustained = {
      hold <- max(1L, round(settings$hold_s))
      inside <- as.numeric(dz <= tolerance)
      run_length <- as.numeric(stats::filter(inside, rep(1, hold), sides = 1))
      held <- !is.na(run_length) & run_length >= hold
      first_true(c(held[-seq_len(hold - 1L)], rep(FALSE, hold - 1L)))     # index of the first sample of the held window
    },
    first = first_true(dz <= tolerance),
    frac90 = first_true(dz <= 0.1 * setup$A_end),
    expfit = {
      elapsed <- time - setup$t_end
      decay_amplitude <- function(tau) { decay <- exp(-elapsed / tau); sum(dz * decay) / sum(decay^2) }
      sse <- function(log_tau) {
        tau <- exp(log_tau)
        sum((dz - decay_amplitude(tau) * exp(-elapsed / tau))^2)
      }
      tau <- exp(optimize(sse, interval = log(c(2, 5000)))$minimum)
      if (decay_amplitude(tau) <= 0 || tau > 4000) censored
      else list(time = tau, censored = FALSE, status = "recovered", cens_time = censor_time)
    },
    stop("Unknown recovery definition: ", definition, call. = FALSE)
  )
}

#' Recovery metrics of one state series under every definition. `rest_window` = c(start, end) of the
#' rest period (s). NULL when the rest period is too short or the end of exercise is unknown.
recovery_metrics <- function(t, x, rest_window, t_end, settings, definitions = RECOVERY_DEFINITIONS) {
  grid <- to_grid(t, x, smooth_s = settings$smooth_s)
  rest <- grid[t >= rest_window[1] & t < rest_window[2]]
  if (nrow(rest) < 15 || !is.finite(t_end)) return(NULL)
  setup <- recovery_setup(grid, mean(rest$x), sd(rest$x), t_end, settings)
  if (is.null(setup)) return(NULL)

  rbindlist(lapply(definitions, function(definition) {
    result <- tryCatch(recovery_time(setup, definition, settings),
                       error = function(e) list(time = NA_real_, censored = NA, status = "error", cens_time = NA_real_))
    data.table(definition = definition, time = result$time, censored = result$censored, status = result$status,
               cens_time = result$cens_time, A_end = setup$A_end, rest_sd = setup$sd0)
  }))
}


# ---- Phase-plane geometry -------------------------------------------------------------------------------------------------
#' Geometry of the (x_p, x_s) path in units of the baseline SDs (energy coordinates): path length,
#' absolute shoelace area of the path closed on its start, largest distance, energy U = (u^2 + v^2)/2
#' and tortuosity. NULL for paths shorter than 5 samples.
path_geometry <- function(t, x_p, x_s, sd_p, sd_s, smooth_s = 5) {
  grid_p <- to_grid(t, x_p, smooth_s = smooth_s)
  grid_s <- to_grid(t, x_s, smooth_s = smooth_s)
  n <- min(nrow(grid_p), nrow(grid_s))
  if (n < 5) return(NULL)
  u <- grid_p$x[1:n] / sd_p
  v <- grid_s$x[1:n] / sd_s
  step_length <- sqrt(diff(u)^2 + diff(v)^2)
  area <- abs(0.5 * sum(u * c(v[-1], v[1]) - c(u[-1], u[1]) * v))
  energy <- 0.5 * (u^2 + v^2)
  data.table(path_length = sum(step_length), loop_area = area, max_dist = max(sqrt(u^2 + v^2)),
             U_max = max(energy), U_mean = mean(energy),
             tortuosity = sum(step_length) / max(sqrt((u - u[1])^2 + (v - v[1])^2)))
}
