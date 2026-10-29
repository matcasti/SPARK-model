# ==============================================================================
# R/signals/physiology.R -- analyses that use the optional covariates (tilt angle, blood pressure,
# respiration). All functions take plain vectors (time in seconds), so they work for any
# dataset; scripts skip them for a subject that lacks the signal.
#
#   tilt_events()  angle_levels()  plateau_means()    tilt angle: every tilt (onset / return), plateaus
#   xcorr_surrogate()                                  lagged correlation with surrogate p-value
#   resp_frequency()  resp_coupling()                  respiration and its coupling with a state
#   brs_alpha()  bp_summary()                          baroreflex alpha index, blood-pressure summary
#   resample_series()                                  beat-wise series -> regular grid
# ==============================================================================


# ---- Resampling and spectra ----------------------------------------------------------------------------------
#' Beat-wise series onto a regular grid (cubic spline). NULL with fewer than 20 finite points.
resample_series <- function(t, x, fs = 4) {
  finite <- is.finite(t) & is.finite(x)
  t <- t[finite]
  x <- x[finite]
  if (length(x) < 20) return(NULL)
  grid <- seq(ceiling(t[1] * fs) / fs, floor(t[length(t)] * fs) / fs, by = 1 / fs)
  data.table(t = grid, x = spline(t, x, xout = grid, method = "natural")$y)
}

#' Odd smoothing span (in periodogram bins) of about `spans_hz` for a series of n samples at fs Hz; NULL = none.
spectrum_spans <- function(n, fs, spans_hz = 0.01) {
  bins <- max(1L, round(spans_hz / (fs / n)))
  if (bins >= 2) 2 * (bins %/% 2) + 1 else NULL
}

#' Tapered, detrended, smoothed periodogram of one series (R's spec.pgram; two-sided density).
smoothed_spectrum <- function(x, fs, spans_hz = 0.01) {
  spec.pgram(ts(x, frequency = fs), spans = spectrum_spans(length(x), fs, spans_hz),
             taper = 0.1, detrend = TRUE, plot = FALSE)
}


# ---- Tilt angle ------------------------------------------------------------------------------------------------
#' Every tilt of a tilt-angle signal. Baseline = median angle in the first `baseline_s` seconds; an onset is
#' the first time the (median-filtered) angle exceeds baseline + threshold for `hold_s`; its return is the first
#' later time within the threshold of baseline for `hold_s`; the search for the next onset resumes at that
#' return. A final tilt without a return is kept with `t_return = NA`. Times are relative to the first sample.
#' @return list(detected, n_events, events = data.table(event, t_onset, t_return, peak_angle), base_angle,
#'         peak_angle, t_onset, t_return); `t_onset` / `t_return` repeat the first event (single-tilt use)
tilt_events <- function(time, angle, tilt_cfg) {
  t <- time - time[1]
  finite <- is.finite(angle)
  no_event <- list(detected = FALSE, n_events = 0L,
                   events = data.table(event = integer(), t_onset = numeric(), t_return = numeric(), peak_angle = numeric()),
                   base_angle = NA_real_, peak_angle = NA_real_, t_onset = NA_real_, t_return = NA_real_)
  if (sum(finite) < 30) return(no_event)

  filled <- angle
  filled[!finite] <- approx(t[finite], angle[finite], xout = t[!finite], rule = 2)$y
  smoothed <- as.numeric(runmed(filled, k = min(5L, 2L * ((length(filled) - 1L) %/% 2L) + 1L)))
  baseline <- median(smoothed[t < tilt_cfg$baseline_s])
  above <- (smoothed - baseline) > tilt_cfg$threshold_deg

  hold <- max(1L, round(tilt_cfg$hold_s / median(diff(t))))
  held_start <- function(flag) {                         # first index of every window of `hold` samples with `flag` set
    run_length <- as.numeric(stats::filter(as.numeric(flag), rep(1, hold), sides = 1))
    which(!is.na(run_length) & run_length >= hold) - hold + 1L
  }
  onset_starts <- held_start(above)
  return_starts <- held_start(!above)

  rows <- list()
  resume <- 1L
  repeat {
    onset_index <- onset_starts[onset_starts >= resume][1]
    if (is.na(onset_index)) break
    return_index <- return_starts[return_starts > onset_index][1]
    last_index <- if (is.na(return_index)) length(t) else return_index
    rows[[length(rows) + 1L]] <- data.table(
      event = length(rows) + 1L, t_onset = t[onset_index],
      t_return = if (is.na(return_index)) NA_real_ else t[return_index],
      peak_angle = max(smoothed[onset_index:last_index]))
    if (is.na(return_index)) break
    resume <- return_index
  }
  if (!length(rows)) return(modifyList(no_event, list(base_angle = baseline)))

  events <- rbindlist(rows)
  list(detected = TRUE, n_events = nrow(events), events = events, base_angle = baseline,
       peak_angle = max(events$peak_angle), t_onset = events$t_onset[1], t_return = events$t_return[1])
}

#' Angle plateaus: contiguous runs at the same level (angle rounded to `level_bin_deg`) lasting at
#' least `plateau_min_s`. @return data.table(level, t_start, t_end)
angle_levels <- function(time, angle, tilt_cfg) {
  t <- time - time[1]
  finite <- is.finite(angle)
  if (sum(finite) < 30) return(data.table(level = numeric(), t_start = numeric(), t_end = numeric()))
  smoothed <- as.numeric(runmed(approx(t[finite], angle[finite], xout = t, rule = 2)$y, k = 5))
  level <- round(smoothed / tilt_cfg$level_bin_deg) * tilt_cfg$level_bin_deg
  runs <- rle(level)
  last <- cumsum(runs$lengths)
  first <- last - runs$lengths + 1L
  plateaus <- data.table(level = runs$values, t_start = t[first], t_end = t[last])
  plateaus[t_end - t_start >= tilt_cfg$plateau_min_s][]
}

#' Mean of each variable in `variables` (a named list of beat-wise vectors) over the second half of
#' every angle plateau. NULL when there is no plateau.
plateau_means <- function(time, angle, variables, tilt_cfg) {
  plateaus <- angle_levels(time, angle, tilt_cfg)
  if (!nrow(plateaus)) return(NULL)
  t <- time - time[1]
  rbindlist(lapply(seq_len(nrow(plateaus)), function(i) {
    second_half <- t >= (plateaus$t_start[i] + plateaus$t_end[i]) / 2 & t <= plateaus$t_end[i]
    as.data.table(c(list(level = plateaus$level[i], t_start = plateaus$t_start[i]),
                    lapply(variables, function(v) mean(v[second_half], na.rm = TRUE))))
  }))
}


# ---- Lagged correlation with surrogates -----------------------------------------------------------------------------
#' Lagged correlation of two beat-wise series after removing slow trends (centered moving mean of
#' `hp_s` seconds, at 1 Hz). Positive lag: y lags x. Circular-shift surrogates preserve each series'
#' autocorrelation and give a null for the peak |r|.
#' @return list(peak_r, lag_s, p, r0, n); NULL for series shorter than 4 x max_lag_s
xcorr_surrogate <- function(t, x, y, max_lag_s = 30, hp_s = 60, n_surr = 200L, seed = 1L) {
  grid_x <- to_grid(t, x, smooth_s = 0)
  grid_y <- to_grid(t, y, smooth_s = 0)
  n <- min(nrow(grid_x), nrow(grid_y))
  if (n < 4 * max_lag_s) return(NULL)

  remove_slow_trend <- function(v) {
    window <- max(3L, 2L * round(hp_s / 2) + 1L)
    detrended <- v - as.numeric(stats::filter(v, rep(1 / window, window), sides = 2))
    detrended[is.na(detrended)] <- 0
    detrended
  }
  inner <- seq(round(hp_s / 2) + 1L, n - round(hp_s / 2))              # drop the edges affected by the filter
  a <- remove_slow_trend(grid_x$x[1:n])[inner]
  b <- remove_slow_trend(grid_y$x[1:n])[inner]

  cross_correlation <- function(u, v) {
    r <- ccf(u, v, lag.max = max_lag_s, plot = FALSE, na.action = na.pass)
    list(lag = r$lag[, 1, 1], r = r$acf[, 1, 1])
  }
  observed <- cross_correlation(a, b)
  peak <- which.max(abs(observed$r))
  surrogate_peaks <- with_seed(seed, vapply(seq_len(n_surr), function(i) {
    shift <- sample(seq(2 * max_lag_s, length(b) - 2 * max_lag_s), 1)
    max(abs(cross_correlation(a, b[c((shift + 1):length(b), 1:shift)])$r))
  }, numeric(1)))

  # ccf(u, v) at lag k pairs u[t + k] with v[t]; the reported lag is that of y relative to x
  list(peak_r = observed$r[peak], lag_s = -observed$lag[peak],
       p = (1 + sum(surrogate_peaks >= abs(observed$r[peak]))) / (n_surr + 1),
       r0 = observed$r[observed$lag == 0], n = length(a))
}


# ---- Respiration -------------------------------------------------------------------------------------------------------
#' Dominant respiratory frequency (Hz) of a beat-wise respiration signal and the share of the
#' in-band power within +/- 0.03 Hz of the peak.
resp_frequency <- function(time, resp, band = c(0.10, 0.50), fs = 4) {
  missing <- c(f_resp = NA_real_, peak_frac = NA_real_)
  series <- resample_series(time, resp, fs)
  if (is.null(series) || nrow(series) < 16 * fs) return(missing)
  spectrum <- smoothed_spectrum(series$x, fs)
  in_band <- spectrum$freq >= band[1] & spectrum$freq <= band[2]
  if (!any(in_band)) return(missing)
  peak <- which(in_band)[which.max(spectrum$spec[in_band])]
  c(f_resp = spectrum$freq[peak],
    peak_frac = sum(spectrum$spec[abs(spectrum$freq - spectrum$freq[peak]) <= 0.03]) / sum(spectrum$spec[in_band]))
}

#' Coupling of a beat-wise series with respiration at the respiratory frequency `f_resp`:
#'   coh          magnitude-squared coherence over f_resp +/- halfwidth, weighted by the respiratory
#'                power (amplitude-blind);
#'   band_sd      SD of the series within that band (units of x, e.g. Hz);
#'   share_total  fraction of the whole (detrended) variance within the band.
resp_coupling <- function(time, x, resp, f_resp, halfwidth = 0.04, fs = 4) {
  missing <- c(coh = NA_real_, band_sd = NA_real_, share_total = NA_real_)
  series <- resample_series(time, x, fs)
  respiration <- resample_series(time, resp, fs)
  if (is.null(series) || is.null(respiration) || !is.finite(f_resp)) return(missing)
  n <- min(nrow(series), nrow(respiration))
  if (n < 16 * fs) return(missing)

  spectra <- spec.pgram(ts(cbind(series$x[1:n], respiration$x[1:n]), frequency = fs),
                        spans = spectrum_spans(n, fs), taper = 0.1, detrend = TRUE, plot = FALSE)
  in_band <- abs(spectra$freq - f_resp) <= halfwidth
  power <- 2 * spectra$spec[, 1] * (spectra$freq[2] - spectra$freq[1])      # one-sided power per bin
  weight <- spectra$spec[in_band, 2]                                        # respiratory power in the band
  c(coh = sum(spectra$coh[in_band] * weight) / sum(weight),
    band_sd = sqrt(sum(power[in_band])), share_total = sum(power[in_band]) / sum(power))
}


# ---- Blood pressure -------------------------------------------------------------------------------------------------------
#' Baroreflex alpha index: sqrt(P_RR / P_SBP), averaged over the frequencies of each band at which
#' the squared coherence between SBP and RR is at least `min_coh` (ms / mmHg). RR in seconds.
#' @return named vector alpha_<band>, coh_<band> (NA where it cannot be computed)
brs_alpha <- function(time, sbp, rr, bands, min_coh = 0.5, fs = 4) {
  result <- setNames(rep(NA_real_, 2 * length(bands)), c(paste0("alpha_", names(bands)), paste0("coh_", names(bands))))
  systolic <- resample_series(time, sbp, fs)
  interval <- resample_series(time, 1000 * rr, fs)
  if (is.null(systolic) || is.null(interval)) return(result)
  n <- min(nrow(systolic), nrow(interval))
  if (n < 32 * fs) return(result)

  spectra <- spec.pgram(ts(cbind(systolic$x[1:n], interval$x[1:n]), frequency = fs),
                        spans = spectrum_spans(n, fs), taper = 0.1, detrend = TRUE, plot = FALSE)
  for (name in names(bands)) {
    in_band <- spectra$freq >= bands[[name]][1] & spectra$freq < bands[[name]][2]
    result[paste0("coh_", name)] <- mean(spectra$coh[in_band])
    coherent <- in_band & spectra$coh >= min_coh
    if (any(coherent)) result[paste0("alpha_", name)] <- mean(sqrt(spectra$spec[coherent, 2] / spectra$spec[coherent, 1]))
  }
  result
}

#' Blood-pressure summary (mmHg): mean systolic, diastolic, MAP, pulse pressure and SD of MAP.
bp_summary <- function(systolic, diastolic, map = NULL) {
  mean_arterial <- if (is.null(map)) (systolic + 2 * diastolic) / 3 else map
  data.table(sbp = mean(systolic, na.rm = TRUE), dbp = mean(diastolic, na.rm = TRUE),
             map = mean(mean_arterial, na.rm = TRUE), pulse_pressure = mean(systolic - diastolic, na.rm = TRUE),
             sd_map = sd(mean_arterial, na.rm = TRUE))
}
