# ==============================================================================
# R/signals/hrv.R -- conventional HRV comparators and the CWT band power used for temporal resolution.
#   hrv_metrics()        SDNN, RMSSD, mean HR, LF, HF and LF/HF of an RR series (s)
#   band_power()         spectral band power (ms^2) of the RR tachogram
#   cwt_band_series()    time-resolved LF / HF power of heart rate from a Morlet CWT
#   first_departure()    first time a series leaves its baseline band
# ==============================================================================

#' Evenly sampled RR tachogram (ms) at `fs` Hz by natural cubic-spline interpolation.
rr_tachogram <- function(rr, fs = 4) {
  beat_time <- cumsum(rr)
  grid <- seq(beat_time[1], beat_time[length(beat_time)], by = 1 / fs)
  list(t = grid, x = 1000 * spline(beat_time, rr, xout = grid, method = "natural")$y)
}

#' Band power (ms^2): integral of the tapered, detrended periodogram of the tachogram over `band`.
band_power <- function(rr, band, fs = 4) {
  tachogram <- rr_tachogram(rr, fs)
  if (length(tachogram$x) < 8 * fs) return(NA_real_)
  spectrum <- spec.pgram(ts(tachogram$x, frequency = fs), taper = 0.1, detrend = TRUE, plot = FALSE)
  in_band <- spectrum$freq >= band[1] & spectrum$freq < band[2]
  if (!any(in_band)) return(NA_real_)
  2 * sum(spectrum$spec[in_band]) * (spectrum$freq[2] - spectrum$freq[1])   # one-sided: R reports a two-sided density
}

#' SDNN and RMSSD (ms), mean HR (bpm), LF and HF power (ms^2) and LF/HF. The spectral entries are NA
#' for windows shorter than hrv_cfg$min_dur_lfhf_s.
hrv_metrics <- function(rr, hrv_cfg) {
  rr <- rr[is.finite(rr)]
  duration <- sum(rr)
  lf <- hf <- NA_real_
  if (duration >= hrv_cfg$min_dur_lfhf_s) {
    lf <- band_power(rr, hrv_cfg$lf, hrv_cfg$fs)
    hf <- band_power(rr, hrv_cfg$hf, hrv_cfg$fs)
  }
  data.table(duration_s = duration, mean_hr = 60 / mean(rr), SDNN = 1000 * sd(rr),
             RMSSD = 1000 * sqrt(mean(diff(rr)^2)), LF = lf, HF = hf, LF_HF = lf / hf)
}

#' Time-resolved LF and HF power of heart rate (Hz^2) from a Morlet CWT (the transform family used by
#' the engine's spectral anchor). Values inside the cone of influence are NA.
#' @param tt beat times (s); default the cumulative RR
#' @return data.table(t, LF, HF) on the CWT grid (hrv_cfg$fs Hz)
cwt_band_series <- function(rr, hrv_cfg, tt = cumsum(rr)) {
  fs <- hrv_cfg$fs
  grid <- seq(min(tt), max(tt), by = 1 / fs)
  heart_rate <- splinefun(tt, 1 / rr, method = "monoH.FC")(grid)
  invisible(capture.output(wavelet <- biwavelet::wt(cbind(grid, heart_rate - mean(heart_rate)), do.sig = FALSE)))

  power <- wavelet$power
  frequency <- 1 / wavelet$period
  for (j in seq_len(ncol(power))) power[wavelet$period > wavelet$coi[j], j] <- NA     # outside the cone of influence
  band_mean <- function(band) {
    rows <- which(frequency >= band[1] & frequency < band[2])
    if (!length(rows)) return(rep(NA_real_, ncol(power)))
    colMeans(power[rows, , drop = FALSE], na.rm = TRUE)
  }
  series <- data.table(t = grid, LF = band_mean(hrv_cfg$lf), HF = band_mean(hrv_cfg$hf))
  series[, `:=`(LF = ifelse(is.nan(LF), NA, LF), HF = ifelse(is.nan(HF), NA, HF))][]
}

#' First time a series leaves its baseline band (mean +/- k SD over `baseline_window`) and stays
#' outside for `hold_s` seconds, looking no further than `t_max`; NA if never (or if the baseline has
#' fewer than 10 samples).
first_departure <- function(t, x, baseline_window, k = 2, hold_s = 5, t_max = Inf) {
  finite <- is.finite(x)
  baseline <- x[t >= baseline_window[1] & t < baseline_window[2] & finite]
  if (length(baseline) < 10) return(NA_real_)
  lower <- mean(baseline) - k * sd(baseline)
  upper <- mean(baseline) + k * sd(baseline)
  outside <- finite & (x < lower | x > upper) & t >= baseline_window[2] - 60 & t <= t_max
  hold <- max(1L, round(hold_s / median(diff(t))))
  run_length <- as.numeric(stats::filter(as.numeric(outside), rep(1, hold), sides = 1))
  first <- which(!is.na(run_length) & run_length >= hold)[1]
  if (is.na(first)) NA_real_ else t[first - hold + 1L]
}
