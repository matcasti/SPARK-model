# ==============================================================================
# R/steps/blockade.R -- pharmacological blockade and autonomic maneuvers (supplementary "Pharmacological blockade and muscle sympathetic nerve activity").
# Called by analysis/datasets/15_blockade.R. Data: the public Lancaster "Blocker data" (7 healthy volunteers, saline, atropine, propranolol; doi
# 10.17635/lancaster/researchdata/280), prepared by tools/prep/blockade.R into data/blockade/ (segment table, R-peak times, maneuver times, 1-s nerve signal).
#
#   blockade_available()        are the prepared files there?
#   blockade_peaks()            R-peak times of one session, double detections removed
#   blockade_segment_fit()      M3 fitted to one 5-min segment: amplitudes, time constants, anchor, heart rate
#   blockade_ratios()           ratio of each condition to saline of the same day, per subject and breathing mode
#   blockade_ratio_summary()    geometric mean over subjects with a bootstrap 95% interval
#   blockade_maneuver_session() states, heart rate and nerve signal during hyperventilation and breath-hold (saline period)
# ==============================================================================
BLOCKADE_CONDITIONS <- c("atropine_high", "double", "propranolol")       # low-dose atropine is paradoxical (vagal stimulation) and is not analysed
BLOCKADE_BREATHING <- c("spontaneous", "paced_0.25")                   # the two modes recorded in every condition (0.1-Hz pacing exists in one session only)

blockade_available <- function(cfg) {
  d <- cfg$paths$blockade
  file.exists(file.path(d, "protocol.csv")) && file.exists(file.path(d, "maneuvers.csv")) && length(list.files(file.path(d, "rpeaks"), pattern = "\\.csv$")) > 0
}

#' A peak that makes an RR shorter than `frac` of the running median of its neighbours is the same beat detected twice and is dropped.
clean_double_detections <- function(tk, frac = 0.45) {
  repeat {
    rr <- diff(tk); if (length(rr) < 12) break
    bad <- which(rr < frac * stats::runmed(rr, 11, endrule = "median")); if (!length(bad)) break
    tk <- tk[-(bad[1] + 1)]
  }
  tk
}

blockade_peaks <- function(cfg, subject, day) {
  clean_double_detections(scan(file.path(cfg$paths$blockade, "rpeaks", sprintf("rpeaks_Subject%d_day%d.csv", subject, day)), quiet = TRUE))
}

#' M3 (production specification, own spectral anchor) fitted to the beats of one segment, skipping its first `trim_s` s. NULL when the segment has too few
#' plausible intervals or the fit fails. SD_p = sqrt(lambda_P sigma^2) and SD_s = sqrt(sigma^2) are the stationary SDs of the two drives (Hz).
blockade_segment_fit <- function(segment, tk, settings) {
  k <- tk >= segment$t0 + settings$trim_s & tk <= segment$t1; rr <- diff(tk[k])
  if (length(rr) < settings$min_beats || mean(rr > 0.25 & rr < 2.5) < 0.97) return(NULL)
  rr <- rr[rr > 0.25 & rr < 2.5]; series <- prepare_series(rr); fit <- fit_variant(series, model_spec("M3"))
  if (!isTRUE(fit$ok)) return(NULL)
  lp <- series$priors$energy_ratio
  data.table(subject = segment$subject, day = segment$day, cond = segment$cond, breathing = segment$breathing, n = length(rr), hr = 60 / mean(rr),
             sd_rr_ms = 1000 * sd(rr), tau_p = fit$tau_p, tau_s = fit$tau_s, lambda_P = lp, SD_p = sqrt(lp * fit$sig2), SD_s = sqrt(fit$sig2), sd_xp = sd(fit$X_p), sd_xs = sd(fit$X_s))
}

#' Ratio of each condition to saline of the same day (day 1 for atropine and double blockade, day 2 for propranolol alone), per subject and breathing mode.
blockade_ratios <- function(fits) {
  rbindlist(lapply(BLOCKADE_CONDITIONS, function(cn) {
    base_day <- if (cn == "propranolol") 2L else 1L
    a <- fits[cond == cn & breathing %in% BLOCKADE_BREATHING]
    b <- fits[cond == "saline" & day == base_day, .(subject, breathing, SD_p0 = SD_p, SD_s0 = SD_s, hr0 = hr, lambda0 = lambda_P, sd_rr0 = sd_rr_ms)]
    x <- merge(a, b, by = c("subject", "breathing"))
    x[, .(subject, breathing, cond, SD_p = SD_p / SD_p0, SD_s = SD_s / SD_s0, rr_var = (sd_rr_ms / sd_rr0)^2, lambda_ratio = lambda_P / lambda0, d_hr = hr - hr0)]
  }))
}

#' Geometric mean over subjects and percentile bootstrap 95% interval (resampling subjects).
geometric_mean_ci <- function(x, n_boot, seed) {
  gm <- function(v) exp(mean(log(v)))
  b <- with_seed(seed, replicate(n_boot, gm(sample(x, replace = TRUE))))
  c(gm = gm(x), lo = unname(quantile(b, 0.025)), hi = unname(quantile(b, 0.975)))
}

blockade_ratio_summary <- function(ratios, settings) {
  ratios[, {
    p <- geometric_mean_ci(SD_p, settings$n_boot, settings$seed); s <- geometric_mean_ci(SD_s, settings$n_boot, settings$seed + 1L)
    .(n = .N, SD_p = p[["gm"]], SD_p_lo = p[["lo"]], SD_p_hi = p[["hi"]], SD_s = s[["gm"]], SD_s_lo = s[["lo"]], SD_s_hi = s[["hi"]],
      rr_var = exp(mean(log(rr_var))), n_SD_p_below_0.3 = sum(SD_p < 0.3), n_SD_p_at_least_0.8 = sum(SD_p >= 0.8), n_SD_s_below_third = sum(SD_s < 1 / 3), d_hr = median(d_hr))
  }, by = .(breathing, cond)]
}

#' One session (subject x day), saline period only: M3 fitted to the RR from the spontaneous baseline to the first drug dose; filtered states, heart rate and the
#' 10-s smoothed nerve signal as z-scores against the baseline segment; their mean change during hyperventilation and during the second half of the breath-hold;
#' lagged correlation of each signal with the nerve signal (positive lag: nerve signal later).
blockade_maneuver_session <- function(m, protocol, cfg, lags = seq(-40, 40, by = 5)) {
  base <- protocol[subject == m$subject & day == m$day & cond == "saline" & breathing == "spontaneous"][1]
  if (!nrow(base) || anyNA(c(m$hv0, m$bh0, m$bh1, m$drug0))) return(NULL)
  tk <- blockade_peaks(cfg, m$subject, m$day); k <- tk >= base$t0 & tk < m$drug0; tb <- tk[k][-1]; rr <- diff(tk[k]); ok <- rr > 0.25 & rr < 2.5; tb <- tb[ok]; rr <- rr[ok]
  fit <- fit_variant(prepare_series(rr), model_spec("M3")); if (!isTRUE(fit$ok)) return(NULL)
  nerve <- fread(file.path(cfg$paths$blockade, "msna", sprintf("msna_Subject%d_day%d.csv", m$subject, m$day)))
  grid <- seq(ceiling(min(tb)), floor(max(tb)))
  msna <- as.numeric(stats::filter(approx(nerve$t, nerve$mean, xout = grid)$y, rep(1 / 10, 10), sides = 2))
  in_base <- grid >= base$t0 + 20 & grid <= base$t1
  z <- function(x) (x - mean(x[in_base], na.rm = TRUE)) / sd(x[in_base], na.rm = TRUE)
  Z <- list(MSNA = z(msna), x_s = z(approx(tb, fit$X_s, xout = grid)$y), x_p = z(approx(tb, fit$X_p, xout = grid)$y), HR = z(approx(tb, 60 / rr, xout = grid)$y))
  windows <- list(hyperventilation = c(m$hv0, m$bh0), breath_hold_late = c(m$bh0 + (m$bh1 - m$bh0) / 2, m$bh1))
  change <- rbindlist(lapply(names(windows), function(w) { a <- grid >= windows[[w]][1] & grid <= windows[[w]][2]
    data.table(subject = m$subject, day = m$day, window = w, signal = names(Z), change = vapply(Z, function(v) mean(v[a], na.rm = TRUE), numeric(1))) }))
  lagged <- rbindlist(lapply(c("x_s", "x_p", "HR"), function(sg) data.table(subject = m$subject, day = m$day, signal = sg, lag = lags,
    r = vapply(lags, function(L) { j <- seq_along(grid); jj <- j + L; v <- jj >= 1 & jj <= length(grid); suppressWarnings(cor(Z[[sg]][j[v]], Z$MSNA[jj[v]], use = "complete.obs")) }, numeric(1)))))
  list(change = change, lagged = lagged)
}
