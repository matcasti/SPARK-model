# ==============================================================================
# experiments/blockade/lib.R -- spectra of RR segments under blockade, empirical shares w(f), the production partition and alternative forms.
# Source from the project root:  source("experiments/blockade/lib.R")
# ==============================================================================
suppressMessages(library(data.table))
EDGES <- c(0.01, 0.02, 0.03, 0.04, 0.06, 0.08, 0.10, 0.125, 0.15, 0.20, 0.25, 0.30, 0.40)   # Hz; log-ish bins, 0.04-0.40 is the production partition range
BINS <- data.table(lo = head(EDGES, -1), hi = tail(EDGES, -1))[, `:=`(f = sqrt(lo * hi), bin = .I)]

#' Power of one segment per frequency bin: RR (s) at the beat times `tk` inside [t0, t1] -> heart rate (Hz) or ln RR, monotone spline on a 4 Hz grid, linear detrend,
#' smoothed periodogram (Daniell); variance in each bin = integral of the PSD. Returns NULL if the segment has too few beats or implausible RR.
segment_power <- function(tk, t0, t1, units = c("rate", "lnrr"), fs = 4, trim = 20) {
  units <- match.arg(units); k <- tk >= t0 + trim & tk <= t1
  rr <- c(NA, diff(tk))[k]; tt <- tk[k]; ok <- is.finite(rr) & rr > 0.25 & rr < 2.5
  if (sum(ok) < 100 || mean(ok) < 0.97) return(NULL)
  tt <- tt[ok]; rr <- rr[ok]; y <- if (units == "rate") 1 / rr else log(rr)
  grid <- seq(min(tt), max(tt), by = 1 / fs); x <- splinefun(tt, y, method = "monoH.FC")(grid)
  sp <- spec.pgram(ts(x, frequency = fs), spans = c(5, 5), taper = 0.1, detrend = TRUE, plot = FALSE)
  df <- diff(sp$freq)[1]; p <- sp$spec
  bp <- vapply(seq_len(nrow(BINS)), function(i) { s <- sp$freq >= BINS$lo[i] & sp$freq < BINS$hi[i]; if (!any(s)) NA_real_ else sum(p[s]) * df }, numeric(1))
  list(power = bp, n = length(tt), mean_rr = mean(rr), sd_rr = sd(rr), total = var(y))
}

#' Remove double detections: a peak that makes an RR below 45% of the median of its neighbours is dropped (the same beat counted twice).
clean_peaks <- function(tk, frac = 0.45) {
  repeat { rr <- diff(tk); if (length(rr) < 12) break
    loc <- stats::runmed(rr, 11, endrule = "median"); bad <- which(rr < frac * loc); if (!length(bad)) break
    tk <- tk[-(bad[1] + 1)] }
  tk
}

#' Segments x bins table: one row per segment of experiments/blockade/_out/protocol.csv, with the power of each bin (`p01`..`p12`), for heart rate (Hz) and ln RR.
build_segments <- function(peak_dir = "experiments/blockade/_out/rpeaks", protocol = "experiments/blockade/_out/protocol.csv") {
  pr <- fread(protocol); out <- list()
  for (i in seq_len(nrow(pr))) { r <- pr[i]; f <- sprintf("%s/rpeaks_Subject%d_day%d.csv", peak_dir, r$subject, r$day)
    tk <- clean_peaks(scan(f, quiet = TRUE)); row <- list()
    for (u in c("rate", "lnrr")) { sp <- segment_power(tk, r$t0, r$t1, units = u); if (is.null(sp)) next
      row[[u]] <- data.table(r[, .(subject, day, cond, breathing, t0)], units = u, n = sp$n, mean_rr = sp$mean_rr, sd_rr = sp$sd_rr, total = sp$total, bin = BINS$bin, f = BINS$f, power = sp$power) }
    out[[i]] <- rbindlist(row) }
  rbindlist(out)
}

# ---- Model spectrum helpers (used by ou_fit.R and extras/figure_ou.R) ----
lor <- function(f, V, tau) 4 * V * tau / (1 + (2 * pi * f * tau)^2)
peak <- function(f, A, fr, w) A * (w / pi) / ((f - fr)^2 + w^2)

#' Periodogram of ln RR (cosine taper, 4 Hz grid) between 0.011 and 0.5 Hz for the beats inside [t0 + trim, t1].
pgram <- function(tk, t0, t1, trim = 20, fs = 4) { k <- tk >= t0 + trim & tk <= t1; rr <- c(NA, diff(tk))[k]; tt <- tk[k]; ok <- is.finite(rr) & rr > 0.25 & rr < 2.5
  if (mean(ok) < 0.97 || sum(ok) < 100) return(NULL); tt <- tt[ok]; y <- log(rr[ok]); grid <- seq(min(tt), max(tt), by = 1 / fs); x <- splinefun(tt, y, method = "monoH.FC")(grid)
  sp <- spec.pgram(ts(x, frequency = fs), taper = 0.1, detrend = TRUE, plot = FALSE, fast = FALSE); s <- sp$freq >= 0.011 & sp$freq <= 0.5; list(f = sp$freq[s], I = sp$spec[s], n = length(x), rr = mean(rr[ok])) }
