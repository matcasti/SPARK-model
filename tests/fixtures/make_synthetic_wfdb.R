# Writes synthetic WFDB records (and TMST text files) in the layout tools/prep-data.R reads, so the
# real preparation script can be run end to end. Needs Python with `wfdb` (via reticulate).
suppressPackageStartupMessages(library(data.table))

pulse_wave <- function(t, onsets, sbp, dbp) {
  k <- findInterval(t, onsets); out <- rep(dbp[1], length(t))
  ok <- k >= 1 & k < length(onsets)
  kk <- k[ok]; u <- (t[ok] - onsets[kk]) / (onsets[kk + 1] - onsets[kk])
  e0 <- exp(-0.82 / 0.3)
  shape <- ifelse(u < 0.18, sin(pi / 2 * u / 0.18), (exp(-(u - 0.18) / 0.3) - e0) / (1 - e0))
  out[ok] <- dbp[kk] + (sbp[kk] - dbp[kk]) * shape
  out
}

#' @param kind "fantasia" or "prcp"; sim from the simulators; returns the truth per beat.
write_wfdb_record <- function(wfdb, np, dir, id, sim, kind, fs = 250L, has_bp = TRUE, has_extra = TRUE, f_resp = 0.25, delay = 0.18, seed = 1, spikes = TRUE) {
  set.seed(seed); t1 <- runif(1, 1, 5)
  beats <- t1 + c(0, sim$t - (sim$t[1] - sim$rr[1]))                 # R-peak times (s); K = length(beats) - 1 intervals
  beats <- round(beats * fs) / fs; K <- length(beats) - 1L; rr <- diff(beats)
  off <- t1 - (sim$t[1] - sim$rr[1])                                 # file clock = sim clock + off
  s <- if (kind == "prcp") sin(((attr(sim, "angle")(beats[-1] - off)) * pi / 180)) else rep(0, K)
  sbp <- 118 + 35 * (rr - mean(rr)) - 6 * s + rnorm(K, 0, 1.5); dbp <- sbp - (45 - 12 * s) + 3 * s
  n <- ceiling((beats[K + 1] + 3) * fs); tt <- (0:(n - 1)) / fs
  onsets <- c(beats + delay, beats[K + 1] + delay + median(rr))
  sbp_all <- c(sbp, sbp[K]); dbp_all <- c(dbp, dbp[K])
  spike_k <- integer()
  if (spikes && K > 60) { spike_k <- c(30L, 45L); sbp_all[spike_k] <- sbp_all[spike_k] + 90 }   # artifacts to be cleaned
  bp <- pulse_wave(tt, onsets, sbp_all, dbp_all)
  if (spikes && K > 100) bp[(80 * fs):(81 * fs)] <- NaN                                           # a 1-s dropout
  chans <- list(ECG = rnorm(n, 0, 0.05)); units <- list("mV")
  extra <- NULL
  if (kind == "fantasia") {
    if (has_extra) { chans$RESP <- sin(2 * pi * f_resp * tt) + rnorm(n, 0, 0.1); units <- c(units, "au") }
    if (has_bp) { chans$BP <- bp; units <- c(units, "mmHg") }
  } else {
    if (has_bp) { chans$ABP <- bp; units <- c(units, "mmHg") }
    if (has_extra) { chans$Angle <- attr(sim, "angle")(tt - off) + rnorm(n, 0, 0.3); units <- c(units, "deg") }
  }
  P <- do.call(cbind, chans)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  wfdb$wrsamp(record_name = id, fs = fs, units = as.list(unlist(units)), sig_name = as.list(names(chans)), p_signal = P,
              fmt = as.list(rep("16", ncol(P))), write_dir = dir)
  ext <- if (kind == "fantasia") "ecg" else "wqrs"
  smp <- as.integer(round(beats * fs)); sym <- rep("N", length(smp))
  extra_s <- as.integer(c(round(0.5 * fs), round(beats[40] * fs) + 60)); extra_sym <- c("+", "~")      # non-beat annotations
  allsmp <- c(smp, extra_s); allsym <- c(sym, extra_sym); o <- order(allsmp)
  py$write_ann(id, ext, as.list(allsmp), as.list(allsym), dir)
  if (kind == "prcp") {
    ab <- as.integer(round((beats + delay) * fs))
    py$write_ann(id, "wabp", as.list(ab), as.list(rep("N", length(ab))), dir)
  }
  data.table(beat = seq_len(K), time = beats[-1], rri = rr, sbp = sbp, dbp = dbp, spike = seq_len(K) %in% spike_k,
             resp = if (kind == "fantasia") sin(2 * pi * f_resp * (beats[-length(beats)] + rr / 2)) else NA_real_,
             angle = if (kind == "prcp") attr(sim, "angle")((beats[-length(beats)] + rr / 2) - off) else NA_real_)
}

write_synthetic_wfdb <- function(root, seed = 5, n_ex = 2) {
  library(reticulate); wfdb <- import("wfdb"); np <- import("numpy")
  # wrann is called from Python (array types are not preserved through reticulate's conversion)
  py_run_string("
import numpy as np, wfdb
def write_ann(rec, ext, samples, symbols, wdir):
    o = np.argsort(np.array(samples, dtype='int64'), kind='stable')
    wfdb.wrann(rec, ext, sample=np.array(samples, dtype='int64')[o], symbol=[symbols[i] for i in o], write_dir=wdir)
")
  truth <- list()
  raw <- file.path(root, "data-raw")
  # Fantasia-like: f1y01 (all signals), f1y02 (no BP), f1o01 (no respiration)
  specs <- list(list("f1y01", "young", TRUE, TRUE), list("f1y02", "young", FALSE, TRUE), list("f1o01", "older", TRUE, FALSE))
  for (i in seq_along(specs)) { sp <- specs[[i]]
    set.seed(seed + i); sim <- sim_rest(seed + i, sp[[2]], 420, m_p = rsa_drive(0.25, 0.04))
    truth[[sp[[1]]]] <- write_wfdb_record(wfdb, np, file.path(raw, "fantasia"), sp[[1]], sim, "fantasia", has_bp = sp[[3]], has_extra = sp[[4]], seed = seed + i) }
  # PRCP-like: prcp01 (BP + angle, step), prcp02 (BP + angle, stairs), prcp03 (no angle channel), prcp04 (no ABP channel)
  pspecs <- list(list("prcp01", "step", TRUE, TRUE), list("prcp02", "stairs", TRUE, TRUE), list("prcp03", "step", TRUE, FALSE), list("prcp04", "step", FALSE, TRUE))
  for (i in seq_along(pspecs)) { sp <- pspecs[[i]]
    sim <- sim_tilt_angle(seed * 10 + i, 420, sp[[2]])
    truth[[sp[[1]]]] <- write_wfdb_record(wfdb, np, file.path(raw, "prcp"), sp[[1]], sim, "prcp", has_bp = sp[[3]], has_extra = sp[[4]], seed = seed + i) }
  # TMST-like: intervals in ms, 36 min long (the script cuts at 30 min)
  dir.create(file.path(raw, "tmst"), recursive = TRUE, showWarnings = FALSE)
  for (i in seq_len(n_ex)) { sim <- sim_exercise(seed * 20 + i, 2200, h = 0.02)
    rr <- sim$rr; rr[c(50, 120)] <- c(rr[50] * 2, rr[120] * 0.4)                       # two artifacts
    writeLines(sprintf("%.1f", 1000 * rr), file.path(raw, "tmst", sprintf("ex%02d.txt", i)))
    truth[[sprintf("ex%02d", i)]] <- data.table(beat = seq_along(rr), rri = sim$rr) }
  file.copy(file.path(raw, "tmst", "ex01.txt"), file.path(raw, "tmst", "ex99.txt"))                # a repeated file: counts once
  invisible(truth)
}
