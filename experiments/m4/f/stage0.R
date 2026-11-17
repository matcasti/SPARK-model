# ==============================================================================
# experiments/m4/f/stage0.R -- EXPERIMENTAL. Stage-0 check of M4f (README.md): does the respiratory-band amplitude of the heart rate (Hz) track the vagal level?
#   Rscript experiments/m4/f/stage0.R
# Development data only: the M4e exercise dev sample (seed 21) and the Lancaster blockade segments (calibration).
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")

#' Respiratory-band envelope of the instantaneous heart rate (Hz): 4 Hz resampling, 20 s detrend, FFT band-pass, Hilbert transform.
hf_env <- function(rr, band = c(0.15, 0.8), fs = 4) {
  tt <- cumsum(rr); g <- seq(ceiling(tt[1]), floor(tt[length(tt)]), by = 1 / fs); f <- approx(tt, 1 / rr, g)$y
  f <- f - as.numeric(stats::filter(f, rep(1 / (20 * fs), 20 * fs), sides = 2)); ok <- !is.na(f); g <- g[ok]; f <- f[ok]
  n <- length(f); X <- fft(f); fr <- (0:(n - 1)) * fs / n; fr <- pmin(fr, fs - fr); X[fr < band[1] | fr > band[2]] <- 0
  h <- rep(0, n); h[1] <- 1; if (n %% 2 == 0) { h[2:(n / 2)] <- 2; h[n / 2 + 1] <- 1 } else h[2:((n + 1) / 2)] <- 2
  list(t = g, env = Mod(fft(X * h, inverse = TRUE) / n), rate = approx(tt, 1 / rr, g)$y) }

# ---- Step-test exercise, M4e development sample ----
d <- load_dataset("exercise", cfg); set.seed(21); dev <- sort(sample(names(d$rr), 60))
run_ex <- function(band) rbindlist(lapply(dev, function(id) { rr <- d$rr[[id]]; ev <- tryCatch(events_from_rr(prepare_series(rr)$dy, "exercise", cfg$segmentation), error = function(e) NULL)
  if (is.null(ev) || !nrow(ev)) return(NULL); h <- hf_env(rr, band); t_on <- ev$t0[1]; t_end <- ev$t1[1]
  ph <- list(rest = h$t >= 20 & h$t < t_on - 5, early = h$t >= t_on + 10 & h$t < t_on + 40, late = h$t >= t_end - 30 & h$t <= t_end)
  if (any(sapply(ph, sum) < 40)) return(NULL)
  pk <- function(k) { x <- h$rate[k] - mean(h$rate[k]); n <- length(x); P <- Mod(fft(x))^2; fr <- (0:(n - 1)) * 4 / n; ok <- fr >= band[1] & fr <= band[2]; fr[ok][which.max(P[ok])] }
  data.table(id = id, phase = names(ph), amp = sapply(ph, function(k) mean(h$env[k])), rate = sapply(ph, function(k) mean(h$rate[k])), fpk = sapply(ph, pk)) }))
for (b in list(c(0.15, 0.5), c(0.15, 0.8))) { W <- dcast(run_ex(b), id ~ phase, value.var = c("amp", "rate", "fpk"))
  cat(sprintf("Exercise, band %.2f-%.2f Hz, n = %d, rest HR %.0f bpm\n", b[1], b[2], nrow(W), 60 * median(W$rate_rest)))
  for (p in c("early", "late")) { r <- W[[paste0("amp_", p)]] / W$amp_rest; dr <- W[[paste0("rate_", p)]] - W$rate_rest
    cat(sprintf("  %-5s amp ratio %.2f [IQR %.2f-%.2f], share < 0.5 %.2f, dHR %.0f bpm, Spearman(log ratio, dHR) %.2f, peak %.2f -> %.2f Hz\n", p, median(r), quantile(r, .25), quantile(r, .75),
                mean(r < 0.5), 60 * median(dr), cor(log(r), dr, method = "spearman"), median(W$fpk_rest), median(W[[paste0("fpk_", p)]]))) } }

# ---- Blockade: ratio to saline of the same day and breathing mode ----
protocol <- fread(file.path(cfg$paths$blockade, "protocol.csv")); st <- cfg$blockade
seg_rr <- function(g) { tk <- blockade_peaks(cfg, g$subject, g$day); k <- tk >= g$t0 + st$trim_s & tk <= g$t1; rr <- diff(tk[k]); if (length(rr) < st$min_beats) NULL else rr[rr > 0.25 & rr < 2.5] }
BL <- rbindlist(lapply(which(protocol$breathing != "paced_0.1"), function(j) { g <- protocol[j]; rr <- seg_rr(g); if (is.null(rr)) return(NULL); h <- hf_env(rr, c(0.15, 0.5))
  data.table(subject = g$subject, day = g$day, cond = g$cond, breathing = g$breathing, amp = mean(h$env), rate = mean(h$rate)) }))
R <- merge(BL[cond != "saline"], BL[cond == "saline", .(subject, day, breathing, amp0 = amp, rate0 = rate)], by = c("subject", "day", "breathing"))
S <- R[, .(n = .N, amp_ratio = round(exp(mean(log(amp / amp0))), 2), drate_Hz = round(mean(rate - rate0), 3)), by = .(breathing, cond)][order(breathing, cond)]
print(S, row.names = FALSE)
cat(sprintf("Implied gamma (atropine high, all of the shift vagal): %.1f per Hz\n", S[cond == "atropine_high", -mean(log(amp_ratio)) / mean(drate_Hz)]))
