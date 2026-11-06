# EXPERIMENTAL. Oscillator frequency range: 0.15-0.40 Hz (used so far) vs 0.15-0.70 Hz (breathing during exercise can exceed 0.4 Hz), on the exercise recordings.
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); d <- load_dataset("exercise", cfg); jobs <- lapply(names(d$rr), function(id) list(id = id, rr = d$rr[[id]]))
diag_of <- function(f) { r <- f$run(full = TRUE); k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < f$ctl$diag_gate_max); z <- r$innovations_diagnostic[k]; n <- length(z)
  c(ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = 40, plot = FALSE)$acf)) > qnorm(.975) / sqrt(n)), LL = r$LL_star) }
one <- function(j) { series <- prepare_series(j$rr); rows <- list()
  for (g in list(narrow = c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4), wide = c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7))) { assign("OSC_GRID_HZ", g, envir = globalenv()); f <- fit_ext(series, list(osc = TRUE)); dg <- diag_of(f)
    rows[[length(rows) + 1]] <- data.table(id = j$id, grid = names(which(sapply(list(narrow = 6, wide = 9), function(n) n == length(g)))), hz = f$omega / 2 / pi, alpha = f$p$osc$alpha, ks_ok = dg[["ks_ok"]], pacf = dg[["pacf"]], LL = dg[["LL"]], tau_p = f$tau_p) }
  rbindlist(rows) }
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording")); saveRDS(res, "experiments/m3plus/_out/osc_range.rds")
print(res[, .(n = .N, ks_pass = round(100 * mean(ks_ok), 1), pacf = round(mean(pacf), 1), LL_med = round(median(LL), 1), f_hz_med = round(median(hz), 2), pct_above_0.4 = round(100 * mean(hz > 0.4), 1), tau_p_med = round(median(tau_p), 1)), by = grid], row.names = FALSE)
w <- dcast(res, id ~ grid, value.var = "LL"); cat("\nLog-likelihood gain of the wide grid over the narrow one: median", round(median(w$wide - w$narrow), 1), "; recordings with gain > 2:", round(100 * mean(w$wide - w$narrow > 2), 1), "%\n")
