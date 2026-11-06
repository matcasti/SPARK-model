# EXPERIMENTAL. Real recordings: in-sample fit (log-likelihood, AIC-type penalty, KS, PACF, tau_p) and the HELD-OUT predictive score of each variant:
# fit on the first 70% of the beats (anchor from that part only), then filter the whole series with the fitted parameters fixed (nu0, sigma^2 included)
# and average the one-step predictive log-density over the last 30% (beats with the artifact multiplier below the likelihood gate).
# Variants: base (M3), osc (damped oscillator), het (noise scale rate^-gamma), osc+het. 'robust' is only compared in-sample (its weights use the observation itself).
# Run from the project root:  Rscript experiments/m3plus/real_yardstick.R [quick]
quick <- length(commandArgs(TRUE)) > 0
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
VARIANTS <- list(base = list(), osc = list(osc = TRUE), het = list(het = TRUE), osc_het = list(osc = TRUE, het = TRUE), robust = list(robust = TRUE))
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr)
  if (quick) { set.seed(5); ids <- sort(sample(ids, min(length(ids), c(prcp = 3, fantasia = 6, exercise = 8)[[ds]]))) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
diag_of <- function(f) { r <- f$run(full = TRUE); keep <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < f$ctl$diag_gate_max); z <- r$innovations_diagnostic[keep]; n <- length(z)
  c(LL = r$LL_star, ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf_viol = sum(abs(as.numeric(pacf(z, lag.max = min(f$ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(0.975) / sqrt(n))) }
tail_score <- function(fh, dy_full, n_head, mean_rate_head) {
  ctl <- fh$ctl; p <- fh$p; n <- length(dy_full); gm <- gate_mult(dy_full, ctl); qs <- NULL
  if (isTRUE(fh$spec$het)) { r <- rate_local(dy_full) * mean(1 / dy_full) / mean_rate_head; qs <- pmin(pmax(r^(-p$gamma), 0.1), 10) }
  h <- fh$run(); res <- ext_run(dy_full, ctl, p$kp, p$ks, fh$lp, p$lR, osc = p$osc, qs = qs, gm = gm, full = TRUE, beta_fixed = h$nu0, sig2_fixed = h$sig2)
  if (!is.finite(res$LL_star)) return(NA_real_); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]; mean(res$logdens[idx])
}
one <- function(job) {
  rows <- list(); dy <- job$rr; n <- length(dy); n_head <- floor(0.7 * n)
  series <- prepare_series(dy); series_h <- prepare_series(dy[seq_len(n_head)]); mr_head <- mean(1 / series_h$dy)
  for (v in names(VARIANTS)) {
    f <- tryCatch(fit_ext(series, VARIANTS[[v]]), error = function(e) list(ok = FALSE)); fh <- if (v == "robust") list(ok = FALSE) else tryCatch(fit_ext(series_h, VARIANTS[[v]]), error = function(e) list(ok = FALSE))
    row <- data.table(ds = job$ds, id = job$id, variant = v, n_par = if (isTRUE(f$ok)) f$n_par else NA_integer_)
    if (isTRUE(f$ok)) { dg <- diag_of(f); row[, `:=`(LL = dg[["LL"]], cost = f$cost, ks_ok = dg[["ks_ok"]], pacf_viol = dg[["pacf_viol"]], tau_p = f$tau_p, tau_s = f$tau_s,
        osc_hz = if (!is.null(f$p$osc)) f$omega / 2 / pi else NA_real_, osc_alpha = if (!is.null(f$p$osc)) f$p$osc$alpha else NA_real_, osc_lo = if (!is.null(f$p$osc)) f$p$osc$lo else NA_real_, gamma = f$p$gamma)] }
    row[, tail_logdens := if (isTRUE(fh$ok)) tail_score(fh, series$dy, n_head, mr_head) else NA_real_]
    rows[[v]] <- row }
  rbindlist(rows, fill = TRUE)
}
log_msg("Held-out yardstick | %d recordings | %d variants", length(jobs), length(VARIANTS))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording"), fill = TRUE)
dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(res, sprintf("experiments/m3plus/_out/real_%s.rds", if (quick) "quick" else "full"))
b <- res[variant == "base", .(ds, id, LL0 = LL, T0 = tail_logdens, pv0 = pacf_viol)]; r <- merge(res, b, by = c("ds", "id")); r[, `:=`(dLL = LL - LL0, dT = tail_logdens - T0, aic_gain = 2 * (LL - LL0) - 2 * (n_par - 3))]
cat("\n== In-sample (all recordings; LL not comparable for 'robust') ==\n")
print(r[, .(n = .N, extra_par = median(n_par - 3, na.rm = TRUE), dLL_med = round(median(dLL, na.rm = TRUE), 1), pct_AIC_better = round(100 * mean(aic_gain > 0, na.rm = TRUE), 1), ks_pass = round(100 * mean(ks_ok, na.rm = TRUE), 1), pacf_viol = round(mean(pacf_viol, na.rm = TRUE), 1), tau_p_med = round(median(tau_p, na.rm = TRUE), 1), pct_tau_p_in_prior = round(100 * mean(tau_p >= 3 & tau_p <= 6, na.rm = TRUE), 1)), by = variant][order(variant)], row.names = FALSE)
cat("\n== Held-out predictive score of the last 30% (nats per beat; difference to base, paired over recordings) ==\n")
print(r[variant != "base" & variant != "robust", .(n = sum(is.finite(dT)), mean_gain = round(mean(dT, na.rm = TRUE), 3), median_gain = round(median(dT, na.rm = TRUE), 3), pct_better = round(100 * mean(dT > 0, na.rm = TRUE), 1), p_wilcoxon = signif(tryCatch(wilcox.test(dT, mu = 0)$p.value, error = function(e) NA), 2)), by = variant], row.names = FALSE)
cat("\n== Held-out gain by dataset ==\n"); print(r[variant != "base" & variant != "robust", .(mean_gain = round(mean(dT, na.rm = TRUE), 3), pct_better = round(100 * mean(dT > 0, na.rm = TRUE), 1)), by = .(variant, ds)][order(variant, ds)], row.names = FALSE)
cat("\n== Oscillator: fitted frequency (Hz) and damping alpha ==\n"); print(r[variant == "osc", .(f_hz_med = round(median(osc_hz, na.rm = TRUE), 2), f_hz_q = paste(round(quantile(osc_hz, c(.1, .9), na.rm = TRUE), 2), collapse = "-"), alpha_med = round(median(osc_alpha, na.rm = TRUE), 2), lo_med = round(median(osc_lo, na.rm = TRUE), 2)), by = ds], row.names = FALSE)
cat("\nHeteroskedastic gamma (median, IQR):\n"); print(r[variant == "het", .(gamma_med = round(median(gamma, na.rm = TRUE), 2), q25 = round(quantile(gamma, .25, na.rm = TRUE), 2), q75 = round(quantile(gamma, .75, na.rm = TRUE), 2)), by = ds], row.names = FALSE)
