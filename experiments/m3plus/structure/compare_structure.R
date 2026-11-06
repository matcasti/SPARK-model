# EXPERIMENTAL. Which minimal addition explains the unmodeled structure (KS, PACF) while keeping the reading "fast vagal step, slow sympathetic ramp"?
# Variants: M3; M3O (oscillator); M3S_s / M3S (stochastic set-point of the sympathetic / of both drives); M3OS (oscillator + both set-points); M3A (adaptive noise scale);
# M3OSA (all three). Measures: KS pass, PACF violations, log-likelihood, held-out forward score (non-adaptive variants), time constants and, for exercise, whether the
# vagal state still reacts fast and the sympathetic state slowly (time to 50% of the rest-to-exercise change) and with the expected signs.
# Run from the project root:  Rscript experiments/m3plus/structure/compare_structure.R [quick]
arg <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""; quick <- arg == "quick"; sub100 <- arg == "sub"
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
vset <- if (length(commandArgs(TRUE)) >= 2) commandArgs(TRUE)[2] else "all"
VAR <- if (vset == "adapt") list(M3 = list(), M3A10 = list(adapt = 10), M3A25 = list(adapt = 25), M3A60 = list(adapt = 60), M3OA = list(osc = TRUE, adapt = 25)) else list(M3 = list(), M3O = list(osc = TRUE), M3S_s = list(sp = "s"), M3S = list(sp = c("p", "s")), M3OS = list(osc = TRUE, sp = c("p", "s")), M3A = list(adapt = 25), M3OSA = list(osc = TRUE, sp = c("p", "s"), adapt = 25))
jobs <- list(); for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (quick) { set.seed(5); ids <- sort(sample(ids, min(length(ids), c(prcp = 3, fantasia = 5, exercise = 8)[[ds]]))) }; if (sub100 && ds == "exercise") { set.seed(5); ids <- sort(sample(ids, 100)) }; for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
diag_of <- function(f) { r <- f$run(full = TRUE); k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < f$ctl$diag_gate_max); z <- r$innovations_diagnostic[k]; n <- length(z)
  list(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(f$ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n)), LL = r$LL_star, r = r) }
tail_score <- function(fh, dy_full, n_head) { ctl <- fh$ctl; p <- fh$p; gm <- gate_mult(dy_full, ctl); h <- fh$run()
  res <- ext2_run(dy_full, ctl, p$kp, p$ks, fh$lp, p$lR, osc = p$osc, q = p$q, gm = gm, full = TRUE, beta_fixed = h$nu0, sig2_fixed = h$sig2); if (!is.finite(res$LL_star)) return(NA_real_)
  idx <- (n_head + 1):length(dy_full); idx <- idx[gm[idx] < ctl$lik_gate_max]; mean(res$logdens[idx]) }
t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2; k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
one <- function(job) {
  dy <- job$rr; series <- prepare_series(dy); dy <- series$dy; n <- length(dy); n_head <- floor(0.7 * n); series_h <- prepare_series(dy[seq_len(n_head)]); rows <- list()
  ev <- if (job$ds == "exercise") events_from_rr(dy, "exercise", cfg$segmentation) else NULL; tt <- cumsum(dy)
  for (v in names(VAR)) { sp <- VAR[[v]]; f <- tryCatch(fit_ext2(series, sp), error = function(e) list(ok = FALSE)); if (!isTRUE(f$ok)) next
    dg <- diag_of(f); fh <- if (is.null(sp$adapt)) tryCatch(fit_ext2(series_h, sp), error = function(e) list(ok = FALSE)) else list(ok = FALSE)
    row <- data.table(ds = job$ds, id = job$id, variant = v, n_par = f$n_par, LL = dg$LL, ks_ok = dg$ks, pacf = dg$pacf, tau_p = f$tau_p, tau_s = f$tau_s, q_p = f$p$q[1], q_s = f$p$q[2], tail = if (isTRUE(fh$ok)) tail_score(fh, dy, n_head) else NA_real_)
    if (!is.null(ev)) { r <- dg$r; t_on <- ev$t0[1]; t_end <- ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) { xp <- r$X_p + if (!is.null(f$p$osc)) r$X_o else 0; xs <- r$X_s
        row[, `:=`(t50_vagal = t50(xp, tt, t_on, mean(xp[rest]), mean(xp[endw])), t50_symp = t50(xs, tt, t_on, mean(xs[rest]), mean(xs[endw])), d_vagal = mean(xp[endw]) - mean(xp[rest]), d_symp = mean(xs[endw]) - mean(xs[rest]))] } }
    rows[[v]] <- row }
  rbindlist(rows, fill = TRUE) }
log_msg("Structure comparison | %d recordings | %d variants", length(jobs), length(VAR)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording"), fill = TRUE)
dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(res, sprintf("experiments/m3plus/_out/structure_%s%s.rds", if (quick) "quick" else if (sub100) "sub" else "full", if (vset == "all") "" else paste0("_", vset)))
b <- res[variant == "M3", .(ds, id, LL0 = LL, T0 = tail, par0 = n_par)]; r <- merge(res, b, by = c("ds", "id")); r[, `:=`(dLL = LL - LL0, dT = tail - T0, aic = 2 * (LL - LL0) - 2 * (n_par - par0))]
cat("\n== All recordings ==\n"); print(r[, .(n = .N, extra_par = median(n_par - par0), dLL = round(median(dLL), 1), pct_AIC_better = round(100 * mean(aic > 0), 1), ks_pass = round(100 * mean(ks_ok), 1), pacf = round(mean(pacf), 1), heldout_gain = round(mean(dT, na.rm = TRUE), 3), pct_heldout_better = round(100 * mean(dT > 0, na.rm = TRUE), 1), tau_p = round(median(tau_p), 1), pct_tau_p_in_prior = round(100 * mean(tau_p >= 3 & tau_p <= 6), 1)), by = variant][order(match(variant, names(VAR)))], row.names = FALSE)
cat("\n== KS pass (%) by dataset ==\n"); print(dcast(r[, .(ks = round(100 * mean(ks_ok), 1)), by = .(ds, variant)], variant ~ ds, value.var = "ks")[order(match(variant, names(VAR)))], row.names = FALSE)
cat("\n== Exercise: does the vagal state react fast and the sympathetic state slowly? (time to 50% of the rest-to-end-of-exercise change, s; sign check) ==\n")
print(r[ds == "exercise" & is.finite(t50_vagal), .(n = .N, t50_vagal_med = round(median(t50_vagal), 1), t50_symp_med = round(median(t50_symp, na.rm = TRUE), 1), pct_vagal_faster = round(100 * mean(t50_vagal < t50_symp, na.rm = TRUE), 1), pct_expected_signs = round(100 * mean(d_vagal < 0 & d_symp > 0), 1)), by = variant][order(match(variant, names(VAR)))], row.names = FALSE)
