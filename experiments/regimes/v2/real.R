# EXPERIMENTAL (round 2). Real recordings (exercise 60, Fantasia 20, PRCP 10): M3, M3O (oscillator), HMM (regularised 2 regimes) and HMMO (HMM with the oscillator of the M3O fit).
# Models: M3, M3O, HMM0 (round-1 agnostic prior), HMMr (reciprocity prior), HMMOr (reciprocity + oscillator). Per model: explicit log-likelihood vs M3, KS, PACF, tau_p, exercise realism (t50 vagal/sympathetic, signs), PRCP hits within 30 s of the annotated onsets, and the held-out forward score
# (fit first 70%, score the last 30%; regimes inferred online) against M3 and M3O. Run from the project root:  Rscript experiments/regimes/v2/real.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); ev_ann <- fread("data/prcp_events.csv"); REG <- list(recip = 20); GATE <- 9   # round-1 HMM = agnostic prior (HMM0); HMMr = reciprocity prior; HMMOr = reciprocity prior + oscillator
jobs <- list(); for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (ds == "exercise") { set.seed(21); ids <- sort(sample(ids, 60)) }; if (ds == "fantasia") { set.seed(22); ids <- sort(sample(ids, 20)) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]], off = if (ds == "prcp") d$rec[[id]]$time[1] - d$rec[[id]]$rr[1] else 0) }
z_diag <- function(z, gm, ctl) { k <- which(is.finite(z) & gm < ctl$diag_gate_max); z <- z[k]; n <- length(z); c(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n))) }
t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2; k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
one <- function(job) {
  series <- prepare_series(job$rr); dy <- series$dy; n <- length(dy); tt <- cumsum(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); n_eff <- sum(gm < ctl$lik_gate_max); expl <- function(l) l - n_eff / 2 * (log(2 * pi) + 1)
  base <- fit_ext(series); fo <- fit_ext(series, list(osc = TRUE)); cp <- cp_search(series, base, min_gain = -Inf); ll0 <- expl(base$run()$LL_star); rows <- list()
  ex_ev <- if (job$ds == "exercise") events_from_rr(dy, "exercise", cfg$segmentation) else NULL
  stats <- function(m, xp, xs, z, ll, tau_p, tau_s, onsets) { dg <- z_diag(z, gm, ctl); row <- data.table(ds = job$ds, id = job$id, model = m, ks_ok = dg[["ks"]], pacf = dg[["pacf"]], dLL = ll - ll0, tau_p = tau_p, tau_s = tau_s, gain = if (is.null(cp)) NA_real_ else cp$gain, n_onsets = length(onsets))
    if (!is.null(ex_ev)) { t_on <- ex_ev$t0[1]; t_end <- ex_ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) row[, `:=`(t50_vagal = t50(xp, tt, t_on, mean(xp[rest]), mean(xp[endw])), t50_symp = t50(xs, tt, t_on, mean(xs[rest]), mean(xs[endw])), d_vagal = mean(xp[endw]) - mean(xp[rest]), d_symp = mean(xs[endw]) - mean(xs[rest]))] }
    if (job$ds == "prcp") { a1 <- ev_ann[as.character(id) == job$id & is.finite(onset_s)]; an <- a1$onset_s - job$off; ret <- a1[is.finite(return_s), return_s] - job$off; every <- c(an, ret)
      row[, `:=`(n_annot = length(an), hits = if (length(onsets)) sum(vapply(an, function(a) min(abs(onsets - a)) <= 30, TRUE)) else 0L, false_alarms = if (length(onsets)) sum(vapply(onsets, function(o) all(abs(every - o) > 30), TRUE)) else 0L)] }
    rows[[m]] <<- row }
  e <- state_estimates(base, FALSE); stats("M3", e$xp, e$xs, e$r$std_innov, expl(e$r$LL_star), base$tau_p, base$tau_s, numeric(0))
  if (isTRUE(fo$ok)) { e <- state_estimates(fo, FALSE); stats("M3O", e$xp, e$xs, e$r$std_innov, expl(e$r$LL_star), fo$tau_p, fo$tau_s, numeric(0)) }
  if (!is.null(cp)) { h0 <- fit_regime(series, 2, cp, base, 7, 0.7); if (isTRUE(h0$ok)) { r <- h0$run(); stats("HMM0", r$X_p, r$X_s, r$z, r$LL, h0$tau_p, h0$tau_s, regime_onsets2(r$W, tt)) }
    h <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), REG)); if (isTRUE(h$ok)) { r <- h$run(); stats("HMMr", r$X_p, r$X_s, r$z, r$LL, h$tau_p, h$tau_s, regime_onsets2(r$W, tt)) }
    if (isTRUE(fo$ok)) { h2 <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), c(REG, list(osc = fo$p$osc)))); if (isTRUE(h2$ok)) { r <- h2$run(); stats("HMMOr", r$X_p + r$X_o, r$X_s, r$z, r$LL, h2$tau_p, h2$tau_s, regime_onsets2(r$W, tt)) } } }
  # held-out forward score
  n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); bh <- fit_ext(sh); foh <- fit_ext(sh, list(osc = TRUE)); cph <- cp_search(sh, bh, min_gain = -Inf); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]
  hr <- bh$run(); r3 <- ext_run(dy, ctl, bh$p$kp, bh$p$ks, sh$priors$energy_ratio, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2); m3 <- mean(r3$logdens[idx]); mo <- NA_real_
  if (isTRUE(foh$ok)) { hro <- foh$run(); ro <- ext_run(dy, ctl, foh$p$kp, foh$p$ks, sh$priors$energy_ratio, foh$p$lR, osc = foh$p$osc, gm = gm, full = TRUE, beta_fixed = hro$nu0, sig2_fixed = hro$sig2); mo <- mean(ro$logdens[idx]) }
  hm0 <- hm <- hmo <- NA_real_; if (!is.null(cph)) { hh0 <- fit_regime(sh, 2, cph, bh, 7, 0.7); if (isTRUE(hh0$ok)) hm0 <- mean(hh0$run(dy, gm)$logdens[idx])
    hh <- do.call(fit_regime, c(list(sh, 2, cph, bh, 7, 0.7), REG)); if (isTRUE(hh$ok)) hm <- mean(hh$run(dy, gm)$logdens[idx])
    if (isTRUE(foh$ok)) { hh2 <- do.call(fit_regime, c(list(sh, 2, cph, bh, 7, 0.7), c(REG, list(osc = foh$p$osc)))); if (isTRUE(hh2$ok)) hmo <- mean(hh2$run(dy, gm)$logdens[idx]) } }
  held <- data.table(ds = job$ds, id = job$id, m3 = m3, m3o = mo, hmm0 = hm0, hmmr = hm, hmmor = hmo, gain_head = if (is.null(cph)) NA_real_ else cph$gain, tail_bad = mean(dy[idx] > 2))
  list(rows = rbindlist(rows, fill = TRUE), held = held) }
log_msg("Round-2 real data | %d recordings", length(jobs)); out <- drop_failed(map_recordings(jobs, one, n_cores), "recording")
res <- rbindlist(lapply(out, `[[`, "rows"), fill = TRUE); held <- rbindlist(lapply(out, `[[`, "held")); saveRDS(list(res = res, held = held, gate = GATE), "experiments/regimes/v2/_out/real.rds"); cat("Saved", nrow(res), "rows for", uniqueN(res$id), "recordings\n")
