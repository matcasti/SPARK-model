# EXPERIMENTAL. Regime models on real recordings: exercise (100 random recordings), Fantasia (40, rest: no regime expected), PRCP (10, with the annotated tilt and stand-up onsets as an independent check).
# Models: M3; M3E_det (events from the HR segmentation); CP (one rectangular regime estimated with the SPARK likelihood); HMM2 (2 latent regimes). The gated versions (regimes only if the change-point
# gain >= gate) are formed afterwards. Measures: KS pass, PACF, log-likelihood (explicit, AIC-penalised), held-out forward score (fit first 70%, regime inference continues online in the last 30%),
# time constants, realism of the exercise response (vagal faster than sympathetic, expected signs) and, in PRCP, hits and false alarms against the annotated onsets.
# Run from the project root:  Rscript experiments/regimes/run_real.R [stay_prior] [gate]
args <- commandArgs(TRUE); STAY <- if (length(args) >= 1) as.numeric(args[1]) else 7; GATE <- if (length(args) >= 2) as.numeric(args[2]) else 10
source("analysis/00_setup.R"); source("experiments/regimes/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); ev_ann <- fread("data/prcp_events.csv")
jobs <- list(); for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (ds == "exercise") { set.seed(5); ids <- sort(sample(ids, 100)) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]], off = if (ds == "prcp") d$rec[[id]]$time[1] - d$rec[[id]]$rr[1] else 0) }
z_diag <- function(z, gm, ctl) { k <- which(is.finite(z) & gm < ctl$diag_gate_max); z <- z[k]; n <- length(z); c(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n))) }
t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2; k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
one <- function(job) {
  dy <- job$rr; series <- prepare_series(dy); dy <- series$dy; n <- length(dy); tt <- cumsum(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); n_eff <- sum(gm < ctl$lik_gate_max); expl <- function(llstar) llstar - n_eff / 2 * (log(2 * pi) + 1)
  mode <- c(exercise = "exercise", prcp = "tilt")[job$ds]; ev_det <- if (!is.na(mode)) events_from_rr(dy, mode, cfg$segmentation) else NULL
  base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); cp_ok <- !is.null(cp) && cp$gain >= GATE; h <- if (!is.null(cp)) fit_regime(series, 2, cp, base, stay_prior = STAY, stay_sd = 0.7) else list(ok = FALSE)
  fits <- list(M3 = base, M3E_det = if (is.null(ev_det)) base else fit_ext(series, list(inputs = list(events = ev_det))), CP = if (cp_ok) fit_ext(series, list(inputs = list(events = data.frame(t0 = cp$t0, t1 = cp$t1)))) else base)
  ex_ev <- if (job$ds == "exercise") events_from_rr(dy, "exercise", cfg$segmentation) else NULL; rows <- list()
  stats <- function(m, xp, xs, z, ll, n_par, tau_p, tau_s, onsets) { dg <- z_diag(z, gm, ctl); row <- data.table(ds = job$ds, id = job$id, model = m, ks_ok = dg[["ks"]], pacf = dg[["pacf"]], dLL = ll - expl(base$run()$LL_star), n_par = n_par, tau_p = tau_p, tau_s = tau_s, gain = if (is.null(cp)) NA_real_ else cp$gain, n_onsets = length(onsets))
    if (!is.null(ex_ev)) { t_on <- ex_ev$t0[1]; t_end <- ex_ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) row[, `:=`(t50_vagal = t50(xp, tt, t_on, mean(xp[rest]), mean(xp[endw])), t50_symp = t50(xs, tt, t_on, mean(xs[rest]), mean(xs[endw])), d_vagal = mean(xp[endw]) - mean(xp[rest]), d_symp = mean(xs[endw]) - mean(xs[rest]))] }
    if (job$ds == "prcp") { a1 <- ev_ann[as.character(id) == job$id & is.finite(onset_s)]; an <- a1$onset_s - job$off; ret <- a1[is.finite(return_s), return_s] - job$off; every <- c(an, ret)
      row[, `:=`(n_annot = length(an), n_changes = length(every), hits = if (length(onsets)) sum(vapply(an, function(a) min(abs(onsets - a)) <= 30, TRUE)) else 0L, hits_any = if (length(onsets)) sum(vapply(every, function(a) min(abs(onsets - a)) <= 30, TRUE)) else 0L, false_alarms = if (length(onsets)) sum(vapply(onsets, function(o) all(abs(every - o) > 30), TRUE)) else 0L)] }
    rows[[m]] <<- row }
  for (m in names(fits)) { f <- fits[[m]]; if (!isTRUE(f$ok)) next; e <- state_estimates(f, FALSE); r <- e$r; on <- if (m == "M3E_det" && !is.null(ev_det)) ev_det$t0 else if (m == "CP" && cp_ok) cp$t0 else numeric(0)
    npar <- c(M3 = 5, M3E_det = 7, CP = 9)[[m]]; if (m == "CP" && !cp_ok) npar <- 5; if (m == "M3E_det" && is.null(ev_det)) npar <- 5
    stats(m, e$xp, e$xs, r$std_innov, expl(r$LL_star), npar, f$tau_p, f$tau_s, on) }
  if (isTRUE(h$ok)) { r <- h$run(); stats("HMM2", r$X_p, r$X_s, r$z, r$LL, h$n_par, h$tau_p, h$tau_s, regime_onsets2(r$W, tt)) }
  # held-out forward score (first 70% -> last 30%): M3 and HMM2 (online regime inference in the tail)
  n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); hh <- if (!is.null(cph)) fit_regime(sh, 2, cph, bh, stay_prior = STAY, stay_sd = 0.7) else list(ok = FALSE); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]
  hr <- bh$run(); r3 <- ext_run(dy, ctl, bh$p$kp, bh$p$ks, sh$priors$energy_ratio, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2)
  held <- data.table(ds = job$ds, id = job$id, m3_tail = mean(r3$logdens[idx]), hmm_tail = if (isTRUE(hh$ok)) mean(hh$run(dy, gm)$logdens[idx]) else NA_real_, gain_head = if (is.null(cph)) NA_real_ else cph$gain)
  list(rows = rbindlist(rows, fill = TRUE), held = held) }
log_msg("Regime experiments, real data | %d recordings | stay prior %g | gate %g", length(jobs), STAY, GATE); out <- drop_failed(map_recordings(jobs, one, n_cores), "recording")
res <- rbindlist(lapply(out, `[[`, "rows"), fill = TRUE); held <- rbindlist(lapply(out, `[[`, "held")); dir.create("experiments/regimes/_out", showWarnings = FALSE); saveRDS(list(res = res, held = held, stay = STAY, gate = GATE), "experiments/regimes/_out/real.rds")
cat("Saved", nrow(res), "rows for", uniqueN(res$id), "recordings\n")
