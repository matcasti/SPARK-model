# EXPERIMENTAL (round 2). Screening of regularised regime fits on a subset: exercise (40, including the six worst held-out cases of round 1), PRCP (10), Fantasia (20).
# Per recording and configuration: in-sample (dLL vs M3, KS, PACF, tau_p, mu) and held-out forward score (fit first 70%, last 30% scored), plus the change-point gain.
# Run from the project root:  Rscript experiments/regimes/v2/screen.R
source("analysis/00_setup.R"); source("experiments/regimes/lib.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
CONF <- list(base = list(), reg = list(mu_sd = 0.4, mu_max = 1.5), anc = list(anchor_sd = 0.5), reganc = list(mu_sd = 0.4, mu_max = 1.5, anchor_sd = 0.5))
jobs <- list(); hard <- c("104", "099", "211", "208", "207", "138")
for (ds in c("exercise", "prcp", "fantasia")) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (ds == "exercise") { set.seed(11); ids <- sort(unique(c(hard, sample(setdiff(ids, hard), 34)))) }; if (ds == "fantasia") { set.seed(12); ids <- sort(sample(ids, 20)) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
z_diag <- function(z, gm, ctl) { k <- which(is.finite(z) & gm < ctl$diag_gate_max); z <- z[k]; n <- length(z); c(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n))) }
one <- function(job) {
  series <- prepare_series(job$rr); dy <- series$dy; n <- length(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); n_eff <- sum(gm < ctl$lik_gate_max); expl <- function(l) l - n_eff / 2 * (log(2 * pi) + 1)
  base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); ll0 <- expl(base$run()$LL_star); eM3 <- state_estimates(base, FALSE); dgM <- z_diag(eM3$r$std_innov, gm, ctl)
  n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]
  hr <- bh$run(); r3 <- ext_run(dy, ctl, bh$p$kp, bh$p$ks, sh$priors$energy_ratio, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2); rows <- list()
  for (cn in names(CONF)) { a <- CONF[[cn]]; h <- if (!is.null(cp)) do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), a)) else list(ok = FALSE); hh <- if (!is.null(cph)) do.call(fit_regime, c(list(sh, 2, cph, bh, 7, 0.7), a)) else list(ok = FALSE)
    row <- data.table(ds = job$ds, id = job$id, conf = cn, gain = if (is.null(cp)) NA_real_ else cp$gain, gain_head = if (is.null(cph)) NA_real_ else cph$gain, ks_m3 = dgM[["ks"]], pacf_m3 = dgM[["pacf"]])
    if (isTRUE(h$ok)) { r <- h$run(); dg <- z_diag(r$z, gm, ctl); row[, `:=`(dLL = expl(r$LL) - ll0, ks = dg[["ks"]], pacf = dg[["pacf"]], tau_p = h$tau_p, tau_s = h$tau_s, mu_p = h$p$mu[2, 1], mu_s = h$p$mu[2, 2], sig2 = h$p$sig2, lR = h$p$lR)] }
    if (isTRUE(hh$ok)) { ld <- hh$run(dy, gm)$logdens; row[, `:=`(held = mean(ld[idx]) - mean(r3$logdens[idx]), held_worst = min(ld[idx] - r3$logdens[idx]))] }
    rows[[cn]] <- row }
  rbindlist(rows, fill = TRUE) }
log_msg("Round-2 screening of regularised regime fits | %d recordings x %d configurations", length(jobs), length(CONF)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording"), fill = TRUE)
saveRDS(res, "experiments/regimes/v2/_out/screen.rds"); cat("Saved", nrow(res), "rows\n")
