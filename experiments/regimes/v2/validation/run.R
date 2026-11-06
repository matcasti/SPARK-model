# Frozen validation run (see protocol.md). Run from the project root:  Rscript experiments/regimes/v2/validation/run.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); GATE <- 9; ids <- readLines("experiments/regimes/v2/validation/ids.txt"); d <- load_dataset("exercise", cfg)
jobs <- lapply(ids, function(id) list(id = id, rr = d$rr[[id]]))
z_diag <- function(z, gm, ctl) { k <- which(is.finite(z) & gm < ctl$diag_gate_max); z <- z[k]; n <- length(z); c(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n))) }
t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2; k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
one <- function(job) {
  t_start <- Sys.time(); series <- prepare_series(job$rr); dy <- series$dy; n <- length(dy); tt <- cumsum(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); n_eff <- sum(gm < ctl$lik_gate_max); expl <- function(l) l - n_eff / 2 * (log(2 * pi) + 1)
  base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); open <- !is.null(cp) && cp$gain >= GATE; ll0 <- expl(base$run()$LL_star); ex_ev <- events_from_rr(dy, "exercise", cfg$segmentation); rows <- list()
  stats <- function(m, xp, xs, z, ll, tau_p, fb) { dg <- z_diag(z, gm, ctl); row <- data.table(id = job$id, model = m, fallback = fb, ks_ok = dg[["ks"]], pacf = dg[["pacf"]], dLL = ll - ll0, tau_p = tau_p)
    if (!is.null(ex_ev)) { t_on <- ex_ev$t0[1]; t_end <- ex_ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) row[, `:=`(t50_vagal = t50(xp, tt, t_on, mean(xp[rest]), mean(xp[endw])), t50_symp = t50(xs, tt, t_on, mean(xs[rest]), mean(xs[endw])), d_vagal = mean(xp[endw]) - mean(xp[rest]), d_symp = mean(xs[endw]) - mean(xs[rest]))] }
    rows[[m]] <<- row }
  e <- state_estimates(base, FALSE); stats("M3", e$xp, e$xs, e$r$std_innov, ll0, base$tau_p, FALSE); m3 <- function(m, fb = TRUE) stats(m, e$xp, e$xs, e$r$std_innov, ll0, base$tau_p, fb)
  ev_det <- events_from_rr(dy, "exercise", cfg$segmentation); fd <- if (is.null(ev_det)) NULL else fit_ext(series, list(inputs = list(events = ev_det))); if (!is.null(fd) && isTRUE(fd$ok)) { ed <- state_estimates(fd, FALSE); stats("M3E_det", ed$xp, ed$xs, ed$r$std_innov, expl(ed$r$LL_star), fd$tau_p, FALSE) } else m3("M3E_det")
  ok_h <- FALSE; ok_y <- FALSE
  if (open) { h <- fit_regime(series, 2, cp, base, 7, 0.7, recip = 20); if (isTRUE(h$ok)) { ok_h <- TRUE; r <- h$run(smooth = TRUE); stats("HMMg", r$X_p, r$X_s, r$z, r$LL, h$tau_p, FALSE)
      Ws <- regime_smooth(r$W, r$A, h$p$Pi); ev_h <- regime_episodes(Ws, tt, max_ev = 8); fy <- if (is.null(ev_h)) NULL else tryCatch(fit_ext(series, list(inputs = list(events = ev_h))), error = function(e) NULL)
      if (!is.null(fy) && isTRUE(fy$ok)) { ok_y <- TRUE; ey <- state_estimates(fy, FALSE); stats("HYB", ey$xp, ey$xs, ey$r$std_innov, expl(ey$r$LL_star), fy$tau_p, FALSE) } } }
  if (!ok_h) m3("HMMg"); if (!ok_y) m3("HYB")
  # forward score (M3 and HMMg)
  n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]; hr <- bh$run()
  r3 <- ext_run(dy, ctl, bh$p$kp, bh$p$ks, sh$priors$energy_ratio, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2); m3tail <- mean(r3$logdens[idx]); hmm <- m3tail; open_h <- !is.null(cph) && cph$gain >= GATE
  if (open_h) { hh <- fit_regime(sh, 2, cph, bh, 7, 0.7, recip = 20); if (isTRUE(hh$ok)) hmm <- mean(hh$run(dy, gm)$logdens[idx]) }
  held <- data.table(id = job$id, m3 = m3tail, hmmg = hmm, open_head = open_h, tail_bad = mean(dy[idx] > 2), n_tail = length(idx))
  list(rows = rbindlist(rows, fill = TRUE)[, `:=`(open = open, gain = if (is.null(cp)) NA_real_ else cp$gain, secs = as.numeric(Sys.time() - t_start, units = "secs"))], held = held) }
log_msg("Frozen validation | %d untouched exercise recordings", length(jobs)); out <- drop_failed(map_recordings(jobs, one, n_cores), "recording")
saveRDS(list(res = rbindlist(lapply(out, `[[`, "rows"), fill = TRUE), held = rbindlist(lapply(out, `[[`, "held")), ids = ids, n_failed = length(jobs) - length(out)), "experiments/regimes/v2/_out/validation.rds"); cat("Saved;", length(out), "of", length(jobs), "recordings\n")
