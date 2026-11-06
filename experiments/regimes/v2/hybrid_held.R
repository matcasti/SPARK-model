# EXPERIMENTAL (round 2b). Forward (held-out) score of the hybrid: everything is detected and fitted on the first 70% of the record; the last 30% is scored with the head's parameters and regressor
# amplitudes. Events found on the head; an episode still active at the head's end is assumed to persist to the end of the record (what a forecaster would assume).
# Models: M3, M3E_det (pipeline detector on the head), HYB (episodes of the smoothed regime probability on the head), HMMr (regime inferred online in the tail). Exercise 60 and PRCP 10.
# Run from the project root:  Rscript experiments/regimes/v2/hybrid_held.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); GATE <- 9
jobs <- list(); for (ds in c("exercise", "prcp")) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (ds == "exercise") { set.seed(21); ids <- sort(sample(ids, 60)) }; for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
persist <- function(ev, t_head, t_end) { if (is.null(ev) || !nrow(ev)) return(NULL); open <- ev$t1 >= t_head - 5; ev$t1[open] <- t_end + 1000; ev }
one <- function(job) {
  series <- prepare_series(job$rr); dy <- series$dy; n <- length(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); tt <- cumsum(dy); n_head <- floor(0.7 * n); t_head <- tt[n_head]; t_end <- tt[n]
  sh <- prepare_series(dy[seq_len(n_head)]); dyh <- sh$dy; tth <- cumsum(dyh); bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); open <- !is.null(cph) && cph$gain >= GATE; lp <- sh$priors$energy_ratio
  idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]; hr <- bh$run(); rows <- list(); tail_bad <- mean(dy[idx] > 2)
  score_ev <- function(ev) { f <- fit_ext(sh, list(inputs = list(events = ev))); if (!isTRUE(f$ok)) return(NA_real_); rr <- f$run(); X <- input_design(dy, f$p$kp, f$p$ks, ev)
    r <- ext_run(dy, ctl, f$p$kp, f$p$ks, lp, f$p$lR, Xextra = X, gm = gm, full = TRUE, beta_fixed = rr$beta, sig2_fixed = rr$sig2); mean(r$logdens[idx]) }
  r3 <- ext_run(dy, ctl, bh$p$kp, bh$p$ks, lp, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2); m3 <- mean(r3$logdens[idx])
  mode <- c(exercise = "exercise", prcp = "tilt")[[job$ds]]; ev_det <- persist(events_from_rr(dyh, mode, cfg$segmentation), t_head, t_end); det <- if (is.null(ev_det)) m3 else score_ev(ev_det)
  hmm <- hyb <- m3; if (open) { h <- fit_regime(sh, 2, cph, bh, 7, 0.7, recip = 20); if (isTRUE(h$ok)) { hmm <- mean(h$run(dy, gm)$logdens[idx]); r <- h$run(smooth = TRUE); Ws <- regime_smooth(r$W, r$A, h$p$Pi); ev_h <- persist(regime_episodes(Ws, tth, max_ev = 8), t_head, t_end)
      if (!is.null(ev_h)) { v <- tryCatch(score_ev(ev_h), error = function(e) NA_real_); if (is.finite(v)) hyb <- v } } }
  data.table(ds = job$ds, id = job$id, m3 = m3, det = det, hmm = hmm, hyb = hyb, open = open, gain_head = if (is.null(cph)) NA_real_ else cph$gain, n_det = if (is.null(ev_det)) 0L else nrow(ev_det), tail_bad = tail_bad) }
log_msg("Hybrid held-out | %d recordings", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording"), fill = TRUE); saveRDS(res, "experiments/regimes/v2/_out/hybrid_held.rds"); cat("Saved", nrow(res), "rows\n")
x <- res[tail_bad < 0.1]; l <- melt(x, id.vars = c("ds", "id", "m3", "open"), measure.vars = c("det", "hmm", "hyb"), variable.name = "model", value.name = "ll"); l[, d := ll - m3]
cat("\nHeld-out forward score minus M3 (nats/beat); degenerate tails excluded; HMM and HYB equal M3 where the head gate is closed\n")
print(l[, .(n = .N, open_pct = round(100 * mean(open)), mean = round(mean(d, na.rm = TRUE), 3), median = round(median(d, na.rm = TRUE), 3), trimmed = round(mean(d, trim = 0.1, na.rm = TRUE), 3), better_pct = round(100 * mean(d > 0, na.rm = TRUE)), worst = round(min(d, na.rm = TRUE), 2)), by = .(ds, model)][order(ds, match(model, c("det", "hmm", "hyb")))])
