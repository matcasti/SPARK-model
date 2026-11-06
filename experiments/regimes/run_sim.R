# EXPERIMENTAL. Regime models vs M3 on simulated recordings with known truth. Models: M3; M3E_true (true event times: an upper bound for event-based models); M3E_det (events from the HR segmentation of
# the pipeline); CP (one rectangular regime with onset/end estimated by the SPARK likelihood); HMM2 and HMM3 (2 or 3 latent regimes). Scenarios: regime (the model is correct), exercise
# (the pipeline's exercise profile: misspecified), two (two regimes, tilt-like) and rest (no regime: false alarms and cost of the extra flexibility).
# Run from the project root:  Rscript experiments/regimes/run_sim.R
source("analysis/00_setup.R"); source("experiments/regimes/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); N <- 60L
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cor = cor(est, truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed); s <- z$s; series <- prepare_series(s$rr); dy <- series$dy; tt <- cumsum(dy); off <- z$off; ctl <- series$control; n_eff <- sum(gate_mult(dy, ctl) < ctl$lik_gate_max)
  ev_true <- if (!is.null(z$ev)) data.frame(t0 = z$ev$t0 - off, t1 = z$ev$t1 - off) else NULL; mode <- if (job$scenario == "two") "tilt" else "exercise"
  base <- fit_ext(series); ev_det <- events_from_rr(dy, mode, cfg$segmentation); cp <- cp_search(series, base)
  fits <- list(M3 = base, M3E_true = if (!is.null(ev_true)) fit_ext(series, list(inputs = list(events = ev_true))), M3E_det = if (is.null(ev_det)) base else fit_ext(series, list(inputs = list(events = ev_det))),
               CP = if (is.null(cp)) base else fit_ext(series, list(inputs = list(events = data.frame(t0 = cp$t0, t1 = cp$t1)))))
  hm <- list(HMM2 = fit_regime(series, 2, cp, base), HMM3 = fit_regime(series, 3, cp, base))
  rest <- if (!is.null(ev_true)) tt >= cfg$exercise$rest_skip_s & tt < ev_true$t0[1] - 5 else rep(FALSE, length(tt)); endw <- if (!is.null(ev_true)) tt >= ev_true$t1[1] - 15 & tt <= ev_true$t1[1] else rep(FALSE, length(tt))
  rows <- list(); ll_m3 <- base$run()$LL_star - n_eff / 2 * (log(2 * pi) + 1)
  add <- function(m, xp, xpv, xs, xsv, tau_p, tau_s, onsets, ll, n_par) { a <- score(xp, xpv, s$x_p); b <- score(xs, xsv, s$x_s)
    chg <- function(est, tr) if (any(rest) && any(endw)) (mean(est[endw]) - mean(est[rest])) / (mean(tr[endw]) - mean(tr[rest])) else NA_real_
    rows[[m]] <<- data.table(model = m, rmse_p = a[["rmse"]], cor_p = a[["cor"]], cover_p = a[["cover"]], rmse_s = b[["rmse"]], cor_s = b[["cor"]], cover_s = b[["cover"]], chg_p = chg(xp, s$x_p), chg_s = chg(xs, s$x_s),
      tau_p_err = log(tau_p / z$tau_p), tau_s_err = log(tau_s / z$tau_s), n_onsets = length(onsets), onset_err = if (length(onsets) && !is.null(ev_true)) onsets[which.min(abs(onsets - ev_true$t0[1]))] - ev_true$t0[1] else NA_real_, dLL = ll - ll_m3, n_par = n_par) }
  for (m in names(fits)) { f <- fits[[m]]; if (is.null(f) || !isTRUE(f$ok)) next; e <- state_estimates(f, FALSE); ll <- e$r$LL_star - n_eff / 2 * (log(2 * pi) + 1)
    on <- if (m == "M3E_true") ev_true$t0 else if (m == "M3E_det") (if (is.null(ev_det)) numeric(0) else ev_det$t0) else if (m == "CP") (if (is.null(cp)) numeric(0) else cp$t0) else numeric(0)
    add(m, e$xp, e$xpv, e$xs, e$xsv, f$tau_p, f$tau_s, on, ll, f$n_par + 5 - length(f$theta) + length(f$theta) + if (!is.null(f$spec$inputs)) 4 else 2) }
  for (m in names(hm)) { h <- hm[[m]]; if (!isTRUE(h$ok)) next; r <- h$run(); add(m, r$X_p, r$X_p_var, r$X_s, r$X_s_var, h$tau_p, h$tau_s, regime_onsets(r$W, tt), r$LL, h$n_par) }
  out <- rbindlist(rows); out[, `:=`(scenario = job$scenario, sim = job$sim, cp_found = !is.null(cp))]; out }
jobs <- unlist(lapply(c("regime", "exercise", "two", "rest"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Regime experiments | %d simulated recordings", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"))
dir.create("experiments/regimes/_out", showWarnings = FALSE); saveRDS(res, "experiments/regimes/_out/sim.rds")
ord <- c("M3", "M3E_true", "M3E_det", "CP", "HMM2", "HMM3"); res[, model := factor(model, ord)]
for (sc in c("regime", "exercise", "two", "rest")) { cat(sprintf("\n===== %s =====\n", sc)); d <- res[scenario == sc]
  print(d[, .(n = .N, rmse_vagal = round(mean(rmse_p), 3), rmse_symp = round(mean(rmse_s), 3), cor_symp = round(mean(cor_s), 3), cover_vagal = round(mean(cover_p), 2), cover_symp = round(mean(cover_s), 2), chg_symp = round(median(chg_s, na.rm = TRUE), 2), chg_vagal = round(median(chg_p, na.rm = TRUE), 2),
            tau_p_abs_err = round(median(abs(tau_p_err)), 2), onset_err_med = round(median(onset_err, na.rm = TRUE), 1), onsets = round(mean(n_onsets), 2), dLL = round(median(dLL), 1)), by = model][order(model)], row.names = FALSE) }
