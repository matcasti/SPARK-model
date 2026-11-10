# ==============================================================================
# experiments/hmm_states/sim.R -- the second-stage HMM on simulated recordings with known truth.
#   Rscript experiments/hmm_states/sim.R test        recovery of regime means and boundaries on synthetic states
#   Rscript experiments/hmm_states/sim.R             null (no regimes) for the gate, then exercise-like and tilt-like recordings
# ==============================================================================
source("experiments/hmm_states/lib.R")
args <- commandArgs(TRUE)
if (length(args) && args[1] == "test") { set.seed(4); T <- 200; V <- c(1e-3, 4e-4); truth <- c(rep(1, 80), rep(2, 60), rep(1, 60)); mu <- rbind(c(0, 0), c(-0.04, 0.05))
  vt <- cbind(runif(T, 2e-4, 8e-4), runif(T, 1e-4, 4e-4)); x <- t(sapply(seq_len(T), function(i) mu[truth[i], ] + rnorm(2, 0, sqrt(c(V[1], V[2])))))   # truth around the regime means, then measured with the known error
  y <- x + cbind(rnorm(T, 0, sqrt(vt[, 1])), rnorm(T, 0, sqrt(vt[, 2]))); Vtot <- V + 0 * V
  f <- hmm_fit(y, vt, V - 0 * V, 2, "aware"); acc <- max(mean(f$path == truth), mean(f$path == 3 - truth)); cat(sprintf("synthetic: path accuracy %.2f; means\n", acc)); print(round(f$mu, 3)); cat("true:\n"); print(mu)
  # null: no regime; the log-likelihood gain of K = 2 over K = 1 under the aware emission
  lrs <- replicate(60, { x0 <- t(sapply(1:T, function(i) rnorm(2, 0, sqrt(V)))); y0 <- x0 + cbind(rnorm(T, 0, sqrt(vt[, 1])), rnorm(T, 0, sqrt(vt[, 2]))); b <- data.table(t = 1:T * 5, y1 = y0[, 1], y2 = y0[, 2], v1 = vt[, 1], v2 = vt[, 2]); detect_lr(b, V)$lr })
  cat(sprintf("null LR (independent bins): median %.1f, 95th percentile %.1f\n", median(lrs), quantile(lrs, .95))); quit(save = "no") }

source("analysis/00_setup.R"); source("experiments/m3plus/lib.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); DUR <- 600; set.seed(1)
states_from_rr <- function(rr, t = cumsum(rr)) { series <- prepare_series(rr); f <- fit_variant(series, model_spec("M3")); if (!isTRUE(f$ok)) return(NULL)
  list(b = bin_states(t, f$X_p, f$X_s, f$X_p_var, f$X_s_var), V = c(series$priors$energy_ratio * f$sig2, f$sig2)) }
stats_one <- function(st, gate = NULL) { b <- st$b; out <- list()
  for (m in c("aware", "naive")) { d <- detect_lr(b, st$V, m); e <- regime_events(d$fit2, b$t, n0 = 15L); f3 <- hmm_fit(as.matrix(b[, .(y1, y2)]), as.matrix(b[, .(v1, v2)]), st$V, 3, m)
    bic <- c(-2 * d$fit1$ll + d$fit1$npar * log(nrow(b)), -2 * d$fit2$ll + d$fit2$npar * log(nrow(b)), -2 * f3$ll + f3$npar * log(nrow(b)))
    out[[m]] <- data.table(mode = m, lr = d$lr, K_bic = which.min(bic), onset = e$onset, end = e$end, away = e$away, shift1 = diff(d$fit2$mu[order(match(seq_len(2), c(e$base, 3 - e$base))), 1]), shift2 = diff(d$fit2$mu[order(match(seq_len(2), c(e$base, 3 - e$base))), 2])) }
  rbindlist(out) }
# ---- null: no regimes -------------------------------------------------------------------------------------------------------------
null_job <- function(job) { truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]()); sim <- simulate_rr(DUR, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  st <- states_from_rr(sim$rr); if (is.null(st)) return(NULL); cbind(scenario = job$scenario, sim = job$sim, stats_one(st)) }
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(1:25, function(i) list(scenario = sc, sim = i, seed = id_seed(paste("hmmnull", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Null: %d stationary recordings", length(jobs)); nul <- rbindlist(drop_failed(map_recordings(jobs, null_job, n_cores), "null")); saveRDS(nul, "experiments/hmm_states/_out/null.rds")
gate <- nul[, .(gate95 = quantile(lr, .95), median = median(lr), K3_share = mean(K_bic == 3), K2_share = mean(K_bic == 2)), by = mode]; cat("\n== Null (stationary M3 drives, 100 recordings): LR of K = 2 over K = 1 ==\n"); print(gate, row.names = FALSE)
g <- setNames(gate$gate95, gate$mode)
# ---- known perturbations -----------------------------------------------------------------------------------------------------------
osc_draw <- function() list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1))))
build <- function(sc, seed) { set.seed(seed); d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))), osc = osc_draw())
  mk <- function(mp, ms) sim_ext(DUR, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = d$osc, m_p_fun = mp, m_s_fun = ms)
  if (sc == "exercise") { p <- draw_exercise_pars(seed, DUR); return(c(d, list(true = c(p$t_on, p$t_off), s = mk(prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s))))) }
  t_on <- runif(1, 0.28, 0.45) * DUR; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
  c(d, list(true = c(t_on, t_off), s = mk(prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))), prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100))))))) }
ev_job <- function(job) { z <- build(job$scenario, job$seed); s <- z$s; st <- states_from_rr(s$rr, t = s$t); if (is.null(st)) return(NULL); o <- stats_one(st); tp <- s$x_p + s$z1; ts <- s$x_s; tt <- s$t
  chg <- function(x) mean(x[tt >= z$true[2] - 30 & tt < z$true[2]]) - mean(x[tt >= z$true[1] - 90 & tt < z$true[1] - 10])
  cbind(scenario = job$scenario, sim = job$sim, o, t_on = z$true[1], t_off = z$true[2], true_shift1 = chg(tp), true_shift2 = chg(ts)) }
jobs <- unlist(lapply(c("exercise", "tilt"), function(sc) lapply(1:60, function(i) list(scenario = sc, sim = i, seed = id_seed(paste("hmmev", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Perturbations: %d recordings", length(jobs)); ev <- rbindlist(drop_failed(map_recordings(jobs, ev_job, n_cores), "perturbation")); saveRDS(ev, "experiments/hmm_states/_out/events.rds")
ev[, `:=`(detected = lr > g[mode], onset_err = onset - t_on, end_err = end - t_off, ratio1 = shift1 / true_shift1, ratio2 = shift2 / true_shift2)]
cat("\n== Known perturbations (truth: onset and end of the profile), K = 2, detection at the null 95th percentile ==\n"); options(width = 200)
print(ev[, .(n = .N, detected = round(mean(detected), 2), onset_hit_30s = round(mean(detected & abs(onset_err) <= 30, na.rm = TRUE), 2), onset_err_median_s = round(median(onset_err[detected], na.rm = TRUE), 1), onset_err_IQR = sprintf("%.0f..%.0f", quantile(onset_err[detected], .25, na.rm = TRUE), quantile(onset_err[detected], .75, na.rm = TRUE)),
              end_err_median_s = round(median(end_err[detected], na.rm = TRUE), 1), shift_x_p_ratio = round(median(ratio1[detected], na.rm = TRUE), 2), shift_x_s_ratio = round(median(ratio2[detected], na.rm = TRUE), 2), K_bic_is_2 = round(mean(K_bic == 2), 2)), by = .(scenario, mode)][order(scenario, mode)], row.names = FALSE)
