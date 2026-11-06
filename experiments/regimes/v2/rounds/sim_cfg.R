# EXPERIMENTAL (round 2). Which regularisation helps or hurts state reconstruction? HMM only, scenarios exercise and regime, with the oscillator in the simulation. Run: Rscript experiments/regimes/v2/sim_cfg.R [N]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 30L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
CONF <- list(base = list(), mu = list(mu_sd = 0.4, mu_max = 1.5), mu_wide = list(mu_sd = 0.8, mu_max = 3), anc = list(anchor_sd = 0.5), anc_wide = list(anchor_sd = 1.5))
score <- function(est, truth) sqrt(mean((est - truth)^2)) / sd(truth)
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); e <- state_estimates(base, FALSE)
  rows <- list(data.table(conf = "M3", rmse_p = score(e$xp, s$x_p + s$z1), rmse_s = score(e$xs, s$x_s), tau_p = base$tau_p, gain = if (is.null(cp)) NA_real_ else cp$gain, true_tau_p = z$tau_p))
  if (!is.null(cp)) for (cn in names(CONF)) { h <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), CONF[[cn]])); if (!isTRUE(h$ok)) next; r <- h$run()
    rows[[length(rows) + 1]] <- data.table(conf = cn, rmse_p = score(r$X_p, s$x_p + s$z1), rmse_s = score(r$X_s, s$x_s), tau_p = h$tau_p, gain = cp$gain, true_tau_p = z$tau_p) }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim)] }
jobs <- unlist(lapply(c("exercise", "regime"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/regimes/v2/_out/sim_cfg.rds")
m <- res[conf == "M3", .(scenario, sim, p0 = rmse_p, s0 = rmse_s)]; x <- merge(res[conf != "M3"], m, by = c("scenario", "sim"))
print(x[, .(n = .N, vagal = round(100 * (mean(rmse_p) / mean(p0) - 1)), symp = round(100 * (mean(rmse_s) / mean(s0) - 1)), tau_p_med = round(median(tau_p), 1), true_tau_p = round(median(true_tau_p), 1)), by = .(scenario, conf)][order(scenario, conf)])
