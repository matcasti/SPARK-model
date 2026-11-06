# EXPERIMENTAL (round 2). Reciprocity prior with one start versus three starts (cp targets; tau_p = 5 s with typical targets; tau_p = 5 s with cp targets), best objective kept. Run: Rscript experiments/regimes/v2/sim_multi.R [N]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 30L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); score <- function(est, truth) sqrt(mean((est - truth)^2)) / sd(truth)
CONF <- list(recip1 = list(recip = 20), recip3 = list(recip = 20, multi = TRUE))
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = 9); e <- state_estimates(base, FALSE); tp <- s$x_p + s$z1
  rows <- list(data.table(conf = "M3", rmse_p = score(e$xp, tp), rmse_s = score(e$xs, s$x_s), tau_p = base$tau_p, secs = 0))
  if (!is.null(cp)) for (cn in names(CONF)) { t0 <- Sys.time(); h <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), CONF[[cn]])); if (!isTRUE(h$ok)) next; r <- h$run(); rows[[length(rows) + 1]] <- data.table(conf = cn, rmse_p = score(r$X_p, tp), rmse_s = score(r$X_s, s$x_s), tau_p = h$tau_p, secs = as.numeric(Sys.time() - t0, units = "secs")) }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim, true_tau_p = z$tau_p)] }
jobs <- unlist(lapply(c("exercise", "regime", "two"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/regimes/v2/_out/sim_multi.rds")
m <- res[conf == "M3", .(scenario, sim, p0 = rmse_p, s0 = rmse_s)]; x <- merge(res[conf != "M3"], m, by = c("scenario", "sim"))
print(x[, .(n = .N, vagal_mean = round(100 * (mean(rmse_p) / mean(p0) - 1)), symp_mean = round(100 * (mean(rmse_s) / mean(s0) - 1)), vagal_med = round(median(rmse_p / p0), 2), worse_vagal_pct = round(100 * mean(rmse_p > p0)), worst = round(max(rmse_p / p0), 1), tau_p_med = round(median(tau_p), 1), secs = round(mean(secs), 1)), by = .(scenario, conf)][order(scenario, conf)])
