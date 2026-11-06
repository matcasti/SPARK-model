# EXPERIMENTAL (round 2). Reciprocity prior (hinge penalty on mu_p * mu_s > 0, i.e. both drives moving the rate in opposite directions) versus the agnostic prior; HMM only. Run: Rscript experiments/regimes/v2/sim_recip.R [N]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 30L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); score <- function(est, truth) sqrt(mean((est - truth)^2)) / sd(truth)
CONF <- list(agnostic = list(), recip20 = list(recip = 20), recip100 = list(recip = 100))
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = 9); e <- state_estimates(base, FALSE); tp <- s$x_p + s$z1
  rows <- list(data.table(conf = "M3", rmse_p = score(e$xp, tp), rmse_s = score(e$xs, s$x_s), open = !is.null(cp)))
  if (!is.null(cp)) for (cn in names(CONF)) { h <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), CONF[[cn]])); if (!isTRUE(h$ok)) next; r <- h$run(); rows[[length(rows) + 1]] <- data.table(conf = cn, rmse_p = score(r$X_p, tp), rmse_s = score(r$X_s, s$x_s), open = TRUE) }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim)] }
jobs <- unlist(lapply(c("exercise", "regime", "two"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/regimes/v2/_out/sim_recip.rds")
m <- res[conf == "M3", .(scenario, sim, p0 = rmse_p, s0 = rmse_s)]; x <- merge(res[conf != "M3"], m, by = c("scenario", "sim"))
print(x[, .(n = .N, vagal_mean = round(100 * (mean(rmse_p) / mean(p0) - 1)), symp_mean = round(100 * (mean(rmse_s) / mean(s0) - 1)), vagal_median_ratio = round(median(rmse_p / p0), 2), symp_median_ratio = round(median(rmse_s / s0), 2), worse_vagal_pct = round(100 * mean(rmse_p > p0)), worse_symp_pct = round(100 * mean(rmse_s > s0)), worst_ratio = round(max(rmse_p / p0), 1)), by = .(scenario, conf)][order(scenario, conf)])
