# EXPERIMENTAL (round 2). Graded 3-state chain versus the 2-regime HMM and M3 in the simulated scenarios. Run: Rscript experiments/regimes/v2/sim_chain.R [N]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 30L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
score <- function(est, var, truth) c(sqrt(mean((est - truth)^2)) / sd(truth), mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); tp <- s$x_p + s$z1
  add <- function(m, xp, xpv, xs, xsv, extra = NA_real_) { a <- score(xp, xpv, tp); b <- score(xs, xsv, s$x_s); data.table(model = m, rmse_p = a[1], cover_p = a[2], rmse_s = b[1], cover_s = b[2], extra = extra) }
  e <- state_estimates(base, FALSE); rows <- list(add("M3", e$xp, e$xpv, e$xs, e$xsv)); g <- if (is.null(cp)) NA_real_ else cp$gain
  if (!is.null(cp)) { h <- fit_regime(series, 2, cp, base, 7, 0.7); if (isTRUE(h$ok)) { r <- h$run(); rows[[2]] <- add("HMM2", r$X_p, r$X_p_var, r$X_s, r$X_s_var, rev(range(r$W[, 2]))[1]) }
    h3 <- fit_regime3c(series, cp, base); if (isTRUE(h3$ok)) { r <- h3$run(); rows[[3]] <- add("CHAIN3", r$X_p, r$X_p_var, r$X_s, r$X_s_var, h3$p$f) } }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim, gain = g)] }
jobs <- unlist(lapply(c("exercise", "regime", "two", "rest"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/regimes/v2/_out/sim_chain.rds")
m <- res[model == "M3", .(scenario, sim, p0 = rmse_p, s0 = rmse_s)]; x <- merge(res[model != "M3"], m, by = c("scenario", "sim"))
print(x[, .(n = .N, vagal = round(100 * (mean(rmse_p) / mean(p0) - 1)), symp = round(100 * (mean(rmse_s) / mean(s0) - 1)), cover_v = round(mean(cover_p), 2), cover_s = round(mean(cover_s), 2), f_median = round(median(extra), 2)), by = .(scenario, model)][order(scenario, model)])
