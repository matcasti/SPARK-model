# EXPERIMENTAL (round 2). Smoothing of the regime model: filtered states, GPB1 + RTS, and Kim's regime smoother + soft-input RTS, against M3 filtered and M3 RTS. Run: Rscript experiments/regimes/v2/sim_smooth.R [N]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 24L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); f <- function(x, t) sqrt(mean((x - t)^2)) / sd(t)
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = 9); if (is.null(cp)) return(NULL)
  h <- fit_regime(series, 2, cp, base, 7, 0.7); if (!isTRUE(h$ok)) return(NULL); tp <- s$x_p + s$z1; e <- state_estimates(base, FALSE); es <- state_estimates(base, TRUE); r <- h$run(smooth = TRUE); sm <- regime_states_smooth(h)
  data.table(scenario = job$scenario, sim = job$sim, M3_filt_p = f(e$xp, tp), M3_rts_p = f(es$xp, tp), HMM_filt_p = f(r$X_p, tp), HMM_rts_p = f(r$Xs_p, tp), HMM_kim_p = f(sm$X_p, tp),
             M3_filt_s = f(e$xs, s$x_s), M3_rts_s = f(es$xs, s$x_s), HMM_filt_s = f(r$X_s, s$x_s), HMM_rts_s = f(r$Xs_s, s$x_s), HMM_kim_s = f(sm$X_s, s$x_s)) }
jobs <- unlist(lapply(c("exercise", "regime", "two"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/regimes/v2/_out/sim_smooth.rds"); print(res[, c(list(n = .N), lapply(.SD, function(x) round(mean(x), 3))), by = scenario, .SDcols = patterns("_p$|_s$")])
