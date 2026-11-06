# EXPERIMENTAL (round 2). Simulations with truth (scenarios of round 1, now with a respiratory oscillator in the vagal pathway): state RMSE and 95% coverage of the filtered and smoothed
# estimates for M3, M3O, the regularised 2-regime HMM (HMM) and the HMM with the oscillator (HMMO); gated versions (regimes only if the change-point gain >= 9) are formed in the summary.
# Run from the project root:  Rscript experiments/regimes/v2/sim.R [N per scenario]
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 30L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
REG <- list(recip = 20)
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
one <- function(job) {
  z <- sim_regime(job$scenario, job$seed, osc = list(alpha = 0.2, f_hz = 0.25, lo = 0.3)); s <- z$s; series <- prepare_series(s$rr); dy <- series$dy; tt <- cumsum(dy); truth_p <- s$x_p + s$z1; truth_s <- s$x_s
  base <- fit_ext(series); fo <- fit_ext(series, list(osc = TRUE)); cp <- cp_search(series, base, min_gain = -Inf); rows <- list(); dy <- series$dy; n_eff <- sum(gate_mult(dy, series$control) < series$control$lik_gate_max); expl <- function(l) l - n_eff / 2 * (log(2 * pi) + 1); ll_m3 <- expl(base$run()$LL_star); ll_m3o <- if (isTRUE(fo$ok)) expl(fo$run()$LL_star) else NA_real_
  add <- function(model, mode, xp, xpv, xs, xsv, gain, dll = NA_real_) { a <- score(xp, xpv, truth_p); b <- score(xs, xsv, truth_s); rows[[length(rows) + 1]] <<- data.table(model = model, mode = mode, rmse_p = a[["rmse"]], rmse_s = b[["rmse"]], cover_p = a[["cover"]], cover_s = b[["cover"]], gain = gain, dll = dll) }
  g <- if (is.null(cp)) NA_real_ else cp$gain
  for (sm in c(FALSE, TRUE)) { e <- state_estimates(base, sm); add("M3", if (sm) "smooth" else "filter", e$xp, e$xpv, e$xs, e$xsv, g)
    if (isTRUE(fo$ok)) { e <- state_estimates(fo, sm); add("M3O", if (sm) "smooth" else "filter", e$xp, e$xpv, e$xs, e$xsv, g) } }
  if (!is.null(cp)) { h <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), REG)); if (isTRUE(h$ok)) { r <- h$run(smooth = TRUE); add("HMM", "filter", r$X_p, r$X_p_var, r$X_s, r$X_s_var, g, r$LL - ll_m3); add("HMM", "smooth", r$Xs_p, r$Xs_p_var, r$Xs_s, r$Xs_s_var, g, r$LL - ll_m3) }
    if (isTRUE(fo$ok)) { h2 <- do.call(fit_regime, c(list(series, 2, cp, base, 7, 0.7), c(REG, list(osc = fo$p$osc)))); if (isTRUE(h2$ok)) { r <- h2$run(smooth = TRUE); add("HMMO", "filter", r$X_p + r$X_o, r$X_v_var, r$X_s, r$X_s_var, g, r$LL - ll_m3o); add("HMMO", "smooth", r$Xs_p + r$Xs_o, r$Xs_v_var, r$Xs_s, r$Xs_s_var, g, r$LL - ll_m3o) } } }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim)] }
jobs <- unlist(lapply(c("regime", "exercise", "two", "rest"), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes2", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Round-2 simulations | %d jobs", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"), fill = TRUE); saveRDS(res, "experiments/regimes/v2/_out/sim.rds"); cat("Saved", nrow(res), "rows\n")
