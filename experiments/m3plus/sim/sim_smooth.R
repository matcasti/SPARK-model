# EXPERIMENTAL. Ideas 2 and 3: RTS smoother and propagation of parameter uncertainty (Laplace samples) for the states of M3, scored against the truth on the
# pipeline's simulated scenarios (same seeds). Run from the project root:  Rscript experiments/m3plus/sim_smooth.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); N_SAMPLES <- 40L
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cor = cor(est, truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
one <- function(job) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  series <- prepare_series(sim$rr); ctl <- series$control; pri <- series$priors; dy <- series$dy
  fit <- fit_variant(series, model_spec("M3"), full = FALSE); if (!isTRUE(fit$ok)) return(NULL)
  par <- fit$par; f <- ext_run(dy, ctl, par$kp, par$ks, par$lp, par$lR, full = TRUE, smooth = TRUE)
  rows <- list()
  add <- function(method, xp, xpv, xs, xsv) rows[[length(rows) + 1]] <<- data.table(method = method, state = c("x_p", "x_s"),
        rbind(score(xp, xpv, sim$x_p), score(xs, xsv, sim$x_s)))
  add("filtered (MAP)", f$X_p, f$X_p_var, f$X_s, f$X_s_var)
  add("smoothed (MAP)", f$X_p_smooth, f$X_p_var_smooth, f$X_s_smooth, f$X_s_var_smooth)
  # Laplace samples of theta = (ln ks, ln(kp-ks), ln lR) from the numerical Hessian of the M3 objective
  spec <- model_spec("M3"); H <- tryCatch(numDeriv::hessian(function(th) variant_objective(th, dy, spec, pri, ctl), fit$theta), error = function(e) NULL)
  if (!is.null(H) && all(is.finite(H))) {
    ev <- eigen((H + t(H)) / 2, symmetric = TRUE); vals <- pmax(ev$values, 1e-3)                       # ridge: keep the covariance finite
    Sigma <- ev$vectors %*% diag(1 / vals) %*% t(ev$vectors); L <- tryCatch(chol(Sigma + 1e-10 * diag(3)), error = function(e) NULL)
    if (!is.null(L)) {
      th <- t(fit$theta + t(matrix(rnorm(N_SAMPLES * 3), N_SAMPLES) %*% L))
      runs <- lapply(seq_len(N_SAMPLES), function(i) { p <- unpack_theta(th[i, ], spec, pri$energy_ratio); ext_run(dy, ctl, p$kp, p$ks, p$lp, p$lR, full = TRUE, smooth = TRUE) })
      runs <- runs[vapply(runs, function(r) is.finite(r$LL_star), TRUE)]
      mix <- function(field, vfield) { M <- sapply(runs, `[[`, field); V <- sapply(runs, `[[`, vfield); list(mean = rowMeans(M), var = rowMeans(V) + apply(M, 1, var)) }
      a <- mix("X_p", "X_p_var"); b <- mix("X_s", "X_s_var"); add("filtered + Laplace", a$mean, a$var, b$mean, b$var)
      a <- mix("X_p_smooth", "X_p_var_smooth"); b <- mix("X_s_smooth", "X_s_var_smooth"); add("smoothed + Laplace", a$mean, a$var, b$mean, b$var)
    }
  }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim)]
}
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Smoother and Laplace on %d simulated recordings", length(jobs))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); res[, method := factor(method, c("filtered (MAP)", "smoothed (MAP)", "filtered + Laplace", "smoothed + Laplace"))]
dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(res, "experiments/m3plus/_out/sim_smooth.rds")
cat("\n== All scenarios (mean over simulations) ==\n"); print(res[, .(n = .N, rmse = round(mean(rmse), 3), cor = round(mean(cor), 3), cover95 = round(mean(cover), 3)), by = .(state, method)][order(state, method)], row.names = FALSE)
cat("\n== Coverage by scenario ==\n"); print(dcast(res, scenario + state ~ method, value.var = "cover", fun.aggregate = function(x) round(mean(x), 2)), row.names = FALSE)
