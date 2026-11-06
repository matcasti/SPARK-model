# EXPERIMENTAL. Idea 3 (variant): intervals that include the uncertainty of the anchor. lambda_P ~ lognormal(ln anchor, s^2) around the spectral anchor,
# kinetic parameters at the MAP; smoothed states; mixture of the per-sample Gaussians. Run from the project root:  Rscript experiments/m3plus/sim_anchor_unc.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); S_VALUES <- c(0.5, 1, 1.5); N_SAMPLES <- 40L
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))), width = mean(2 * 1.96 * sqrt(pmax(var, 0))) / sd(truth))
one <- function(job) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  series <- prepare_series(sim$rr); ctl <- series$control; pri <- series$priors; dy <- series$dy
  fit <- fit_variant(series, model_spec("M3"), full = FALSE); if (!isTRUE(fit$ok)) return(NULL)
  par <- fit$par; rows <- list()
  add <- function(method, f, field = "_smooth") rows[[length(rows) + 1]] <<- data.table(method = method, state = c("x_p", "x_s"), rbind(score(f[[paste0("xp", field)]], f[[paste0("xpv", field)]], sim$x_p), score(f[[paste0("xs", field)]], f[[paste0("xsv", field)]], sim$x_s)))
  base <- ext_run(dy, ctl, par$kp, par$ks, par$lp, par$lR, full = TRUE, smooth = TRUE)
  add("smoothed (MAP)", list(xp_smooth = base$X_p_smooth, xpv_smooth = base$X_p_var_smooth, xs_smooth = base$X_s_smooth, xsv_smooth = base$X_s_var_smooth))
  for (s in S_VALUES) {
    runs <- lapply(seq_len(N_SAMPLES), function(i) ext_run(dy, ctl, par$kp, par$ks, par$lp * exp(s * rnorm(1)), par$lR, full = TRUE, smooth = TRUE)); runs <- runs[vapply(runs, function(r) is.finite(r$LL_star), TRUE)]
    mix <- function(a, b) { M <- sapply(runs, `[[`, a); V <- sapply(runs, `[[`, b); list(m = rowMeans(M), v = rowMeans(V) + apply(M, 1, var)) }
    p <- mix("X_p_smooth", "X_p_var_smooth"); q <- mix("X_s_smooth", "X_s_var_smooth")
    add(sprintf("smoothed + anchor uncertainty (s=%.1f)", s), list(xp_smooth = p$m, xpv_smooth = p$v, xs_smooth = q$m, xsv_smooth = q$v))
  }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim)]
}
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"))
saveRDS(res, "experiments/m3plus/_out/sim_anchor_unc.rds")
cat("\n== Coverage / width (95% bands) by method and state, all scenarios ==\n"); print(res[, .(rmse = round(mean(rmse), 3), cover95 = round(mean(cover), 3), width_sd = round(mean(width), 2)), by = .(state, method)][order(state, method)], row.names = FALSE)
cat("\n== Coverage by scenario ==\n"); print(dcast(res, scenario + state ~ method, value.var = "cover", fun.aggregate = function(x) round(mean(x), 2)), row.names = FALSE)
