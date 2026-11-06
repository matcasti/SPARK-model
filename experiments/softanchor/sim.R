# ==============================================================================
# experiments/softanchor/sim.R -- EXPERIMENTAL. Which anchor variant recovers the TRUE latent states best? Same simulated scenarios and seeds as the
# pipeline (02_simulation_recovery.R). Variants: hard anchor (M3), soft anchors (s = 0.15 ... 1.2), no anchor (s = Inf) and an ORACLE hard anchor
# (lambda_P fixed at its true value: the best any anchor could do). Run from the project root:  Rscript experiments/softanchor/sim.R
# ==============================================================================
source("analysis/00_setup.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
spec_free <- drop_term(model_spec("M3"), "anchor")
VARIANTS <- list(hard = 0, soft0.15 = 0.15, soft0.3 = 0.3, soft0.6 = 0.6, soft1.2 = 1.2, free = Inf, oracle = NA)

fit_variant_s <- function(series, s, lp_true) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy
  if (is.na(s)) { pri$energy_ratio <- lp_true; series$priors <- pri; s <- 0 }                      # oracle: the anchor is the truth
  ln_anchor <- log(pri$energy_ratio)
  if (s == 0) { spec <- model_spec("M3"); o <- optimize_spec(series, spec); par <- unpack_theta(o$par, spec, pri$energy_ratio)
  } else {
    cost <- function(th) { v <- variant_objective(th, dy, spec_free, pri, ctl); if (v >= 1e11) return(1e12); if (is.finite(s)) v + 0.5 * ((th[4] - ln_anchor) / s)^2 else v }
    runs <- lapply(.m3_starts(pri, ctl), function(st) tryCatch(optim(c(unname(.par_to_theta(st)), ln_anchor), cost, method = "BFGS", control = list(maxit = ctl$maxit)), error = function(e) NULL))
    vals <- vapply(runs, function(o) if (is.null(o) || !is.finite(o$value)) Inf else o$value, 0); k <- which(vals <= min(vals) + 0.01)[1]
    par <- unpack_theta(runs[[k]]$par, spec_free, pri$energy_ratio)
  }
  f <- filter_eval(dy, par$kp, par$ks, par$lp, par$lR, ctl, full = TRUE)
  list(ok = TRUE, spec = "M3", conv = 0, tau_p = 1 / par$kp, tau_s = 1 / par$ks, par = par, X_p = f$X_p, X_s = f$X_s, X_p_var = f$X_p_var, X_s_var = f$X_s_var,
       share_p = var(f$X_p) / (var(f$X_p) + var(f$X_s)), rmse_ms = 1000 * sqrt(mean((dy - f$implied_rr)^2)), LL_star = f$LL_star)
}
one <- function(job) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  series <- prepare_series(sim$rr); true_share <- var(sim$x_p) / (var(sim$x_p) + var(sim$x_s))
  rbindlist(lapply(names(VARIANTS), function(v) { r <- recovery_row(fit_variant_s(series, VARIANTS[[v]], truth$lambda_p), truth, sim, true_share); r[, `:=`(variant = v, anchor_err = log(series$priors$energy_ratio / truth$lambda_p))] }))[, `:=`(scenario = job$scenario, sim = job$sim)]
}
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Anchor variants on %d simulated recordings x %d variants", length(jobs), length(VARIANTS))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"))
dir.create("experiments/softanchor/_out", showWarnings = FALSE); saveRDS(res, "experiments/softanchor/_out/sim.rds")
res[, variant := factor(variant, names(VARIANTS))]
tab <- res[, .(n = .N, abs_ln_err_lp = round(median(abs(ln_err_lp)), 2), cor_xp = round(median(cor_x_p), 3), cor_xs = round(median(cor_x_s), 3), rmse_xp = round(median(rmse_x_p), 3), rmse_xs = round(median(rmse_x_s), 3),
               abs_ln_err_tau_p = round(median(abs(ln_err_tau_p)), 2), abs_ln_err_tau_s = round(median(abs(ln_err_tau_s)), 2), cover_xp = round(mean(cover_x_p), 2), cover_xs = round(mean(cover_x_s), 2)), by = .(scenario, variant)]
for (sc in names(RECOVERY_SCENARIOS)) { cat("\n==", sc, "==\n"); print(tab[scenario == sc][order(variant)][, !"scenario"], row.names = FALSE) }
cat("\n== All scenarios ==\n"); print(res[, .(abs_ln_err_lp = round(median(abs(ln_err_lp)), 2), rmse_xp = round(mean(rmse_x_p), 3), rmse_xs = round(mean(rmse_x_s), 3), cor_xp = round(mean(cor_x_p), 3), cor_xs = round(mean(cor_x_s), 3)), by = variant][order(variant)], row.names = FALSE)
