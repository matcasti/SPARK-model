# ==============================================================================
# R/steps/simulation_recovery.R -- recovery of the latent parameters and states under known truth.
# Called by analysis/02_simulation_recovery.R. Each simulated recording is fitted with every
# specification (M0, M1, M2, M3) and compared with the generating values.
# ==============================================================================

# True time scales (s) and lambda_P drawn for each scenario
RECOVERY_SCENARIOS <- list(
  in_prior  = function() list(tau_p = runif(1, 3.5, 5.5), tau_s = runif(1, 120, 250),
                              lambda_p = exp(runif(1, log(1.5), log(4)))),      # inside the prior ranges
  off_prior = function() list(tau_p = runif(1, 8, 12), tau_s = runif(1, 40, 70),
                              lambda_p = exp(runif(1, log(1.5), log(4)))),      # outside them
  vagal_dom = function() list(tau_p = runif(1, 3.5, 5.5), tau_s = runif(1, 120, 250),
                              lambda_p = exp(runif(1, log(6), log(12)))),       # vagal-dominant variance
  symp_dom  = function() list(tau_p = runif(1, 3.5, 5.5), tau_s = runif(1, 120, 250),
                              lambda_p = exp(runif(1, log(0.2), log(0.5))))     # sympathetic-dominant variance
)

#' Errors of one fit relative to the generating values.
recovery_row <- function(fit, truth, simulated, true_share) {
  if (!isTRUE(fit$ok)) return(data.table(spec = fit$spec, ok = FALSE))
  data.table(
    spec = fit$spec, ok = TRUE, conv = fit$conv,
    ln_err_tau_p = log(fit$tau_p / truth$tau_p), ln_err_tau_s = log(fit$tau_s / truth$tau_s),
    ln_err_lp = log(fit$par$lp / truth$lambda_p),
    cor_x_p = cor(fit$X_p, simulated$x_p), cor_x_s = cor(fit$X_s, simulated$x_s),
    rmse_x_p = sqrt(mean((fit$X_p - simulated$x_p)^2)) / sd(simulated$x_p),
    rmse_x_s = sqrt(mean((fit$X_s - simulated$x_s)^2)) / sd(simulated$x_s),
    share_err = fit$share_p - true_share, true_share = true_share, est_share = fit$share_p,
    rmse_rr_ms = fit$rmse_ms, LL_star = fit$LL_star,
    tau_p_true = truth$tau_p, tau_s_true = truth$tau_s, lp_true = truth$lambda_p
  )
}

#' Simulate one recording of `job$scenario` and fit all specifications.
#' @param job list(scenario, sim, seed)
recovery_one_simulation <- function(job, cfg) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  simulated <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s,
                           lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  fits <- fit_specs(prepare_series(simulated$rr))
  true_share <- var(simulated$x_p) / (var(simulated$x_p) + var(simulated$x_s))
  rows <- rbindlist(lapply(fits, recovery_row, truth = truth, simulated = simulated, true_share = true_share))
  rows[, `:=`(scenario = job$scenario, sim = job$sim)]
}

#' Median errors by scenario and specification.
recovery_summarise <- function(table) {
  table[ok == TRUE, .(n = .N,
                      med_abs_ln_err_tau_p = median(abs(ln_err_tau_p)), med_abs_ln_err_tau_s = median(abs(ln_err_tau_s)),
                      med_cor_x_p = median(cor_x_p), med_cor_x_s = median(cor_x_s),
                      med_abs_share_err = median(abs(share_err)), med_rmse_rr_ms = median(rmse_rr_ms)),
        by = .(scenario, spec)]
}
