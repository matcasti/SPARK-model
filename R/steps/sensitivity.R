# ==============================================================================
# R/steps/sensitivity.R -- how much do the reconstructions depend on the physiological
# hyperparameters? Called by analysis/04_sensitivity.R. Reference model: M3;
# every setting is refitted from the engine's start values.
#
#   sensitivity_hyperparameters()  one hyperparameter changed at a time
#   sensitivity_leave_one_out()    one prior component removed at a time
#   sensitivity_prior_posterior()  where the estimates lie relative to the priors; contraction
#   sensitivity_analyse_recording() chains the three analyses for one recording
# ==============================================================================

# Hyperparameter variations (arguments of spark_priors()); the reference values are the engine defaults
SENSITIVITY_SETTINGS <- list(
  tau_p_bounds = list(narrow = c(3.5, 5), wide = c(2, 9), low = c(1.5, 3), high = c(6, 12)),
  tau_s_bounds = list(narrow = c(150, 250), wide = c(60, 500), low = c(50, 150), high = c(200, 600)),
  topo_target_ratio = list(`10` = 10, `40` = 40),
  topo_prior_weight = list(`5` = 5, `20` = 20),
  jitter_bounds_ms = list(narrow = c(1, 10), wide = c(0.2, 50)),
  crossover_hz = list(`0.08` = 0.08, `0.12` = 0.12),
  sns_leakage_at_boundary = list(`0.02` = 0.02, `0.10` = 0.10),
  lower_lf_vagal_fraction = list(`0.1` = 0.1, `0.5` = 0.5)
)

#' The engine priors with one hyperparameter changed.
sensitivity_priors <- function(param, value) {
  do.call(spark_priors, setNames(list(value), param))
}

#' How a fit differs from the reference fit: log-ratios of the time constants, change in
#' RR error, correlation / normalized RMS difference of the states, share and likelihood.
compare_to_reference <- function(fit, reference) {
  data.table(
    ln_ratio_tau_p = log(fit$tau_p / reference$tau_p), ln_ratio_tau_s = log(fit$tau_s / reference$tau_s),
    ln_ratio_lp = log(fit$par$lp / reference$par$lp), d_rmse_ms = fit$rmse_ms - reference$rmse_ms,
    cor_x_p = suppressWarnings(cor(fit$X_p, reference$X_p)), cor_x_s = suppressWarnings(cor(fit$X_s, reference$X_s)),
    rel_rmsd_x_p = sqrt(mean((fit$X_p - reference$X_p)^2)) / sd(reference$X_p),
    rel_rmsd_x_s = sqrt(mean((fit$X_s - reference$X_s)^2)) / sd(reference$X_s),
    d_share_p = fit$share_p - reference$share_p, d_LL = fit$LL_star - reference$LL_star
  )
}


# ---- The four analyses -------------------------------------------------------------------------------
sensitivity_hyperparameters <- function(rr, reference, spec) {
  rows <- list()
  for (param in names(SENSITIVITY_SETTINGS)) {
    for (label in names(SENSITIVITY_SETTINGS[[param]])) {
      series <- prepare_series(rr, sensitivity_priors(param, SENSITIVITY_SETTINGS[[param]][[label]]))
      fit <- fit_variant(series, spec)
      if (!isTRUE(fit$ok)) next
      rows[[length(rows) + 1L]] <- cbind(
        data.table(kind = "hyperparameter", setting = paste0(param, "=", label), param = param),
        compare_to_reference(fit, reference))
    }
  }
  rows
}

sensitivity_leave_one_out <- function(base_series, reference, spec) {
  rows <- list()
  for (term in c("tau", "topo", "jitter", "anchor")) {
    fit <- fit_variant(base_series, drop_term(spec, term))
    if (!isTRUE(fit$ok)) next
    rows[[length(rows) + 1L]] <- cbind(
      data.table(kind = "leave_one_out", setting = paste0("M3 without ", term), param = term),
      compare_to_reference(fit, reference))
  }
  rows
}

#' Where the M3 estimates lie relative to the priors (standardized distance, inside the 95% interval),
#' posterior contraction of ln(tau) from the Laplace approximation, and the M3 / M0 ratio.
sensitivity_prior_posterior <- function(base_series, reference, spec) {
  priors <- base_series$priors
  curvature <- variant_hessian(reference, base_series, spec)
  data_only <- fit_variant(base_series, model_spec("M0"), full = FALSE)       # no priors
  inside <- function(value, bounds) value >= bounds[1] & value <= bounds[2]
  data.table(
    tau_p = reference$tau_p, tau_s = reference$tau_s,
    z_tau_p = (log(reference$tau_p) - priors$log_tau_p_mu) / priors$sigma_log_tau_p,
    z_tau_s = (log(reference$tau_s) - priors$log_tau_s_mu) / priors$sigma_log_tau_s,
    in_prior_tau_p = inside(reference$tau_p, priors$tau_p_bounds),
    in_prior_tau_s = inside(reference$tau_s, priors$tau_s_bounds),
    contraction_tau_p = 1 - (curvature$sd_ln_tau_p / priors$sigma_log_tau_p)^2,
    contraction_tau_s = 1 - (curvature$sd_ln_tau_s / priors$sigma_log_tau_s)^2,
    ln_ratio_M3_M0_tau_p = if (isTRUE(data_only$ok)) log(reference$tau_p / data_only$tau_p) else NA_real_,
    ln_ratio_M3_M0_tau_s = if (isTRUE(data_only$ok)) log(reference$tau_s / data_only$tau_s) else NA_real_,
    jitter_ms = reference$rr_jitter_ms, lambda_P = reference$par$lp, lambda_R = reference$par$lR
  )
}


# ---- One recording ----------------------------------------------------------------------------------------
#' @param job list(ds = dataset, id = recording id, rr = RR vector)
sensitivity_analyse_recording <- function(job) {
  spec <- model_spec("M3")
  base_series <- prepare_series(job$rr)
  reference <- fit_variant(base_series, spec)
  if (!isTRUE(reference$ok)) stop("reference fit failed")

  rows <- c(sensitivity_hyperparameters(job$rr, reference, spec),
            sensitivity_leave_one_out(base_series, reference, spec))
  label <- function(table) table[, `:=`(Dataset = job$ds, SubjectID = job$id)]

  list(sens = label(rbindlist(rows)),
       prior_post = label(sensitivity_prior_posterior(base_series, reference, spec)))
}

#' Median across recordings of each comparison, by setting.
sensitivity_summarise <- function(sens) {
  sens[, .(n = .N,
           med_abs_ln_tau_p = median(abs(ln_ratio_tau_p)), med_abs_ln_tau_s = median(abs(ln_ratio_tau_s)),
           med_cor_x_p = median(cor_x_p), med_cor_x_s = median(cor_x_s),
           med_rel_rmsd_x_p = median(rel_rmsd_x_p), med_rel_rmsd_x_s = median(rel_rmsd_x_s),
           med_d_rmse_ms = median(d_rmse_ms)),
       by = .(kind, setting)]
}
