# ==============================================================================
# 11_supplement_checks.R -- numerical verification of statements in the supplementary material.
#
# Reads   results/fit_<dataset>.rds   (for the fitted regime; optional)
# Writes  results/supp_checks.rds   bound, bound_mc, fitted_regime, ou_block, cpp_vs_r, calibration
# Checks  (a) bound on the terminal-threshold approximation of the first-passage CDF (+ Monte Carlo)
#         (b) fitted regime of eps = relative timing jitter in the cohorts
#         (c) closed-form OU covariance (incl. Taylor branch) vs Van Loan numerical integration
#         (d) compiled filter vs the R engine
#         (e) calibration of the innovation diagnostics under a correctly specified model
#         The checks are in R/steps/checks.R.
# ==============================================================================
source("analysis/00_setup.R")

# ---- (a) bound on the terminal-threshold approximation ---------------------------------------------------------
bound <- first_passage_bound_table()
bound_mc <- rbindlist(lapply(c(0.1, 0.3), first_passage_monte_carlo))
log_msg("(a) bound: max ratio sup_R / bound = %.3f (must be <= 1)", max(bound$ratio))

# ---- (b) fitted regime ---------------------------------------------------------------------------------------------
fitted_regime <- fitted_regime_table(cfg)
if (nrow(fitted_regime)) {
  log_msg("(b) fitted eps: median %.4f, max %.4f; max CDF bound %.4f",
          median(fitted_regime$eps), max(fitted_regime$eps), max(fitted_regime$bound_cdf))
}

# ---- (c) OU covariance ---------------------------------------------------------------------------------------------------
ou_block <- ou_block_check_table()
log_msg("(c) OU block vs Van Loan: max rel err %.2e (closed form) / %.2e (Taylor branch)",
        max(ou_block[taylor == FALSE, max_rel_err]), max(ou_block[taylor == TRUE, max_rel_err]))

# ---- (d) compiled filter vs engine --------------------------------------------------------------------------------------------
cpp_vs_r <- cpp_vs_r_table(cfg$seed)
if (!is.null(cpp_vs_r)) log_msg("(d) C++ vs R filter: max discrepancy %.1e", max(unlist(cpp_vs_r[, .(d_LL, d_x_p, d_x_s, d_rr)])))

# ---- (e) calibration of the innovation diagnostics ---------------------------------------------------------------------------------
calibration <- calibration_table(cfg)
if (nrow(calibration)) {
  log_msg("(e) calibration (n=%d): KS D > crit in %.0f%% (nominal 5%%); mean PACF violations %.2f (nominal 2.0); mean |z| sd %.3f",
          nrow(calibration), 100 * mean(calibration$KS_Stat > calibration$crit),
          mean(calibration$pacf_viol), mean(calibration$z_sd))
}

# ---- Save ---------------------------------------------------------------------------------------------------------------------------
write_result(list(bound = bound, bound_mc = bound_mc, fitted_regime = fitted_regime, ou_block = ou_block,
                  cpp_vs_r = cpp_vs_r, calibration = calibration), "supp_checks", cfg)
print(bound)
print(bound_mc)
print(ou_block[, .(max_rel_err = max(max_rel_err)), by = taylor])
