# ==============================================================================
# experiments/ou_logrr/compare/lib.R -- the candidate (and its sub-variants) against production M3: fits, unit conversion, scoring.
# Load after analysis/00_setup.R:  source("experiments/ou_logrr/compare/lib.R")
# ==============================================================================
source("experiments/ou_logrr/lib.R"); source("experiments/m3plus/lib.R")   # candidate core; M3+ helpers (RTS smoothed M3, simulators)

#' Sub-variants of the candidate (roadmap steps 3-5, plus two ablations of mine).
#'   F1 spec + anchor      tau priors of the sheet (tau_p 3-6, tau_s 60-100 s), V_p/V_s fixed by the production CWT anchor   (step 3)
#'   F2 spec free V        same priors, V_p and V_s free                                                                   (step 4)
#'   F3 effector + anchor  tau_p 0.5-2 s, tau_s 10-20 s, anchor                                                            (step 5)
#'   F4 wide tau_p anchor  tau_p 3-60 s (my addition, motivated by the exercise recovery), anchor
#'   F5 M3 priors anchor   the production tau priors (tau_p 3-6, tau_s 100-300 s) with the candidate's observation model: isolates the observation layer
CAND_VARIANTS <- list(
  F1 = list(tau_prior = "spark", anchor = TRUE), F2 = list(tau_prior = "spark", anchor = FALSE), F3 = list(tau_prior = "effector", anchor = TRUE),
  F4 = list(tau_prior = "spark", bounds = list(p = c(3, 60)), anchor = TRUE), F5 = list(tau_prior = "spark", bounds = list(s = c(100, 300)), anchor = TRUE))

#' Fit one candidate variant to a prepared series (the anchor, when used, is the production one).
fit_cand_variant <- function(series, v, sigma_t = 0.001) {
  pri <- ou_priors(v$tau_prior, sigma_t = sigma_t, lambda_P = if (isTRUE(v$anchor)) series$priors$energy_ratio else NULL, bounds = v$bounds)
  tryCatch(fit_candidate(series$dy, pri, hessian = FALSE), error = function(e) NULL)
}

#' Candidate states in the units of the other model: Hz (x_Hz = nu x, nu = exp(-c): first-order map between ln RR and rate) or ln RR.
cand_states <- function(fit, smooth = FALSE, units = c("hz", "ln")) {
  units <- match.arg(units); k <- if (units == "hz") exp(-fit$c) else 1
  x <- if (smooth) fit$x_smooth else fit$x_filt; P <- if (smooth) fit$P_smooth_c else fit$P_filt_c      # variances include the uncertainty of the profiled level, as M3 does for nu0
  list(xp = k * x[, 1], xs = k * x[, 2], xpv = k^2 * P[, 1], xsv = k^2 * P[, 2])
}

#' M3 (production) states: filtered from the pipeline fit; smoothed from the M3+ RTS (same model). Hz, or converted to ln RR with nu0.
m3_states <- function(series, smooth = FALSE, units = c("hz", "ln")) {
  units <- match.arg(units)
  if (!smooth) { f <- fit_variant(series, model_spec("M3")); if (!isTRUE(f$ok)) return(NULL); k <- if (units == "hz") 1 else 1 / f$nu0
                 return(list(xp = k * f$X_p, xs = k * f$X_s, xpv = k^2 * f$X_p_var, xsv = k^2 * f$X_s_var, fit = f)) }
  f <- fit_ext(series); if (!isTRUE(f$ok)) return(NULL); e <- state_estimates(f, TRUE); k <- if (units == "hz") 1 else 1 / e$r$nu0
  list(xp = k * e$xp, xs = k * e$xs, xpv = k^2 * e$xpv, xsv = k^2 * e$xsv, fit = f)
}

#' Standardised RMSE (RMSE / SD of the truth), correlation and 95% coverage.
score <- function(est, var, truth) { k <- is.finite(est) & is.finite(truth) & is.finite(var)
  c(rmse = sqrt(mean((est[k] - truth[k])^2)) / sd(truth[k]), cor = cor(est[k], truth[k]), cover = mean(abs(est[k] - truth[k]) <= 1.96 * sqrt(pmax(var[k], 0)))) }

#' Share of the true change of a state recovered: mean over the last 15 s of the event minus the mean before it (est / true).
change_ratio <- function(est, truth, t, t0, t1) { pre <- t >= t0 - 120 & t < t0 - 10; en <- t >= t1 - 15 & t < t1
  if (sum(pre) < 5 || sum(en) < 5) return(NA_real_); (mean(est[en]) - mean(est[pre])) / (mean(truth[en]) - mean(truth[pre])) }
