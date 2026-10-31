# ==============================================================================
# SPARK engine (M3): phase-domain state-space model of cardiac autonomic
# dynamics.
#
# Model
#   x_p, x_s: independent OU drives observed through the integrated IPFM phase.
#   GLS profiles out nu0 and sigma^2; MAP estimation runs on
#   theta = (log ks, log(kp - ks), log lR).
#   lambda_P = V_p / V_s is fixed before the optimization from the logistic
#   partition of the wavelet energy of the heart rate (spark_priors).
#   Priors: log-normal on tau_p, tau_s and the RR timing jitter; Beta on the
#   separation U = (kp - ks) / (kp + ks).
#
# Usage
#   pri <- spark_priors(); ctl <- spark_control()
#   fit <- fit_model(dy, pri, ctl)
#   report_model(fit); visualize_model(fit); diagnose_model(fit)
#   batch <- batch_process_hrv(rr_list, pri, ctl); report_batch_metrics(batch)
#
# Source one engine file per R session (function names are shared).
# ==============================================================================

library(ggplot2)
library(biwavelet)
library(patchwork)

# Configuration ----------------------------------------------------------------
# spark_priors(): prior hyperparameters. spark_control(): filter, optimizer and
# reporting settings. Both are stored in the fit, so downstream functions reuse
# exactly the settings that produced it.

# Validation helpers ----
#' Assert `x` is an increasing positive interval c(lower, upper).
.check_positive_interval <- function(x, name) {
  if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) || any(x <= 0) || x[1] >= x[2])
    stop(sprintf("`%s` must be a length-2 increasing positive interval, e.g. c(lower, upper).", name), call. = FALSE)
}
#' Assert `x` is a finite scalar within the given open/closed bounds.
.check_scalar <- function(x, name, lower = -Inf, upper = Inf, lower_open = FALSE, upper_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) &&
    (if (lower_open) x > lower else x >= lower) &&
    (if (upper_open) x < upper else x <= upper)
  if (!ok) stop(sprintf("`%s` must be a single finite number in %s%g, %g%s.", name,
                        if (lower_open) "(" else "[", lower, upper, if (upper_open) ")" else "]"), call. = FALSE)
}

#' Prior hyperparameters (single source of truth for all priors).
#'
#' @param tau_p_bounds,tau_s_bounds Central `ci_coverage` interval (s) of the
#'   log-normal priors on the vagal / sympathetic time constants.
#' @param ci_coverage Coverage of the central intervals above and in `jitter_bounds_ms`.
#' @param topo_target_ratio,topo_prior_weight Beta prior on U = (kp - ks) / (kp + ks);
#'   its mode implies kp/ks = `topo_target_ratio`, `topo_prior_weight` is its strength.
#' @param jitter_bounds_ms Central interval (ms) of the log-normal prior on RR timing
#'   jitter (ECG ~1-5 ms, PPG up to ~20 ms); data-independent.
#' @param lf_lower_hz,lf_upper_hz,hf_upper_hz Band edges (Hz); `lf_upper_hz` is the
#'   LF/HF boundary.
#' @param crossover_hz,sns_leakage_at_boundary Logistic SNS/vagal weighting: weights
#'   cross at `crossover_hz`, SNS weight is `sns_leakage_at_boundary` at `lf_upper_hz`.
#' @param lower_lf_vagal_fraction Floor of the vagal weight at low frequencies: the vagal
#'   weight is floor + (1 - floor) * L(f) and the sympathetic weight (1 - floor) * (1 - L(f)),
#'   with L the logistic transition at `crossover_hz`.
#' @param cwt_fs CWT resampling rate (Hz).
#' @return Validated `spark_priors` list with derived constants.
spark_priors <- function(
    tau_p_bounds = c(3, 6),
    tau_s_bounds = c(100, 300),
    ci_coverage  = 0.95,
    topo_target_ratio = 20,
    topo_prior_weight = 10,
    jitter_bounds_ms = c(0.5, 20),
    lf_lower_hz             = 0.04,
    lf_upper_hz             = 0.15,
    hf_upper_hz             = 0.40,
    crossover_hz            = 0.10,
    sns_leakage_at_boundary = 0.05,
    lower_lf_vagal_fraction = 0.25,
    cwt_fs                  = 4
) {
  .check_positive_interval(tau_p_bounds, "tau_p_bounds")
  .check_positive_interval(tau_s_bounds, "tau_s_bounds")
  .check_positive_interval(jitter_bounds_ms, "jitter_bounds_ms")
  .check_scalar(ci_coverage, "ci_coverage", 0, 1, TRUE, TRUE)
  .check_scalar(topo_target_ratio, "topo_target_ratio", 1, lower_open = TRUE)
  .check_scalar(topo_prior_weight, "topo_prior_weight", 0, lower_open = TRUE)
  .check_scalar(sns_leakage_at_boundary, "sns_leakage_at_boundary", 0, 0.5, TRUE, TRUE)
  .check_scalar(lower_lf_vagal_fraction, "lower_lf_vagal_fraction", 0, 1, FALSE, TRUE)
  .check_scalar(cwt_fs, "cwt_fs", 0, lower_open = TRUE)
  if (!(lf_lower_hz > 0 && lf_lower_hz < crossover_hz && crossover_hz < lf_upper_hz && lf_upper_hz < hf_upper_hz))
    stop("Spectral bands must satisfy 0 < lf_lower_hz < crossover_hz < lf_upper_hz < hf_upper_hz.", call. = FALSE)

  # Derived constants.
  # Log-normal from a central interval: mu = mean(log bounds), sigma = log-range / (2 z).
  z_score <- qnorm((ci_coverage + 1.0) / 2.0)

  log_tau_p_mu    <- mean(log(tau_p_bounds))
  sigma_log_tau_p <- diff(log(tau_p_bounds)) / (2 * z_score)
  log_tau_s_mu    <- mean(log(tau_s_bounds))
  sigma_log_tau_s <- diff(log(tau_s_bounds)) / (2 * z_score)

  jitter_bounds_s <- jitter_bounds_ms / 1000
  log_jit_mu      <- mean(log(jitter_bounds_s))
  sigma_log_jit   <- diff(log(jitter_bounds_s)) / (2 * z_score)

  # Beta mode (a - 1) / (a + b - 2) = U_mode, so a - 1 = U_mode * w and b - 1 = (1 - U_mode) * w.
  U_mode        <- (topo_target_ratio - 1.0) / (topo_target_ratio + 1.0)
  alpha_minus_1 <- U_mode * topo_prior_weight
  beta_minus_1  <- (1.0 - U_mode) * topo_prior_weight

  spectral_slope <- -log(sns_leakage_at_boundary / (1 - sns_leakage_at_boundary)) /
    (lf_upper_hz - crossover_hz)                                  # logistic steepness: SNS weight = leakage at lf_upper_hz

  structure(list(
    tau_p_bounds = tau_p_bounds, tau_s_bounds = tau_s_bounds, ci_coverage = ci_coverage,
    topo_target_ratio = topo_target_ratio, topo_prior_weight = topo_prior_weight,
    jitter_bounds_ms = jitter_bounds_ms,
    lf_lower_hz = lf_lower_hz, lf_upper_hz = lf_upper_hz, hf_upper_hz = hf_upper_hz,
    crossover_hz = crossover_hz, sns_leakage_at_boundary = sns_leakage_at_boundary,
    lower_lf_vagal_fraction = lower_lf_vagal_fraction, cwt_fs = cwt_fs,
    ks_mode = 1 / exp(log_tau_s_mu), kp_mode = 1 / exp(log_tau_p_mu),
    log_tau_p_mu = log_tau_p_mu, sigma_log_tau_p = sigma_log_tau_p,
    log_tau_s_mu = log_tau_s_mu, sigma_log_tau_s = sigma_log_tau_s,
    log_jit_mu = log_jit_mu, sigma_log_jit = sigma_log_jit,
    alpha_minus_1 = alpha_minus_1, beta_minus_1 = beta_minus_1,
    spectral_slope = spectral_slope
  ), class = "spark_priors")
}

#' Print the prior specification (incl. data-derived V_p/V_s once available).
print.spark_priors <- function(x, ...) {
  pct <- 100 * x$ci_coverage
  cat("SPARK prior specification (spark_priors)\n")
  cat(sprintf("  tau_p   ~ log-normal, %g%% interval [%g, %g] s\n", pct, x$tau_p_bounds[1], x$tau_p_bounds[2]))
  cat(sprintf("  tau_s   ~ log-normal, %g%% interval [%g, %g] s\n", pct, x$tau_s_bounds[1], x$tau_s_bounds[2]))
  cat(sprintf("  U       ~ Beta (topology): target kp/ks = %g, weight = %g\n", x$topo_target_ratio, x$topo_prior_weight))
  cat(sprintf("  jitter  ~ log-normal, %g%% interval [%g, %g] ms\n", pct, x$jitter_bounds_ms[1], x$jitter_bounds_ms[2]))
  cat(sprintf("  V_p     : CWT prior, LF [%g, %g] Hz, HF <= %g Hz, crossover %g Hz, SNS leakage %g, lower-LF vagal fraction %g, fs %g Hz\n",
              x$lf_lower_hz, x$lf_upper_hz, x$hf_upper_hz, x$crossover_hz, x$sns_leakage_at_boundary,
              x$lower_lf_vagal_fraction, x$cwt_fs))
  if (!is.null(x$energy_ratio))
    cat(sprintf("  data-derived: V_p/V_s = %.5f, nu0_init = %.4f Hz\n", x$energy_ratio, x$nu0_init))
  invisible(x)
}

#' Non-prior settings: artifact multiplier, optimizer, diagnostics.
#'
#' @param jump_threshold,jump_power Artifact multiplier: jitter variance is multiplied by
#'   1 + (J / jump_threshold)^jump_power, J = largest relative RR jump at the beat.
#' @param lik_gate_max Beats with multiplier >= this are excluded from sigma^2 and LL.
#' @param lR_start BFGS start value of lambda_R.
#' @param kp_gap_floor Floor on (kp_mode - ks_mode) for the start value of log(kp - ks).
#' @param maxit BFGS iteration cap.
#' @param diag_gate_max Beats with multiplier >= this are dropped from diagnostics.
#' @param pacf_lag_max Maximum PACF lag.
#' @return Validated `spark_control` list.
spark_control <- function(
    jump_threshold = 0.25,
    jump_power     = 10,
    lik_gate_max   = 100,
    lR_start     = 0.001,
    kp_gap_floor = 0.05,
    maxit        = 1000,
    diag_gate_max      = 1.5,
    pacf_lag_max       = 40
) {
  .check_scalar(jump_threshold, "jump_threshold", 0, lower_open = TRUE)
  .check_scalar(jump_power, "jump_power", 0, lower_open = TRUE)
  .check_scalar(lik_gate_max, "lik_gate_max", 1)
  .check_scalar(lR_start, "lR_start", 0, lower_open = TRUE)
  .check_scalar(kp_gap_floor, "kp_gap_floor", 0, lower_open = TRUE)
  .check_scalar(maxit, "maxit", 1)
  .check_scalar(diag_gate_max, "diag_gate_max", 1)
  .check_scalar(pacf_lag_max, "pacf_lag_max", 1)
  structure(list(
    jump_threshold = jump_threshold, jump_power = jump_power, lik_gate_max = lik_gate_max,
    lR_start = lR_start, kp_gap_floor = kp_gap_floor, maxit = maxit,
    diag_gate_max = diag_gate_max, pacf_lag_max = as.integer(pacf_lag_max)
  ), class = "spark_control")
}

#' Print control settings.
print.spark_control <- function(x, ...) {
  cat("SPARK control settings (spark_control)\n")
  cat(sprintf("  artifact multiplier : jump_threshold = %g, jump_power = %g, lik_gate_max = %g\n",
              x$jump_threshold, x$jump_power, x$lik_gate_max))
  cat(sprintf("  optimizer           : lR_start = %g, kp_gap_floor = %g, maxit = %d\n",
              x$lR_start, x$kp_gap_floor, as.integer(x$maxit)))
  cat(sprintf("  diagnostics         : innovations kept if artifact multiplier < %g; PACF lags = %d\n",
              x$diag_gate_max, x$pacf_lag_max))
  invisible(x)
}

# Spectral prior ----------------------------------------------------------------
#' Data-derived spectral prior for V_p / V_s.
#'
#' Resamples HR (Hz) to `cwt_fs`, computes a Morlet CWT, averages power over time
#' inside the cone of influence, and splits it into vagal / sympathetic energy with
#' logistic frequency weights (the vagal branch keeps a floor fraction at low frequency).
#'
#' @param dy RR intervals (s).
#' @param priors `spark_priors` object.
#' @param verbose Print progress.
#' @return `priors` plus `energy_ratio` (V_p/V_s). `nu0_init`, `p_lf` and `p_hf`
#'   are diagnostics and are not used in the fit.
extract_spectral_priors <- function(dy, priors = spark_priors(), verbose = TRUE) {
  stopifnot(inherits(priors, "spark_priors"))
  if (verbose) cat("Extracting Amplitude, Bandwidth, and MAP priors...\n")
  time_cum <- cumsum(dy)
  hr_hz <- 1 / dy

  fs <- priors$cwt_fs
  t_grid <- seq(min(time_cum), max(time_cum), by = 1/fs)

  # Monotone spline: no ringing between beats
  hr_interp <- splinefun(time_cum, hr_hz, method = "monoH.FC")(t_grid)
  hr_centered <- hr_interp - mean(hr_interp)

  data_mat <- cbind(t_grid, hr_centered)
  invisible(capture.output(wt_res <- biwavelet::wt(data_mat, do.sig = FALSE)))

  # Mask cone-of-influence edge artifacts
  power_masked <- wt_res$power
  for (i in seq_len(ncol(power_masked))) {
    power_masked[wt_res$period > wt_res$coi[i], i] <- NA
  }
  global_power <- rowMeans(power_masked, na.rm = TRUE)
  freqs_raw <- 1 / wt_res$period

  # Ascending frequencies for trapezoidal integration
  ord <- order(freqs_raw)
  freqs <- freqs_raw[ord]
  power <- global_power[ord]

  # Spectral partition: logistic SNS -> vagal crossover (vagal keeps a floor fraction)
  L <- 1 / (1 + exp(-priors$spectral_slope * (freqs - priors$crossover_hz)))
  w_vagal <- priors$lower_lf_vagal_fraction + ((1.0 - priors$lower_lf_vagal_fraction) * L)
  w_symp  <- (1.0 - priors$lower_lf_vagal_fraction) * (1 - L)

  valid_idx <- which(freqs >= priors$lf_lower_hz & freqs <= priors$hf_upper_hz)

  # Trapezoid rule: log-spaced scales give a non-uniform grid
  trapz <- function(x, y) sum(diff(x) * (y[-1] + y[-length(y)]) / 2)

  vagal_energy <- trapz(freqs[valid_idx], (w_vagal * power)[valid_idx])
  symp_energy  <- trapz(freqs[valid_idx], (w_symp * power)[valid_idx])

  vp_vs_ratio <- vagal_energy / symp_energy
  if (verbose) cat("Wavelet-derived V_p/V_s ratio:", round(vp_vs_ratio, 5), "\n")

  # Unweighted LF / HF band power (diagnostic only)
  lf_idx <- which(freqs >= priors$lf_lower_hz & freqs < priors$lf_upper_hz)
  hf_idx <- which(freqs >= priors$lf_upper_hz & freqs <= priors$hf_upper_hz)

  priors$nu0_init     <- mean(hr_hz)
  priors$energy_ratio <- vp_vs_ratio
  priors$p_lf         <- trapz(freqs[lf_idx], power[lf_idx])
  priors$p_hf         <- trapz(freqs[hf_idx], power[hf_idx])
  priors
}

# OU integrals -----------------------------------------------------------------
#' Exact OU integrals over one interval.
#'
#' For dx = -k x dt + sqrt(2kV) dW and I = integral of x, returns the transition and
#' covariance entries the filter needs. Parameterizing by the stationary variance V
#' keeps amplitude independent of kinetics (sigma^2 = 2kV). For k * dt < 1e-3 a
#' Taylor series replaces the closed form to avoid cancellation in 1 - exp(-x).
#'
#' @param k Clearance rate (1/s).
#' @param dt Interval (s).
#' @param V Stationary variance.
#' @return list(phi_int = x -> I transition (independent of V), q11 = Var(x),
#'   q13 = Cov(x, I), q33 = Var(I)).
.compute_ou_block <- function(k, dt, V) {
  x <- k * dt

  # sigma^2 = 2kV fixes the stationary variance at V
  sigma2 <- 2 * k * V

  # Taylor series avoids cancellation in 1 - exp(-x) for small k * dt
  if (x < 1e-3) {
    phi_int <- dt - (k * dt^2)/2 + (k^2 * dt^3)/6
    q11 <- sigma2 * (dt - k * dt^2 + (2 * k^2 * dt^3)/3)
    q13 <- sigma2 * ((dt^2)/2 - (k * dt^3)/2 + (7 * k^2 * dt^4)/24)
    q33 <- sigma2 * ((dt^3)/3 - (k * dt^4)/4 + (7 * k^2 * dt^5)/60)
  } else {
    phi_int <- (1 - exp(-x)) / k
    q11 <- sigma2 * (1 - exp(-2 * x)) / (2 * k)
    q13 <- sigma2 * (1 - 2*exp(-x) + exp(-2*x)) / (2 * k^2)
    q33 <- sigma2 * (2*x - 3 + 4*exp(-x) - exp(-2*x)) / (2 * k^3)
  }

  list(phi_int = phi_int, q11 = q11, q13 = q13, q33 = q33)
}

# Filter ------------------------------------------------------------------------
#' Causal dual-filter GLS.
#'
#' Two Kalman recursions share one gain: X0 filters y = 1 (nu0 = 0) and Xv filters
#' y = dt (unit-nu0 response). nu0 is profiled by GLS, sigma^2 by its concentrated ML
#' estimate, and states are X0 - nu0 * Xv.
#'
#' @param dy RR intervals (s).
#' @param kp,ks Vagal / sympathetic clearance rates (1/s), kp > ks.
#' @param lp Stationary variance of x_p relative to x_s (V_s = 1).
#' @param lR Threshold-jitter variance rate, relative to sigma^2.
#' @param control `spark_control` object.
#' @return list(LL_star, nu0, sig2, states, std_innov, innovations_diagnostic,
#'   time_residuals, implied_rr, gate_multiplier). A non-positive nu0 returns
#'   LL_star = -Inf and NA fields.
run_gls_filter <- function(dy, kp, ks, lp, lR, control = spark_control()) {
  N <- length(dy)
  jump_threshold <- control$jump_threshold
  jump_power <- control$jump_power
  lik_gate_max <- control$lik_gate_max

  # Artifact multiplier: inflate jitter variance at large relative RR jumps
  dy_safe <- pmax(dy[-N], 1e-6)
  jump_back <- c(0, abs(diff(dy)) / dy_safe)
  jump_fwd  <- c(abs(diff(dy)) / dy_safe, 0)
  J_k <- pmax(jump_back, jump_fwd, na.rm = TRUE)
  gate_multiplier <- pmin(1 + (J_k / jump_threshold)^jump_power, 1e10)
  gate_multiplier[!is.finite(gate_multiplier)] <- 1e10

  # One causal pass over the beats
  .gls_pass <- function() {
    X0 <- matrix(0, 4, 1); Xv <- matrix(0, 4, 1)
    P  <- diag(c(lp, 1, 0, 0))                 # stationary init: V_p = lp, V_s = 1, phases = 0
    H  <- matrix(c(0, 0, -1, 1), 1, 4)

    S_vec <- numeric(N); v0_vec <- numeric(N); vv_vec <- numeric(N)
    S_plain_vec <- numeric(N)
    X0_store <- matrix(0, N, 4); Xv_store <- matrix(0, N, 4)
    P_store  <- array(0, dim = c(4, 4, N))
    Phi <- matrix(0, 4, 4); Q <- matrix(0, 4, 4)

    for (k in 1:N) {
      dt <- dy[k]
      dynamic_lR <- min(lR * gate_multiplier[k], 1e12, na.rm = TRUE)

      blk_s <- .compute_ou_block(ks, dt, 1.0)
      blk_p <- .compute_ou_block(kp, dt, lp)

      Phi[1,1] <- exp(-kp * dt); Phi[3,1] <- .compute_ou_block(kp, dt, 1.0)$phi_int
      Phi[2,2] <- exp(-ks * dt); Phi[4,2] <- blk_s$phi_int

      X0_pred <- Phi %*% X0; Xv_pred <- Phi %*% Xv

      Q[1,1] <- blk_p$q11; Q[1,3] <- blk_p$q13; Q[3,1] <- blk_p$q13; Q[3,3] <- blk_p$q33
      Q[2,2] <- blk_s$q11; Q[2,4] <- blk_s$q13; Q[4,2] <- blk_s$q13; Q[4,4] <- blk_s$q33

      P_pred <- Phi %*% P %*% t(Phi) + Q

      # Innovation variance with artifact-inflated jitter
      S <- as.numeric(H %*% P_pred %*% t(H) + dynamic_lR * dt)
      S <- max(S, 1e-12)
      K <- P_pred %*% t(H) / S

      # Same without the artifact multiplier; standardizes the diagnostics
      S_true <- as.numeric(H %*% P_pred %*% t(H) + lR * dt)
      S_plain_vec[k] <- max(S_true, 1e-12)

      v0 <- 1.0 - as.numeric(H %*% X0_pred)
      vv <- dt  - as.numeric(H %*% Xv_pred)

      X0 <- X0_pred + K * v0; Xv <- Xv_pred + K * vv
      I_KH <- diag(4) - K %*% H

      P <- I_KH %*% P_pred %*% t(I_KH) + K %*% (dynamic_lR * dt) %*% t(K)
      P <- (P + t(P)) / 2

      X0_store[k, ] <- as.numeric(X0); Xv_store[k, ] <- as.numeric(Xv)
      P_store[,,k]  <- P

      # IPFM reset of phase states at the beat boundary
      X0[3:4,] <- 0; Xv[3:4,] <- 0
      P[3:4, ] <- 0; P[, 3:4] <- 0

      S_vec[k] <- S; v0_vec[k] <- v0; vv_vec[k] <- vv
    }

    den_nu  <- sum((vv_vec^2) / S_vec, na.rm = TRUE)
    nu0_hat <- if (is.finite(den_nu) && den_nu > 1e-12) sum(v0_vec * vv_vec / S_vec, na.rm = TRUE) / den_nu else 0

    list(X0_store = X0_store, Xv_store = Xv_store, P_store = P_store,
         S_vec = S_vec, S_plain_vec = S_plain_vec,
         v0_vec = v0_vec, vv_vec = vv_vec, nu0_hat = nu0_hat)
  }

  res <- .gls_pass()
  res$gate_multiplier <- gate_multiplier

  # Guard: non-positive nu0 would give infinite implied RR
  if (!is.finite(res$nu0_hat) || res$nu0_hat <= 0) {
    na_states <- matrix(NA_real_, N, 7)
    colnames(na_states) <- c("X_p", "X_s", "Phase_p", "Phase_s", "X_p_var", "X_s_var", "rr_var")
    return(list(LL_star = -Inf, nu0 = NA_real_, sig2 = NA_real_,
                states = na_states, std_innov = rep(NA_real_, N),
                innovations_diagnostic = rep(NA_real_, N),
                time_residuals = rep(NA_real_, N), implied_rr = rep(NA_real_, N),
                gate_multiplier = gate_multiplier))
  }

  # Outputs
  v_final  <- res$v0_vec - res$nu0_hat * res$vv_vec
  ok       <- res$gate_multiplier < lik_gate_max
  if (sum(ok) < 2) ok <- rep(TRUE, N)
  n_eff    <- sum(ok)
  sig2_hat <- max(sum(v_final[ok]^2 / res$S_vec[ok]) / n_eff, 1e-12)
  LL_star  <- -0.5 * sum(log(res$S_vec[ok])) - (n_eff / 2) * log(sig2_hat)

  states_final <- (res$X0_store - res$nu0_hat * res$Xv_store)

  # Var(nu0) = sigma^2 / den_nu propagates into the state variances
  den_nu  <- sum((res$vv_vec^2) / res$S_vec, na.rm = TRUE)
  inv_den <- if (is.finite(den_nu) && den_nu > 1e-12) 1 / den_nu else NA_real_

  X_p_var_fwd <- (res$P_store[1, 1, ] + (res$Xv_store[, 1]^2) * inv_den) * sig2_hat
  X_s_var_fwd <- (res$P_store[2, 2, ] + (res$Xv_store[, 2]^2) * inv_den) * sig2_hat
  rr_var_fwd  <- (res$P_store[3, 3, ] + res$P_store[4, 4, ] - 2 * res$P_store[3, 4, ]) *
    sig2_hat / (res$nu0_hat^2)

  states_fwd_out <- cbind(states_final, X_p_var_fwd, X_s_var_fwd, rr_var_fwd)
  colnames(states_fwd_out) <- c("X_p", "X_s", "Phase_p", "Phase_s", "X_p_var", "X_s_var", "rr_var")

  std_innov            <- v_final / sqrt(res$S_vec * sig2_hat)
  std_innov_diagnostic <- v_final / sqrt(res$S_plain_vec * sig2_hat)
  implied_rr           <- (1.0 + states_final[,3] - states_final[,4]) / res$nu0_hat
  time_residuals       <- dy - implied_rr

  list(LL_star = LL_star, nu0 = res$nu0_hat, sig2 = sig2_hat,
       states = states_fwd_out, std_innov = std_innov,
       innovations_diagnostic = std_innov_diagnostic,
       time_residuals = time_residuals, implied_rr = implied_rr,
       gate_multiplier = res$gate_multiplier)
}


# MAP objective -----------------------------------------------------------------
#' Map unconstrained theta to physical (ks, kp, lR); exp links enforce
#' positivity and kp > ks.
.theta_to_par <- function(theta) {
  ks <- exp(theta[1])
  c(ks = ks, kp = ks + exp(theta[2]), lR = exp(theta[3]))
}
#' Inverse of .theta_to_par.
.par_to_theta <- function(par) {
  c(log(par[["ks"]]), log(par[["kp"]] - par[["ks"]]), log(par[["lR"]]))
}

#' Log-prior terms (up to constants) at a filter result `res`.
#'
#' Log-normal on tau_s, tau_p and RR jitter; Beta on topology.
#' @return Named list with one term per prior.
.log_prior_terms <- function(par, res, dy, priors) {
  ks <- par[["ks"]]; kp <- par[["kp"]]; lR <- par[["lR"]]

  tau_s <- 1.0 / ks
  tau_p <- 1.0 / kp
  tau_s_term <- -0.5 * ((log(tau_s) - priors$log_tau_s_mu) / priors$sigma_log_tau_s)^2
  tau_p_term <- -0.5 * ((log(tau_p) - priors$log_tau_p_mu) / priors$sigma_log_tau_p)^2

  U     <- (kp - ks) / (kp + ks)
  One_U <- (2.0 * ks) / (kp + ks)
  topology_term <- (priors$alpha_minus_1 * log(U)) + (priors$beta_minus_1 * log(One_U))

  # RR timing jitter (s) at the median RR; log-normal prior in physical units
  s_rr  <- sqrt(res$sig2 * lR * median(dy)) / res$nu0
  lR_term <- -0.5 * ((log(s_rr) - priors$log_jit_mu) / priors$sigma_log_jit)^2

  list(tau_s = tau_s_term, tau_p = tau_p_term, topology = topology_term, lR = lR_term)
}

#' Negative log-posterior in theta; non-finite values return 1e12 so BFGS backs off.
objective_map <- function(theta, dy, priors, control) {
  par <- .theta_to_par(theta)
  res <- run_gls_filter(dy, kp = par[["kp"]], ks = par[["ks"]], lp = priors$energy_ratio,
                        lR = par[["lR"]], control = control)
  pen <- .log_prior_terms(par, res, dy, priors)

  val <- -(res$LL_star + pen$tau_s + pen$tau_p + pen$topology + pen$lR)
  if (!is.finite(val)) return(1e12)
  return(val)
}

# Fit ---------------------------------------------------------------------------
#' Fit SPARK by MAP (BFGS, single start).
#'
#' @param dy RR intervals (s); NAs are dropped.
#' @param priors `spark_priors` object (V_p/V_s is added from the CWT).
#' @param control `spark_control` object.
#' @param verbose Print progress.
#' @return `spark_fit`: params, states, implied_rr, innovations, priors and
#'   control (everything report_model, visualize_model and diagnose_model need).
fit_model <- function(dy, priors = spark_priors(), control = spark_control(), verbose = TRUE) {
  stopifnot(inherits(priors, "spark_priors"), inherits(control, "spark_control"))
  dy <- as.numeric(dy[!is.na(dy)])
  priors <- extract_spectral_priors(dy, priors, verbose = verbose)

  theta0 <- c(log(priors$ks_mode),
              log(max(priors$kp_mode - priors$ks_mode, control$kp_gap_floor)),
              log(control$lR_start))

  if (verbose) cat("\nRunning 3D BFGS (M3) from a single start...\n")
  opt_err <- NULL
  opt_res <- tryCatch(
    optim(par = theta0, fn = function(th) objective_map(th, dy, priors, control),
          method = "BFGS", control = list(maxit = control$maxit)),
    error = function(e) { opt_err <<- conditionMessage(e); NULL }
  )
  if (is.null(opt_res) || !is.finite(opt_res$value))
    stop("Optimization failed", if (!is.null(opt_err)) paste0(": ", opt_err) else ".", call. = FALSE)
  if (verbose) cat(sprintf("  objective %.3f (conv %d)\n", opt_res$value, opt_res$convergence))

  par <- .theta_to_par(opt_res$par)
  Vp  <- priors$energy_ratio
  final_res <- run_gls_filter(dy, kp = par[["kp"]], ks = par[["ks"]], lp = Vp,
                              lR = par[["lR"]], control = control)

  structure(list(
    dy = dy, time = cumsum(dy), model_type = "M3",
    priors = priors, control = control,
    params = list(
      ks = par[["ks"]], kp = par[["kp"]], lR = par[["lR"]], lp = Vp,
      nu0 = final_res$nu0, sig2 = final_res$sig2,
      sig_p = sqrt(2 * par[["kp"]] * Vp * final_res$sig2),
      sig_s = sqrt(2 * par[["ks"]] * 1.0 * final_res$sig2),
      rr_jitter_ms = 1000 * sqrt(final_res$sig2 * par[["lR"]] * median(dy)) / final_res$nu0
    ),
    LL_star = final_res$LL_star,
    states = final_res$states,
    implied_rr = final_res$implied_rr,
    innovations = final_res$std_innov,
    innovations_diagnostic = final_res$innovations_diagnostic,
    time_residuals = final_res$time_residuals,
    gate_multiplier = final_res$gate_multiplier,
    opt_raw = opt_res
  ), class = "spark_fit")
}

#' Assert `fit_obj` is a `spark_fit`.
.check_fit <- function(fit_obj) {
  if (!inherits(fit_obj, "spark_fit"))
    stop("Expected a 'spark_fit' object returned by fit_model().", call. = FALSE)
  invisible(fit_obj)
}

# Report ------------------------------------------------------------------------
#' Parameter table with delta-method standard errors.
#'
#' SEs come from the numerical Hessian of the MAP objective (pseudo-inverse tolerates
#' singularity). Profiled, fixed and derived parameters have no SE.
#' @param fit_obj `spark_fit`.
#' @return The parameter table (invisibly).
report_model <- function(fit_obj) {
  .check_fit(fit_obj)
  if (!requireNamespace("numDeriv", quietly = TRUE)) stop("Package 'numDeriv' is required but not installed. Please install it manually.")
  if (!requireNamespace("MASS", quietly = TRUE))     stop("Package 'MASS' is required but not installed. Please install it manually.")

  dy  <- fit_obj$dy
  p   <- fit_obj$params
  pri <- fit_obj$priors
  ctl <- fit_obj$control

  theta_opt <- .par_to_theta(c(ks = p$ks, kp = p$kp, lR = p$lR))

  cat("\nCalculating standard errors via exact Delta Method...\n")
  hess_precise <- numDeriv::hessian(func = function(th) objective_map(th, dy, pri, ctl), x = theta_opt)
  inv_hess <- tryCatch(MASS::ginv(hess_precise), error = function(e) matrix(NA, 3, 3))

  # Jacobian: theta -> physical (ks, kp, lR)
  J <- matrix(0, 3, 3)
  J[1, 1] <- p$ks
  J[2, 1] <- p$ks; J[2, 2] <- p$kp - p$ks
  J[3, 3] <- p$lR
  cov_physical <- J %*% inv_hess %*% t(J)

  se_ks <- sqrt(max(cov_physical[1, 1], 0)); se_kp <- sqrt(max(cov_physical[2, 2], 0))
  se_lR <- sqrt(max(cov_physical[3, 3], 0))
  # Var(1/k) ~ Var(k) / k^4 (delta method)
  se_taus <- sqrt(max(cov_physical[1, 1] / (p$ks^4), 0))
  se_taup <- sqrt(max(cov_physical[2, 2] / (p$kp^4), 0))

  params_list <- list(
    list(name = "Nu_0", group = "Baseline", interp = "Intrinsic SA node pacing rate absent autonomic tone (Hz)", status = "Profiled (GLS)", est = p$nu0, se = NA),
    list(name = "Sigma_sys^2", group = "Scale", interp = "Total stochastic system variance multiplier", status = "Profiled (GLS)", est = p$sig2, se = NA),
    list(name = "Kappa_S", group = "Kinetics", interp = "Clearance rate of sympathetic neurotransmitters (Hz)", status = "Estimated (MAP)", est = p$ks, se = se_ks),
    list(name = "Kappa_P", group = "Kinetics", interp = "Clearance rate of vagal neurotransmitters (Hz)", status = "Estimated (MAP)", est = p$kp, se = se_kp),
    list(name = "Tau_S", group = "Kinetics", interp = "E-folding time 1/kappa of sympathetic response (s)", status = "Derived (MAP)", est = 1/p$ks, se = se_taus),
    list(name = "Tau_P", group = "Kinetics", interp = "E-folding time 1/kappa of parasympathetic response (s)", status = "Derived (MAP)", est = 1/p$kp, se = se_taup),
    list(name = "Sigma_S", group = "Volatility", interp = "Absolute amplitude of sympathetic neural drive", status = "Profiled (GLS)", est = p$sig_s, se = NA),
    list(name = "Sigma_P", group = "Volatility", interp = "Amplitude of parasympathetic neural drive", status = "Derived (Spectral)", est = p$sig_p, se = NA),
    list(name = "Lambda_P", group = "Ratios", interp = "Vagal-to-Sympathetic energy balance (Parseval anchor)", status = "Profiled (Spectral)", est = p$lp, se = NA),
    list(name = "Lambda_R", group = "Ratios", interp = "Noise-to-signal ratio from SA node threshold jitter", status = "Estimated (MAP)", est = p$lR, se = se_lR),
    list(name = "Jitter_RR", group = "Ratios", interp = "Equivalent RR timing jitter SD at median RR (ms); penalized by log-normal prior", status = "Derived (MAP)", est = p$rr_jitter_ms, se = NA),
    list(name = "Theta_Jump", group = "Robustness", interp = "Relative fractional jump threshold for outlier rejection (spark_control)", status = "Fixed (Hyper)", est = ctl$jump_threshold, se = NA),
    list(name = "Power_Jump", group = "Robustness", interp = "Polynomial steepness of asymptotic null-gain gate (spark_control)", status = "Fixed (Hyper)", est = ctl$jump_power, se = NA)
  )

  df <- do.call(rbind, lapply(params_list, function(x) {
    if (is.na(x$se)) { ci_str <- "-"; se_str <- "-" } else {
      ci_str <- sprintf("[%.4f, %.4f]", x$est - 1.96 * x$se, x$est + 1.96 * x$se)
      se_str <- sprintf("%.5f", x$se)
    }
    data.frame(Parameter = x$name, Estimate = sprintf("%.5f", x$est), SE = se_str, `CI_95%` = ci_str,
               Status = x$status, Interpretation = x$interp, stringsAsFactors = FALSE, check.names = FALSE)
  }))

  cat("\n=========================================================================================================\n")
  cat("                   HRV PHASE-DOMAIN SDE-IPFM FITTING REPORT  (M3)\n")
  cat("=========================================================================================================\n")
  print(df, row.names = FALSE, right = FALSE)
  cat("=========================================================================================================\n")
  cat(sprintf("* Log-likelihood (profiled) = %.2f\n", fit_obj$LL_star))
  cat("* Note: 'Profiled' parameters are solved exactly (GLS/concentrated likelihood) or by deterministic wavelet\n")
  cat("  mapping, so they carry no numerical MLE standard errors.\n\n")
  print(pri)
  cat("\n")
  invisible(df)
}

# Visualization -----------------------------------------------------------------
#' Filtered RR with 95% CI, latent drives and phase-space energy map.
#' @param fit_obj `spark_fit`.
#' @return patchwork object.
visualize_model <- function(fit_obj) {
  .check_fit(fit_obj)
  dy <- fit_obj$dy; N <- length(dy); time_cum <- fit_obj$time; p <- fit_obj$params
  ctl <- fit_obj$control

  X_p <- fit_obj$states[, "X_p"]; X_s <- fit_obj$states[, "X_s"]
  X_p_se <- sqrt(pmax(fit_obj$states[, "X_p_var"], 0)); X_s_se <- sqrt(pmax(fit_obj$states[, "X_s_var"], 0))
  rr_se  <- sqrt(pmax(fit_obj$states[, "rr_var"], 0))

  implied_rr <- fit_obj$implied_rr
  df_fit <- data.frame(Time = time_cum, Observed = dy, Implied = implied_rr,
                       Lower = implied_rr - qnorm(0.975) * rr_se, Upper = implied_rr + qnorm(0.975) * rr_se)

  pA <- ggplot(df_fit, aes(x = Time)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "darkred", alpha = 0.15) +
    geom_line(aes(y = Observed, color = "Observed RR"), linewidth = 1/3, alpha = 0.4) +
    geom_line(aes(y = Implied, color = "Filtered RR"), linewidth = 0.5) +
    scale_color_manual(values = c("Observed RR" = "gray50", "Filtered RR" = "darkred")) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "Phase-Domain Filtering & Predictive Fit (with 95% CI)", y = "RR Interval (s)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  df_X <- data.frame(Time = rep(time_cum, 2), Value = c(X_p, X_s), SE = c(X_p_se, X_s_se),
                     State = factor(rep(c("Parasympathetic Drive", "Sympathetic Drive"), each = N),
                                    levels = c("Parasympathetic Drive", "Sympathetic Drive")))
  df_X$Lower <- df_X$Value - qnorm(0.975) * df_X$SE; df_X$Upper <- df_X$Value + qnorm(0.975) * df_X$SE

  pB <- ggplot(df_X, aes(x = Time, y = Value, color = State, fill = State)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.3) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = c("Parasympathetic Drive" = "#4DB0D0", "Sympathetic Drive" = "#DC5050")) +
    scale_fill_manual(values = c("Parasympathetic Drive" = "#4DB0D0", "Sympathetic Drive" = "#DC5050")) +
    scale_x_continuous(expand = c(0,0)) +
    labs(title = "Native Autonomic Drivers (95% CI)", y = "Amplitude (Hz)", x = "") +
    theme_classic() + theme(legend.title = element_blank(), legend.position = "top")

  # Energy map: U = 0.5 x' Sigma^-1 x with the stationary covariance; clipped at the 95th percentile
  Sigma_stat <- diag(c(p$sig_p^2 / (2 * p$kp), p$sig_s^2 / (2 * p$ks)))
  inv_Sigma  <- solve(Sigma_stat)
  p_range <- range(X_p); s_range <- range(X_s)
  margin_p <- diff(p_range) * 0.2; if (margin_p == 0) margin_p <- 0.01
  margin_s <- diff(s_range) * 0.2; if (margin_s == 0) margin_s <- 0.01
  grid <- expand.grid(X_p = seq(p_range[1] - margin_p, p_range[2] + margin_p, length.out = 120),
                      X_s = seq(s_range[1] - margin_s, s_range[2] + margin_s, length.out = 120))
  grid$Energy <- apply(grid, 1, function(v) { 0.5 * as.numeric(t(v) %*% inv_Sigma %*% v) })
  grid$Energy <- pmin(grid$Energy, quantile(grid$Energy, 0.95))
  df_traj <- data.frame(X_p = X_p, X_s = X_s, Time = time_cum)

  pC <- ggplot() +
    geom_raster(data = grid, aes(x = X_p, y = X_s, fill = Energy)) +
    geom_contour(data = grid, aes(x = X_p, y = X_s, z = Energy), color = "white", alpha = 0.15, bins = 15) +
    geom_path(data = df_traj, aes(x = X_p, y = X_s, color = Time, alpha = Time),
              arrow = ggplot2::arrow(type = "closed", length = unit(0.05, "inches")), show.legend = FALSE) +
    scale_fill_viridis_c(option = "mako", name = "Energy (U)") +
    scale_color_viridis_c(option = "plasma", guide = "none") +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    labs(title = "Autonomic Phase Space Topology", x = "Parasympathetic Drive (Hz)", y = "Sympathetic Drive (Hz)") +
    theme_classic() + theme(legend.position = "top")

  ((pA / pB) | pC)
}


# Diagnostics -------------------------------------------------------------------
#' Single definition of the innovation diagnostics (used by plots and the batch table).
#'
#' Keeps beats with finite innovations and artifact multiplier < `diag_gate_max`.
#' KS tests U = pnorm(z) against Uniform(0, 1) (time-rescaling); PACF violations are
#' lags outside +/- 1.96 / sqrt(N). Innovations are standardized with S_ungated
#' (no artifact multiplier).
#' @param fit_obj `spark_fit`.
#' @return list(z, U, keep, N, N_total, ks, ks_bound, pacf, pacf_band, pacf_violations).
.innovation_diagnostics <- function(fit_obj) {
  .check_fit(fit_obj)
  ctl   <- fit_obj$control
  z_all <- fit_obj$innovations_diagnostic
  keep  <- which(is.finite(z_all) & fit_obj$gate_multiplier < ctl$diag_gate_max)
  z     <- z_all[keep]
  N     <- length(z)
  if (N < 10) stop("Fewer than 10 valid innovations available for diagnostics.", call. = FALSE)

  U_k <- pnorm(z)
  ks_res <- ks.test(U_k, "punif")

  lag_max   <- min(ctl$pacf_lag_max, max(1, N - 5))
  pacf_vals <- as.numeric(pacf(z, lag.max = lag_max, plot = FALSE)$acf)
  pacf_band <- qnorm(0.975) / sqrt(N)

  list(z = z, U = U_k, keep = keep, N = N, N_total = length(z_all),
       ks = ks_res, ks_bound = 1.36 / sqrt(N),
       pacf = pacf_vals, pacf_band = pacf_band,
       pacf_violations = sum(abs(pacf_vals) > pacf_band))
}

#' Q-Q plot, KS plot, PACF and histogram of the standardized innovations.
#' @param fit_obj `spark_fit`.
#' @return patchwork object.
diagnose_model <- function(fit_obj) {
  dg <- .innovation_diagnostics(fit_obj)
  z <- dg$z; U_k <- dg$U

  print(dg$ks)

  df_qq <- data.frame(z = z)
  pA <- ggplot(df_qq, aes(sample = z)) +
    stat_qq(color = "darkred", pch = 16) +
    stat_qq_line(linetype = "dashed") +
    labs(title = "Q-Q Plot of Phase Innovations",
         x = "Theoretical Normal",
         y = "Empirical") +
    theme_classic()

  empirical_cdf <- ecdf(U_k)
  df_cdf <- data.frame(U = sort(U_k), CDF = empirical_cdf(sort(U_k)))

  pB <- ggplot(df_cdf, aes(x = U, y = CDF)) +
    geom_step(color = "darkred", linewidth = 1) +
    geom_abline(slope = 1, intercept = 0) +
    geom_abline(slope = 1, intercept = dg$ks_bound, color = "gray") +
    geom_abline(slope = 1, intercept = -dg$ks_bound, color = "gray") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
    labs(title = "Time-Rescaling KS-Plot",
         x = "Theoretical Uniform(0,1)",
         y = "Empirical CDF") +
    theme_classic()

  df_acf <- data.frame(lag = seq_along(dg$pacf), acf = dg$pacf)

  pC <- ggplot(df_acf, aes(x = lag, y = acf)) +
    geom_segment(aes(xend = lag, yend = 0),
                 linewidth = 1) +
    geom_hline(yintercept = c(-dg$pacf_band, dg$pacf_band),
               color = "gray") +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits = c(NA, 1)) +
    labs(title = "Autocorrelation of Innovations",
         x = "Lag (Heartbeats)",
         y = "Partial Autocorrelation Function") +
    theme_classic()

  pD <- ggplot(df_qq, aes(x = z)) +
    geom_histogram(aes(y = after_stat(density)),
                   fill = "#7E6148",
                   alpha = 0.7, bins = 30, color = "white") +
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = 1),
                  color = "black",
                  linetype = "dashed") +
    labs(title = "Exact Gaussian Innovations",
         x = "Standardized Residuals",
         y = "Density") +
    theme_classic()

  (pA / pB / pC / pD) +
    plot_annotation(caption = sprintf("%d of %d beats retained (artifact multiplier < %g); KS D = %.3f, PACF violations = %d",
                                      dg$N, dg$N_total, fit_obj$control$diag_gate_max,
                                      as.numeric(dg$ks$statistic), dg$pacf_violations))
}

# Batch processing --------------------------------------------------------------
#' One fit -> summary row, per-beat state table and padded PACF vector.
.summarise_fit <- function(fit_res, s_name, dataset_tag) {
  p   <- fit_res$params
  ctl <- fit_res$control
  dy  <- fit_res$dy                              # NA-cleaned series actually fitted
  implied_rr <- fit_res$implied_rr
  dg  <- .innovation_diagnostics(fit_res)

  rmse_ms  <- sqrt(mean((dy - implied_rr)^2)) * 1000
  mape_pct <- mean(abs((dy - implied_rr) / dy)) * 100

  pacf_vals <- dg$pacf
  if (length(pacf_vals) < ctl$pacf_lag_max) pacf_vals <- c(pacf_vals, rep(NA, ctl$pacf_lag_max - length(pacf_vals)))

  N_beats <- length(dy)
  summary_row <- data.frame(
    Dataset = dataset_tag, SubjectID = s_name,
    Converged = (fit_res$opt_raw$convergence == 0), Iterations = fit_res$opt_raw$counts[1],
    NumBeats = N_beats, NumValidBeats = dg$N, Nu_0 = p$nu0, Sigma_sys2 = p$sig2,
    Kappa_S = p$ks, Kappa_P = p$kp, Tau_S = 1 / p$ks, Tau_P = 1 / p$kp, Lambda_R = p$lR,
    LL_star = fit_res$LL_star,
    RMSE_ms = rmse_ms, MAPE_pct = mape_pct,
    KS_Stat = as.numeric(dg$ks$statistic), KS_p_value = as.numeric(dg$ks$p.value),
    PACF_Violations = dg$pacf_violations, stringsAsFactors = FALSE)

  states_row <- data.frame(
    Dataset = dataset_tag, SubjectID = s_name, Beat = seq_len(N_beats), Time = fit_res$time,
    RR_Interval = dy, X_p = fit_res$states[, "X_p"], X_s = fit_res$states[, "X_s"],
    Phase_p = fit_res$states[, "Phase_p"], Phase_s = fit_res$states[, "Phase_s"],
    Innovation = fit_res$innovations,
    Valid = seq_len(N_beats) %in% dg$keep, stringsAsFactors = FALSE)

  list(summary = summary_row, states = states_row, pacf = pacf_vals)
}

#' Fit a list of RR series; failed recordings are skipped with a message.
#' @param rr_list Named list of RR vectors (s).
#' @param priors,control `spark_priors` / `spark_control` objects.
#' @param dataset_tag Label stored in the output tables.
#' @return list(summary, states, pacf_matrix, models, priors, control).
batch_process_hrv <- function(rr_list, priors = spark_priors(), control = spark_control(),
                              dataset_tag = "Unknown") {
  stopifnot(inherits(priors, "spark_priors"), inherits(control, "spark_control"))
  N_batch <- length(rr_list)
  subj_names <- names(rr_list)
  if (is.null(subj_names)) subj_names <- paste0("Subject_", seq_len(N_batch))

  summary_rows <- list(); states_list <- list(); raw_fits <- list(); pacf_list <- list()
  cat(sprintf("\nInitializing Batch Processing for Dataset: %s (%d recordings)\n", dataset_tag, N_batch))

  for (i in seq_len(N_batch)) {
    s_name <- subj_names[i]
    cat(sprintf(" Processing [%d/%d]: %s... ", i, N_batch, s_name))

    out <- tryCatch({
      fit_res <- suppressMessages(fit_model(rr_list[[i]], priors = priors, control = control, verbose = FALSE))
      list(fit = fit_res, sm = .summarise_fit(fit_res, s_name, dataset_tag))
    }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); NULL })
    if (is.null(out)) next
    cat("SUCCESS\n")

    raw_fits[[s_name]]     <- out$fit
    summary_rows[[s_name]] <- out$sm$summary
    states_list[[s_name]]  <- out$sm$states
    pacf_list[[s_name]]    <- out$sm$pacf
  }
  if (length(summary_rows) == 0) stop("All recordings failed to fit.", call. = FALSE)

  summary_df <- do.call(rbind, summary_rows); states_df <- do.call(rbind, states_list)
  pacf_mat   <- do.call(rbind, pacf_list)
  rownames(summary_df) <- NULL; rownames(states_df) <- NULL
  list(summary = summary_df, states = states_df, pacf_matrix = pacf_mat, models = raw_fits,
       priors = priors, control = control)
}

#' Print cohort-level convergence, goodness-of-fit and parameter summaries.
#' @param batch_obj Output of batch_process_hrv().
report_batch_metrics <- function(batch_obj) {
  df <- batch_obj$summary; N <- nrow(df); ctl <- batch_obj$control
  cat("\n=========================================================================\n")
  cat(sprintf("               BATCH VALIDATION REPORT (Dataset: %s)\n", df$Dataset[1]))
  cat("=========================================================================\n")
  cat(sprintf("Total Managed Series      : %d\n", N))
  cat(sprintf("Optimization Convergence  : %.1f%%\n", (sum(df$Converged) / N) * 100))
  cat(sprintf("Average Sequence Length   : %.1f beats (%.1f valid for diagnostics)\n", mean(df$NumBeats), mean(df$NumValidBeats)))
  cat("\n--- Statistical Diagnostics (Point-Process Validity) ---\n")
  cat(sprintf("Time-Rescaling KS Pass Rate (p > 0.05)   : %.1f%% (Mean Stat: %.4f)\n", (sum(df$KS_p_value > 0.05) / N) * 100, mean(df$KS_Stat)))
  cat(sprintf("Mean Whiteness Violations (out of %d)    : %.2f lags\n", ctl$pacf_lag_max, mean(df$PACF_Violations)))
  cat("\n--- Kinetic Parameter Estimations (Mean +/- SD) ---\n")
  cat(sprintf("Baseline Pacing (Nu_0)         : %.3f +/- %.3f Hz\n", mean(df$Nu_0), sd(df$Nu_0)))
  cat(sprintf("Sympathetic Clearance (Kappa_S): %.4f +/- %.4f Hz (Tau: %.2fs)\n", mean(df$Kappa_S), sd(df$Kappa_S), mean(df$Tau_S)))
  cat(sprintf("Parasympathetic Clear. (Kappa_P): %.4f +/- %.4f Hz (Tau: %.2fs)\n", mean(df$Kappa_P), sd(df$Kappa_P), mean(df$Tau_P)))
  cat(sprintf("Observation Jitter (Lambda_R)  : %.5f +/- %.5f\n", mean(df$Lambda_R), sd(df$Lambda_R)))
  cat("=========================================================================\n\n")
}
