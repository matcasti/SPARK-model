# ==============================================================================
# R/model/model_variants.R -- the model specifications M0-M3, built on the unmodified engine.
#
# Only the objective assembly differs from the engine (its likelihood, priors and spectral
# anchor are reused), so the engine API and its filter are untouched (the compiled filter is
# checked against the engine filter in tests/test_pipeline.R).
#
#   M0   structural constraints only (positivity, kp > ks); lambda_P free
#   M1   + time-scale (log-normal) and separation (Beta) priors; lambda_P free
#   M2   + spectral anchor (lambda_P fixed from the wavelet partition) only
#   M3   M1 + M2 + timing-jitter prior (the reported model)
#
# Parameter vector theta = (ln ks, ln(kp - ks), ln lambda_R, [ln lambda_P if free]).
#
#   model_spec()       the flags of a specification        drop_term()      leave one prior out
#   prepare_series()   RR + priors + control, prepared once
#   fit_variant()      BFGS fit of one specification       fit_specs()      several specifications
#   variant_hessian()  curvature at the optimum            filter_eval()    filter (C++ or R)
# ==============================================================================


# ---- Filter evaluation (compiled port when available) ---------------------------------------------
load_fast_filter <- function(cfg) {
  compiled <- tryCatch({
    Rcpp::sourceCpp(cfg$paths$fast, env = globalenv())
    TRUE
  }, error = function(e) {
    log_msg("C++ filter unavailable (%s); using the slower R engine filter.", conditionMessage(e))
    FALSE
  })
  invisible(compiled)
}

has_cpp <- function() exists("gls_filter_cpp", mode = "function", envir = globalenv())

#' Filter evaluation: the C++ port when available, otherwise the engine's R filter.
#' Returns list(LL_star, nu0, sig2) and, with `full`, the states and innovations.
filter_eval <- function(dy, kp, ks, lp, lR, control, full = FALSE) {
  if (has_cpp()) {
    return(gls_filter_cpp(dy, kp, ks, lp, lR, control$jump_threshold, control$jump_power, control$lik_gate_max, full))
  }
  engine <- run_gls_filter(dy, kp, ks, lp, lR, control)
  out <- list(LL_star = engine$LL_star, nu0 = engine$nu0, sig2 = engine$sig2)
  if (full && is.finite(engine$LL_star)) {
    out <- c(out, list(
      X_p = engine$states[, "X_p"], X_s = engine$states[, "X_s"],
      X_p_var = engine$states[, "X_p_var"], X_s_var = engine$states[, "X_s_var"],
      implied_rr = engine$implied_rr, std_innov = engine$std_innov,
      innovations_diagnostic = engine$innovations_diagnostic,
      gate_multiplier = engine$gate_multiplier
    ))
  }
  out
}


# ---- Specifications ---------------------------------------------------------------------------------------
#' Which prior components a specification includes (tau, topo, jitter, anchor).
model_spec <- function(name) {
  none <- list(name = name, tau = FALSE, topo = FALSE, jitter = FALSE, anchor = FALSE)
  switch(name,
    M0  = none,
    M1  = modifyList(none, list(tau = TRUE, topo = TRUE)),
    M2  = modifyList(none, list(anchor = TRUE)),
    M3  = modifyList(none, list(tau = TRUE, topo = TRUE, jitter = TRUE, anchor = TRUE)),
    stop("Unknown model spec: ", name, call. = FALSE)
  )
}

#' Leave-one-out specification: `term` is one of "tau", "topo", "jitter", "anchor".
drop_term <- function(spec, term) {
  stopifnot(term %in% c("tau", "topo", "jitter", "anchor"))
  spec[[term]] <- FALSE
  spec$name <- paste0(spec$name, "-", term)
  spec
}

#' A series prepared once: NA-free RR intervals, the priors with their spectral anchor, and the control settings.
prepare_series <- function(dy, priors = spark_priors(), control = spark_control()) {
  dy <- as.numeric(dy[!is.na(dy)])
  list(dy = dy, priors = extract_spectral_priors(dy, priors, verbose = FALSE), control = control)
}


# ---- Parameters and objective -----------------------------------------------------------------------------------
#' theta -> physical parameters (lambda_P is fixed by the anchor, or free).
unpack_theta <- function(theta, spec, lambda_p_fixed) {
  ks <- exp(theta[1])
  kp <- ks + exp(theta[2])
  lR <- exp(theta[3])
  next_index <- 4L
  lp <- lambda_p_fixed
  if (!spec$anchor) {
    lp <- exp(theta[next_index])
    next_index <- next_index + 1L
  }
  list(ks = ks, kp = kp, lR = lR, lp = lp)
}

#' The engine's start values for a specification (a free lambda_P starts at the anchor).
default_theta <- function(spec, series) {
  priors <- series$priors
  control <- series$control
  theta <- c(log(priors$ks_mode),
             log(max(priors$kp_mode - priors$ks_mode, control$kp_gap_floor)),
             log(control$lR_start))
  if (!spec$anchor) theta <- c(theta, log(priors$energy_ratio))
  theta
}

#' Negative log-posterior of a specification (1e12 for non-finite values, as in the engine).
#' `engine_control` is spark_control(); it is not called `control` because optim() has its own.
variant_objective <- function(theta, dy, spec, priors, engine_control) {
  par <- unpack_theta(theta, spec, priors$energy_ratio)
  filtered <- filter_eval(dy, par$kp, par$ks, par$lp, par$lR, engine_control)
  if (!is.finite(filtered$LL_star)) return(1e12)

  penalties <- .log_prior_terms(c(ks = par$ks, kp = par$kp, lR = par$lR), filtered, dy, priors)
  log_prior <- (if (spec$tau) penalties$tau_s + penalties$tau_p else 0) +
               (if (spec$topo) penalties$topology else 0) +
               (if (spec$jitter) penalties$lR else 0)
  cost <- -(filtered$LL_star + log_prior)
  if (is.finite(cost)) cost else 1e12
}


# ---- Fitting ----------------------------------------------------------------------------------------------------------
#' Fit a specification by BFGS from `theta0` (default: the engine start). `ok = FALSE` on failure.
#' @return list(ok, spec, theta, par, tau_p, tau_s, objective, LL_star, nu0, sig2, conv, iters,
#'         rr_jitter_ms and, with `full`, the states, innovations, rmse_ms and share_p)
fit_variant <- function(series, spec, theta0 = NULL, full = TRUE) {
  start <- if (is.null(theta0)) default_theta(spec, series) else theta0
  optimum <- tryCatch(
    optim(start, variant_objective, dy = series$dy, spec = spec, priors = series$priors,
          engine_control = series$control, method = "BFGS", control = list(maxit = series$control$maxit)),
    error = function(e) NULL
  )
  if (is.null(optimum) || !is.finite(optimum$value) || optimum$value >= 1e11) {
    return(list(ok = FALSE, spec = spec$name))
  }
  par <- unpack_theta(optimum$par, spec, series$priors$energy_ratio)
  filtered <- filter_eval(series$dy, par$kp, par$ks, par$lp, par$lR, series$control, full = full)
  if (!is.finite(filtered$LL_star)) return(list(ok = FALSE, spec = spec$name))

  fit <- list(
    ok = TRUE, spec = spec$name, theta = optimum$par, par = par, tau_p = 1 / par$kp, tau_s = 1 / par$ks,
    objective = optimum$value, LL_star = filtered$LL_star, nu0 = filtered$nu0, sig2 = filtered$sig2,
    conv = optimum$convergence, iters = unname(optimum$counts[1]),
    rr_jitter_ms = 1000 * sqrt(filtered$sig2 * par$lR * median(series$dy)) / filtered$nu0
  )
  if (full) {
    fit$X_p <- filtered$X_p; fit$X_s <- filtered$X_s
    fit$X_p_var <- filtered$X_p_var; fit$X_s_var <- filtered$X_s_var
    fit$implied_rr <- filtered$implied_rr; fit$std_innov <- filtered$std_innov
    fit$innovations_diagnostic <- filtered$innovations_diagnostic
    fit$gate_multiplier <- filtered$gate_multiplier
    fit$rmse_ms <- 1000 * sqrt(mean((series$dy - filtered$implied_rr)^2))
    var_p <- var(filtered$X_p); var_s <- var(filtered$X_s)
    fit$share_p <- var_p / (var_p + var_s)                   # share of x_p in the variance of the states
  }
  fit
}

#' Fit several specifications of one prepared series from the default starts.
fit_specs <- function(series, names = c("M0", "M1", "M2", "M3"), full = TRUE) {
  setNames(lapply(names, function(name) fit_variant(series, model_spec(name), full = full)), names)
}

#' Numerical Hessian of the objective at the optimum: eigenvalues, condition number and the
#' Laplace SDs of ln(tau_p), ln(tau_s) by the delta method (pseudo-inverse, as in the engine's report).
variant_hessian <- function(fit, series, spec) {
  objective <- function(theta) variant_objective(theta, series$dy, spec, series$priors, series$control)
  hessian <- numDeriv::hessian(objective, fit$theta)
  eigenvalues <- eigen((hessian + t(hessian)) / 2, symmetric = TRUE, only.values = TRUE)$values
  covariance <- tryCatch(MASS::ginv(hessian), error = function(e) matrix(NA_real_, nrow(hessian), ncol(hessian)))

  ks <- fit$par$ks
  kp <- fit$par$kp
  jacobian <- matrix(0, 2, length(fit$theta))                  # rows: ln tau_s, ln tau_p
  jacobian[1, 1] <- -1                                          # d ln tau_s / d theta1
  jacobian[2, 1] <- -ks / kp                                    # d ln tau_p / d theta1
  jacobian[2, 2] <- -(kp - ks) / kp                             # d ln tau_p / d theta2
  ln_tau_covariance <- jacobian %*% covariance %*% t(jacobian)

  list(eigen = eigenvalues, min_eig = min(eigenvalues),
       cond = if (min(eigenvalues) > 0) max(eigenvalues) / min(eigenvalues) else Inf,
       sd_ln_tau_p = sqrt(max(ln_tau_covariance[2, 2], 0)), sd_ln_tau_s = sqrt(max(ln_tau_covariance[1, 1], 0)))
}
