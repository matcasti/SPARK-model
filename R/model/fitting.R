# ==============================================================================
# R/model/fitting.R -- the reported fits: one recording or a whole dataset, with the engine.
#
#   fit_recording()   fit one recording (M3) and add the file clock and the covariates
#   fit_dataset()     all recordings of a dataset in parallel (a failed recording is reported)
#   fit_model_fast()  the M3 fit: BFGS on the engine objective, compiled filter when available
#   as_rec()          wrap an RR vector as a recording (for simulations)
# ==============================================================================

#' M3 fit: BFGS on the negative log-posterior (compiled filter inside the optimizer when available;
#' the final pass is run by the R engine). Returns a `spark_fit`.
fit_model_fast <- function(dy, priors = spark_priors(), control = spark_control(), verbose = FALSE) {
  dy <- as.numeric(dy[!is.na(dy)])
  priors <- extract_spectral_priors(dy, priors, verbose = verbose)
  series <- list(dy = dy, priors = priors, control = control)
  spec <- model_spec("M3")

  optimum <- optim(default_theta(spec, series), variant_objective, dy = dy, spec = spec, priors = priors,
                   engine_control = control, method = "BFGS", control = list(maxit = control$maxit))
  if (!is.finite(optimum$value) || optimum$value >= 1e11) stop("Optimization failed.", call. = FALSE)

  par <- unpack_theta(optimum$par, spec, priors$energy_ratio)
  lambda_p <- priors$energy_ratio
  filtered <- run_gls_filter(dy, kp = par[["kp"]], ks = par[["ks"]], lp = lambda_p, lR = par[["lR"]], control = control)
  params <- list(
    ks = par[["ks"]], kp = par[["kp"]], lR = par[["lR"]], lp = lambda_p, nu0 = filtered$nu0, sig2 = filtered$sig2,
    sig_p = sqrt(2 * par[["kp"]] * lambda_p * filtered$sig2), sig_s = sqrt(2 * par[["ks"]] * filtered$sig2),
    rr_jitter_ms = 1000 * sqrt(filtered$sig2 * par[["lR"]] * median(dy)) / filtered$nu0
  )
  structure(list(
    dy = dy, time = cumsum(dy), model_type = "M3", priors = priors, control = control, params = params,
    LL_star = filtered$LL_star, states = filtered$states, implied_rr = filtered$implied_rr,
    innovations = filtered$std_innov, innovations_diagnostic = filtered$innovations_diagnostic,
    time_residuals = filtered$time_residuals, gate_multiplier = filtered$gate_multiplier, opt_raw = optimum),
    class = "spark_fit")
}

#' Wrap an RR vector (s) as a recording (clock = cumulative RR, no covariates).
as_rec <- function(rr, id = "rec") list(id = id, time = cumsum(rr), rr = rr, cov = data.table())

#' Fit one recording (M3) and summarize it with the engine's own .summarise_fit().
#' `rec` comes from load_dataset(): rr (s), time (the file's `time`) and covariates aligned by row.
#' The returned states gain Clock (the file's `time`), Time (seconds from the first beat of the
#' analysis window), X_p_var, X_s_var, Implied_RR and the covariate columns.
fit_recording <- function(rec, id, tag, priors, control) {
  fit <- fit_model_fast(rec$rr, priors, control)
  summarised <- .summarise_fit(fit, id, tag)
  states <- as.data.table(summarised$states)
  stopifnot(nrow(states) == length(rec$time))
  states[, `:=`(X_p_var = fit$states[, "X_p_var"], X_s_var = fit$states[, "X_s_var"], Implied_RR = fit$implied_rr,
                Clock = rec$time, Time = rec$time - rec$time[1])]
  if (ncol(rec$cov)) states <- cbind(states, rec$cov)
  list(summary = as.data.table(summarised$summary), states = states, pacf = summarised$pacf,
       params = fit$params, lp = fit$priors$energy_ratio, p_lf = fit$priors$p_lf, p_hf = fit$priors$p_hf)
}

#' Fit all recordings of a dataset (in parallel across recordings); failures are reported and skipped.
#' @return list(summary, states, pacf_matrix, params, spectral, failed)
fit_dataset <- function(rec_list, tag, cfg, priors = spark_priors(), control = spark_control()) {
  log_msg("Fitting %s: %d recordings (M3)", tag, length(rec_list))
  ids <- names(rec_list)
  fits <- map_recordings(seq_along(rec_list), function(i) {
    fit_recording(rec_list[[i]], ids[i], tag, priors, control)
  }, cfg$cores)
  names(fits) <- ids

  failed <- vapply(fits, is_err, logical(1))
  for (id in ids[failed]) log_msg("  FAILED %s: %s", id, as.character(fits[[id]]))
  fits <- fits[!failed]
  if (!length(fits)) stop("All recordings failed to fit for ", tag, call. = FALSE)

  spectral <- rbindlist(lapply(names(fits), function(id) {
    data.table(SubjectID = id, lambda_P = fits[[id]]$lp, p_lf = fits[[id]]$p_lf, p_hf = fits[[id]]$p_hf)
  }))
  list(summary = rbindlist(lapply(fits, `[[`, "summary")),
       states = rbindlist(lapply(fits, `[[`, "states"), fill = TRUE),
       pacf_matrix = do.call(rbind, lapply(fits, `[[`, "pacf")),
       params = lapply(fits, `[[`, "params"), spectral = spectral, failed = ids[failed])
}
