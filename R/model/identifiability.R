# ==============================================================================
# R/model/identifiability.R -- is the latent decomposition identified by RR alone?
# Dispersed multi-start fits of a specification, ambiguity metrics, and profile surfaces
# over (tau_p, tau_s). Called by analysis/03_identifiability.R.
#
#   draw_start_params()   dispersed starting values, shared across specifications
#   run_multistart()      fit one specification from every start
#   ident_metrics()       distinct solutions, near-equivalent fits, collapse, trajectory dispersion
#   profile_surface()     objective over a (tau_p, tau_s) grid, other parameters optimized
#   identifiability_analyse_recording()   all specifications for one recording
# ==============================================================================


# ---- Starting points --------------------------------------------------------------------------------
#' Dispersed start values (log-uniform), shared across specifications so that models are compared
#' from the same points. Row 1 is the engine default (its values are NA).
draw_start_params <- function(n, seed, tau_s_range = c(20, 600), tau_p_range = c(0.5, 30),
                              lR_range = c(1e-4, 1e-1), lp_range = c(0.05, 20)) {
  log_uniform <- function(range) exp(runif(1, log(range[1]), log(range[2])))
  with_seed(seed, {
    starts <- data.table(start = seq_len(n), tau_s = NA_real_, tau_p = NA_real_, lR = NA_real_, lp = NA_real_)
    for (i in seq_len(n)[-1]) {
      repeat {                                          # vagal faster than sympathetic: tau_p < 0.8 tau_s
        tau_s <- log_uniform(tau_s_range)
        tau_p <- log_uniform(tau_p_range)
        if (tau_p < 0.8 * tau_s) break
      }
      set(starts, i, c("tau_s", "tau_p", "lR", "lp"),
          list(tau_s, tau_p, log_uniform(lR_range), log_uniform(lp_range)))
    }
    starts
  })
}

#' theta for one start row and specification (NA row = the engine default start).
theta_from_start <- function(start, spec, series) {
  if (is.na(start$tau_s)) return(default_theta(spec, series))
  ks <- 1 / start$tau_s
  kp <- 1 / start$tau_p
  theta <- c(log(ks), log(kp - ks), log(start$lR))
  if (!spec$anchor) theta <- c(theta, log(start$lp))
  theta
}


# ---- Multi-start fits ----------------------------------------------------------------------------------
#' Fit `spec` from every start. @return list(table = one row per start, traj = x_p/x_s of the successful fits)
run_multistart <- function(series, spec, starts) {
  rows <- vector("list", nrow(starts))
  trajectories <- vector("list", nrow(starts))
  for (i in seq_len(nrow(starts))) {
    fit <- fit_variant(series, spec, theta0 = theta_from_start(starts[i], spec, series), full = TRUE)
    if (!isTRUE(fit$ok)) {
      rows[[i]] <- data.table(start = starts$start[i], ok = FALSE)
      next
    }
    rows[[i]] <- data.table(
      start = starts$start[i], ok = TRUE, conv = fit$conv, objective = fit$objective, LL_star = fit$LL_star,
      tau_p = fit$tau_p, tau_s = fit$tau_s, lR = fit$par$lR, lp = fit$par$lp, nu0 = fit$nu0, sig2 = fit$sig2,
      jitter_ms = fit$rr_jitter_ms, rmse_ms = fit$rmse_ms, share_p = fit$share_p
    )
    trajectories[[i]] <- cbind(x_p = fit$X_p, x_s = fit$X_s)
  }
  list(table = rbindlist(rows, fill = TRUE), traj = trajectories)
}


# ---- Ambiguity metrics ---------------------------------------------------------------------------------
#' Greedy count of distinct solutions: points closer than `radius` in (ln tau_p, ln tau_s)
#' belong to the same solution.
count_solutions <- function(ln_tau_p, ln_tau_s, radius) {
  centers <- matrix(numeric(), 0, 2)
  for (i in seq_along(ln_tau_p)) {
    point <- c(ln_tau_p[i], ln_tau_s[i])
    distances <- sqrt(rowSums((centers - matrix(point, nrow(centers), 2, byrow = TRUE))^2))
    if (!nrow(centers) || min(distances) > radius) centers <- rbind(centers, point)
  }
  nrow(centers)
}

#' Across-start trajectory dispersion (RMS across-start SD / typical within-trajectory SD) and the
#' median pairwise correlation, for a list of equally long trajectories (column `column`).
trajectory_dispersion <- function(trajectories, column) {
  if (length(trajectories) < 2) return(c(disp = NA_real_, cor = NA_real_))
  matrix_ <- do.call(cbind, lapply(trajectories, function(m) m[, column]))
  typical_sd <- sqrt(mean(apply(matrix_, 2, var)))
  dispersion <- if (typical_sd > 0) sqrt(mean(apply(matrix_, 1, var))) / typical_sd else NA_real_
  correlations <- suppressWarnings(cor(matrix_))
  c(disp = dispersion, cor = median(correlations[upper.tri(correlations)], na.rm = TRUE))
}

#' Ambiguity metrics of one (recording, specification) multi-start result.
ident_metrics <- function(multistart, priors, cfg) {
  settings <- cfg$identifiability
  fits <- multistart$table[ok == TRUE]
  counts <- data.table(n_starts = nrow(multistart$table), n_ok = nrow(fits))
  if (!nrow(fits)) return(cbind(counts, data.table(frac_converged = NA_real_)))

  # Near-equivalent fits: profiled log-likelihood within `ll_tol` of the best
  equivalent <- fits[LL_star >= max(LL_star) - settings$ll_tol]
  with_trajectory <- which(!vapply(multistart$traj, is.null, logical(1)))
  equivalent_trajectories <- multistart$traj[with_trajectory][match(equivalent$start, multistart$table$start[with_trajectory])]
  dispersion_p <- trajectory_dispersion(equivalent_trajectories, "x_p")
  dispersion_s <- trajectory_dispersion(equivalent_trajectories, "x_s")

  sd_or_na <- function(x) if (length(x) > 1) sd(x) else NA_real_
  collapsed <- function(share) mean(share < settings$collapse_share[1] | share > settings$collapse_share[2])
  within <- function(x, bounds) x >= bounds[1] & x <= bounds[2]

  cbind(counts, data.table(
    frac_converged = mean(fits$conv == 0),
    best_LL = max(fits$LL_star),
    best_rmse_ms = fits$rmse_ms[which.max(fits$LL_star)],
    median_rmse_ms = median(fits$rmse_ms),
    n_solutions = count_solutions(log(fits$tau_p), log(fits$tau_s), settings$min_dispersion_log),
    n_equiv = nrow(equivalent),
    n_solutions_equiv = count_solutions(log(equivalent$tau_p), log(equivalent$tau_s), settings$min_dispersion_log),
    sd_ln_tau_p = sd_or_na(log(fits$tau_p)), sd_ln_tau_s = sd_or_na(log(fits$tau_s)),
    sd_ln_tau_p_equiv = sd_or_na(log(equivalent$tau_p)), sd_ln_tau_s_equiv = sd_or_na(log(equivalent$tau_s)),
    sd_ln_lp = sd_or_na(log(fits$lp)),
    share_p_range_equiv = diff(range(equivalent$share_p)),
    collapse_frac = collapsed(fits$share_p), collapse_frac_equiv = collapsed(equivalent$share_p),
    traj_disp_p_equiv = dispersion_p[["disp"]], traj_disp_s_equiv = dispersion_s[["disp"]],
    traj_cor_p_equiv = dispersion_p[["cor"]], traj_cor_s_equiv = dispersion_s[["cor"]],
    frac_in_phys_range = mean(within(fits$tau_p, priors$tau_p_bounds) & within(fits$tau_s, priors$tau_s_bounds)),
    frac_tau_p_in_range = mean(within(fits$tau_p, priors$tau_p_bounds)),
    frac_tau_s_in_range = mean(within(fits$tau_s, priors$tau_s_bounds)),
    median_jitter_ms = median(fits$jitter_ms)
  ))
}


# ---- Profile surface -----------------------------------------------------------------------------------
#' Objective over a log-spaced (tau_p, tau_s) grid; at each point the remaining free parameters
#' (lambda_R and, if free, lambda_P) are optimized. Points with tau_p >= tau_s are skipped.
#' @return one row per grid point with obj, LL_star and the loss relative to the best point (dobj, dLL)
profile_surface <- function(series, spec, n = 21L, tau_s_range = c(20, 600), tau_p_range = c(0.5, 30)) {
  tau_p_grid <- exp(seq(log(tau_p_range[1]), log(tau_p_range[2]), length.out = n))
  tau_s_grid <- exp(seq(log(tau_s_range[1]), log(tau_s_range[2]), length.out = n))
  start_theta <- default_theta(spec, series)
  free_start <- start_theta[3:length(start_theta)]             # parameters other than the two time scales

  skipped <- function(tau_p, tau_s) data.table(tau_p = tau_p, tau_s = tau_s, obj = NA_real_, LL_star = NA_real_)
  rows <- vector("list", n * n)
  k <- 0L
  for (tau_p in tau_p_grid) for (tau_s in tau_s_grid) {
    k <- k + 1L
    if (tau_p >= tau_s) { rows[[k]] <- skipped(tau_p, tau_s); next }
    ks <- 1 / tau_s
    kp <- 1 / tau_p
    objective <- function(free) {
      variant_objective(c(log(ks), log(kp - ks), free), series$dy, spec, series$priors, series$control)
    }
    optimum <- tryCatch({
      if (length(free_start) == 1) {
        r <- optimize(function(x) objective(x), interval = log(c(1e-6, 10)))
        list(value = r$objective, par = r$minimum)
      } else {
        optim(free_start, objective, method = "Nelder-Mead", control = list(maxit = 200))
      }
    }, error = function(e) NULL)
    if (is.null(optimum) || optimum$value >= 1e11) { rows[[k]] <- skipped(tau_p, tau_s); next }

    par <- unpack_theta(c(log(ks), log(kp - ks), optimum$par), spec, series$priors$energy_ratio)
    filtered <- filter_eval(series$dy, par$kp, par$ks, par$lp, par$lR, series$control)
    rows[[k]] <- data.table(tau_p = tau_p, tau_s = tau_s, obj = optimum$value, LL_star = filtered$LL_star)
  }
  surface <- rbindlist(rows)
  surface[, dobj := obj - min(obj, na.rm = TRUE)][, dLL := max(LL_star, na.rm = TRUE) - LL_star][]
}


# ---- One recording: all specifications ----------------------------------------------------------------------
#' Multi-start fits of every specification for one recording, with the ambiguity metrics, the
#' Hessian diagnostics of M0 and M3 and, for the first recordings, profile surfaces and trajectories.
#' @param job list(ds, id, rr, surface, keep_traj)
IDENTIFIABILITY_SPECS <- c("M0", "M1", "M2", "M3")

identifiability_analyse_recording <- function(job, cfg) {
  settings <- cfg$identifiability
  series <- prepare_series(job$rr)
  starts <- draw_start_params(settings$n_starts, id_seed(job$id, cfg$seed))
  label <- function(table, spec_name) cbind(data.table(Dataset = job$ds, SubjectID = job$id, spec = spec_name), table)
  out <- list(metrics = list(), starts = list(), hess = list(), surf = list(), traj = list())

  for (spec_name in IDENTIFIABILITY_SPECS) {
    spec <- model_spec(spec_name)
    multistart <- run_multistart(series, spec, starts)
    out$metrics[[spec_name]] <- label(ident_metrics(multistart, series$priors, cfg), spec_name)
    out$starts[[spec_name]]  <- label(multistart$table, spec_name)
    if (job$keep_traj) {
      out$traj[[spec_name]] <- list(t = cumsum(series$dy), table = multistart$table, mats = multistart$traj)
    }
    if (spec_name %in% c("M0", "M3")) {                       # local curvature and profile surface
      best <- fit_variant(series, spec, full = FALSE)
      if (isTRUE(best$ok)) {
        curvature <- variant_hessian(best, series, spec)[c("min_eig", "cond", "sd_ln_tau_p", "sd_ln_tau_s")]
        out$hess[[spec_name]] <- label(as.data.table(curvature), spec_name)
      }
      if (job$surface) {
        out$surf[[spec_name]] <- label(profile_surface(series, spec, settings$surface_grid), spec_name)
      }
    }
  }
  out
}

#' Median of every ambiguity metric by dataset and specification.
identifiability_summarise <- function(metrics) {
  metric_columns <- setdiff(names(metrics), c("Dataset", "SubjectID", "spec"))
  metrics[, lapply(.SD, function(x) median(as.numeric(x), na.rm = TRUE)),
          by = .(Dataset, spec), .SDcols = metric_columns]
}
