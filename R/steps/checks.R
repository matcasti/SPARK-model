# ==============================================================================
# R/steps/checks.R -- numerical verification of statements in the supplementary material, and
# checks of the data-preparation choices. Called by analysis/validation/11_supplement_checks.R and
# analysis/validation/12_preprocessing_checks.R.
#
# Supplement checks
#   first_passage_bound_table()   bound on the terminal-threshold approximation (+ Monte Carlo)
#   fitted_regime_table()         relative timing jitter eps in the fitted recordings
#   ou_block_check_table()        closed-form OU covariance vs Van Loan numerical integration
#   cpp_vs_r_table()              compiled filter vs the R engine
#   calibration_table()           innovation diagnostics under a correctly specified model
# Preprocessing checks
#   prep_check_one_simulation()   earlier vs current file format, rri vs reconstructed RR, cleaning gate
# ==============================================================================


# ---- (a) Terminal-threshold approximation of the first-passage CDF ------------------------------------------------
#' For eps = sigma / sqrt(nu): sup of the neglected reflection term R(t) against its bound
#' eps / (2 sqrt(2 pi)); must be <= 1. (nu = 1: the problem is scale-free in eps.)
first_passage_bound_table <- function(eps_grid = c(5e-4, 1e-3, 5e-3, 1e-2, 2.5e-2, 5e-2, 0.1, 0.2, 0.4)) {
  times <- exp(seq(log(0.05), log(20), length.out = 4000))
  rbindlist(lapply(eps_grid, function(eps) {
    x <- (times + 1) / (eps * sqrt(times))
    reflection <- exp(2 / eps^2 + pnorm(-x, log.p = TRUE))
    bound <- eps / (2 * sqrt(2 * pi))
    data.table(eps = eps, sup_R = max(reflection), bound = bound, ratio = max(reflection) / bound,
               t_at_sup = times[which.max(reflection)])
  }))
}

#' Monte Carlo of the drift-diffusion first passage (Euler steps) against the closed-form CDF
#' and the terminal-threshold CDF. A coarse check of the closed form, not of the bound itself.
first_passage_monte_carlo <- function(eps, n_paths = 4000L, step = 2e-3, t_max = 3, seed = 1L) {
  with_seed(seed, {
    phase <- numeric(n_paths)
    hit_time <- rep(Inf, n_paths)
    alive <- rep(TRUE, n_paths)
    for (k in seq_len(ceiling(t_max / step))) {
      phase[alive] <- phase[alive] + step + eps * sqrt(step) * rnorm(sum(alive))
      crossed <- alive & phase >= 1
      hit_time[crossed] <- k * step
      alive[crossed] <- FALSE
    }
    grid <- seq(0.3, t_max, by = 0.01)
    empirical <- vapply(grid, function(x) mean(hit_time <= x), numeric(1))
    terminal <- pnorm((grid - 1) / (eps * sqrt(grid)))
    closed_form <- terminal + exp(2 / eps^2 + pnorm(-(grid + 1) / (eps * sqrt(grid)), log.p = TRUE))
    data.table(eps = eps, sup_emp_vs_closed = max(abs(empirical - closed_form)),
               sup_emp_vs_terminal = max(abs(empirical - terminal)),
               sup_closed_vs_terminal = max(abs(closed_form - terminal)), mc_se = 0.5 / sqrt(n_paths))
  })
}


# ---- (b) Fitted regime ---------------------------------------------------------------------------------------------------
#' Relative timing jitter eps = jitter / median RR and the implied CDF bound, per fitted recording.
fitted_regime_table <- function(cfg) {
  rbindlist(lapply(cfg$datasets, function(dataset) {
    fits <- read_result(paste0("fit_", dataset), cfg, FALSE)
    if (is.null(fits)) return(NULL)
    median_rr <- fits$states[, .(median_rr = median(RR_Interval)), by = SubjectID]
    jitter <- data.table(SubjectID = names(fits$params),
                         jitter_ms = vapply(fits$params, function(p) p$rr_jitter_ms, numeric(1)))
    table <- merge(jitter, median_rr, by = "SubjectID")
    table[, `:=`(Dataset = dataset, eps = jitter_ms / 1000 / median_rr,
                 bound_cdf = (jitter_ms / 1000 / median_rr) / (2 * sqrt(2 * pi)))]
  }))
}


# ---- (c) OU covariance: closed form vs Van Loan -----------------------------------------------------------------------------
#' Phase-integral transition and process covariance by numerical (Van Loan) integration.
van_loan_ou <- function(k, dt, variance) {
  drift <- matrix(c(-k, 0, 1, 0), 2, byrow = TRUE)
  diffusion <- matrix(c(2 * k * variance, 0, 0, 0), 2)
  van_loan <- rbind(cbind(-drift, diffusion), cbind(matrix(0, 2, 2), t(drift))) * dt
  expm <- as.matrix(Matrix::expm(van_loan))
  phi <- t(expm[3:4, 3:4])
  list(Phi = phi, Q = phi %*% expm[1:2, 3:4])
}

ou_block_check_table <- function() {
  cases <- expand.grid(k = c(0.01, 0.2, 1.5), dt = c(1e-4, 3e-3, 0.05, 0.4, 1.2, 3), V = c(0.5, 2))
  rbindlist(lapply(seq_len(nrow(cases)), function(i) {
    k <- cases$k[i]; dt <- cases$dt[i]; V <- cases$V[i]
    closed <- .compute_ou_block(k, dt, V)
    numeric <- van_loan_ou(k, dt, V)
    q_numeric <- c(numeric$Q[1, 1], numeric$Q[1, 2], numeric$Q[2, 2])
    q_closed <- c(closed$q11, closed$q13, closed$q33)
    data.table(k = k, dt = dt, V = V, kdt = k * dt, taylor = k * dt < 1e-3,
               max_rel_err = max(abs(q_closed - q_numeric) / pmax(abs(q_numeric), 1e-300)),
               phi_err = abs(closed$phi_int - numeric$Phi[2, 1]))
  }))
}


# ---- (d) Compiled filter vs the R engine ------------------------------------------------------------------------------------------
cpp_vs_r_table <- function(seed) {
  if (!has_cpp()) return(NULL)
  set.seed(seed)
  n <- 500
  clock <- cumsum(rep(0.85, n))
  rr <- 0.85 + 0.03 * sin(2 * pi * clock / 9) + 0.02 * sin(2 * pi * clock / 60) + rnorm(n, 0, 0.01)
  rr[c(100, 250)] <- c(1.6, 0.45)                                           # artifacts exercise the multiplier
  control <- spark_control()
  parameter_sets <- list(c(0.3, 0.006, 2.7, 0.1), c(0.5, 0.02, 1.2, 0.01), c(0.2, 1e-4, 5, 0.3))
  rbindlist(lapply(parameter_sets, function(p) {
    engine <- run_gls_filter(rr, p[1], p[2], p[3], p[4], control)
    compiled <- filter_eval(rr, p[1], p[2], p[3], p[4], control, full = TRUE)
    data.table(kp = p[1], ks = p[2], lp = p[3], lR = p[4], d_LL = abs(engine$LL_star - compiled$LL_star),
               d_x_p = max(abs(engine$states[, "X_p"] - compiled$X_p)),
               d_x_s = max(abs(engine$states[, "X_s"] - compiled$X_s)),
               d_rr = max(abs(engine$implied_rr - compiled$implied_rr)))
  }))
}


# ---- (e) Calibration of the innovation diagnostics ---------------------------------------------------------------------------------
#' One correctly specified simulation (time scales inside the prior ranges), fitted with the engine
#' specification; KS distance, PACF violations and innovation SD.
calibration_one_simulation <- function(i, cfg) {
  seed <- id_seed(paste("calib", i), cfg$seed)
  truth <- with_seed(seed, list(tau_p = runif(1, 3.5, 5.5), tau_s = runif(1, 120, 250)))
  simulated <- simulate_rr(cfg$calibration$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s,
                           lR = 0.02, seed = seed, h = cfg$simulation$h)
  id <- paste0("sim", i)
  fit <- fit_recording(as_rec(simulated$rr, id), id, "calibration", spark_priors(), spark_control())
  summary <- fit$summary
  data.table(sim = i, N = summary$NumValidBeats, KS_Stat = summary$KS_Stat, crit = 1.36 / sqrt(summary$NumValidBeats),
             pacf_viol = summary$PACF_Violations, z_mean = mean(fit$states[Valid == TRUE, Innovation]),
             z_sd = sd(fit$states[Valid == TRUE, Innovation]))
}

calibration_table <- function(cfg) {
  results <- map_recordings(seq_len(cfg$calibration$n_sims), calibration_one_simulation, cfg$cores, cfg = cfg)
  rbindlist(Filter(Negate(is_err), results))
}


# ---- Preprocessing checks ------------------------------------------------------------------------------------------------------------
#' Gate settings compared on a simulated artifact-free series: beats changed, RMSSD and SDNN ratios.
cleaning_gate_rows <- function(raw_rr, simulation_id, group) {
  gates <- list(c(2, 5), c(3, 5), c(3.5, 2), c(4, 1))
  rbindlist(lapply(gates, function(gate) {
    cleaned <- clean_sample(raw_rr, n_mad = gate[1], n_iter = gate[2])
    data.table(sim = simulation_id, group = group, n_mad = gate[1], n_iter = gate[2],
               frac_changed = mean(abs(cleaned - raw_rr) > 0.001),
               rmssd_ratio = sqrt(mean(diff(cleaned)^2)) / sqrt(mean(diff(raw_rr)^2)),
               sdnn_ratio = sd(cleaned) / sd(raw_rr))
  }))
}

#' One simulated Fantasia-like and PRCP-like recording written by the earlier ("legacy") and the
#' current preparation script; fits from the rri column (and, for legacy files, from diff(time)).
prep_check_one_simulation <- function(i, cfg, directory) {
  group <- if (i %% 2) "young" else "older"
  fantasia_sim <- sim_rest(900 + i, group, duration_s = cfg$simulation$duration_s,
                           m_p = rsa_drive(0.25, if (group == "young") 0.05 else 0.025))
  prcp_sim <- sim_tilt_angle(950 + i, cfg$simulation$duration_s, "step")
  sims <- list(fantasia = fantasia_sim, prcp = prcp_sim)

  for (script in c("legacy", "current")) {
    folder <- file.path(directory, script)
    dir.create(folder, showWarnings = FALSE)
    is_legacy <- script == "legacy"
    write_fantasia_csv(fantasia_sim, sprintf("f%02d", i), folder, has_bp = FALSE, has_resp = FALSE, seed = i, legacy = is_legacy)
    write_prcp_csv(prcp_sim, sprintf("p%02d", i), folder, has_bp = FALSE, has_angle = FALSE, seed = i, legacy = is_legacy)
  }

  rows <- list()
  for (kind in c("fantasia", "prcp")) for (script in c("legacy", "current")) {
    truth <- sims[[kind]]
    file <- file.path(directory, script, sprintf("%s%02d.csv", substr(kind, 1, 1), i))
    for (source in if (script == "legacy") c("rri", "dtime") else "rri") {
      rec <- read_recording(file, cfg, rr_source = source)
      fit <- fit_variant(prepare_series(rec$rr), model_spec("M3"))
      rows[[length(rows) + 1]] <- data.table(
        kind = kind, script = script, src = source, sim = i, n = length(rec$rr),
        consistency = rec$qc$time_rri_consistency, dup_frac = mean(diff(rec$rr) == 0),
        sd_rr_ms = 1000 * sd(rec$rr), rmssd_ms = 1000 * sqrt(mean(diff(rec$rr)^2)),
        jitter_ms = fit$rr_jitter_ms, tau_p = fit$tau_p, tau_s = fit$tau_s,
        true_tau_p = 1 / attr(truth, "truth")$kp, true_sd_rr_ms = 1000 * sd(truth$rr),
        true_rmssd_ms = 1000 * sqrt(mean(diff(truth$rr)^2)),
        true_dup_frac = mean(diff(round(truth$rr / 0.004) * 0.004) == 0)
      )
    }
  }
  list(src = rbindlist(rows), clean = cleaning_gate_rows(fantasia_sim$rr, i, group))
}
