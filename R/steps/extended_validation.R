# ==============================================================================
# R/steps/extended_validation.R -- three supplementary analyses that probe the identification structure from outside the main results.
# Called by analysis/validation/13_identifiability_boundary.R, analysis/datasets/14_fantasia_pooling.R and analysis/validation/16_state_descriptors.R.
#
#   Identifiability boundary   simulated recordings over the time-scale ratio, the variance ratio and the length; M0-M3 fitted to each
#     boundary_jobs(), boundary_job()
#   Partial pooling            Fantasia windows A and B refitted with the production priors, weak priors and empirical-Bayes priors from the other half
#     pooling_fit(), pooling_hyperparameters(), pooling_reliability()
#   Perturbation descriptors   STATE descriptors (R/signals/state_descriptors.R) on simulated recordings with known states and on the real recordings
#     descriptor_simulation(), descriptor_sim_job(), descriptor_real_job()
# ==============================================================================

# ---- Identifiability boundary ---------------------------------------------------------------------------------------------
#' One job per simulated recording: axis "time-scale ratio" (tau_s / tau_p with lambda_P fixed) and axis "variance ratio" (lambda_P with tau_s / tau_p fixed), each at every length.
boundary_jobs <- function(settings, seed) {
  cells <- rbind(CJ(axis = "time-scale ratio", ratio = settings$ratios, lambda = settings$lambda_fixed, duration = settings$durations),
                 CJ(axis = "variance ratio", ratio = settings$ratio_fixed, lambda = settings$lambdas, duration = settings$durations))
  unlist(lapply(seq_len(nrow(cells)), function(i) lapply(seq_len(settings$n_per_cell), function(j)
    c(as.list(cells[i]), list(sim = j, seed = id_seed(paste("ident", i, j), seed))))), recursive = FALSE)
}

#' Fit every specification to one simulated recording and compare with the truth (recovery_row(), R/steps/simulation_recovery.R).
boundary_job <- function(job, settings, h) {
  tau_p <- settings$tau_p; tau_s <- job$ratio * tau_p
  sim <- simulate_rr(job$duration, kp = 1 / tau_p, ks = 1 / tau_s, lp = job$lambda, sig2 = 9e-4, lR = 0.02, h = h, seed = job$seed)
  fits <- fit_specs(prepare_series(sim$rr)); true_share <- var(sim$x_p) / (var(sim$x_p) + var(sim$x_s))
  rows <- rbindlist(lapply(fits, recovery_row, truth = list(tau_p = tau_p, tau_s = tau_s, lambda_p = job$lambda), simulated = sim, true_share = true_share), fill = TRUE)
  rows[, `:=`(axis = job$axis, ratio = job$ratio, lambda = job$lambda, duration = job$duration, sim = job$sim)]
}


# ---- Partial pooling of the time constants --------------------------------------------------------------------------------------
#' M3 with the time-scale priors replaced by 95% bounds `bounds_p`, `bounds_s` (s); amplitudes and time constants of the fit, NULL on failure.
pooling_fit <- function(rr, bounds_p, bounds_s) {
  fit <- tryCatch(fit_model_fast(rr, spark_priors(tau_p_bounds = bounds_p, tau_s_bounds = bounds_s), spark_control()), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  p <- fit$params; data.table(tau_p = 1 / p$kp, tau_s = 1 / p$ks, SD_p = sqrt(p$lp * p$sig2), SD_s = sqrt(p$sig2))
}

#' Hyperparameters of ln tau from the weak-prior fits of the subjects in `ids`: the mean over both windows and the between-subject SD from the covariance of the two
#' windows (test-retest), floored at 0.02.
pooling_hyperparameters <- function(weak, ids) {
  a <- weak[win == "A" & id %in% ids][order(id)]; b <- weak[win == "B" & id %in% ids][order(id)]; keep <- intersect(a$id, b$id); a <- a[id %in% keep]; b <- b[id %in% keep]
  one <- function(v) { la <- log(a[[v]]); lb <- log(b[[v]]); c(mu = mean(c(la, lb)), sd_b = sqrt(max(cov(la, lb), 0.02^2))) }
  list(p = one("tau_p"), s = one("tau_s"))
}

#' ICC(3,1) and Spearman correlation of the log values between windows A and B, young-vs-older AUC (window A) and the share of tau inside the production prior.
pooling_reliability <- function(fits, groups) {
  rbindlist(lapply(unique(fits$setting), function(st) {
    d <- fits[setting == st]; a <- d[win == "A"][order(id)]; b <- d[win == "B"][order(id)]; keep <- intersect(a$id, b$id); a <- a[id %in% keep]; b <- b[id %in% keep]
    g <- groups[a$id]
    rbindlist(lapply(c("tau_p", "tau_s", "SD_p", "SD_s"), function(v) {
      la <- log(a[[v]]); lb <- log(b[[v]]); r <- rank(a[[v]]); n1 <- sum(g == "young"); n2 <- sum(g == "older")
      data.table(setting = st, variable = v, n = length(keep), icc = icc_consistency(la, lb), spearman = cor(la, lb, method = "spearman"),
                 auc_young_older = (sum(r[g == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2), sd_log = sd(c(la, lb)),
                 in_production_prior = switch(v, tau_p = mean(a[[v]] >= 3 & a[[v]] <= 6), tau_s = mean(a[[v]] >= 100 & a[[v]] <= 300), NA_real_))
    }))
  }))
}


# ---- Perturbation-response descriptors -------------------------------------------------------------------------------------------
#' One simulated recording with known states (simulate_rr_osc(): a respiratory oscillator the model does not represent). Scenarios: "exercise" (rest -> exercise ->
#' recovery, draw_exercise_pars()), "tilt" (smaller, tilt-like step) and "shift" (the exercise displacement persists: a true regime shift, no recovery).
descriptor_simulation <- function(scenario, seed, duration_s) {
  set.seed(seed)
  d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))),
            osc = list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1)))))
  make <- function(mp, ms) simulate_rr_osc(duration_s, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = d$osc, m_p_fun = mp, m_s_fun = ms)
  if (scenario %in% c("exercise", "shift")) {
    p <- draw_exercise_pars(seed, duration_s); if (scenario == "shift") p$toff_p <- p$toff_s <- 1e6
    return(c(d, list(pert = c(p$t_on, p$t_off), sim = make(prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s)))))
  }
  t_on <- runif(1, 0.28, 0.45) * duration_s; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
  c(d, list(pert = c(t_on, t_off), sim = make(prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))),
                                              prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100)))))))
}

#' Descriptors of the TRUE path (whitened by the true stationary SDs) against those of the M3 reconstruction (posterior draws, scale from M3 refitted on the
#' resting beats), plus the variance inflation of the full-recording fit relative to the truth and Ptransition.
descriptor_sim_job <- function(job, settings) {
  z <- descriptor_simulation(job$scenario, job$seed, settings$duration_s); s <- z$sim; H <- settings$horizon[[job$scenario]]; rest <- c(20, z$pert[1] - 5)
  series <- prepare_series(s$rr); f <- fit_variant(series, model_spec("M3")); if (!isTRUE(f$ok)) return(NULL)
  rest_series <- prepare_series(s$rr[s$t < z$pert[1] - 5]); fr <- fit_variant(rest_series, model_spec("M3")); if (!isTRUE(fr$ok)) return(NULL)
  V_true <- c(z$lp * 9e-4, 9e-4); V_full <- c(series$priors$energy_ratio * f$sig2, f$sig2); V_rest <- c(rest_series$priors$energy_ratio * fr$sig2, fr$sig2)
  gt <- states_on_grid(s$t, s$x_p + s$z1, s$x_s, 0 * s$t, 0 * s$t); r <- gt$t >= rest[1] & gt$t < rest[2]
  truth <- distance_descriptors(gt$t, whitened_distance((gt$xp - mean(gt$xp[r])) / sqrt(V_true[1]), (gt$xs - mean(gt$xs[r])) / sqrt(V_true[2])), z$pert, H, settings$hold_s, rest_mask = r)
  ge <- states_on_grid(s$t, f$X_p, f$X_s, f$X_p_var, f$X_s_var)
  est <- summarise_draws(descriptor_draws(ge, rest, z$pert, V_rest, c(fr$tau_p, fr$tau_s), settings$n_draw, H, settings$hold_s, seed = job$seed), settings$horizons)
  if (is.null(est)) return(NULL)
  pt <- regime_transition(ge, rest, z$pert, V_rest, H)
  data.table(scenario = job$scenario, sim = job$sim, H = H, Trec_true = truth[["Trec"]], Dmax_true = truth[["Dmax"]], Preturn_true = truth[["Preturn"]],
             Dmax = est[["Dmax"]], Dmax_lo = est[["Dmax_lo"]], Dmax_hi = est[["Dmax_hi"]], Trec = est[["Trec"]], Trec_lo = est[["Trec_lo"]], Trec_hi = est[["Trec_hi"]],
             Preturn = est[["Preturn"]], rest_exceed = est[["rest_exceed"]], K = pt[["K"]], transition = pt[["trans"]],
             V_ratio_full_p = V_full[1] / V_true[1], V_ratio_full_s = V_full[2] / V_true[2], V_ratio_rest_p = V_rest[1] / V_true[1], V_ratio_rest_s = V_rest[2] / V_true[2])
}

#' One real perturbation: `job` has the beat-wise states of the recording (t, x_p, x_s, vp, vs, rr), the rest and perturbation windows, the horizon, an optional end of
#' follow-up (t_max) and whether the HMM sees the whole recording (exercise) or only the event window (tilts, whose recordings hold several events).
descriptor_real_job <- function(job, settings) {
  s <- job$s; g <- states_on_grid(s$t, s$x_p, s$x_s, s$vp, s$vs); if (is.finite(job$t_max)) g <- g[t <= job$t_max]
  k <- s$t >= job$rest[1] & s$t < job$rest[2]; if (sum(k) < 60) return(NULL)
  rest_series <- prepare_series(s$rr[k]); fr <- fit_variant(rest_series, model_spec("M3")); if (!isTRUE(fr$ok)) return(NULL)
  V <- c(rest_series$priors$energy_ratio * fr$sig2, fr$sig2)
  est <- summarise_draws(descriptor_draws(g, job$rest, job$pert, V, c(fr$tau_p, fr$tau_s), settings$n_draw, job$H, settings$hold_s, seed = job$seed), settings$horizons)
  if (is.null(est)) return(NULL)
  pt <- regime_transition(if (job$whole) g else g[t >= job$rest[1]], job$rest, job$pert, V, job$H)
  rest_rr <- s$rr[k]
  data.table(id = job$id, subject = job$subject, kind = job$kind, H = job$H, t(est), K = pt[["K"]], transition = pt[["trans"]],
             rest_hr = 60 / mean(rest_rr), rest_sdnn_ms = 1000 * sd(rest_rr), rest_rmssd_ms = 1000 * sqrt(mean(diff(rest_rr)^2)))
}
