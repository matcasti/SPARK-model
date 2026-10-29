# ==============================================================================
# R/steps/benchmarks.R -- methodological benchmarks on simulated recordings with known truth.
#   Recovery-time definition and exercise-end rule  (analysis/07_recovery_benchmark.R)
#   RR-based event detection                        (analysis/08_event_detection_validation.R)
# ==============================================================================

# ---- Recovery definition and exercise-end rule -----------------------------------------------------------
BENCHMARK_DURATION_S <- 420
EXERCISE_END_RULES <- c("peak", "cp", "profile")

#' Spearman correlation; NA when fewer than 4 finite pairs or a constant input.
rank_correlation <- function(a, b) {
  finite <- is.finite(a) & is.finite(b)
  if (sum(finite) < 4 || sd(a[finite]) == 0 || sd(b[finite]) == 0) return(NA_real_)
  cor(a[finite], b[finite], method = "spearman")
}

#' Recovery metrics of both states for one realization and one exercise-end rule.
benchmark_rule_rows <- function(simulated_time, fit, events, rule, truth, realization, seed, cfg) {
  t_end <- switch(rule, peak = events$t_end_peak, cp = events$t_end_cp, profile = events$t_end_profile)
  t_on <- if (is.finite(events$t_onset_profile)) events$t_onset_profile else events$t_onset
  if (!is.finite(t_end) || !is.finite(t_on)) return(NULL)
  rows <- lapply(c("p", "s"), function(state) {
    metrics <- recovery_metrics(simulated_time, fit[[paste0("X_", state)]],
                                c(cfg$exercise$rest_skip_s, t_on - 5), t_end, cfg$exercise)
    if (is.null(metrics)) return(NULL)
    cbind(data.table(seed = seed, real = realization, state = state, end_method = rule,
                     true_tau = if (state == "p") truth$toff_p else truth$toff_s,
                     end_err = t_end - truth$t_off), metrics)
  })
  rbindlist(rows)
}

#' One true profile -> two noise realizations, each fitted and segmented once; every
#' exercise-end rule is evaluated on the same fits.
benchmark_one_profile <- function(seed, cfg) {
  truth <- draw_exercise_pars(seed, BENCHMARK_DURATION_S)
  rows <- list()
  for (realization in 1:2) {
    simulated <- sim_exercise(seed, BENCHMARK_DURATION_S, pars = truth,
                              noise_seed = id_seed(paste(seed, realization), cfg$seed), h = cfg$simulation$h)
    fit <- fit_variant(prepare_series(simulated$rr), model_spec("M3"))
    if (!isTRUE(fit$ok)) next
    events <- detect_events(simulated$rr, "exercise", cfg$segmentation)
    for (rule in EXERCISE_END_RULES) {
      rows[[length(rows) + 1L]] <- benchmark_rule_rows(cumsum(simulated$rr), fit, events, rule, truth,
                                                       realization, seed, cfg)
    }
  }
  rbindlist(rows)
}

#' Median error of the end of exercise by rule; the rule with the smallest median absolute error wins.
benchmark_end_rule_table <- function(table) {
  unique(table[, .(seed, real, end_method, end_err)])[
    , .(median_abs_err_s = median(abs(end_err), na.rm = TRUE), median_err_s = median(end_err, na.rm = TRUE), n = .N),
    by = end_method]
}

#' Stops when two rows share a seed, realization, state, definition and end rule: the simulated profiles must be distinct (an earlier seed
#' function repeated profiles, which silently turned the test-retest correlation into a correlation of counts).
assert_unique_realizations <- function(table) {
  keys <- table[, .(seed, real, state, end_method, definition)]
  if (anyDuplicated(keys)) stop("The benchmark has repeated (seed, realization) pairs: the simulated profiles must have distinct seeds.", call. = FALSE)
  invisible(TRUE)
}

#' Success rate, rank correlation with the true recovery constant and test-retest correlation,
#' per definition and state (for the chosen end rule).
benchmark_definition_performance <- function(table, end_rule) {
  assert_unique_realizations(table)
  table[end_method == end_rule][, {
    recovered <- status == "recovered"
    by_realization <- dcast(.SD[, .(seed, real, time)], seed ~ real, value.var = "time")
    .(n = .N, success = mean(recovered),
      rho_truth = rank_correlation(time[recovered], true_tau[recovered]),
      rho_retest = if (all(c("1", "2") %in% names(by_realization))) rank_correlation(by_realization[["1"]], by_realization[["2"]]) else NA_real_,
      median_rel_err = if (definition[1] == "expfit" && any(recovered)) median((time[recovered] - true_tau[recovered]) / true_tau[recovered]) else NA_real_)
  }, by = .(definition, state)]
}

#' Score = success * (truth rank correlation + retest correlation) / 2, averaged over the states.
#' @return list(definition, end_method, performance, aggregate, end_rule, n_profiles)
choose_recovery_definition <- function(table, n_profiles) {
  end_rules <- benchmark_end_rule_table(table)
  best_end <- end_rules[which.min(median_abs_err_s), end_method]
  performance <- benchmark_definition_performance(table, best_end)
  aggregate <- performance[, .(success = mean(success), rho_truth = mean(rho_truth, na.rm = TRUE),
                               rho_retest = mean(rho_retest, na.rm = TRUE)), by = definition]
  aggregate[, score := success * (rho_truth + rho_retest) / 2]
  if (all(is.na(aggregate$score))) {
    log_msg("No definition could be scored; falling back to 'sustained'.")
    aggregate[definition == "sustained", score := 0]
  }
  list(definition = aggregate[which.max(score), definition], end_method = best_end, performance = performance,
       aggregate = aggregate[order(-score)], end_rule = end_rules, n_profiles = n_profiles)
}


# ---- RR-based event detection --------------------------------------------------------------------------------
MIN_SHIFT_CANDIDATES_HZ <- c(0.05, 0.08, 0.10, 0.12)

#' Detection of a simulated tilt: detected?, onset error, number of segments.
detection_tilt_row <- function(seed, cfg) {
  simulated <- sim_tilt(seed, h = cfg$simulation$h)
  truth <- attr(simulated, "events")
  detected <- detect_events(simulated$rr, "tilt", cfg$segmentation)
  data.table(seed = seed, detected = detected$detected, onset_err = detected$t_onset - truth$t_onset,
             K = detected$K, total_shift_hz = truth$amp_s - truth$amp_p)
}

#' Detection of a simulated exercise bout: onset and end errors for each rule.
detection_exercise_row <- function(seed, cfg) {
  simulated <- sim_exercise(seed, h = cfg$simulation$h)
  truth <- attr(simulated, "events")
  detected <- detect_events(simulated$rr, "exercise", cfg$segmentation)
  data.table(seed = seed, detected = detected$detected,
             onset_err_cp = detected$t_onset - truth$t_onset, onset_err_profile = detected$t_onset_profile - truth$t_onset,
             end_err_peak = detected$t_end_peak - truth$t_end, end_err_cp = detected$t_end_cp - truth$t_end,
             end_err_profile = detected$t_end_profile - truth$t_end)
}

#' Detection of repeated simulated tilts: each true tilt is matched to the nearest detected onset
#' (within `cfg$tilt$match_tol_s`); `n_extra` counts detected events that match no tilt.
detection_repeated_row <- function(seed, cfg, n_tilts = 3L) {
  simulated <- sim_tilt_angle(seed, 150 + n_tilts * 420, "repeated", n_tilts = n_tilts, h = cfg$simulation$h)
  truth <- attr(simulated, "events")$tilts
  detected <- detect_events(simulated$rr, "tilt", cfg$segmentation)
  matched <- match_events(truth$t_onset, detected$events$t_onset, cfg$tilt$match_tol_s)
  found <- !is.na(matched$detected)
  data.table(seed = seed, tilt = seq_len(nrow(truth)), matched = found, onset_err = matched$error,
             return_err = ifelse(found, detected$events$t_recovery[matched$detected], NA_real_) - truth$t_return,
             n_extra = detected$n_events - sum(found))
}

#' False detections in a stationary recording for each minimum-shift candidate.
detection_rest_rows <- function(seed, cfg) {
  simulated <- sim_rest(seed, if (seed %% 2) "young" else "older", duration_s = cfg$simulation$duration_s, h = cfg$simulation$h)
  rbindlist(lapply(MIN_SHIFT_CANDIDATES_HZ, function(shift) {
    segmentation <- modifyList(cfg$segmentation, list(min_shift_hz = shift))
    data.table(seed = seed, min_shift_hz = shift, detected = detect_events(simulated$rr, "tilt", segmentation)$detected)
  }))
}

#' Detection rate of simulated tilts for one minimum-shift candidate.
detection_rate_for_shift <- function(shift, seeds, cfg) {
  segmentation <- modifyList(cfg$segmentation, list(min_shift_hz = shift))
  detected <- map_recordings(seeds, function(seed) {
    detect_events(sim_tilt(seed, h = cfg$simulation$h)$rr, "tilt", segmentation)$detected
  }, cfg$cores)
  data.table(min_shift_hz = shift, detection_rate = mean(unlist(Filter(Negate(is_err), detected))))
}
