# ==============================================================================
# R/model/synthetic_data.R -- CSV writers for simulated recordings with known ground truth, in two formats:
#   * current tools/prep-data.R format (default): one row per ECG beat, time = end of the interval,
#     rri = diff(time), covariates from the same beat;
#   * legacy = TRUE: the earlier script (rolling join of rri onto the beat times of a second signal,
#     2 MAD x 5 cleaning), kept to document and test why it was replaced.
# Used by the test suite and analysis/validation/12_preprocessing_checks.R. The outputs are synthetic.
#
#   write_synthetic_data()   the three dataset folders (prcp, fantasia, tmst) + ground truth
#   write_fantasia_csv(), write_prcp_csv(), write_tmst_csv()   one recording each
# (The order of random draws is part of the seeded output: do not reorder them.)
# ==============================================================================
suppressPackageStartupMessages(library(data.table))


# ---- Beat table helpers -----------------------------------------------------------------------------------
#' Offset between the file clock and the simulation clock: file time = sim time + offset.
.offset <- function(sim, t1) t1 - (sim$t[1] - sim$rr[1])

#' ECG peak times (s, 250 Hz quantization); the first peak precedes the first simulated interval.
.peaks <- function(sim, t1) round((t1 + c(0, sim$t - (sim$t[1] - sim$rr[1]))) * 250) / 250

#' Legacy join of the earlier script: rri (on the cumulative-RR clock) rolled onto the beat times of a
#' second signal (record clock), after the earlier 2 MAD x 5 cleaning of the RR intervals.
.prep_join <- function(peaks, signal_rows) {
  rr_interval <- clean_sample(diff(peaks), n_mad = 2, n_iter = 5)
  rr_table <- data.table(time = cumsum(rr_interval), rri = rr_interval)
  setkey(rr_table, time)
  rr_table[signal_rows, roll = Inf]
}

#' Covariates per beat of a simulation: SBP / DBP follow RR (baroreflex-like) and the tilt angle;
#' resp is a sinusoid at f_resp. Plausible shapes, not a physiological model.
.beat_covariates <- function(time, rr, angle = NULL, f_resp = NULL, seed = 1) {
  set.seed(seed)
  n <- length(rr)
  tilt_effect <- sin((if (is.null(angle)) 0 else angle) * pi / 180)
  slow_variation <- as.numeric(stats::filter(rnorm(n, 0, 1.2), 0.9, "recursive"))
  systolic <- 118 + 35 * (rr - mean(rr)) - 6 * tilt_effect + slow_variation + rnorm(n, 0, 0.8)
  pulse_pressure <- 45 - 12 * tilt_effect + rnorm(n, 0, 1)
  diastolic <- systolic - pulse_pressure + 3 * tilt_effect
  list(systolic = systolic, diastolic = diastolic,
       resp = if (is.null(f_resp)) NULL else sin(2 * pi * f_resp * time) + rnorm(n, 0, 0.15))
}

#' Write a recording table with MAP added; returns it invisibly.
.write_recording <- function(table, id, directory) {
  table[, map := (systolic + 2 * diastolic) / 3]
  fwrite(table, file.path(directory, paste0(id, ".csv")))
  invisible(table)
}

#' A covariate cleaned with the gate of the chosen format, or NA when the signal is absent.
.cleaned_or_na <- function(values, present, gate) {
  if (present) do.call(clean_sample, c(list(values), gate)) else NA_real_
}


# ---- One recording ------------------------------------------------------------------------------------------
write_fantasia_csv <- function(sim, id, dir, has_bp = TRUE, has_resp = TRUE, f_resp = 0.25, seed = 1, legacy = FALSE) {
  set.seed(seed)
  t1 <- runif(1, 0.9, 6)
  peaks <- .peaks(sim, t1)
  beat_times <- peaks[-length(peaks)]
  rr <- diff(peaks)
  covariates <- .beat_covariates(beat_times, rr, NULL, f_resp, seed)
  gate <- if (legacy) list(n_mad = 2, n_iter = 5) else list(n_mad = 3.5, n_iter = 2)

  signals <- data.table(systolic = .cleaned_or_na(covariates$systolic, has_bp, gate),
                        diastolic = .cleaned_or_na(covariates$diastolic, has_bp, gate),
                        resp = .cleaned_or_na(covariates$resp, has_resp, gate))
  table <- if (legacy) {
    .prep_join(peaks, cbind(data.table(time = beat_times), signals))
  } else {
    cbind(data.table(time = round(peaks[-1], 4), rri = do.call(clean_sample, c(list(rr), gate))), signals)
  }
  .write_recording(table, id, dir)
}

write_prcp_csv <- function(sim, id, dir, has_bp = TRUE, has_angle = TRUE, seed = 1, legacy = FALSE) {
  set.seed(seed)
  t1 <- runif(1, 0.9, 6)
  peaks <- .peaks(sim, t1)
  rr <- diff(peaks)
  angle_function <- attr(sim, "angle")
  offset <- .offset(sim, t1)
  tilt_angle <- function(times, noise_sd) {
    if (is.null(angle_function)) rep(0, length(times)) else angle_function(times - offset) + rnorm(length(times), 0, noise_sd)
  }

  if (legacy) {
    # Rows are arterial pulses: ECG peak + pulse transit delay, about 0.5% of the pulses missed
    pulse_onsets <- peaks[-length(peaks)] + 0.18 + rnorm(length(rr), 0, 0.004)
    pulse_onsets <- pulse_onsets[runif(length(pulse_onsets)) > 0.005]
    angle <- tilt_angle(pulse_onsets, 0.6)
    covariates <- .beat_covariates(pulse_onsets, approx(peaks[-1], rr, xout = pulse_onsets, rule = 2)$y, angle, NULL, seed)
    signals <- data.table(time = pulse_onsets,
                          systolic = if (has_bp) clean_sample(covariates$systolic, n_mad = 2, n_iter = 5) else NA_real_,
                          diastolic = if (has_bp) clean_sample(covariates$diastolic, n_mad = 2, n_iter = 5) else NA_real_,
                          angle = if (has_angle) angle else NA_real_)
    table <- .prep_join(peaks, signals)
  } else {
    beat_ends <- peaks[-1]
    beat_middles <- peaks[-length(peaks)] + rr / 2
    angle <- tilt_angle(beat_middles, 0.3)
    covariates <- .beat_covariates(beat_ends, rr, angle, NULL, seed)
    table <- data.table(time = round(beat_ends, 4), rri = clean_sample(rr),
                        systolic = if (has_bp) clean_sample(covariates$systolic) else NA_real_,
                        diastolic = if (has_bp) clean_sample(covariates$diastolic) else NA_real_,
                        angle = if (has_angle) angle else NA_real_)
  }
  .write_recording(table, id, dir)
}

write_tmst_csv <- function(sim, id, dir, exclude_n = 5L, max_min = 30) {
  rr <- round(clean_sample(sim$rr), 3)
  rr <- rr[-c(seq_len(exclude_n), (length(rr) - exclude_n + 1):length(rr))]        # edge beats removed
  time <- round(cumsum(rr), 3)
  keep <- (time / 60) < max_min
  table <- data.table(time = time[keep], rri = rr[keep])
  fwrite(table, file.path(dir, paste0(id, ".csv")))
  invisible(table)
}


# ---- The three datasets ----------------------------------------------------------------------------------------
#' Write the three dataset folders (prcp, fantasia, tmst) and return the ground truth.
#' PRCP: the protocols cycle through `prcp_protocols`; the last subject lacks blood pressure and the one before lacks the angle. Fantasia: the last
#' older subject lacks blood pressure and the last young subject lacks respiration.
write_synthetic_data <- function(root, n_prcp = 4, n_ex = 6, n_fant = 5, seed = 1,
                                 prcp_duration = 600, ex_duration = 900, fant_duration = 480,
                                 prcp_protocols = c("step", "stairs")) {
  dirs <- c(prcp = file.path(root, "prcp"), fantasia = file.path(root, "fantasia"), exercise = file.path(root, "tmst"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  truth <- list(); maneuver_rows <- list()

  for (i in seq_len(n_prcp)) {
    protocol <- prcp_protocols[(i - 1L) %% length(prcp_protocols) + 1L]
    duration <- if (protocol %in% c("repeated", "mixed")) max(prcp_duration, 150 + 3 * 420) else prcp_duration
    sim <- sim_tilt_angle(seed * 100 + i, duration, protocol)
    events <- attr(sim, "events")
    id <- sprintf("prcp%02d", i)
    has_bp <- i != n_prcp
    has_angle <- i != n_prcp - 1
    write_prcp_csv(sim, id, dirs[["prcp"]], has_bp = has_bp, has_angle = has_angle, seed = seed + i)
    truth[[id]] <- data.table(Dataset = "prcp", SubjectID = id, t_onset = events$t_onset, t_end = events$t_return,
                              protocol = protocol, has_bp = has_bp, has_angle = has_angle, has_resp = NA,
                              tau_rec_p = NA_real_, tau_rec_s = NA_real_)
    if (!is.null(events$maneuvers)) maneuver_rows[[id]] <- cbind(id = id, events$maneuvers, label = events$maneuvers$kind)
  }
  if (length(maneuver_rows)) fwrite(rbindlist(maneuver_rows), file.path(root, "prcp_events.csv"))      # protocol annotations, as written by tools/prep-data.R
  for (i in seq_len(n_ex)) {
    sim <- sim_exercise(seed * 200 + i, ex_duration)
    events <- attr(sim, "events")
    id <- sprintf("ex%02d", i)
    write_tmst_csv(sim, id, dirs[["exercise"]])
    truth[[id]] <- data.table(Dataset = "exercise", SubjectID = id, t_onset = events$t_onset, t_end = events$t_end,
                              protocol = NA_character_, has_bp = NA, has_angle = NA, has_resp = NA,
                              tau_rec_p = events$tau_rec_p, tau_rec_s = events$tau_rec_s)
  }
  for (group in c("young", "older")) for (i in seq_len(n_fant)) {
    f_resp <- runif(1, 0.2, 0.32)
    rsa_amplitude <- if (group == "young") 0.05 else 0.025
    sim_seed <- seed * 300 + i + (group == "older") * 50
    set.seed(sim_seed)
    sim <- sim_rest(sim_seed, group, duration_s = fant_duration, m_p = rsa_drive(f_resp, rsa_amplitude))
    id <- sprintf("f1%s%02d", substr(group, 1, 1), i)
    has_bp <- !(group == "older" && i == n_fant)
    has_resp <- !(group == "young" && i == n_fant)
    write_fantasia_csv(sim, id, dirs[["fantasia"]], has_bp = has_bp, has_resp = has_resp, f_resp = f_resp, seed = seed + i)
    truth[[id]] <- data.table(Dataset = "fantasia", SubjectID = id, t_onset = NA_real_, t_end = NA_real_,
                              protocol = NA_character_, has_bp = has_bp, has_angle = NA, has_resp = has_resp,
                              f_resp = f_resp, tau_rec_p = NA_real_, tau_rec_s = NA_real_)
  }
  list(dirs = dirs, truth = rbindlist(truth, fill = TRUE))
}
