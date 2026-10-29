# ==============================================================================
# tools/prep-data.R -- prepare one CSV per recording for the SPARK analysis pipeline.
#
# Output (data/<dataset>/<id>.csv), ONE ROW PER RR INTERVAL (one ECG beat):
#   time       s   END of the interval: time[k] - time[k-1] == rri[k]
#                  (record clock for Fantasia/PRCP; cumulative clock for TMST)
#   rri        s   RR interval ending at `time`, from the ECG annotations
#   systolic, diastolic, map   mmHg  blood pressure of the cardiac cycle of the interval
#   resp       a.u. mean respiration over the interval (Fantasia)
#   angle      deg  mean tilt-table angle over the interval (PRCP)
# Signals a subject does not have are written as all-NA columns.
#
# Differences from the earlier version (all checked in tests/test_prep.R):
#  * no rolling join: RR and covariates share one beat table, so `rri` is exactly
#    diff(time) and no interval is repeated or skipped (the join left 4-11% repeated values);
#  * PRCP: RR comes from the ECG (wqrs) beats and blood pressure is taken from the ABP
#    waveform in the pulse window of each ECG beat (before, rows were ABP pulses);
#  * blood pressure is read in a window shifted by the estimated ECG-to-pulse delay, so
#    each window holds one complete pulse;
#  * non-beat annotation symbols are excluded; implausible blood pressure is set to NA;
#  * the cleaning gate is 3.5 MAD x 2 passes (the old 2 MAD x 5 changed 4-7% of artifact-free
#    beats and lowered RMSSD by 8-11%; 3.5 MAD changes ~0.2% and still corrects ~99.7% of
#    injected artifacts; see tests/benchmark_cleaning_gate.R);
#  * vectorized (no rbind in a loop); records are processed independently and failures are
#    reported instead of stopping the run; a QC table is written (data/prep_qc.csv).
#
#  * PRCP: the protocol annotations (.anI) are parsed into data/prcp_events.csv, one row per maneuver (kind: slow tilt, rapid tilt,
#    stand-up; onset / return times). The analysis uses them to label the tilts, to add the stand-ups (which the tilt angle cannot
#    show) and to space baselines after every maneuver; the cleaning gate leaves beats around maneuver boundaries untouched.
#
# Usage: Rscript tools/prep-data.R   (from the project root; needs data-raw/ with the WFDB files)
# Override before sourcing: RAW_DIR, OUT_DIR, PREP_DATASETS, PREP_SETTINGS, PREP_AUTORUN.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ---- Settings -----------------------------------------------------------------------------------------
if (!exists("RAW_DIR")) RAW_DIR <- "data-raw"
if (!exists("OUT_DIR")) OUT_DIR <- "data"
if (!exists("PREP_DATASETS")) PREP_DATASETS <- c("fantasia", "prcp", "tmst")
if (!exists("PREP_AUTORUN")) PREP_AUTORUN <- TRUE

PREP_DEFAULTS <- list(
  # Adaptive outlier gate for rri / blood pressure / respiration
  clean = list(window = 50L, n_mad = 3.5, n_iter = 2L),
  # WFDB annotation symbols that are heartbeats (others, e.g. rhythm or noise marks, are ignored)
  beat_symbols = c("N", "L", "R", "B", "A", "a", "J", "S", "V", "r", "F", "e", "j", "n", "E", "/", "f", "Q", "?"),
  bp_range = list(sbp = c(50, 300), dbp = c(20, 200), min_pulse_pressure = 10),
  min_valid_frac = 0.8,                    # a signal with fewer valid values counts as absent
  gap_s = 2.5,                             # intervals longer than this are reported as gaps
  delay_beats = 600L, delay_max_s = 0.6,   # ECG-to-pulse delay estimation
  tmst = list(exclude_n = 5L, max_minutes = 30),   # edge beats removed; recording cut at 30 min
  # PRCP protocol annotations (extension .anI: start / end of every tilt, stand-up, return to supine). They are written to
  # prcp_events.csv, and the RR cleaning gate does not replace physiological beats within [-before, +after] seconds of a
  # maneuver boundary, where genuine fast transitions (a return to supine can change RR by 0.4 s in a few beats) look like outliers.
  events = list(extension = "anI", protect_before_s = 5, protect_after_s = 30, rr_limits_s = c(0.3, 2.0))
)
PREP_SETTINGS <- if (exists("PREP_SETTINGS")) modifyList(PREP_DEFAULTS, PREP_SETTINGS) else PREP_DEFAULTS


# ---- Cleaning ---------------------------------------------------------------------------------------------
#' Adaptive outlier cleaning: running-median / running-MAD gate, rejected samples replaced by linear
#' interpolation, repeated `n_iter` times. Length is always preserved.
#' @return numeric vector; `attr(, "changed")` flags the replaced samples.
clean_sample <- function(x, window = PREP_SETTINGS$clean$window, n_mad = PREP_SETTINGS$clean$n_mad,
                         n_iter = PREP_SETTINGS$clean$n_iter) {
  original <- x
  one_pass <- function(x) {
    was_na <- is.na(x)
    if (all(was_na)) return(x)
    if (any(was_na)) {
      i <- seq_along(x)
      x[was_na] <- approx(i[!was_na], x[!was_na], xout = i[was_na], rule = 2)$y
    }
    n <- length(x)
    k <- min(window, if (n %% 2 == 1) n else n - 1L)
    if (k %% 2 == 0) k <- k + 1L
    if (k < 3L) return(x)
    running_median <- runmed(x, k = k)
    abs_dev <- abs(x - running_median)
    running_mad <- pmax(1.4826 * runmed(abs_dev, k = k), median(abs_dev) * 0.1)   # floor: no zero gate on flat stretches
    inside <- x >= running_median - n_mad * running_mad & x <= running_median + n_mad * running_mad & !was_na
    if (sum(inside) < 2) {
      warning("Not enough valid samples to interpolate.")
      return(rep(NA_real_, n))
    }
    i <- seq_len(n)
    approx(i[inside], x[inside], xout = i, rule = 2)$y
  }
  cleaned <- one_pass(x)
  if (n_iter >= 2) for (pass in seq_len(n_iter - 1L)) cleaned <- one_pass(cleaned)
  attr(cleaned, "changed") <- is.finite(original) & is.finite(cleaned) & abs(cleaned - original) > 1e-9
  cleaned
}


# ---- Beat table -------------------------------------------------------------------------------------------------
#' Beat annotation sample numbers (0-based), restricted to beat symbols.
read_beats <- function(wfdb, record_path, extension) {
  annotation <- wfdb$rdann(record_name = record_path, extension = extension)
  samples <- as.numeric(annotation$sample)
  symbols <- as.character(unlist(annotation$symbol))
  is_beat <- if (length(symbols) == length(samples)) symbols %in% PREP_SETTINGS$beat_symbols else rep(TRUE, length(samples))
  if (!any(is_beat)) {
    warning(sprintf("%s.%s: no standard beat symbols; using all annotations.", basename(record_path), extension))
    is_beat[] <- TRUE
  }
  sort(unique(samples[is_beat]))
}

#' Per-beat summaries of a continuous signal over the windows [start_k, start_{k+1}) (0-based sample
#' indices), vectorized. @return data.table(beat, max, min, mean) for beats 1..K (NA where empty).
window_stats <- function(signal, starts) {
  K <- length(starts) - 1L
  window_of_sample <- findInterval(seq_along(signal) - 1, starts)
  values <- data.table(id = window_of_sample, v = signal)[id >= 1L & id <= K & is.finite(v)]
  stats <- values[, .(max = max(v), min = min(v), mean = mean(v)), by = id]
  out <- data.table(beat = seq_len(K), max = NA_real_, min = NA_real_, mean = NA_real_)
  out[stats$id, `:=`(max = stats$max, min = stats$min, mean = stats$mean)]
  out
}

#' Delay (samples) from the ECG R peak to the foot of the arterial pulse: ensemble-average the R-to-R
#' segments and take the minimum of the smoothed average (the diastolic foot).
estimate_pulse_delay <- function(bp, peaks, fs, settings = PREP_SETTINGS) {
  n_beats <- length(peaks) - 1L
  segment_length <- as.integer(round(median(diff(peaks))))
  chosen <- unique(round(seq(1, n_beats, length.out = min(settings$delay_beats, n_beats))))
  segment_of <- function(k) {
    samples <- peaks[k] + seq_len(segment_length)
    if (max(samples) > length(bp)) rep(NA_real_, segment_length) else bp[samples]
  }
  ensemble <- rowMeans(vapply(chosen, segment_of, numeric(segment_length)), na.rm = TRUE)
  if (!all(is.finite(ensemble))) return(0L)
  k <- max(3L, as.integer(round(0.05 * fs)))
  smoothed <- as.numeric(stats::filter(ensemble, rep(1 / k, k), sides = 2))
  smoothed[is.na(smoothed)] <- ensemble[is.na(smoothed)]
  search_limit <- as.integer(round(settings$delay_max_s * fs))
  which.min(smoothed[seq_len(min(search_limit, segment_length))]) - 1L     # position 1 = the R-peak sample itself
}


# ---- One WFDB record (Fantasia, PRCP) -------------------------------------------------------------------------------
#' Blood-pressure columns (cleaned systolic / diastolic per beat) of a record, or NULL when absent or
#' mostly implausible. Updates `qc` with the pulse delay and the invalid share.
blood_pressure_columns <- function(bp, peaks, fs, qc, settings) {
  if (is.null(bp) || sum(is.finite(bp)) == 0) return(list(columns = NULL, qc = qc))
  delay <- estimate_pulse_delay(bp, peaks, fs, settings)
  qc$delay_s <- delay / fs
  per_beat <- window_stats(bp, peaks + delay)
  systolic <- per_beat$max
  diastolic <- per_beat$min
  limits <- settings$bp_range
  implausible <- !is.finite(systolic) | !is.finite(diastolic) |
    systolic < limits$sbp[1] | systolic > limits$sbp[2] |
    diastolic < limits$dbp[1] | diastolic > limits$dbp[2] |
    (systolic - diastolic) < limits$min_pulse_pressure
  qc$bp_invalid_frac <- mean(implausible)
  if (mean(!implausible) < settings$min_valid_frac) return(list(columns = NULL, qc = qc))
  systolic[implausible] <- NA
  diastolic[implausible] <- NA
  qc$has_bp <- TRUE
  list(columns = list(systolic = as.numeric(clean_sample(systolic)), diastolic = as.numeric(clean_sample(diastolic))), qc = qc)
}

# ---- Protocol annotations (PRCP .anI) ---------------------------------------------------------------------------------------
#' Free-text protocol annotations of one record: data.table(time_s, label); NULL when the file is absent or empty.
read_event_annotations <- function(wfdb, record_path, extension, fs) {
  if (is.null(extension) || !file.exists(paste0(record_path, ".", extension))) return(NULL)
  annotation <- tryCatch(wfdb$rdann(record_name = record_path, extension = extension), error = function(e) NULL)
  if (is.null(annotation)) return(NULL)
  samples <- as.numeric(annotation$sample)
  labels <- as.character(unlist(annotation$aux_note))
  if (!length(samples) || length(labels) != length(samples)) return(NULL)
  data.table(time_s = samples / fs, label = labels)
}

#' Maneuvers of one record from its protocol labels (PRCP .anI). The wording differs between records, and all of these are recognized
#' (case, spaces, verbs and typos tolerated):
#'   tilts:      "Start|Initiation of|Initiate  slow|fast|rapid tilt up|down", "End|Conclusion of|Conclude  ... tilt up|down|back down"
#'   stand-ups:  "Stand-up from supine", "Stand up"  ...  "Back to supine", "Transition back to supine"
#'   Finapres:   "Blood pressure signal fades", "Switched off Finapres" ... "Switched Finapres on" (arterial pressure unusable in between);
#'               "Blood pressure artifacts" is only counted. Comment lines ("## ...") are ignored.
#' Labels are read in time order; a tilt opens at its start of ascent and closes at the end of its descent (or at the next ascent, when an
#' annotation is missing); a stand-up opens at "stand" and closes at "back to supine". Anything else is returned in attr(, "unparsed").
#' @return data.table(id, kind, onset_s, onset_end_s, return_s, return_end_s, label); times in s of the record clock (NA when not annotated);
#'         attributes: unparsed (labels), bp_gaps (data.table(from_s, to_s) without usable arterial pressure), n_bp_flags
parse_maneuvers <- function(id, time_s, label) {
  empty <- data.table(id = character(), kind = character(), onset_s = numeric(), onset_end_s = numeric(),
                      return_s = numeric(), return_end_s = numeric(), label = character())
  no_gaps <- data.table(from_s = numeric(), to_s = numeric())
  if (!length(time_s)) return(structure(empty, unparsed = character(), bp_gaps = no_gaps, n_bp_flags = 0L))
  o <- order(time_s)
  time_s <- time_s[o]; label <- label[o]
  text <- tolower(trimws(gsub("[[:space:]]+", " ", label)))
  is_comment <- !nzchar(text) | grepl("^#", text)
  is_bp <- !is_comment & grepl("finapres|blood pressure|arterial pressure|\\bbp\\b", text)
  bp_on <- is_bp & grepl("\\bon\\b|restor|resum", text) & !grepl("\\boff\\b", text)
  bp_off <- is_bp & !bp_on & grepl("\\boff\\b|fade|lost|loss|drop", text)
  is_tilt <- !is_comment & !is_bp & grepl("tilt", text)
  is_stand <- !is_comment & !is_bp & !is_tilt & grepl("stand", text) & !grepl("back", text)
  is_back <- !is_comment & !is_bp & !is_tilt & grepl("back to supin|return to supin|supin", text)
  speed <- ifelse(grepl("slow", text), "slow", ifelse(grepl("fast|rapid", text), "rapid", NA_character_))
  is_start <- grepl("^(start|initiat|begin|onset|commenc)", text)
  is_end <- grepl("^(end|conclu|complet|finish|terminat|stop)", text)
  goes_down <- grepl("dow|back", text)
  rows <- list(); open_tilt <- NULL; open_stand <- NULL; unparsed <- character()
  close_tilt <- function() { if (!is.null(open_tilt)) rows[[length(rows) + 1L]] <<- open_tilt; open_tilt <<- NULL }
  for (k in seq_along(text)) {
    if (is_comment[k] || is_bp[k]) next
    if (is_tilt[k] && is_start[k] && !goes_down[k]) {                      # start of the ascent: a new tilt
      close_tilt()
      kind <- if (identical(speed[k], "rapid")) "rapid tilt" else if (identical(speed[k], "slow")) "slow tilt" else "tilt"
      open_tilt <- data.table(id = id, kind = kind, onset_s = time_s[k], onset_end_s = NA_real_, return_s = NA_real_, return_end_s = NA_real_, label = trimws(label[k]))
    } else if (is_tilt[k] && is_end[k] && !goes_down[k] && !is.null(open_tilt)) {
      open_tilt$onset_end_s <- time_s[k]
    } else if (is_tilt[k] && is_start[k] && goes_down[k] && !is.null(open_tilt)) {
      open_tilt$return_s <- time_s[k]
    } else if (is_tilt[k] && is_end[k] && goes_down[k] && !is.null(open_tilt)) {
      open_tilt$return_end_s <- time_s[k]; close_tilt()
    } else if (is_stand[k]) {
      rows[[length(rows) + 1L]] <- data.table(id = id, kind = "stand-up", onset_s = time_s[k], onset_end_s = NA_real_, return_s = NA_real_, return_end_s = NA_real_, label = trimws(label[k]))
      open_stand <- length(rows)
    } else if (is_back[k] && !is.null(open_stand) && is.na(rows[[open_stand]]$return_s)) {
      rows[[open_stand]]$return_s <- time_s[k]; open_stand <- NULL
    } else unparsed <- c(unparsed, trimws(label[k]))
  }
  close_tilt()
  out <- if (length(rows)) rbindlist(rows) else empty
  setorder(out, onset_s)
  # intervals without usable arterial pressure: from the first sign of loss to the next switch-on (open-ended when there is none)
  gaps <- list(); from <- NA_real_
  for (k in which(bp_off | bp_on)) {
    if (bp_off[k]) { if (is.na(from)) from <- time_s[k] }
    else if (!is.na(from)) { gaps[[length(gaps) + 1L]] <- data.table(from_s = from, to_s = time_s[k]); from <- NA_real_ }
  }
  if (!is.na(from)) gaps[[length(gaps) + 1L]] <- data.table(from_s = from, to_s = Inf)
  structure(out, unparsed = unparsed, bp_gaps = if (length(gaps)) rbindlist(gaps) else no_gaps, n_bp_flags = sum(is_bp & !bp_on & !bp_off))
}

#' Beats that lie within [-before, +after] seconds of a maneuver boundary.
protected_beats <- function(time, maneuvers, before, after) {
  boundaries <- sort(unique(stats::na.omit(unlist(maneuvers[, .(onset_s, onset_end_s, return_s, return_end_s)]))))
  if (!length(boundaries)) return(rep(FALSE, length(time)))
  i <- findInterval(time, boundaries)                                      # last boundary at or before each beat
  before_gap <- ifelse(i >= 1L, time - boundaries[pmax(i, 1L)], Inf)       # seconds since the last boundary
  after_gap <- ifelse(i < length(boundaries), boundaries[pmin(i + 1L, length(boundaries))] - time, Inf)   # seconds to the next boundary
  before_gap <= after | after_gap <= before
}


#' Beat table of an ECG-annotated WFDB record with optional BP / respiration / angle channels.
#' @param bp_name,resp_name,angle_name channel names (NULL to skip the signal)
#' @return list(table = data.table, qc = one-row data.table)
prep_wfdb_record <- function(wfdb, record_path, beat_ext, bp_name = NULL, resp_name = NULL, angle_name = NULL,
                             settings = PREP_SETTINGS, event_ext = NULL) {
  record <- wfdb$rdrecord(record_name = record_path)
  fs <- as.numeric(record$fs)
  channel_names <- as.character(unlist(record$sig_name))
  channel <- function(name) {
    column <- if (is.null(name)) integer() else which(channel_names == name)
    if (length(column)) as.numeric(record$p_signal[, column[1]]) else NULL
  }

  peaks <- read_beats(wfdb, record_path, beat_ext)
  if (length(peaks) < 10L) stop("fewer than 10 beat annotations")
  n_beats <- length(peaks) - 1L
  rri_raw <- diff(peaks) / fs
  table <- data.table(time = round(peaks[-1] / fs, 4), rri = rri_raw)      # time = END of the interval
  qc <- data.table(n_beats = n_beats, duration_s = round(sum(rri_raw), 1), fs = fs,
                   n_gaps = sum(rri_raw > settings$gap_s), longest_gap_s = round(max(rri_raw), 2),
                   delay_s = NA_real_, bp_invalid_frac = NA_real_,
                   has_bp = FALSE, has_resp = FALSE, has_angle = FALSE)

  # Blood pressure: one pulse window per ECG beat, shifted by the ECG-to-pulse delay
  table[, `:=`(systolic = NA_real_, diastolic = NA_real_)]
  pressure <- blood_pressure_columns(channel(bp_name), peaks, fs, qc, settings)
  qc <- pressure$qc
  if (!is.null(pressure$columns)) table[, `:=`(systolic = pressure$columns$systolic, diastolic = pressure$columns$diastolic)]

  # Respiration or tilt angle: mean over the interval
  if (!is.null(resp_name)) {
    table[, resp := NA_real_]
    resp <- channel(resp_name)
    if (!is.null(resp) && sum(is.finite(resp)) > 0) {
      per_beat <- window_stats(resp, peaks)$mean
      if (mean(is.finite(per_beat)) >= settings$min_valid_frac) {
        table[, resp := as.numeric(clean_sample(per_beat))]
        qc$has_resp <- TRUE
      }
    }
  }
  if (!is.null(angle_name)) {
    table[, angle := NA_real_]
    angle <- channel(angle_name)
    if (!is.null(angle) && sum(is.finite(angle)) > 0) {
      per_beat <- window_stats(angle, peaks)$mean
      if (mean(is.finite(per_beat)) >= settings$min_valid_frac) {
        table[, angle := per_beat]
        qc$has_angle <- TRUE
      }
    }
  }

  cleaned_rr <- clean_sample(rri_raw)
  changed <- attr(cleaned_rr, "changed")
  maneuvers <- NULL
  qc$n_maneuvers <- 0L; qc$n_protected <- 0L; qc$n_unparsed <- 0L; qc$n_bp_masked <- 0L; qc$n_bp_flags <- 0L
  qc$beat_coverage <- round(max(table$time) / (nrow(record$p_signal) / fs), 3)       # beats should span the whole record
  if (qc$beat_coverage < 0.9) message(sprintf("[%s] the beat annotations cover only %.0f%% of the record", gsub("-", "", basename(record_path)), 100 * qc$beat_coverage))
  if (!is.null(event_ext)) {                                               # protocol annotations: maneuvers, and beats that the gate must not replace
    id <- gsub("-", "", basename(record_path))
    annotations <- read_event_annotations(wfdb, record_path, event_ext, fs)
    if (!is.null(annotations)) {
      maneuvers <- parse_maneuvers(id, annotations$time_s, annotations$label)
      qc$n_maneuvers <- nrow(maneuvers); qc$n_unparsed <- length(attr(maneuvers, "unparsed"))
      if (qc$n_unparsed) message(sprintf("[%s] annotations not recognised: %s", id, paste(unique(attr(maneuvers, "unparsed")), collapse = "; ")))
      gaps <- attr(maneuvers, "bp_gaps")                                    # arterial pressure lost (Finapres switched off) between annotated times
      if (!is.null(gaps) && nrow(gaps)) {
        masked <- vapply(table$time, function(t) any(t >= gaps$from_s & t <= gaps$to_s), logical(1))
        table[masked, `:=`(systolic = NA_real_, diastolic = NA_real_)]
        qc$n_bp_masked <- sum(masked)
      }
      qc$n_bp_flags <- attr(maneuvers, "n_bp_flags")
      keep <- protected_beats(table$time, maneuvers, settings$events$protect_before_s, settings$events$protect_after_s) &
        is.finite(rri_raw) & rri_raw >= settings$events$rr_limits_s[1] & rri_raw <= settings$events$rr_limits_s[2]
      cleaned_rr <- ifelse(keep, rri_raw, cleaned_rr)
      changed <- changed & !keep
      qc$n_protected <- sum(keep)
    }
  }
  table[, rri := as.numeric(cleaned_rr)]
  qc$rr_changed_frac <- mean(changed)
  table[, map := (systolic + 2 * diastolic) / 3]
  setcolorder(table, intersect(c("time", "rri", "systolic", "diastolic", "resp", "angle", "map"), names(table)))
  list(table = table, qc = qc, events = maneuvers)
}


# ---- One TMST text file -----------------------------------------------------------------------------------------------
#' TMST: RR intervals (ms) from a text file; edge beats removed; cut at max_minutes.
prep_tmst_record <- function(path, settings = PREP_SETTINGS) {
  intervals_ms <- suppressWarnings(as.numeric(readLines(path, warn = FALSE)))
  intervals_ms <- intervals_ms[is.finite(intervals_ms)]
  if (length(intervals_ms) < 3 * settings$tmst$exclude_n) stop("too few intervals")
  cleaned <- clean_sample(intervals_ms)
  changed_frac <- mean(attr(cleaned, "changed"))
  rr <- round(as.numeric(cleaned) / 1000, 3)
  n <- length(rr)
  rr <- rr[-c(seq_len(settings$tmst$exclude_n), (n - settings$tmst$exclude_n + 1L):n)]
  time <- round(cumsum(rr), 3)
  keep <- (time / 60) < settings$tmst$max_minutes
  qc <- data.table(n_beats = sum(keep), duration_s = round(sum(rr[keep]), 1), fs = NA_real_,
                   n_gaps = sum(rr[keep] > settings$gap_s), longest_gap_s = max(rr[keep]),
                   delay_s = NA_real_, bp_invalid_frac = NA_real_,
                   has_bp = FALSE, has_resp = FALSE, has_angle = FALSE, rr_changed_frac = changed_frac)
  list(table = data.table(time = time[keep], rri = rr[keep]), qc = qc)
}


# ---- Dataset drivers -------------------------------------------------------------------------------------------------------
#' Record names of a WFDB folder (from the .hea files).
wfdb_records <- function(dir) sub("\\.hea$", "", list.files(dir, pattern = "\\.hea$"))

#' Prepare every record of one dataset; failures are reported and skipped. Output ids drop '-'.
prep_dataset <- function(dataset, wfdb) {
  in_dir <- file.path(RAW_DIR, dataset)
  out_dir <- file.path(OUT_DIR, dataset)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(in_dir)) {
    message(sprintf("[%s] folder not found: %s (skipped)", dataset, in_dir))
    return(NULL)
  }
  records <- if (dataset == "tmst") sub("\\.txt$", "", list.files(in_dir, pattern = "\\.txt$")) else wfdb_records(in_dir)

  prepare_record <- function(record) {
    path <- file.path(in_dir, record)
    switch(dataset,
      fantasia = prep_wfdb_record(wfdb, path, "ecg", bp_name = "BP", resp_name = "RESP"),
      prcp     = prep_wfdb_record(wfdb, path, "wqrs", bp_name = "ABP", angle_name = "Angle", event_ext = PREP_SETTINGS$events$extension),
      tmst     = prep_tmst_record(paste0(path, ".txt")))
  }
  qc_rows <- list(); event_rows <- list()
  for (record in records) {
    id <- gsub("-", "", record)
    prepared <- tryCatch(prepare_record(record), error = function(e) {
      message(sprintf("[%s] %s FAILED: %s", dataset, record, conditionMessage(e)))
      NULL
    })
    if (is.null(prepared)) next
    fwrite(prepared$table, file.path(out_dir, paste0(id, ".csv")))
    qc_rows[[id]] <- cbind(data.table(dataset = dataset, id = id), prepared$qc)
    if (!is.null(prepared$events) && nrow(prepared$events)) event_rows[[id]] <- prepared$events
    message(sprintf("[%s] %-10s %5d beats | RR changed %.2f%% | BP %s resp %s angle %s", dataset, id,
                    prepared$qc$n_beats, 100 * prepared$qc$rr_changed_frac,
                    prepared$qc$has_bp, prepared$qc$has_resp, prepared$qc$has_angle))
  }
  if (length(event_rows)) {                                                # one table for the dataset: data/<dataset>_events.csv
    fwrite(rbindlist(event_rows), file.path(OUT_DIR, paste0(dataset, "_events.csv")))
    message(sprintf("[%s] %d maneuvers of %d records written to %s_events.csv", dataset, sum(vapply(event_rows, nrow, 1L)), length(event_rows), dataset))
  }
  rbindlist(qc_rows, fill = TRUE)
}

prep_all <- function() {
  wfdb <- NULL
  if (any(c("fantasia", "prcp") %in% PREP_DATASETS)) {
    library(reticulate)
    # reticulate >= 1.41 installs wfdb on demand; with older versions install wfdb in the active Python
    if (utils::packageVersion("reticulate") >= "1.41") py_require("wfdb")
    wfdb <- import("wfdb")
  }
  qc <- rbindlist(lapply(PREP_DATASETS, prep_dataset, wfdb = wfdb), fill = TRUE)
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (nrow(qc)) fwrite(qc, file.path(OUT_DIR, "prep_qc.csv"))
  invisible(qc)
}

if (PREP_AUTORUN) qc <- prep_all()
