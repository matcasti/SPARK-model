# ==============================================================================
# R/core/io.R -- reading the CSV recordings written by tools/prep-data.R.
#
#   load_dataset()     a whole dataset folder  -> index (one row per recording) + recordings
#   read_recording()   one CSV                 -> list(id, time, rr, cov, qc)
#   window_rec()       cut a recording to a time window
#   clean_sample()     outlier gate (same as tools/prep-data.R); used to rebuild RR from `time`
#
# A recording is: time (s), rr (s), cov (data.table of the covariates with enough valid
# data: systolic, diastolic, map, resp, angle) and qc (one-row data.table of file checks).
# ==============================================================================


# ---- Checks and small helpers -----------------------------------------------------------------
#' Stop with a clear message if a configured data path is still a placeholder or missing.
check_data_path <- function(path, label, cfg) {
  if (cfg$is_placeholder(path)) {
    stop(sprintf("Set cfg$paths$%s in config.R (still '%s').", label, path), call. = FALSE)
  }
  if (!dir.exists(path)) stop(sprintf("Folder not found for '%s': %s", label, path), call. = FALSE)
  invisible(path)
}

#' "young" / "older" from the record name (cfg$fantasia regexes); NA if ambiguous.
fantasia_group <- function(ids, cfg) {
  is_young <- grepl(cfg$fantasia$young_regex, ids, ignore.case = TRUE)
  is_older <- grepl(cfg$fantasia$old_regex,   ids, ignore.case = TRUE)
  group <- ifelse(is_young & !is_older, "young", ifelse(is_older & !is_young, "older", NA_character_))
  factor(group, levels = c("young", "older"))
}


# ---- Outlier gate ----------------------------------------------------------------------------------
#' Adaptive outlier cleaning (port of clean_sample() in tools/prep-data.R): values beyond
#' `n_mad` running median absolute deviations are replaced by linear interpolation; repeated
#' `n_iter` times. Defaults are the gate of the current tools/prep-data.R (the earlier script
#' used n_mad = 2, n_iter = 5). Length is preserved.
clean_sample <- function(x, window = 50, n_mad = 3.5, n_iter = 2) {
  one_pass <- function(x) {
    was_na <- is.na(x)
    if (any(was_na)) {
      i <- seq_along(x)
      x[was_na] <- approx(i[!was_na], x[!was_na], xout = i[was_na], rule = 2)$y
    }
    k <- if (window %% 2 == 0) window + 1 else window
    k <- min(k, if (length(x) %% 2 == 1) length(x) else length(x) - 1L)
    running_median <- runmed(x, k = k)
    abs_dev <- abs(x - running_median)
    running_mad <- pmax(1.4826 * runmed(abs_dev, k = k), median(abs_dev) * 0.1)
    inside <- x >= running_median - n_mad * running_mad &
              x <= running_median + n_mad * running_mad & !was_na
    if (sum(inside) < 2) return(rep(NA_real_, length(x)))
    i <- seq_along(x)
    approx(i[inside], x[inside], xout = i, rule = 2)$y
  }
  y <- one_pass(x)
  for (pass in seq_len(max(n_iter - 1, 0))) y <- one_pass(y)
  y
}


# ---- One recording ----------------------------------------------------------------------------------
#' Read the CSV and keep the rows with a finite RR interval and time.
#' (Leading NA rows exist in files from the earlier script, whose rolling join had no RR yet.)
read_recording_rows <- function(file, io) {
  raw <- tryCatch(
    data.table::fread(file, data.table = TRUE, na.strings = c("NA", "NaN", "")),
    error = function(e) stop(sprintf("Cannot read %s: %s", basename(file), conditionMessage(e)), call. = FALSE)
  )
  if (!io$rr_col %in% names(raw)) {
    stop(sprintf("%s: column '%s' not found (columns: %s).",
                 basename(file), io$rr_col, paste(names(raw), collapse = ", ")), call. = FALSE)
  }
  rri <- suppressWarnings(as.numeric(raw[[io$rr_col]]))
  has_time <- io$time_col %in% names(raw)
  time <- if (has_time) suppressWarnings(as.numeric(raw[[io$time_col]])) else cumsum(replace(rri, is.na(rri), 0))
  if (has_time && any(diff(time[is.finite(time)]) < 0)) {
    stop(sprintf("%s: `time` is not monotonic.", basename(file)), call. = FALSE)
  }
  keep <- is.finite(rri) & is.finite(time)
  if (sum(keep) < io$min_beats) {
    stop(sprintf("%s: only %d valid beats (< %d).", basename(file), sum(keep), io$min_beats), call. = FALSE)
  }
  median_rr <- median(rri[keep])
  if (median_rr < io$rr_range_s[1] || median_rr > io$rr_range_s[2]) {
    stop(sprintf("%s: median RR = %.3f is outside %s s; RR must be in seconds.",
                 basename(file), median_rr, paste(io$rr_range_s, collapse = "-")), call. = FALSE)
  }
  list(rows = raw[keep], rri = rri[keep], time = time[keep], has_time = has_time,
       n_rows = nrow(raw), n_dropped = sum(!keep))
}

#' How well does `rri` follow `time`?
#'  consistency  share of beats with diff(time)[k] == rri[k + 1] (TMST: ~1; earlier rolling join: low)
#'  drift_s      largest gap between cumsum(rri) and the elapsed `time`
#'  dup_frac     share of consecutive identical rri values
rri_time_alignment <- function(rows, io) {
  if (!rows$has_time) return(list(consistency = 1, drift_s = 0, dup_frac = mean(diff(rows$rri) == 0)))
  list(
    consistency = mean(abs(diff(rows$time) - rows$rri[-1]) <= io$consistency_tol_s),
    drift_s     = max(abs(cumsum(rows$rri) - rows$rri[1] - (rows$time - rows$time[1]))),
    dup_frac    = mean(diff(rows$rri) == 0)
  )
}

#' "rri" or "dtime". "auto" takes rri when it is aligned with `time`, else the differences of `time`.
resolve_rr_source <- function(requested, alignment, has_time, io) {
  source <- if (is.null(requested)) "auto" else requested
  if (source != "auto") return(source)
  if (has_time && alignment$consistency < io$consistency_min) "dtime" else "rri"
}

#' RR series, beat times and the matching CSV rows for the chosen source.
build_rr_series <- function(rows, source, file) {
  if (source == "rri") return(list(rr = rows$rri, time = rows$time, rows = rows$rows))
  if (!rows$has_time) {
    stop(sprintf("%s: rr_source = 'dtime' needs a `time` column.", basename(file)), call. = FALSE)
  }
  rr <- clean_sample(diff(rows$time))                    # interval k = beat k -> beat k + 1
  finite <- is.finite(rr)
  list(rr = rr[finite], time = rows$time[-length(rows$time)][finite],
       rows = rows$rows[-nrow(rows$rows)][finite])
}

#' Covariates with at least `io$min_cov_coverage` valid values, as a data.table aligned with the RR series.
usable_covariates <- function(csv_rows, io) {
  covariates <- data.table(.row = seq_len(nrow(csv_rows)))
  for (name in intersect(io$covariates, names(csv_rows))) {
    values <- suppressWarnings(as.numeric(csv_rows[[name]]))
    if (mean(is.finite(values)) >= io$min_cov_coverage) covariates[, (name) := values]
  }
  covariates[, .row := NULL]
  covariates
}

#' Read one recording and standardize it.
#' @return list(id, time, rr, cov, qc)
read_recording <- function(file, cfg, rr_source = NULL) {
  io <- cfg$io
  rows <- read_recording_rows(file, io)
  alignment <- rri_time_alignment(rows, io)
  source <- resolve_rr_source(rr_source, alignment, rows$has_time, io)
  series <- build_rr_series(rows, source, file)

  positive <- series$rr > 0
  if (!all(positive)) {                                  # kept as before: drop non-positive intervals
    series$rr <- series$rr[positive]
    if (length(series$rr) != length(series$time)) {
      stop(sprintf("%s: non-positive intervals after cleaning.", basename(file)), call. = FALSE)
    }
  }
  covariates <- usable_covariates(series$rows, io)
  have <- names(covariates)

  qc <- data.table(
    SubjectID = tools::file_path_sans_ext(basename(file)),
    n_rows = rows$n_rows, n_na_dropped = rows$n_dropped, rr_source = source,
    drift_s = alignment$drift_s, dup_rri_frac = alignment$dup_frac,
    time_rri_consistency = alignment$consistency, rr_median = median(series$rr),
    has_bp = all(c("systolic", "diastolic") %in% have),
    has_resp = "resp" %in% have, has_angle = "angle" %in% have,
    covariates = paste(have, collapse = ",")
  )
  list(id = qc$SubjectID, time = series$time, rr = series$rr, cov = covariates, qc = qc)
}

#' Rows of a recording with time in (start, end], in seconds of the `time` column.
window_rec <- function(rec, start_s = NULL, duration_s = NULL, end_s = NULL) {
  start <- if (is.null(start_s)) -Inf else start_s
  end <- if (!is.null(end_s)) end_s else if (!is.null(duration_s) && is.finite(start)) start + duration_s else Inf
  in_window <- rec$time > start & rec$time <= end
  list(id = rec$id, time = rec$time[in_window], rr = rec$rr[in_window],
       cov = if (ncol(rec$cov)) rec$cov[in_window] else rec$cov, qc = rec$qc)
}


# ---- One dataset ----------------------------------------------------------------------------------------
#' Is covariate `name` present with enough valid values in this (windowed) recording?
has_covariate <- function(rec, name, min_coverage) {
  name %in% names(rec$cov) && mean(is.finite(rec$cov[[name]])) >= min_coverage
}

#' Load a dataset folder.
#' @return list(index = data.table, one row per recording with file checks and covariate
#'   availability in the analysis window; rec = windowed recordings (named list);
#'   rec_full = whole recordings; rr = named list of the RR vectors of `rec`).
load_dataset <- function(dataset, cfg) {
  dir <- cfg$paths[[dataset]]
  check_data_path(dir, dataset, cfg)
  files <- sort(list.files(dir, pattern = cfg$io$file_pattern, full.names = TRUE))
  if (!length(files)) stop(sprintf("No files matching '%s' in %s", cfg$io$file_pattern, dir), call. = FALSE)

  rr_source <- cfg$io$rr_source[[dataset]]
  if (is.null(rr_source)) rr_source <- "auto"
  rec_full <- lapply(files, read_recording, cfg = cfg, rr_source = rr_source)
  names(rec_full) <- vapply(rec_full, `[[`, character(1), "id")

  window <- cfg$windows[[dataset]]
  rec <- lapply(rec_full, window_rec, start_s = window$start_s, duration_s = window$duration_s)
  enough_beats <- vapply(rec, function(r) length(r$rr), integer(1)) >= cfg$io$min_beats
  if (any(!enough_beats)) {
    log_msg("%s: dropping %d recording(s) with < %d beats in the analysis window: %s", dataset,
            sum(!enough_beats), cfg$io$min_beats, paste(names(rec)[!enough_beats], collapse = ", "))
  }
  rec <- rec[enough_beats]
  rec_full <- rec_full[enough_beats]

  index <- data.table(
    Dataset = dataset, SubjectID = names(rec),
    n_beats = vapply(rec, function(r) length(r$rr), integer(1)),
    duration_s = vapply(rec, function(r) sum(r$rr), numeric(1)),
    n_beats_full = vapply(rec_full, function(r) length(r$rr), integer(1))
  )
  index <- merge(index, rbindlist(lapply(rec_full, `[[`, "qc")), by = "SubjectID", sort = FALSE)
  for (name in c("systolic", "diastolic", "map", "resp", "angle")) {
    index[, (paste0("win_", name)) := vapply(rec[SubjectID], has_covariate, logical(1),
                                              name = name, min_coverage = cfg$io$min_cov_coverage)]
  }
  if (dataset == "fantasia") {
    index[, Group := fantasia_group(SubjectID, cfg)]
    if (anyNA(index$Group)) {
      stop("Fantasia: could not assign young/older from file names for: ",
           paste(index[is.na(Group), SubjectID], collapse = ", "), ". Adjust cfg$fantasia regexes.", call. = FALSE)
    }
  }
  list(index = index, rec = rec, rec_full = rec_full, rr = lapply(rec, `[[`, "rr"))
}

#' Protocol annotations of the PRCP records (written by tools/prep-data.R): data.table(id, kind, onset_s, onset_end_s, return_s,
#' return_end_s, label), one row per maneuver; NULL when the file does not exist (analysis from the tilt angle and heart rate only).
read_prcp_events <- function(cfg) {
  file <- cfg$paths$prcp_events
  if (is.null(file) || !file.exists(file)) return(NULL)
  events <- data.table::fread(file, colClasses = list(character = c("id", "kind", "label")))
  needed <- c("id", "kind", "onset_s", "return_s")
  if (!all(needed %in% names(events))) stop(sprintf("%s lacks the columns: %s", file, paste(setdiff(needed, names(events)), collapse = ", ")), call. = FALSE)
  for (column in c("onset_end_s", "return_end_s")) if (!column %in% names(events)) events[[column]] <- NA_real_
  events[, id := as.character(id)]
  events
}
