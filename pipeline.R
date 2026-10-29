# ==============================================================================
# pipeline.R -- the analysis pipeline as data: which script runs, what it reads, what it writes.
#
# One source of truth for
#   * tools/run_all.R      run order, stale detection, check that each step wrote what it declares,
#   * docs/DATAFLOW.md     the table of scripts (tools/make_docs.R),
#   * tests/test_workflow  the declarations must match the read_result()/write_result() calls in the scripts.
# config.R says WHAT is analysed and how; this file says IN WHICH ORDER and from what to what.
#
# reads   "data"  the CSV recordings (data/);   "name"  results/name.rds is required;
#         "name?" it is used when present;       "*"     every result that exists (the figures)
# writes  results/<name>.rds;   files  other files the step writes ({results} and {figures} are the folders of config.R)
# ==============================================================================
step <- function(id, script, title, reads, writes, used_in = "", files = character()) {
  list(id = id, script = file.path("analysis", script), title = title, reads = reads, writes = writes, used_in = used_in, files = files)
}

PIPELINE <- list(
  step("01", "01_fit_all.R", "Fit every recording with M3",
       reads = "data", writes = c("config_used", "fit_prcp", "fit_fantasia", "fit_exercise", "df_population_summary"),
       files = "{results}/population_summary.csv", used_in = "Fit and innovation diagnostics; time constants (Figure 2)"),
  step("02", "02_simulation_recovery.R", "Parameter and state recovery under known truth",
       reads = character(), writes = "sim_recovery", used_in = "Recovery under known truth (Figure S1)"),
  step("03", "03_identifiability.R", "Multi-start identifiability, M0-M3",
       reads = "data", writes = "identifiability", used_in = "Results 1-2 (Figure 1, Table 4)"),
  step("04", "04_sensitivity.R", "Sensitivity to hyperparameters, prior vs posterior",
       reads = "data", writes = "sensitivity", used_in = "Sensitivity, prior contrast (Figure S2)"),
  step("05", "05_prcp.R", "Orthostatic perturbation (PRCP)",
       reads = c("fit_prcp", "data"), writes = "prcp", used_in = "Result 3: tilts and stand-ups, pressure response, lagged correlation, dose-response, responses by maneuver (Figures 3 and 4; S6 and S7)"),
  step("06", "06_fantasia.R", "Aging, respiration, baroreflex (Fantasia)",
       reads = c("fit_fantasia", "data"), writes = "fantasia", used_in = "Result 4: group contrasts, stability, respiration and baroreflex (Figures 5 and 6)"),
  step("07", "07_recovery_benchmark.R", "Choice of the recovery definition",
       reads = character(), writes = "recovery_choice", used_in = "Recovery definition (Figure S3)"),
  step("08", "08_event_detection_validation.R", "RR-based event detection on simulations",
       reads = character(), writes = "event_detection", used_in = "Event detection (Figure S5)"),
  step("09", "09_exercise.R", "Exercise and recovery (TMST)",
       reads = c("fit_exercise", "recovery_choice?"), writes = "exercise", used_in = "Result 5 (Figure 7; S6)"),
  step("10", "10_hrv_comparators.R", "SDNN / RMSSD / LF/HF comparators",
       reads = c("fit_prcp", "fit_fantasia", "fit_exercise", "prcp?", "exercise?"), writes = "hrv", used_in = "HRV comparators"),
  step("11", "11_supplement_checks.R", "Bound, calibration, equivalence checks",
       reads = c("fit_prcp?", "fit_fantasia?", "fit_exercise?"), writes = "supp_checks", used_in = "Bound, calibration, equivalence checks (Figure S4)"),
  step("12", "12_preprocessing_checks.R", "Data-preparation checks",
       reads = character(), writes = "prep_checks", used_in = "Data-preparation checks (earlier vs current format, cleaning gate)"),
  step("13", "13_figures.R", "All figures",
       reads = "*", writes = "figure_info", files = "{figures}", used_in = "All figures (look and examples: `cfg$figures`)")
)

pipeline_ids <- function() vapply(PIPELINE, `[[`, "", "id")
pipeline_step <- function(id) PIPELINE[[match(id, pipeline_ids())]]

#' Names of the results a step reads (`which`: "all", "required" or "optional"); `*` expands to every result an earlier step writes.
pipeline_reads <- function(step, which = c("all", "required", "optional")) {
  which <- match.arg(which)
  reads <- step$reads
  if (identical(reads, "*")) reads <- paste0(unique(unlist(lapply(PIPELINE[seq_len(match(step$id, pipeline_ids()) - 1L)], `[[`, "writes"))), "?")
  reads <- setdiff(reads, "data")
  optional <- endsWith(reads, "?")
  names <- sub("\\?$", "", reads)
  switch(which, all = names, required = names[!optional], optional = names[optional])
}

#' Files a step produces (project-relative paths; `{results}` and `{figures}` in `files` come from cfg$paths).
pipeline_outputs <- function(step, cfg) {
  files <- sub("{results}", cfg$paths$results, step$files, fixed = TRUE)
  c(file.path(cfg$paths$results, paste0(step$writes, ".rds")), sub("{figures}", cfg$paths$figures, files, fixed = TRUE))
}

#' Files a step depends on, by kind: results of earlier steps, CSV recordings and code.
pipeline_inputs <- function(step, cfg) {
  results <- file.path(cfg$paths$results, paste0(pipeline_reads(step), ".rds"))
  data <- if ("data" %in% step$reads) unlist(lapply(unname(cfg$paths[c("prcp", "fantasia", "exercise")]), list.files, pattern = "\\.csv$", full.names = TRUE))
  code <- c(step$script, "analysis/00_setup.R", "config.R", "pipeline.R", list.files(c("R", "engine"), pattern = "\\.(R|cpp)$", recursive = TRUE, full.names = TRUE))
  list(results = results[file.exists(results)], data = data, code = code)
}

#' Is each step up to date? One row per step: state is "missing" (an output does not exist), "stale" (an input is newer than
#' the oldest output) or "fresh". `because` names the newest input that made it stale.
pipeline_status <- function(cfg) {
  rows <- lapply(PIPELINE, function(step) {
    outputs <- pipeline_outputs(step, cfg)
    exists <- file.exists(outputs) | dir.exists(outputs)
    if (!all(exists)) return(data.frame(id = step$id, state = "missing", because = basename(outputs[!exists][1]), stringsAsFactors = FALSE))
    oldest <- min(file.mtime(outputs[!dir.exists(outputs)]))
    inputs <- unlist(pipeline_inputs(step, cfg), use.names = FALSE)
    newer <- inputs[file.mtime(inputs) > oldest]
    data.frame(id = step$id, state = if (length(newer)) "stale" else "fresh",
               because = if (length(newer)) newer[which.max(file.mtime(newer))] else "", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$script <- vapply(PIPELINE, `[[`, "", "script")
  out
}
