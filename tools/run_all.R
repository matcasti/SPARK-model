# ==============================================================================
# tools/run_all.R -- runs the analysis pipeline (steps and their inputs/outputs: pipeline.R).
#
#   Rscript tools/run_all.R                  every step, in order
#   Rscript tools/run_all.R --list           the steps and what each reads / writes; runs nothing
#   Rscript tools/run_all.R --status         which steps are fresh, stale or missing; runs nothing
#   Rscript tools/run_all.R --stale          only the steps that are stale or missing
#   Rscript tools/run_all.R --only=05,13     only these steps          (--from=05  --to=09  --skip=02)
#   Rscript tools/run_all.R --dry-run        with any of the above: print what would run
#
# After each step the declared outputs are checked (they must exist and be newer than the start of the step).
# A full or partial run writes to results/:
#   run_log.csv     step, status, seconds, time of the last run of each step
#   provenance.txt  git commit and state of the working tree, R and package versions, md5 of config.R, the engine and the data
#   MANIFEST.csv    every file of results/ with size and md5
#   sessionInfo.txt
# Run from the project root.
# ==============================================================================
if (!file.exists("config.R")) stop("Run from the project root (the folder containing config.R).", call. = FALSE)
source("config.R"); source("pipeline.R"); source("tools/check_env.R")

args <- commandArgs(trailingOnly = TRUE)
option <- function(name) { hit <- grep(paste0("^--", name, "="), args, value = TRUE); if (length(hit)) strsplit(sub(".*=", "", hit[1]), ",")[[1]] else NULL }
flag <- function(name) paste0("--", name) %in% args
unknown <- setdiff(sub("=.*", "", args), c("--list", "--status", "--stale", "--dry-run", "--only", "--from", "--to", "--skip"))
if (length(unknown)) stop("Unknown option(s): ", paste(unknown, collapse = " "), call. = FALSE)

check_ids <- function(ids) { bad <- setdiff(ids, pipeline_ids()); if (length(bad)) stop("Unknown step(s): ", paste(bad, collapse = ", "), call. = FALSE); ids }
selected <- pipeline_ids()
if (!is.null(option("only")))  selected <- intersect(selected, check_ids(option("only")))
if (!is.null(option("from")))  selected <- selected[selected >= check_ids(option("from"))]
if (!is.null(option("to")))    selected <- selected[selected <= check_ids(option("to"))]
if (!is.null(option("skip")))  selected <- setdiff(selected, check_ids(option("skip")))

# ---- Listing and status (nothing runs) -------------------------------------------------------------------------
if (flag("list")) {
  for (step in PIPELINE) {
    cat(sprintf("%s  %-34s %s\n", step$id, basename(step$script), step$title))
    cat(sprintf("    reads  %s\n    writes %s\n", paste(step$reads, collapse = ", "), paste(c(step$writes, step$files), collapse = ", ")))
  }
  quit(status = 0)
}
status <- pipeline_status(cfg)
if (flag("status")) {
  print(status[, c("id", "script", "state", "because")], row.names = FALSE)
  quit(status = 0)
}
if (flag("stale")) selected <- intersect(selected, status$id[status$state != "fresh"])
if (!length(selected)) { cat("Nothing to run.\n"); quit(status = 0) }
if (flag("dry-run")) { cat("Would run:", paste(selected, collapse = " "), "\n"); quit(status = 0) }

# ---- Run ---------------------------------------------------------------------------------------------------------------
report <- check_environment(verbose = FALSE)
if (any(!report$ok[report$group == "required"])) stop("Missing required packages: ", paste(report$item[report$group == "required" & !report$ok], collapse = ", "), call. = FALSE)

# A step that reads a result no step in this run produces needs that result to exist already
for (id in selected) {
  step <- pipeline_step(id)
  missing <- setdiff(pipeline_reads(step, "required"), c(unlist(lapply(PIPELINE[match(selected[selected < id], pipeline_ids())], `[[`, "writes")),
                                                          sub("\\.rds$", "", list.files(cfg$paths$results, pattern = "\\.rds$"))))
  if (length(missing)) stop(sprintf("Step %s needs results/%s.rds: run the step that writes it first.", id, paste(missing, collapse = ".rds, results/")), call. = FALSE)
}

log_rows <- list()
for (id in selected) {
  step <- pipeline_step(id)
  cat(sprintf("\n===== %s  %s =====\n", step$script, step$title))
  started <- Sys.time()
  exit_status <- system2("Rscript", shQuote(step$script))
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  outputs <- pipeline_outputs(step, cfg)
  written <- file.exists(outputs) & (dir.exists(outputs) | file.mtime(outputs) >= started - 1)
  state <- if (exit_status != 0) "failed" else if (!all(written)) "outputs missing" else "ok"
  log_rows[[id]] <- data.frame(step = id, script = step$script, status = state, seconds = round(seconds, 1), finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  if (state != "ok") {
    detail <- if (exit_status != 0) "the script stopped with an error" else paste("not written by this run:", paste(outputs[!written], collapse = ", "))
    stop(sprintf("Step %s failed (%s).", step$script, detail), call. = FALSE)
  }
}

# ---- Record the run ------------------------------------------------------------------------------------------------------
ensure <- function(path) dir.create(path, showWarnings = FALSE, recursive = TRUE)
ensure(cfg$paths$results)
md5 <- function(files) unname(tools::md5sum(files))
git <- function(...) tryCatch(system2("git", c(...), stdout = TRUE, stderr = FALSE), error = function(e) character(), warning = function(w) character())
changed <- git("status", "--porcelain", "--", "R", "analysis", "engine", "config.R", "pipeline.R", "tools")
data_files <- sort(list.files(unlist(cfg$paths[c("prcp", "fantasia", "exercise")]), pattern = "\\.csv$", full.names = TRUE))
config_files <- c("config.R", Sys.getenv("SPARK_CONFIG", unset = ""))
config_files <- config_files[nzchar(config_files) & file.exists(config_files)]
packages <- c("data.table", "ggplot2", "patchwork", "ggtext", "biwavelet", "numDeriv", "MASS", "Rcpp", "survival", "ragg")
provenance <- c(
  sprintf("run_finished: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("steps_run: %s", paste(selected, collapse = " ")),
  sprintf("git_commit: %s", paste(git("rev-parse", "HEAD"), collapse = "")),
  sprintf("git_code_files_changed_since_commit: %d", length(changed)),
  sprintf("R: %s | %s", R.version.string, R.version$platform),
  sprintf("packages: %s", paste0(packages, " ", vapply(packages, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "-", ""), collapse = ", ")),
  sprintf("seed: %s", cfg$seed),
  sprintf("config_md5: %s", paste(md5(config_files), collapse = " ")),
  sprintf("engine_md5: %s", md5(cfg$paths$engine)),
  sprintf("data_files: %d | data_md5: %s", length(data_files), if (length(data_files)) { h <- tempfile(); writeLines(md5(data_files), h); md5(h) } else "-"))
writeLines(provenance, file.path(cfg$paths$results, "provenance.txt"))

# run_log.csv keeps the last run of every step: steps not run now keep their previous row
log_file <- file.path(cfg$paths$results, "run_log.csv")
previous <- if (file.exists(log_file)) read.csv(log_file, colClasses = "character") else NULL
previous <- previous[!previous$step %in% selected, , drop = FALSE]
run_log <- rbind(previous, do.call(rbind, lapply(log_rows, function(r) { r$seconds <- as.character(r$seconds); r })))
write.csv(run_log[order(run_log$step), ], log_file, row.names = FALSE)

result_files <- sort(list.files(c(cfg$paths$results, cfg$paths$figures), recursive = TRUE, full.names = TRUE))
result_files <- result_files[!basename(result_files) %in% c("MANIFEST.csv", ".DS_Store")]
write.csv(data.frame(file = result_files, bytes = file.size(result_files), md5 = md5(result_files)), file.path(cfg$paths$results, "MANIFEST.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$results, "sessionInfo.txt"))

cat(sprintf("\nSteps run: %s | %.1f min | provenance, run_log and MANIFEST written to %s/\n", paste(selected, collapse = " "), sum(as.numeric(run_log$seconds[run_log$step %in% selected])) / 60, cfg$paths$results))
cat("All steps completed.\n")
