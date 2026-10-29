# ==============================================================================
# tests/run_tests.R -- runs the test suites and prints one line per suite.
#   Rscript tests/run_tests.R            every suite except the end-to-end run
#   Rscript tests/run_tests.R --quick    structure, text generators and unit checks only (a few minutes)
#   Rscript tests/run_tests.R --e2e      also the whole pipeline on synthetic data (tests/test_pipeline.R with SPARK_E2E=1; ~10 min)
# Exit status 1 when a suite fails; the tail of its output is shown. Run from the project root.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
suites <- data.frame(
  script = c("tests/test_workflow.R", "tests/test_narrative.R", "tests/test_pipeline.R", "tests/test_prep.R", "tests/test_edge_cases.R"),
  what = c("structure, pipeline.R vs scripts, docs", "manuscript text generators", "unit checks (filter, simulator, events, loader, ...)",
           "tools/prep-data.R on synthetic WFDB records", "scripts with absent covariates + manuscript chunks"),
  quick = c(TRUE, TRUE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
if ("--quick" %in% args) suites <- suites[suites$quick, ]
env <- if ("--e2e" %in% args) "SPARK_E2E=1" else character()

results <- lapply(seq_len(nrow(suites)), function(i) {
  log <- tempfile(fileext = ".log"); started <- Sys.time()
  status <- system2("Rscript", shQuote(suites$script[i]), stdout = log, stderr = log, env = env)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  skipped <- status == 0 && any(grepl("skipped", readLines(log, warn = FALSE), ignore.case = TRUE))
  cat(sprintf("%-8s %-26s %5.0f s  %s\n", if (status != 0) "FAILED" else if (skipped) "skipped" else "ok", basename(suites$script[i]), seconds, suites$what[i]))
  if (status != 0) cat(paste0("         ", tail(readLines(log, warn = FALSE), 8)), sep = "\n")
  status
})
failed <- sum(unlist(results) != 0)
cat(sprintf("\n%d of %d suites passed.\n", nrow(suites) - failed, nrow(suites)))
if (failed) quit(status = 1)
