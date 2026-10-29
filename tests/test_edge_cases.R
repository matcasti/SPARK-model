# Runs the analysis scripts on datasets in which covariates are absent everywhere, to check that
# empty tables degrade gracefully (the case that broke script 06: no blood pressure anywhere,
# a single-step tilt protocol). Also renders the R chunks of the manuscript and supplement.
#   Rscript tests/test_edge_cases.R
cases <- c("nobp_step", "noangle", "noevent")
scripts <- c("analysis/01_fit_all.R", "analysis/05_prcp.R", "analysis/06_fantasia.R", "analysis/09_exercise.R", "analysis/10_hrv_comparators.R", "analysis/13_figures.R")
cfg_file <- normalizePath("tests/config_test.R"); fails <- character()
for (cs in cases) {
  root <- file.path(tempfile("edge"), cs); dir.create(root, recursive = TRUE)
  env <- c(SPARK_TEST_ROOT = root, SPARK_CONFIG = cfg_file, EDGE_CASE = cs)
  cat(sprintf("\n=== case '%s' (%s)\n", cs, root))
  st <- system2("Rscript", "tests/make_edge_data.R", env = paste0(names(env), "=", shQuote(env)), stdout = FALSE, stderr = FALSE)
  if (st != 0) { fails <- c(fails, paste(cs, "make_edge_data")); next }
  for (sc in scripts) {
    log <- file.path(root, paste0(basename(sc), ".log"))
    st <- system2("Rscript", sc, env = paste0(names(env), "=", shQuote(env)), stdout = log, stderr = log)
    cat(sprintf("  %-34s %s\n", basename(sc), if (st == 0) "ok" else "FAILED"))
    if (st != 0) { fails <- c(fails, paste(cs, basename(sc))); cat(paste("   ", tail(readLines(log), 4)), sep = "\n") }
  }
  # render the R chunks of the manuscript and the supplement against these results
  file.symlink(normalizePath("R"), file.path(root, "R")); file.symlink(normalizePath("engine"), file.path(root, "engine"))   # the documents source ../R
  md <- file.path(root, "manuscript"); dir.create(md); file.copy(c("manuscript/manuscript.qmd", "manuscript/supplementary.qmd"), md)
  file.copy(file.path(root, "results", "df_population_summary.rds"), md); file.copy("engine/model_engine.R", file.path(md, "model_code.R"))
  for (f in c("manuscript", "supplementary")) {
    r <- tryCatch({ old <- setwd(md); on.exit(setwd(old), add = TRUE)
      suppressWarnings(knitr::knit(paste0(f, ".qmd"), output = paste0(f, ".md"), quiet = TRUE)); setwd(old); "ok" }, error = function(e) { setwd(old); conditionMessage(e) })
    cat(sprintf("  render %-14s %s\n", f, if (identical(r, "ok")) "ok" else paste("FAILED:", substr(r, 1, 120))))
    if (!identical(r, "ok")) fails <- c(fails, paste(cs, "render", f))
  }
}
if (length(fails)) stop("Edge-case failures: ", paste(fails, collapse = "; "), call. = FALSE)
cat("\nAll edge cases passed.\n")
