# Compares result files of two runs numerically (tolerance 1e-8) and reports differences.
#   Rscript tools/checks/compare_results.R <baseline_results_dir> <new_results_dir> [name1 name2 ...]
args <- commandArgs(trailingOnly = TRUE)
base_dir <- args[1]; new_dir <- args[2]
files <- if (length(args) > 2) paste0(args[-(1:2)], ".rds") else sort(list.files(base_dir, pattern = "\\.rds$"))
n_diff <- 0L
for (f in files) {
  a <- file.path(base_dir, f); b <- file.path(new_dir, f)
  if (!file.exists(b)) { cat(sprintf("%-28s MISSING in new run\n", f)); n_diff <- n_diff + 1L; next }
  res <- all.equal(readRDS(a), readRDS(b), tolerance = 1e-8, check.attributes = FALSE)
  ok <- isTRUE(res)
  cat(sprintf("%-28s %s\n", f, if (ok) "identical (tol 1e-8)" else "DIFFERENT"))
  if (!ok) { n_diff <- n_diff + 1L; cat(paste0("    ", head(res, 6)), sep = "\n") }
}
if (n_diff) quit(status = 1)
cat("All compared result files match the baseline.\n")
