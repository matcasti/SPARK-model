# ==============================================================================
# tools/checks/check_env.R -- is this machine able to run the project?
#   Rscript tools/checks/check_env.R          prints the report; exit status 1 when a required package is missing
#   source("tools/checks/check_env.R"); check_environment()   the same, returning the table invisibly
#
#   required     the analysis scripts stop without them
#   recommended  figures and speed: ragg/systemfonts/svglite draw the symbols, a C++ compiler runs the multi-start fits ~10x faster
#   preparation  only tools/prep-data.R (raw WFDB files -> CSV): Python `wfdb` through reticulate
#   rendering    only the manuscript: Quarto
# ==============================================================================
ENV_PACKAGES <- list(
  required    = c("data.table", "ggplot2", "patchwork", "ggtext", "biwavelet", "numDeriv", "MASS", "Rcpp", "survival", "Matrix", "scales", "viridisLite"),
  recommended = c("ragg", "systemfonts", "svglite"),
  preparation = c("reticulate"),
  rendering   = c("knitr", "rmarkdown")
)

check_environment <- function(verbose = TRUE) {
  installed <- rownames(installed.packages())
  rows <- do.call(rbind, lapply(names(ENV_PACKAGES), function(group)
    data.frame(group = group, item = ENV_PACKAGES[[group]], ok = ENV_PACKAGES[[group]] %in% installed, stringsAsFactors = FALSE)))
  rows$version <- vapply(seq_len(nrow(rows)), function(i) if (rows$ok[i]) as.character(packageVersion(rows$item[i])) else "", "")

  has_cpp <- nzchar(Sys.which("clang++")) || nzchar(Sys.which("g++"))
  has_wfdb <- rows$ok[rows$item == "reticulate"] && tryCatch(reticulate::py_module_available("wfdb"), error = function(e) FALSE)
  has_quarto <- nzchar(Sys.which("quarto"))
  tools_rows <- data.frame(group = c("recommended", "preparation", "rendering"), item = c("C++ compiler", "Python wfdb", "quarto"),
                           ok = c(has_cpp, has_wfdb, has_quarto), version = "", stringsAsFactors = FALSE)
  rows <- rbind(rows, tools_rows)

  if (verbose) {
    cat(sprintf("R %s | %s\n", getRversion(), R.version$platform))
    for (group in names(ENV_PACKAGES)) {
      part <- rows[rows$group == group, ]
      cat(sprintf("%-12s %s\n", group, paste0(part$item, if (any(!part$ok)) ifelse(part$ok, "", " (MISSING)") else "", collapse = ", ")))
    }
  }
  invisible(rows)
}

if (sys.nframe() == 0L) {
  report <- check_environment()
  if (any(!report$ok[report$group == "required"])) quit(status = 1)
}
