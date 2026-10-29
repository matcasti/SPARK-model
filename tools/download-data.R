# ==============================================================================
# tools/download-data.R -- fetches the two public datasets from PhysioNet into data-raw/ and verifies every file
#   Rscript tools/download-data.R              # both datasets
#   Rscript tools/download-data.R prcp         # one of them (fantasia | prcp)
#   Rscript tools/download-data.R --check      # only verify what is already in data-raw/
# Files that are already present with the right SHA-256 are kept, so the script can be re-run to resume.
# The exercise cohort (data-raw/tmst/) is not on PhysioNet; it is part of this repository.
# Both datasets are Open Data Commons Attribution (ODC-By 1.0): Fantasia, doi:10.13026/C2RG61; PRCP, doi:10.13026/C2PC7R.
# ==============================================================================
DATASETS <- list(
  fantasia = list(url = "https://physionet.org/files/fantasia/1.0.0/", ext = c("dat", "hea", "ecg")),
  prcp     = list(url = "https://physionet.org/files/prcp/1.0.0/",     ext = c("dat", "hea", "anI", "wabp", "wqrs"))
)
RAW_DIR <- "data-raw"

sha256 <- function(path) {
  f <- get0("sha256sum", envir = asNamespace("tools"), inherits = FALSE)          # R >= 4.5
  if (!is.null(f)) return(unname(f(path)))
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(path, algo = "sha256", file = TRUE))
  out <- system2(if (nzchar(Sys.which("sha256sum"))) "sha256sum" else "shasum", if (nzchar(Sys.which("sha256sum"))) path else c("-a", "256", path), stdout = TRUE)
  strsplit(out, " ")[[1]][1]
}

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
chosen <- setdiff(args, "--check")
if (!length(chosen)) chosen <- names(DATASETS)
bad <- setdiff(chosen, names(DATASETS)); if (length(bad)) stop("unknown dataset: ", paste(bad, collapse = ", "), " (choose from ", paste(names(DATASETS), collapse = ", "), ")")
options(timeout = max(600, getOption("timeout")))

failed <- 0L
for (name in chosen) {
  spec <- DATASETS[[name]]; dir <- file.path(RAW_DIR, name); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  sums_file <- tempfile(fileext = ".txt")
  status <- try(download.file(paste0(spec$url, "SHA256SUMS.txt"), sums_file, mode = "wb", quiet = TRUE), silent = TRUE)
  if (inherits(status, "try-error") || status != 0) { message(name, ": cannot read ", spec$url, "SHA256SUMS.txt"); failed <- failed + 1L; next }
  lines <- trimws(readLines(sums_file, warn = FALSE)); lines <- lines[nzchar(lines)]
  sums <- data.frame(hash = sub(" .*$", "", lines), file = trimws(sub("^[0-9a-f]+ +\\*?", "", lines)), stringsAsFactors = FALSE)
  sums <- sums[!grepl("/", sums$file) & tolower(tools::file_ext(sums$file)) %in% tolower(spec$ext), ]   # top-level WFDB files only
  n_ok <- n_new <- 0L
  for (i in seq_len(nrow(sums))) {
    path <- file.path(dir, sums$file[i])
    if (file.exists(path) && identical(sha256(path), sums$hash[i])) { n_ok <- n_ok + 1L; next }
    if (check_only) { message(name, ": ", sums$file[i], if (file.exists(path)) " does not match its checksum" else " is missing"); failed <- failed + 1L; next }
    ok <- try(download.file(paste0(spec$url, sums$file[i]), path, mode = "wb", quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error") || ok != 0 || !identical(sha256(path), sums$hash[i])) { message(name, ": ", sums$file[i], " failed"); unlink(path); failed <- failed + 1L } else n_new <- n_new + 1L
  }
  message(sprintf("%s: %d files verified, %d downloaded (%d expected) in %s", name, n_ok, n_new, nrow(sums), dir))
}
if (failed) stop(failed, " file(s) missing or corrupt; re-run the script")
