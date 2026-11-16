# ==============================================================================
# tools/prep/blockade.R -- prepares the Lancaster "Blocker data" (doi 10.17635/lancaster/researchdata/280, CC BY-NC) for analysis/datasets/15_blockade.R.
#   Rscript tools/prep/blockade.R            (from the project root; `make prep` runs it too)
# The zip (about 515 MB) must be downloaded by hand into cfg$paths$blockade_zip: the repository's server refuses scripted downloads behind a browser
# challenge, which this project does not try to bypass. Without the zip the script says so and exits without error (the blockade analysis is optional).
# Everything is read straight from the zip with Python's standard zipfile (nothing is extracted or executed): the experiment logs (segment table and
# maneuver times, parsed from their printable text), the ECG (R-peak detection) and the integrated nerve signal (1-s means).
# Writes cfg$paths$blockade: protocol.csv, maneuvers.csv, rpeaks/rpeaks_Subject<s>_day<d>.csv, rpeaks/qc.csv, msna/msna_Subject<s>_day<d>.csv
# Needs Python 3 with numpy, pandas and scipy (the project's .venv, or SPARK_PYTHON).
# ==============================================================================
if (!file.exists("config.R")) stop("Run from the project root (the folder containing config.R).", call. = FALSE)
source("config.R")
zip <- cfg$paths$blockade_zip; out <- cfg$paths$blockade
if (!file.exists(zip)) {
  message("Blockade data not found (", zip, "). Download 'Blocker data' by hand from https://doi.org/10.17635/lancaster/researchdata/280 and save the zip there; skipping.")
  quit(status = 0)
}
# Python: SPARK_PYTHON if set, else the project's .venv, else python3 on PATH (on macOS /usr/bin comes first in PATH for the C++ compiler, and Apple's python3
# has no numpy, hence the explicit choice)
python <- Sys.getenv("SPARK_PYTHON", unset = if (file.exists(".venv/bin/python3")) ".venv/bin/python3" else Sys.which("python3"))
if (!nzchar(python)) stop("python3 not found", call. = FALSE)
ok <- system2(python, c("-I", "-c", shQuote("import numpy, pandas, scipy")), stdout = FALSE, stderr = FALSE) == 0
if (!ok) stop("The Python at ", python, " lacks numpy/pandas/scipy: set SPARK_PYTHON to one that has them.", call. = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
run <- function(script, args) {
  status <- system2(python, c("-I", file.path("tools/prep", script), shQuote(args)))
  if (status != 0) stop(script, " failed", call. = FALSE)
}
run("blockade_protocol.py", c(zip, file.path(out, "protocol.csv")))
run("blockade_maneuvers.py", c(zip, file.path(out, "maneuvers.csv")))
run("blockade_rpeaks.py", c(zip, file.path(out, "rpeaks")))
run("blockade_msna.py", c(zip, file.path(out, "msna")))
n_peaks <- length(list.files(file.path(out, "rpeaks"), pattern = "^rpeaks_.*\\.csv$")); n_msna <- length(list.files(file.path(out, "msna"), pattern = "\\.csv$"))
message(sprintf("Blockade data prepared in %s: %d segments, %d sessions with R peaks, %d with the nerve signal", out, nrow(read.csv(file.path(out, "protocol.csv"))), n_peaks, n_msna))
