# ==============================================================================
# 00_setup.R -- sourced first by every analysis script.
# Loads: packages -> config.R (`cfg`) -> the model engine -> every function in R/ (R/load.R)
#        -> the compiled filter (optional; falls back to the engine's R filter).
# Scripts must be run from the project root (the folder that contains config.R).
# ==============================================================================

if (!file.exists("config.R")) stop("Run from the project root (the folder containing config.R).", call. = FALSE)
source("config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  for (pkg in c("biwavelet", "numDeriv", "MASS", "Rcpp", "survival", "ggtext")) library(pkg, character.only = TRUE)
})
set.seed(cfg$seed)

source(cfg$paths$engine)                                              # estimation engine (unmodified)
source("R/load.R"); load_project_code()                               # every function in R/, layer by layer
invisible(load_fast_filter(cfg))

ensure_dir(cfg$paths$results)
ensure_dir(cfg$paths$figures)
log_msg("Setup done | cores = %d | C++ filter = %s", cfg$cores, has_cpp())
