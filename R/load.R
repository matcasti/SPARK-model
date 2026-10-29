# ==============================================================================
# R/load.R -- loads the project's functions, layer by layer.
#
#   R/core/     helpers every other file may use: logging, results store, CSV loader, statistics
#   R/model/    the model: specifications M0-M3, fits, identifiability, simulator, synthetic recordings
#   R/signals/  what is computed from RR series and states: events, physiology, HRV, summaries, landscapes
#   R/steps/    the per-recording steps of each analysis script (analysis/<group>/NN_*.R only call these)
#   R/report/   notation, plots, figures and manuscript text built from the saved results
#
# A layer may use the layers before it, not the ones after it. Within a layer files load alphabetically,
# and no file runs code at load time beyond defining functions and constants, so the order inside a layer
# does not matter. tests/suites/test_workflow.R checks that every file of R/ sits in a layer and that no function is
# defined twice.
# ==============================================================================
PROJECT_LAYERS <- c("core", "model", "signals", "steps", "report")

#' Source the functions of the given layers (default: all) in layer order.
load_project_code <- function(layers = PROJECT_LAYERS, root = "R") {
  for (layer in layers) {
    files <- sort(list.files(file.path(root, layer), pattern = "\\.R$", full.names = TRUE))
    for (file in files) source(file)
  }
  invisible(NULL)
}
