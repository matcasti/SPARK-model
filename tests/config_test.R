# Test overrides: synthetic data, small workloads, separate output folders.
root <- Sys.getenv("SPARK_TEST_ROOT", unset = "tests/_run")
cfg$paths$prcp <- file.path(root, "data/prcp"); cfg$paths$prcp_events <- file.path(root, "data/prcp_events.csv"); cfg$paths$fantasia <- file.path(root, "data/fantasia")
cfg$paths$exercise <- file.path(root, "data/tmst")
cfg$paths$results <- file.path(root, "results"); cfg$paths$figures <- file.path(root, "figures")
cfg$cores <- 1L
cfg$windows$fantasia <- list(start_s = 60, duration_s = 360, stability_A = c(60, 240), stability_B = c(240, 420))
cfg$identifiability <- modifyList(cfg$identifiability, list(n_starts = 4L, max_recordings = 2L, n_surface_recordings = 1L, surface_grid = 9L))
cfg$sensitivity$max_recordings <- 2L
cfg$simulation <- modifyList(cfg$simulation, list(n_per_scenario = 3L, duration_s = 300, n_detect = 6L))
cfg$calibration <- modifyList(cfg$calibration, list(n_sims = 6L, duration_s = 300))
cfg$stats$n_boot <- 200L
cfg$recovery_benchmark <- list(n = 8L)
