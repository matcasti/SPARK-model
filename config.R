# ==============================================================================
# config.R -- every setting of the analysis, in one place. Sourced by all scripts.
# The defaults reproduce the manuscript; edit the data paths if your folders differ.
#
#   1. General         seed, cores, datasets, which fit path
#   2. Data            folders, CSV layout, analysis windows, Fantasia groups
#   3. Event definitions   tilt angle, heart-rate segmentation
#   4. Analyses        physiology, exercise, HRV, statistics
#   5. Workload        multi-start, sensitivity, simulations (reduce for a quick pass)
#   6. Figures         look, palettes, which recording to showcase, landscape settings
# ==============================================================================
cfg <- list(

  # ---- 1. General ---------------------------------------------------------------------------------
  seed     = 20261004L,
  cores    = max(1L, parallel::detectCores() - 1L),
  datasets = c("prcp", "fantasia", "exercise"),


  # ---- 2. Data ----------------------------------------------------------------------------------------
  paths = list(
    engine   = "engine/model_engine.R",
    fast     = "engine/fast_filter.cpp",
    results  = "results",
    figures  = "manuscript/figures",
    # Folders written by tools/prep-data.R: one .csv per recording
    prcp     = "data/prcp",
    prcp_events = "data/prcp_events.csv",   # protocol annotations of the PRCP records (tools/prep-data.R); optional
    fantasia = "data/fantasia",
    exercise = "data/tmst"          # two-minute step test (rest - exercise - rest)
  ),

  # CSV layout written by tools/prep-data.R: time (s), rri (s) and, when available, systolic,
  # diastolic, map (blood pressure), resp (Fantasia) and angle (PRCP tilt table).
  io = list(
    file_pattern = "\\.csv$",
    time_col     = "time",
    rr_col       = "rri",
    # How the RR series is taken from a file, per dataset:
    #   "rri"   the rri column as written (ECG RR);
    #   "dtime" successive differences of `time`, cleaned with the same gate as tools/prep-data.R;
    #   "auto"  "rri" if rri is the interval ending at `time` (diff(time) matches rri in at least
    #           `consistency_min` of beats; true for files from the current tools/prep-data.R),
    #           else "dtime" (files from the earlier prep script, whose rolling join left
    #           repeated / skipped rri values; there `time` for Fantasia are ECG R peaks).
    # PRCP is fixed to "rri": in files from the earlier script `time` was the ABP pulse onset
    # (diff(time) = pulse interval, not RR), and the current script writes ECG beats anyway.
    rr_source         = list(prcp = "rri", fantasia = "auto", exercise = "rri"),
    consistency_tol_s = 0.0025,
    consistency_min   = 0.90,
    rr_range_s        = c(0.2, 3.0),   # plausibility range of the median RR; guards against ms files
    min_beats         = 60L,
    covariates        = c("systolic", "diastolic", "map", "resp", "angle"),
    min_cov_coverage  = 0.8            # minimum share of valid values for a covariate to be used
  ),

  # Analysis windows, in seconds of the `time` column (NULL = whole recording)
  windows = list(
    prcp     = list(start_s = NULL, duration_s = NULL),
    exercise = list(start_s = NULL, duration_s = NULL),
    # Fantasia records last ~2 h: a bounded window keeps the fits tractable. Windows A and B
    # (non-overlapping) are used for the test-retest stability analysis.
    fantasia = list(start_s = 600, duration_s = 600, stability_A = c(600, 1200), stability_B = c(1200, 1800))
  ),

  # Fantasia group labels come from the file names (young / older only)
  fantasia = list(young_regex = "y[0-9]+", old_regex = "o[0-9]+"),


  # ---- 3. Event definitions ----------------------------------------------------------------------------
  # Tilt angle (PRCP): events from the angle signal when present (the heart-rate-based detection is
  # then validated against it); plateaus for the dose-response over angle.
  # A recording may hold several tilts: every one is an event. Each event is analysed against its own baseline,
  # which starts `washout_s` after the previous return (the states are still recovering before that) and is at
  # most `baseline_max_s` long; an event with less than `baseline_min_s` of baseline is not analysed.
  # The recovery window of an event lasts `washout_s` (or until the next onset). Detected events are matched
  # to the angle-derived ones within `match_tol_s` (validation of the RR-based detection).
  tilt = list(baseline_s = 60, threshold_deg = 10, hold_s = 5, plateau_min_s = 20, level_bin_deg = 5,
              washout_s = 120, baseline_max_s = 180, baseline_min_s = 60, edge_skip_s = 10, match_tol_s = 60,
              # Protocol annotations (data/prcp_events.csv): label every tilt as slow or rapid, add the stand-ups (the table angle does not
              # change) as events, and use every maneuver for the washout. FALSE: tilt angle / heart rate only, as before.
              use_annotations = TRUE, annotation_slack_s = 5),

  # Heart-rate segmentation (tilt without angle, and exercise): no event timings are given.
  # Exercise: one bout, at most `max_changepoints` change points. Tilt: up to `max_tilt_events` perturbations
  # (two change points each), limited by the recording length.
  segmentation = list(bin_s = 5, min_seg_s = 30, max_changepoints = 2L, max_tilt_events = 6L, min_shift_hz = 0.10,
                      smooth_bins = 3L, min_exercise_s = 40, max_exercise_s = 360),


  # ---- 4. Analyses --------------------------------------------------------------------------------------
  # Respiration / blood-pressure analyses (Fantasia, PRCP); subjects without the signal are skipped
  physio = list(resp_band = c(0.10, 0.50), coh_halfwidth_hz = 0.04, n_surrogates = 200L, max_lag_s = 30,
                brs_bands = list(LF = c(0.04, 0.15), HF = c(0.15, 0.40)), brs_min_coh = 0.5),

  # Exercise. recovery_definition / end_method "auto" read the empirical choice made by
  # analysis/validation/07_recovery_benchmark.R (end_method: peak | cp | profile)
  exercise = list(recovery_definition = "auto", end_method = "auto",
                  smooth_s = 10, hold_s = 20, tol_sd = 1.0, rest_skip_s = 20, end_window_s = 15, max_exercise_s = 360),

  hrv   = list(lf = c(0.04, 0.15), hf = c(0.15, 0.40), fs = 4, min_dur_lfhf_s = 120),
  stats = list(n_boot = 10000L, alpha = 0.05),


  # ---- 5. Workload (reduce for a quick pass) --------------------------------------------------------------
  identifiability = list(n_starts = 100, max_recordings = 100, n_surface_recordings = 2L, surface_grid = 100L,
                         min_dispersion_log = 0.10, ll_tol = 2.0, collapse_share = c(0.05, 0.95)),
  sensitivity        = list(max_recordings = 100),
  # n_detect: event-detection validation; n_prep: preprocessing checks (simulated recordings per kind)
  simulation         = list(n_per_scenario = 100L, duration_s = 420, h = 0.01, n_detect = 100L, n_prep = 12L),
  recovery_benchmark = list(n = 100L),
  calibration        = list(n_sims = 100L, duration_s = 420),

  # ---- 6. Figures -----------------------------------------------------------------------------------------
  figures = list(
    format = "png",              # raster copy used by the manuscript: "png", "pdf" or "svg"
    dpi = 300,                   # for png
    vector_format = "pdf",       # vector copy written next to every figure (for submission); NULL for none
    journal_width_in = 7.48,     # width of the journal's double-column figures (190 mm)
    text_width_in = 6.54,        # width of the text block of the Typst PDF (A4, 2.2 cm margins): figures are shown at most this wide
    min_text_pt = 6,             # smallest acceptable axis text (pt) at the narrower of the two widths; smaller text is reported
    base_size = 10,              # base font size (pt) of every figure, at the design width of the figure (axis text is base_size - 1)
    font = "Helvetica",                   # preferred font family; the first installed of font, Arial, Helvetica, Liberation Sans, DejaVu Sans is used

    # Which recording to showcase in the example figures (7 and 9):
    #   "largest"  the recording with the largest standardized response
    #   "median"   the recording whose response is closest to the median
    #   "<id>"     a recording id, e.g. "prcp02" (see results/prcp.rds$response$SubjectID)
    example = list(prcp = "median", exercise = "median"),

    # Colours (any R colour). Edit to restyle every figure.
    palette = list(
      groups   = c(young = "#2A9D8F", older = "#C0392B"),
      datasets = c(prcp = "#3B6FB6", fantasia = "#D98E04", exercise = "#2E8B57"),
      states   = c(x_p = "#2A6EBB", x_s = "#C0392B"),                          # vagal-like, sympathetic-like
      phases   = c(baseline = "#4C78A8", tilt = "#E4572E", rest = "#4C78A8", exercise = "#E4572E",
                   recovery = "#2E8B57"),
      models   = c(M0 = "#7F7F7F", M1 = "#4C78A8", M2 = "#F58518", M3 = "#54A24B"),
      accent   = "#222222"
    ),

    # Density landscapes of the (x_p, x_s) plane (supplementary figure S6)
    landscape = list(
      palette = "viridis",          # viridis option for the potential: "mako", "magma", "viridis", "plasma", "rocket", "cividis"
      n_grid = 90,               # grid resolution of the density
      bandwidth_scale = 1,       # >1 smoother, <1 more detailed
      potential_max = 6,         # cap of U = -ln(p / p_max)
      show_flow = TRUE,          # arrows of the mean velocity (converge on attractors)
      show_basins = TRUE,        # mark the density maxima
      show_contours = TRUE,
      flow_bins = 15, flow_min_n = 10,
      min_rel_density = 0.3,     # a basin must reach this share of the peak density
      limits_quantile = c(0.005, 0.995),
      path_smooth_s = 10         # smoothing of the drawn trajectory (s); the density uses 5 s
    )
  )
)

cfg$is_placeholder <- function(p) grepl("<[^>]+>", p)

# Optional overrides without editing this file: SPARK_CONFIG=path/to/local_config.R (a script that
# modifies `cfg` in place). Used by the test suite.
local_cfg <- Sys.getenv("SPARK_CONFIG", unset = "")
if (nzchar(local_cfg) && file.exists(local_cfg)) source(local_cfg, local = TRUE)
