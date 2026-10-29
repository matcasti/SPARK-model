# ==============================================================================
# R/report/figures.R -- one function per manuscript figure. Each takes `results` (a function returning a saved
# result by name, or NULL when missing) and `cfg`, and returns list(plot, width, height, info).
# analysis/13_figures.R saves them; a figure whose inputs are missing is skipped with a message.
# The look (fonts, palettes, which recording to showcase) is set in config.R, section 6 (cfg$figures);
# the building blocks are in R/report/plots.R.
#
# Main:          fig1 fit diagnostics and time constants | fig4 CWT comparison
#                fig5 non-identifiability | fig7 orthostatic perturbation
#                fig8 Fantasia groups | fig9 exercise and recovery | fig10 Fantasia respiration / baroreflex
# Supplementary: figS1 recovery under known truth | figS2 sensitivity | figS3 recovery benchmark
#                figS4 numerical checks | figS5 event detection | figS6 state-space landscapes | figS7 PRCP responses by maneuver
# Every figure answers one question; panels of a row share their vertical extent and panels that share a
# time axis are facets of one plot (R/report/plots.R).
# ==============================================================================


# ---- Helpers -----------------------------------------------------------------------------------------------
#' Result by name, or a clear error naming what is missing.
need <- function(results, name, hint = NULL) {
  value <- results(name)
  if (is.null(value)) {
    stop(sprintf("needs results/%s.rds%s", name, if (is.null(hint)) "" else paste0(" (", hint, ")")), call. = FALSE)
  }
  value
}

n_label <- function(n, what = "recordings") sprintf("n = %d %s", n, what)

#' Which recording to showcase. `scores`: named numeric vector (response magnitude per recording).
#' `setting`: "largest", "median" or a recording id (cfg$figures$example). An unknown id falls back
#' to "largest" with a warning.
choose_example <- function(scores, setting = "largest") {
  scores <- scores[is.finite(scores)]
  if (setting %in% names(scores)) return(list(id = setting, rule = "chosen by id"))
  if (!setting %in% c("largest", "median")) {
    warning(sprintf("example '%s' not found among: %s; using 'largest'.", setting, paste(names(scores), collapse = ", ")),
            call. = FALSE)
    setting <- "largest"
  }
  id <- if (setting == "largest") names(which.max(scores)) else names(which.min(abs(scores - median(scores))))
  list(id = id, rule = if (setting == "largest") "largest standardized response" else "median standardized response")
}

#' The showcased recording of a perturbation dataset: its series, id and how it was chosen.
#' @param kind "tilt" (results$prcp) or "exercise" (results$exercise)
example_recording <- function(results, kind, cfg) {
  if (kind == "tilt") {
    prcp <- need(results, "prcp", "run analysis/datasets/05_prcp.R")
    if (is.null(prcp$series_all)) stop("results/prcp.rds has no series_all: rerun analysis/datasets/05_prcp.R", call. = FALSE)
    response <- prcp$response
    series <- prcp$series_all
    setting <- cfg$figures$example$prcp
  } else {
    exercise <- need(results, "exercise", "run analysis/datasets/09_exercise.R")
    response <- exercise$rest_response
    series <- exercise$series_all
    setting <- cfg$figures$example$exercise
  }
  scores <- setNames(abs(response$z_d_x_p) + abs(response$z_d_x_s), response$SubjectID)
  scores <- scores[names(scores) %in% names(series)]
  choice <- choose_example(scores, setting)
  list(series = series[[choice$id]], id = choice$id, rule = choice$rule, all_series = series)
}

#' Plain axis labels for log scales (0.01, 1, 100, 10000 instead of 1e-02, ...).
# Labels of log axes: as many decimals as the value needs (1e-4 used to be printed as "0")
log_labels <- function(x) vapply(x, function(v) if (!is.finite(v)) "" else format(v, scientific = FALSE, trim = TRUE, drop0trailing = TRUE, digits = 6), character(1))

#' Boxplot with the recordings overlaid; the median is the box line.
box_points <- function(mapping, width = 0.12, size = 1.6) {
  list(geom_boxplot(mapping, outlier.shape = NA, linewidth = 0.35, alpha = 0.2),
       geom_jitter(mapping, width = width, height = 0, size = size, alpha = 0.85))
}


# ---- Figure 1: fit diagnostics and estimated time constants ----------------------------------------------------
figure_fit_diagnostics <- function(results, cfg) {
  population <- as.data.table(need(results, "df_population_summary"))
  colours <- figure_palette(cfg, "datasets")
  priors <- spark_priors()
  pacf_by_dataset <- rbindlist(lapply(cfg$datasets, function(dataset) {
    fits <- results(paste0("fit_", dataset))
    if (is.null(fits)) return(NULL)
    data.table(Dataset = dataset, lag = seq_len(ncol(fits$pacf_matrix)), mean = colMeans(fits$pacf_matrix, na.rm = TRUE))
  }))
  pacf_by_dataset[, bound := 1.96 / sqrt(population[Dataset == Dataset[1], median(NumValidBeats)]), by = Dataset]
  population[, crit := 1.36 / sqrt(NumValidBeats)]
  medians <- population[, .(RMSE_ms = median(RMSE_ms)), by = Dataset]

  a <- ggplot(population, aes(Dataset, RMSE_ms, fill = Dataset, colour = Dataset)) + box_points(aes(fill = Dataset, colour = Dataset)) +
    scale_x_discrete(labels = setNames(sprintf("%s<br>median %.1f", medians$Dataset, medians$RMSE_ms), medians$Dataset)) +
    scale_fill_manual(values = colours, guide = "none") + scale_colour_manual(values = colours, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(x = NULL, y = "RMSE of filtered RR (ms)", title = "Tracking error") + theme_spark(cfg)
  b <- ggplot(pacf_by_dataset, aes(lag, mean, colour = Dataset)) +
    geom_ribbon(aes(ymin = -bound, ymax = bound, fill = Dataset), alpha = 0.08, colour = NA) +
    geom_hline(yintercept = 0, linewidth = 0.3) + geom_line(linewidth = 0.6) +
    scale_colour_manual(values = colours, name = NULL) + scale_fill_manual(values = colours, guide = "none") +
    labs(x = "Lag (beats)", y = "Mean PACF of innovations", title = "Serial dependence", subtitle = "bands: white-noise bounds") + theme_spark(cfg)
  c <- ggplot(population, aes(crit, KS_Stat, colour = Dataset)) +
    geom_abline(slope = 1, linetype = 2, linewidth = 0.4, colour = "grey40") + geom_point(size = 1.8, alpha = 0.85) +
    scale_colour_manual(values = colours, guide = "none") +
    labs(x = tex("Critical value $1.36/\\sqrt{N}$"), y = tex("KS distance $D$"), title = "Time-rescaling distance", subtitle = tex("below line: $D$ < critical value")) + theme_spark(cfg)

  taus <- melt(population[, .(Dataset, SubjectID, `x_p` = Tau_P, `x_s` = Tau_S)], id.vars = c("Dataset", "SubjectID"))
  prior_bands <- data.table(variable = factor(c("x_p", "x_s")), lo = c(priors$tau_p_bounds[1], priors$tau_s_bounds[1]),
                            hi = c(priors$tau_p_bounds[2], priors$tau_s_bounds[2]), x = 1:2)
  d <- ggplot(taus, aes(variable, value)) + scale_x_discrete(labels = tex_named(c(x_p = "$x_p$", x_s = "$x_s$"))) +
    geom_rect(data = prior_bands, aes(xmin = x - 0.42, xmax = x + 0.42, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey80", alpha = 0.7) +
    geom_jitter(aes(colour = Dataset), width = 0.12, height = 0, size = 1.8, alpha = 0.85) +
    stat_summary(fun = median, geom = "crossbar", width = 0.6, linewidth = 0.35, colour = "black") +
    scale_y_log10(labels = log_labels) + scale_colour_manual(values = colours, guide = "none") +
    labs(x = NULL, y = tex("$\\tau = 1/\\kappa$ (s)"), title = "Estimated time constants", subtitle = "grey: prior; bar: median") + theme_spark(cfg)
  plot <- (figure_grid(list(a, b, c, d), ncol = 2)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  list(plot = plot, width = 9, height = 7.6, info = list())
}


# ---- Figure 5: non-identifiability under RR alone ---------------------------------------------------------------------------
#' Near-equivalent fits (within the likelihood tolerance) of x_p for one recording, M0 vs M3.
near_equivalent_trajectories <- function(trajectories, ll_tolerance) {
  recording <- trajectories[[1]]
  one_spec <- function(spec) {
    fits <- recording[[spec]]
    ok <- which(!vapply(fits$mats, is.null, logical(1)))
    rows <- rbindlist(lapply(ok, function(i) {
      data.table(t = recording$M0$t, x = fits$mats[[i]][, "x_p"], start = i, LL = fits$table$LL_star[i], spec = spec)
    }))
    rows[, dLL := max(LL) - LL][dLL <= ll_tolerance]
  }
  rbind(one_spec("M0"), one_spec("M3"))
}

figure_identifiability <- function(results, cfg) {
  ident <- need(results, "identifiability")
  models <- figure_palette(cfg, "models")
  datasets <- figure_palette(cfg, "datasets")
  priors <- spark_priors()
  tolerance <- cfg$identifiability$ll_tol
  short_labels <- tex_named(c(M0 = "$M_0$: RR only", M3 = "$M_3$: RR + physiological\ninformation"))

  # ---- A: where the best fit puts the variability (one value per recording) -----------------------------------------------------
  collapse_share <- cfg$identifiability$collapse_share
  if (is.null(collapse_share)) collapse_share <- c(0.01, 0.99)
  starts <- ident$starts[ok == TRUE & spec %in% c("M0", "M3") & is.finite(share_p)]
  starts[, dLL := max(LL_star) - LL_star, by = .(SubjectID, spec)]                    # loss relative to the best start of the same specification
  starts[, equivalent := dLL <= tolerance]
  per_recording <- starts[, .(best = share_p[which.max(LL_star)]), by = .(Dataset, SubjectID, spec)]
  per_recording[, collapsed := best < collapse_share[1] | best > collapse_share[2]]
  dataset_levels <- names(datasets)[names(datasets) %in% unique(per_recording$Dataset)]
  per_recording[, `:=`(Dataset = factor(Dataset, dataset_levels), spec_f = factor(spec, names(short_labels), short_labels))]
  percent_of <- function(sp, condition) sprintf("%.0f%%", 100 * mean(per_recording[spec == sp, eval(condition)]))
  dataset_fill <- scale_fill_manual(values = unlist(datasets[dataset_levels]), name = "Dataset", drop = FALSE)
  histogram <- function(column, xlab, title, subtitle) {
    ggplot(per_recording, aes(.data[[column]], fill = Dataset)) +
      geom_histogram(binwidth = 0.05, center = 0, closed = "right", linewidth = 0.2, colour = "white") +
      facet_wrap(~spec_f, ncol = 1, scales = "free_y") +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + coord_cartesian(xlim = c(-0.04, 1.04)) +
      dataset_fill + scale_y_continuous(breaks = function(x) { b <- pretty(x, n = 4); b[b %% 1 == 0] }, expand = expansion(mult = c(0, 0.05))) + guides(fill = guide_legend()) +
      labs(x = tex(xlab), y = "Recordings", title = title, subtitle = tex(subtitle)) + theme_spark(cfg)
  }
  a1 <- histogram("best", "Share of $x_p$ variance at the best fit (0: all in $x_s$; 1: all in $x_p$)", "Variance split at the best fit",
                  sprintf("one process takes all (< %s or > %s):\n$M_0$ %s of recordings, $M_3$ %s", collapse_share[1], collapse_share[2], percent_of("M0", quote(collapsed)), percent_of("M3", quote(collapsed))))

  # ---- B: profile of the objective with the end points of the starts (one recording) -------------------------------------------
  surface <- ident$surfaces[spec %in% c("M0", "M3")]
  b <- NULL
  if (nrow(surface)) {
    sid <- surface$SubjectID[1]
    surface <- surface[SubjectID == sid]
    surface[, spec_f := factor(spec, names(short_labels), short_labels)]
    x_range <- range(surface$tau_s); y_range <- range(surface$tau_p)
    one <- ident$starts[ok == TRUE & SubjectID == sid & spec %in% c("M0", "M3") & is.finite(tau_s) & is.finite(tau_p)]
    one[, `:=`(outside = tau_s < x_range[1] | tau_s > x_range[2] | tau_p < y_range[1] | tau_p > y_range[2],
               tau_s_plot = pmin(pmax(tau_s, x_range[1]), x_range[2]), tau_p_plot = pmin(pmax(tau_p, y_range[1]), y_range[2]),
               spec_f = factor(spec, names(short_labels), short_labels))]
    prior_box <- data.table(xmin = priors$tau_s_bounds[1], xmax = priors$tau_s_bounds[2], ymin = priors$tau_p_bounds[1], ymax = priors$tau_p_bounds[2])
    b <- ggplot(surface, aes(tau_s, tau_p)) +
      geom_contour_filled(aes(z = pmin(dobj, 20)), breaks = c(0, 2, 5, 10, 20), na.rm = TRUE) +
      geom_rect(data = prior_box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = NA, colour = "grey25", linetype = 2, linewidth = 0.35) +
      geom_point(data = one, aes(tau_s_plot, tau_p_plot, shape = outside), colour = "#D95F02", size = 1.7, alpha = 0.9) +
      facet_grid(. ~ spec_f) + coord_cartesian(xlim = x_range, ylim = y_range) +
      scale_x_log10(breaks = c(30, 100, 300), labels = log_labels) +
      scale_y_log10(breaks = c(1, 3, 10, 30), labels = log_labels) +
      scale_fill_brewer(palette = "Blues", direction = -1, name = tex("$\\Delta$ objective")) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "end point of a start", `TRUE` = "beyond the plotted range (at the edge)"), name = NULL) +
      labs(x = tex("$\\tau_s$ (s)"), y = tex("$\\tau_p$ (s)"), title = "Profile of the objective",
           subtitle = "one recording; dashed box: prior") +
      theme_spark(cfg) + theme(legend.box = "vertical", legend.spacing.y = unit(0.05, "cm"))
  }

  # ---- C: reconstructed x_p of the near-equivalent starts (same recording) ----------------------------------------------------
  c <- NULL
  if (length(ident$trajectories)) {
    near_equivalent <- near_equivalent_trajectories(ident$trajectories, tolerance)
    near_equivalent[, spec := factor(spec, names(short_labels), short_labels)]
    c <- ggplot(near_equivalent, aes(t, x, group = start, colour = spec)) +
      geom_line(linewidth = 0.3, alpha = 0.6) + facet_grid(spec ~ ., labeller = labeller(spec = setNames(tex(c("$M_0$", "$M_3$")), short_labels))) +
      scale_colour_manual(values = setNames(unname(models[c("M0", "M3")]), short_labels), guide = "none") +
      labs(x = "Time (s)", y = tex("$x_p$ (Hz)"), title = tex("Reconstructed $x_p$, same recording"),
           subtitle = "all starts within the likelihood tolerance") + theme_spark(cfg)
  }
  upper <- if (!is.null(b)) (tag_panel(a1, "A") | tag_panel(b, "B")) + plot_layout(widths = c(1, 1.1)) else tag_panel(a1, "A")
  plot <- if (is.null(c)) upper else (upper / tag_panel(c, "C")) + plot_layout(heights = c(1.15, 0.8))
  list(plot = plot, width = 9, height = 8, info = list())
}

# ---- Figure 7: orthostatic perturbation --------------------------------------------------------------------------------------
#' Panels of the group-level responses (one point per recording; small points: individual tilts).
prcp_amplitude_panel <- function(prcp, cfg) {
  colour <- figure_palette(cfg, "datasets")[["prcp"]]
  per_recording <- prcp$response
  plot <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#2E8B57", alpha = 0.07) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(data = prcp$response_events, aes(z_d_x_p, z_d_x_s), size = 1.2, alpha = 0.35, colour = colour) +
    geom_point(data = per_recording, aes(z_d_x_p, z_d_x_s), size = 2.7, colour = colour)
  if (nrow(per_recording) <= 12) {
    plot <- plot + if (requireNamespace("ggrepel", quietly = TRUE))   # labels pushed apart; without ggrepel, overlapping labels are dropped
      ggrepel::geom_text_repel(data = per_recording, aes(z_d_x_p, z_d_x_s, label = SubjectID), size = 2.7, colour = "grey30", segment.size = 0.2, min.segment.length = 0.2, max.overlaps = Inf, seed = 1)
    else geom_text(data = per_recording, aes(z_d_x_p, z_d_x_s, label = SubjectID), size = 2.7, hjust = -0.25, colour = "grey30", check_overlap = TRUE)
  }
  plot + labs(x = tex("$\\Delta x_p$ (baseline SD)"), y = tex("$\\Delta x_s$ (baseline SD)"), title = "Response amplitude",
              subtitle = "large: recording median") + theme_spark(cfg)
}

#' Onset latency and step-response time constant of x_p, x_s and MAP across recordings (open: at a search bound).
prcp_timing_panels <- function(prcp, cfg) {
  states <- figure_palette(cfg, "states")
  colour_of <- c(x_p = states[["x_p"]], x_s = states[["x_s"]], MAP = "#6B6B6B")
  response <- prcp$response
  gather <- function(p, s, m) rbind(data.table(v = "x_p", y = p), data.table(v = "x_s", y = s),
                                    if (!is.null(m)) data.table(v = "MAP", y = m))[is.finite(y)]
  latency <- gather(response$lat_p, response$lat_s, prcp$bp$lat_map)[, bound := abs(y) >= 39]
  tau <- gather(response$tau_p_resp, response$tau_s_resp, prcp$bp$tau_map)[, bound := y <= 2.01 | y >= 249]
  panel <- function(data, ylab, title, log = FALSE) {
    p <- ggplot(data, aes(v, y, colour = v)) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.35, fill = NA, width = 0.5) +
      geom_jitter(aes(shape = bound), width = 0.1, height = 0, size = 1.9, alpha = 0.9) +
      scale_colour_manual(values = colour_of, guide = "none") + scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
      scale_x_discrete(labels = tex_named(c(x_p = "$x_p$", x_s = "$x_s$", MAP = "MAP"))) +
      labs(x = NULL, y = ylab, title = title, subtitle = tex("open: at a search bound\n")) + theme_spark(cfg)
    if (log) p + scale_y_log10(labels = log_labels) else p + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3)
  }
  list(latency = panel(latency, "Onset minus tilt onset (s)", "Onset latency"),
       tau = panel(tau, "Time constant (s)", "Response time constant", log = TRUE))
}

#' Peak lagged correlation of each state with MAP (positive lag: MAP lags the state); open: surrogate p >= 0.05.
prcp_lag_panel <- function(prcp, cfg) {
  states <- figure_palette(cfg, "states")
  lagged <- copy(prcp$xcorr)
  lagged[, state := factor(state, c("X_p", "X_s"), c("x_p", "x_s"))]
  ggplot(lagged, aes(lag_s, peak_r, colour = state, fill = ifelse(p < 0.05, as.character(state), NA))) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_point(shape = 21, size = 2.4, stroke = 0.8) +
    scale_colour_manual(values = c(x_p = states[["x_p"]], x_s = states[["x_s"]]), guide = "none") +
    scale_fill_manual(values = c(x_p = states[["x_p"]], x_s = states[["x_s"]]), na.value = "white", guide = "none") +
    labs(x = "Lag of MAP behind the state (s)", y = "Peak detrended correlation", title = "State-MAP coupling",
         subtitle = tex("blue: $x_p$; red: $x_s$")) + theme_spark(cfg)
}

#' Mean state, HR and MAP by tilt-angle plateau, standardized (recordings with several levels).
prcp_dose_panel <- function(prcp, cfg) {
  states <- figure_palette(cfg, "states")
  dose <- melt(prcp$dose_levels, id.vars = "level")
  dose[, z := (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE), by = variable]
  ggplot(dose, aes(level, z, colour = variable)) + geom_line(linewidth = 0.8) + geom_point(size = 2) +
    scale_colour_manual(values = c(x_p = states[["x_p"]], x_s = states[["x_s"]], hr = "#444444", map = "#6B6B6B"), name = NULL) +
    labs(x = "Tilt angle (deg)", y = "Standardized mean", title = "Dose-response", subtitle = "mean over each plateau") + theme_spark(cfg) + theme(legend.position = "none")
}

figure_prcp_group <- function(results, cfg) {
  prcp <- need(results, "prcp")
  example <- example_recording(results, "tilt", cfg)
  locked <- event_locked(prcp$series_all, "tilt", c(-120, 150))
  variables <- intersect(c("x_p", "x_s", "hr", "map"), unique(locked$variable))
  n_tilts <- prcp$summary$n_events

  windows <- example$series$ev$windows
  has_stand <- "kind" %in% names(windows) && any(windows$kind == "stand-up")
  kinds_text <- if (is.null(prcp$by_kind)) "" else {
    k <- prcp$by_kind[kind %in% c("slow tilt", "rapid tilt")]
    if (nrow(k) > 1) paste0(" (", paste(k$n_events, k$kind, collapse = ", "), ")") else ""
  }
  a <- series_stack(example$series, "tilt", cfg, title = "One recording",
                    subtitle = paste0("recording ", example$id))
  b <- plot_event_locked(locked, cfg, variables, "tilt onset", title = "Average response to a tilt",
                         subtitle = sprintf("%s, %d tilts", n_label(uniqueN(locked$id)), n_tilts))

  wide <- dcast(locked[variable %in% c("x_p", "x_s")], id + u ~ variable, value.var = "value")[is.finite(x_p) & is.finite(x_s)]
  wide[, phase := factor(ifelse(u < 0, "baseline", "tilt"), c("baseline", "tilt"))]
  paths <- wide[, .(id, u = x_p, v = x_s, phase)]
  central <- wide[, .(u = median(x_p), v = median(x_s), phase = phase[1]), by = .(time = u)][order(time)][, !"time"]
  c <- plot_state_plane(paths, cfg, central, title = "Path in the state plane",
                        subtitle = "blue: before; orange: during", legend = FALSE)
  timing <- prcp_timing_panels(prcp, cfg)
  bottom <- list(c, prcp_amplitude_panel(prcp, cfg), timing$latency, timing$tau)
  if (!is.null(prcp$xcorr) && nrow(prcp$xcorr)) bottom[[length(bottom) + 1]] <- prcp_lag_panel(prcp, cfg)
  if (!is.null(prcp$dose_levels) && nrow(prcp$dose_levels) >= 3) bottom[[length(bottom) + 1]] <- prcp_dose_panel(prcp, cfg)
  letters <- LETTERS[seq_len(length(bottom) + 2)]
  top <- (tag_panel(a, letters[1]) | tag_panel(b, letters[2])) + plot_layout(widths = c(1.7, 1))
  lower <- wrap_plots(mapply(tag_panel, bottom, letters[-(1:2)], SIMPLIFY = FALSE), ncol = 3)             # three per row: legible at the journal width
  list(plot = top / lower + plot_layout(heights = c(2.2, 1.25 * ceiling(length(bottom) / 3))), width = 9, height = 6 + 2.8 * ceiling(length(bottom) / 3),
       info = list(prcp = example$id, prcp_rule = example$rule))
}


# ---- Figure 4: filtered states vs CWT moving band power --------------------------------------------------------------------
figure_cwt_comparison <- function(results, cfg) {
  prcp <- need(results, "prcp")
  locked <- event_locked(prcp$series_all, "tilt", c(-120, 150), with_cwt = TRUE)
  if (!all(c("LF", "HF") %in% locked$variable)) stop("series carry no CWT band power: rerun analysis/datasets/05_prcp.R", call. = FALSE)
  states <- figure_palette(cfg, "states")
  summary <- locked[, .(median = median(value), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75), n = .N), by = .(variable, u)]
  summary <- summary[n >= max(2, 0.5 * uniqueN(locked$id))]
  kind_colours <- c(`filtered state (causal)` = "#222222", `CWT band power (windowed)` = "#E0A100")
  pair <- function(state, power, colour, title) {
    d <- summary[variable %in% c(state, power)]
    d[, kind := factor(ifelse(variable == state, names(kind_colours)[1], names(kind_colours)[2]), names(kind_colours))]
    ggplot(d, aes(u, median, colour = kind, fill = kind)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) + geom_vline(xintercept = 0, linetype = 2, colour = "grey35", linewidth = 0.3) +
      geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
      scale_colour_manual(values = kind_colours, name = NULL) + scale_fill_manual(values = kind_colours, name = NULL) +
      labs(x = "Time from tilt onset (s)", y = "Change from baseline (SD)", title = title,
           subtitle = "median and quartiles") + theme_spark(cfg)
  }
  a <- pair("x_s", "LF", states[["x_s"]], tex("$x_s$ vs LF power"))
  b <- pair("x_p", "HF", states[["x_p"]], tex("$x_p$ vs HF power"))
  departures <- melt(prcp$lead[, .(SubjectID, `x_s` = lead_x_s, `LF power` = lead_LF, `x_p` = lead_x_p, `HF power` = lead_HF)], id.vars = "SubjectID")
  departures[, kind := factor(ifelse(grepl("power", variable), names(kind_colours)[2], names(kind_colours)[1]), names(kind_colours))]
  c <- ggplot(departures[is.finite(value)], aes(variable, value, colour = kind)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.35, fill = NA, width = 0.5) + geom_jitter(width = 0.1, height = 0, size = 1.7, alpha = 0.85) +
    scale_colour_manual(values = kind_colours, guide = "none") +
    scale_x_discrete(labels = tex_named(c(x_s = "$x_s$", `LF power` = "LF power", x_p = "$x_p$", `HF power` = "HF power"))) +
    labs(x = NULL, y = "First departure from baseline (s)", title = "Timing of the departure",
         subtitle = n_label(uniqueN(departures$SubjectID))) + theme_spark(cfg)
  plot <- (figure_grid(list(a, b, c), ncol = 3, widths = c(1, 1, 0.8))) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  list(plot = plot, width = 9, height = 5, info = list())
}


# ---- Figure 8: Fantasia young vs older ---------------------------------------------------------------------------------------
figure_fantasia_groups <- function(results, cfg) {
  fantasia <- need(results, "fantasia")
  groups <- figure_palette(cfg, "groups")
  families <- c(latent = figure_palette(cfg, "states")[["x_p"]], HRV = "#E45756")
  latent_labels <- c(Tau_P = "$\\tau_p$ (s)", Tau_S = "$\\tau_s$ (s)", lambda_P = "$\\lambda_P$", sd_x_p = "SD of $x_p$ (Hz)",
                     sd_x_s = "SD of $x_s$ (Hz)", Sigma_sys2 = "$\\hat\\sigma^2_{sys}$ (Hz$^2$)")
  latent <- melt(fantasia$table[, c("Group", "SubjectID", names(latent_labels)), with = FALSE], id.vars = c("Group", "SubjectID"))
  latent[, variable := factor(variable, names(latent_labels), tex(unname(latent_labels)))]
  contrasts <- fantasia$contrasts[is.finite(cliff_delta)]
  n_group <- fantasia$table[, .N, by = Group]

  a <- ggplot(latent, aes(Group, value, colour = Group, fill = Group)) + box_points(aes(colour = Group, fill = Group), size = 1.4) +
    facet_wrap(~variable, scales = "free_y", nrow = 1) + scale_colour_manual(values = groups, guide = "none") +
    scale_fill_manual(values = groups, guide = "none") +
    labs(x = NULL, y = NULL, title = "Latent summaries by group",
         subtitle = paste(n_group[order(factor(Group, c("young", "older"))), paste0(Group, " n = ", N)], collapse = ";  ")) + theme_spark(cfg)
  b <- ggplot(contrasts, aes(cliff_delta, reorder(variable, cliff_delta), colour = family)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) + geom_pointrange(aes(xmin = cliff_lo, xmax = cliff_hi), size = 0.4) +
    scale_colour_manual(values = families, name = NULL) + coord_cartesian(xlim = c(-1, 1)) + scale_y_discrete(labels = sym) +
    labs(x = tex("Cliff's $\\delta$, older vs young (95% bootstrap CI)"), y = NULL, title = "Group contrasts",
         subtitle = "latent and HRV summaries") + theme_spark(cfg)
  c <- if (!is.null(fantasia$stability)) {
    ggplot(fantasia$stability, aes(icc_log, reorder(variable, icc_log), fill = family)) + geom_col(width = 0.7) +
      scale_fill_manual(values = families, guide = "none") + coord_cartesian(xlim = c(min(0, fantasia$stability$icc_log, na.rm = TRUE), 1)) +
      scale_y_discrete(labels = sym) +
      labs(x = "Test-retest ICC(3,1) of the log values", y = NULL, title = "Stability across two windows",
           subtitle = "1 = identical ranking") + theme_spark(cfg)
  }
  lower <- if (is.null(c)) tag_panel(b, "B") else (tag_panel(b, "B") | tag_panel(c, "C")) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  list(plot = tag_panel(a, "A") / lower + plot_layout(heights = c(1, 1.25)), width = 9, height = 7.4, info = list())
}


# ---- Figure 9: exercise and recovery ----------------------------------------------------------------------------------------
figure_exercise <- function(results, cfg) {
  exercise <- need(results, "exercise")
  example <- example_recording(results, "exercise", cfg)
  phase_colours <- figure_palette(cfg, "phases")
  locked <- event_locked(exercise$series_all, "exercise", c(-150, 300))
  paths <- rbindlist(lapply(names(exercise$series_all), function(id) {
    cbind(id = id, normalized_path(exercise$series_all[[id]], "exercise", cfg$figures$landscape$path_smooth_s))
  }))[!is.na(phase)]
  recovery <- exercise$recovery[definition == exercise$summary$recovery_definition & status != "no_excursion"]
  recovery[, y := ifelse(censored, cens_time, time)]
  recovery[, state_label := factor(state, c("p", "s"), c("x_p", "x_s"))]
  km_medians <- exercise$km[, .(state_label = factor(state, c("p", "s"), c("x_p", "x_s")), y = median)][is.finite(y)]
  states <- figure_palette(cfg, "states")

  a <- series_stack(example$series, "exercise", cfg, title = "One recording",
                    subtitle = paste0("recording ", example$id))
  b <- plot_event_locked(locked, cfg, intersect(c("x_p", "x_s", "hr"), unique(locked$variable)), "end of exercise", reference = "rest", title = "Recovery after exercise",
                         subtitle = n_label(uniqueN(locked$id)), highlight = example$id)
  c <- ggplot(recovery, aes(state_label, y, colour = state_label)) +
    geom_jitter(aes(shape = censored), width = 0.12, height = 0, size = 2.2) +
    geom_errorbar(data = km_medians, aes(ymin = y, ymax = y), width = 0.5, linewidth = 0.8, colour = "black") +
    scale_colour_manual(values = c(x_p = states[["x_p"]], x_s = states[["x_s"]]), guide = "none") +
    scale_x_discrete(labels = tex_named(c(x_p = "$x_p$", x_s = "$x_s$"))) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("recovered", "censored at recording end"), name = NULL) +
    labs(x = NULL, y = tex(paste0("Recovery time (s):\n", words_for(as.character(exercise$summary$recovery_definition), RECOVERY_DEFINITION_TEXT))), title = "Recovery time",
         subtitle = "bar: Kaplan-Meier median") + theme_spark(cfg)
  d <- plot_state_plane(paths, cfg, NULL, reference = "rest", title = "State-plane path, by phase", highlight = example$id,
                        subtitle = paste0(n_label(uniqueN(paths$id)), "; common axes"))
  top <- (tag_panel(a, "A") | tag_panel(b, "B")) + plot_layout(widths = c(1.4, 1))
  lower <- (tag_panel(c, "C") | tag_panel(d, "D")) + plot_layout(widths = c(0.75, 1.9))
  list(plot = top / lower + plot_layout(heights = c(1.5, 1)), width = 9, height = 9.4,
       info = list(exercise = example$id, exercise_rule = example$rule))
}


# ---- Figure 10: Fantasia respiration and baroreflex -------------------------------------------------------------------------
figure_fantasia_physiology <- function(results, cfg) {
  fantasia <- need(results, "fantasia")
  physiology <- fantasia$physiology
  groups <- figure_palette(cfg, "groups")
  if (is.null(physiology)) stop("needs fantasia physiology", call. = FALSE)
  has_signal <- any(c("resp_sd_x_p", "alpha_LF") %in% names(physiology)) && any(is.finite(c(physiology$resp_sd_x_p, physiology$alpha_LF)))
  if (!has_signal) stop("no respiration or blood-pressure data in any Fantasia recording", call. = FALSE)
  if (!"alpha_LF" %in% names(physiology)) physiology[, alpha_LF := NA_real_]
  if (!"resp_sd_x_p" %in% names(physiology)) physiology[, `:=`(resp_sd_x_p = NA_real_, resp_sd_x_s = NA_real_, f_resp = NA_real_)]

  panels <- list()
  with_resp <- physiology[is.finite(resp_sd_x_p) & is.finite(resp_sd_x_s)]
  if (nrow(with_resp)) {
    band_sd <- melt(with_resp[, .(SubjectID, Group, `x_p` = resp_sd_x_p, `x_s` = resp_sd_x_s)], id.vars = c("SubjectID", "Group"))
    panels$a <- ggplot(band_sd, aes(variable, value)) + geom_line(aes(group = SubjectID), colour = "grey75") +
      scale_x_discrete(labels = tex_named(c(x_p = "$x_p$", x_s = "$x_s$"))) +
      geom_point(aes(colour = Group), size = 2.2) + scale_colour_manual(values = groups, name = NULL) + scale_y_log10(labels = log_labels) +
      labs(x = NULL, y = "SD within the respiratory band (Hz)", title = "Respiratory-band SD",
           subtitle = n_label(nrow(with_resp))) + theme_spark(cfg)
    panels$b <- ggplot(physiology[is.finite(f_resp)], aes(Group, f_resp, colour = Group, fill = Group)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = cfg$hrv$hf[1], ymax = cfg$hrv$hf[2], fill = "grey70", alpha = 0.25) +
      box_points(aes(colour = Group, fill = Group), size = 1.8) +
      scale_colour_manual(values = groups, guide = "none") + scale_fill_manual(values = groups, guide = "none") +
      coord_cartesian(ylim = cfg$physio$resp_band) +
      labs(x = NULL, y = "Respiratory frequency (Hz)", title = "Respiration vs HF band", subtitle = "grey: 0.15-0.40 Hz") + theme_spark(cfg)
  }
  with_bp <- physiology[is.finite(alpha_LF)]
  if (nrow(with_bp)) {
    panels$c <- ggplot(with_bp, aes(Group, alpha_LF, colour = Group, fill = Group)) + box_points(aes(colour = Group, fill = Group), size = 1.8) +
      scale_colour_manual(values = groups, guide = "none") + scale_fill_manual(values = groups, guide = "none") +
      labs(x = NULL, y = tex("Baroreflex $\\alpha$, LF (ms/mmHg)"), title = "Baroreflex sensitivity", subtitle = n_label(nrow(with_bp))) + theme_spark(cfg)
    panels$d <- ggplot(with_bp, aes(sd_x_p, alpha_LF, colour = Group)) + geom_point(size = 2.2) + scale_colour_manual(values = groups, guide = "none") +
      labs(x = tex("SD of $x_p$ (Hz)"), y = tex("Baroreflex $\\alpha$, LF (ms/mmHg)"), title = "Latent state vs baroreflex", subtitle = "one point per recording") + theme_spark(cfg)
  }
  if (!length(panels)) stop("no respiration or blood-pressure data", call. = FALSE)
  list(plot = figure_grid(unname(panels), ncol = 2) + plot_layout(guides = "collect") & theme(legend.position = "bottom"), width = 9, height = 7.6, info = list())
}


# ---- Supplementary: event detection ------------------------------------------------------------------------------------------
figure_event_detection <- function(results, cfg) {
  prcp <- need(results, "prcp")
  detection <- need(results, "event_detection", "run analysis/validation/08_event_detection_validation.R")
  colour <- figure_palette(cfg, "datasets")[["prcp"]]
  validation <- prcp$validation[matched == TRUE]
  real <- melt(validation[, .(SubjectID, event, onset = err_onset, `return` = err_return)], id.vars = c("SubjectID", "event"))[is.finite(value)]
  a <- if (nrow(real)) {
    ggplot(real, aes(variable, value)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
      geom_boxplot(outlier.shape = NA, linewidth = 0.35, fill = NA, width = 0.5, colour = colour) +
      geom_jitter(width = 0.1, height = 0, size = 1.6, alpha = 0.7, colour = colour) +
      labs(x = NULL, y = "Detected minus reference (s)", title = "Real recordings",
           subtitle = sprintf("%d of %d maneuvers matched", nrow(validation), nrow(prcp$validation))) +
      theme_spark(cfg)
  }
  parts <- list(`Tilt onset\n(one tilt)` = detection$tilt$onset_err, `Tilt onset\n(repeated)` = detection$repeated$onset_err,
                `Exercise\nonset\n(profile)` = detection$exercise$onset_err_profile, `Exercise\nend\n(peak)` = detection$exercise$end_err_peak,
                `Exercise\nend\n(change\npoint)` = detection$exercise$end_err_cp, `Exercise\nend\n(profile)` = detection$exercise$end_err_profile)
  simulated <- rbindlist(lapply(names(parts), function(n) data.table(measure = n, error = parts[[n]])))[is.finite(error)]
  simulated[, measure := factor(measure, names(parts))]
  b <- ggplot(simulated, aes(measure, error)) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) + scale_x_discrete(labels = tex) +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.35, fill = "grey93") +
    labs(x = NULL, y = "Detected minus true time (s)", title = "Simulations with known timing",
         subtitle = sprintf("repeated tilts: %.0f%% found", 100 * detection$summary$repeated_recall)) +
    theme_spark(cfg)
  rates <- melt(detection$false_detection[, .(min_shift_hz, `false detections (stationary)` = false_detection_rate, `detected tilts` = detection_rate)],
                id.vars = "min_shift_hz")
  c <- ggplot(rates, aes(min_shift_hz, value, colour = variable)) + geom_vline(xintercept = cfg$segmentation$min_shift_hz, linetype = 2, linewidth = 0.3, colour = "grey40") +
    geom_line(linewidth = 0.7) + geom_point(size = 2) + scale_y_continuous(limits = c(0, 1)) +
    scale_colour_manual(values = c(`false detections (stationary)` = "#C0392B", `detected tilts` = "#2A6EBB"), name = NULL) +
    labs(x = "Minimum shift between segments (Hz)", y = "Share of simulated recordings", title = "Choice of the minimum shift",
         subtitle = "dashed: value used") + theme_spark(cfg) + guides(colour = guide_legend(nrow = 1))
  list(plot = wrap_plots(mapply(tag_panel, Filter(Negate(is.null), list(a, b, c)), LETTERS[1:3], SIMPLIFY = FALSE),
                         design = "AB\nCC", widths = c(1, 2.1), heights = c(1.6, 1)), width = 9, height = 7.4, info = list())
}


# ---- Supplementary: state-space landscapes ------------------------------------------------------------------------------------
#' One dataset's row: a landscape per phase, pooled over recordings (every recording counts equally).
landscape_row <- function(results, kind, cfg) {
  style <- cfg$figures$landscape
  source <- if (kind == "tilt") need(results, "prcp")$series_all else need(results, "exercise")$series_all
  pooled <- pooled_phase_landscapes(source, kind, style)
  n_recordings <- uniqueN(pooled$paths$id)
  panels <- lapply(names(pooled$landscapes), function(current_phase) {
    flow <- flow_field(pooled$paths[phase == current_phase, .(id, t, u, v)], pooled$limits, style$flow_bins, style$flow_min_n)
    n_basins <- nrow(pooled$landscapes[[current_phase]]$basins)
    plot_landscape(pooled$landscapes[[current_phase]], cfg, flow = flow, limits = pooled$limits, show_legend = FALSE,
                   title = sprintf("%s (%d basin%s)", tools::toTitleCase(current_phase), n_basins, if (n_basins != 1) "s" else ""),
                   subtitle = n_label(n_recordings))
  })
  list(panels = panels, modes = pooled$modes)
}

figure_landscapes <- function(results, cfg) {
  tilt <- landscape_row(results, "tilt", cfg)
  exercise <- landscape_row(results, "exercise", cfg)
  style <- cfg$figures$landscape
  plot <- wrap_plots(c(tilt$panels, exercise$panels), ncol = 3, guides = "collect") +
    plot_annotation(title = "Where the reconstructed states spend their time",
                    subtitle = "Top: orthostatic perturbation. Bottom: exercise. Units: baseline SD.",
                    theme = theme(plot.title = element_text(face = "bold", size = cfg$figures$base_size + 2),
                                  plot.subtitle = element_text(colour = "grey35", size = cfg$figures$base_size - 1))) &
    theme(legend.position = "bottom")
  list(plot = plot, width = 9, height = 8, info = list())
}


# ---- Supplementary figures ---------------------------------------------------------------------------------------------------------------
figure_simulation_recovery <- function(results, cfg) {
  recovery <- need(results, "sim_recovery")
  models <- figure_palette(cfg, "models")
  measure_labels <- c(ln_err_tau_p = "Error of ln $\\tau_p$", ln_err_tau_s = "Error of ln $\\tau_s$", cor_x_p = "Correlation of $x_p$\nwith the truth",
                      cor_x_s = "Correlation of $x_s$\nwith the truth", share_err = "Error in the share of\n$x_p$ variance")
  errors <- melt(recovery$table[ok == TRUE, c("scenario", "spec", names(measure_labels)), with = FALSE], id.vars = c("scenario", "spec"))
  errors[, `:=`(variable = factor(variable, names(measure_labels), tex(unname(measure_labels))),
                scenario = factor(scenario, names(SCENARIO_LABELS), tex(unname(SCENARIO_LABELS))))]
  zero_lines <- data.table(variable = factor(tex(unname(measure_labels[c("ln_err_tau_p", "ln_err_tau_s", "share_err")])), levels(errors$variable)), y = 0)
  plot <- ggplot(errors, aes(spec, value, fill = spec, colour = spec)) +
    geom_hline(data = zero_lines, aes(yintercept = y), linewidth = 0.3, linetype = 2, colour = "grey40") +
    geom_boxplot(outlier.size = 0.8, linewidth = 0.35, alpha = 0.25) + facet_grid(variable ~ scenario, scales = "free_y") +
    scale_fill_manual(values = models, guide = "none") + scale_colour_manual(values = models, guide = "none") +
    scale_x_discrete(labels = function(x) tex(spec_tex(x))) +
    labs(x = NULL, y = NULL, title = "Recovery under known truth",
         subtitle = "dashed line: no error") + theme_spark(cfg)
  list(plot = plot, width = 9, height = 8, info = list())
}

figure_sensitivity <- function(results, cfg) {
  sensitivity <- need(results, "sensitivity")
  worst_case <- sensitivity$sensitivity[, .(cor = pmin(cor_x_p, cor_x_s), ln_tau = pmax(abs(ln_ratio_tau_p), abs(ln_ratio_tau_s))),
                                        by = .(kind, setting, SubjectID)]
  worst_case[, `:=`(kind = factor(kind, c("hyperparameter", "leave_one_out"), tex(c("Hyperparameter\nvariations", "One prior\nterm removed"))),
                    setting = factor(setting, unique(setting), tex(setting_axis_label(unique(setting)))))]
  horizontal_boxes <- function(column) {
    ggplot(worst_case, aes(setting, .data[[column]])) + geom_boxplot(outlier.size = 0.8, linewidth = 0.35, fill = "grey93") +
      coord_flip() + facet_grid(kind ~ ., scales = "free_y", space = "free_y")
  }
  a <- horizontal_boxes("ln_tau") + labs(x = NULL, y = tex("Largest $|\\ln(\\mathrm{estimate}/\\mathrm{reference})|$\nof $\\tau_p$, $\\tau_s$"),
                                          title = tex("Change of the\ntime constants"), subtitle = tex("reference: $M_3$")) + theme_spark(cfg)
  b <- horizontal_boxes("cor") + labs(x = NULL, y = tex("Lowest correlation of $x_p$, $x_s$\nwith the reference"),
                                       title = tex("Change of the\nreconstructed states"), subtitle = "1 = identical trajectories") +
    scale_x_discrete(labels = NULL) + theme_spark(cfg) + theme(axis.text.y = element_blank())
  list(plot = (tag_panel(a, "A") | tag_panel(b, "B")) + plot_layout(widths = c(1.9, 1)), width = 9, height = 8, info = list())
}

figure_recovery_benchmark <- function(results, cfg) {
  benchmark <- need(results, "recovery_choice")
  performance <- melt(benchmark$choice$performance, id.vars = c("definition", "state"), measure.vars = c("success", "rho_truth", "rho_retest"))
  performance[, variable := factor(variable, c("success", "rho_truth", "rho_retest"),
                                   c("share with a defined estimate", "rank correlation with the truth", "test-retest rank correlation"))]
  states <- figure_palette(cfg, "states")
  plot <- ggplot(performance, aes(definition, value, fill = state)) + geom_col(position = "dodge", width = 0.75) + facet_wrap(~variable) +
    scale_fill_manual(values = c(p = states[["x_p"]], s = states[["x_s"]]), labels = unname(state_labels()), name = NULL) +
    labs(x = "Recovery definition", y = NULL, title = paste0("Chosen definition: ", benchmark$choice$definition,
                                                            "   |   exercise-end rule: ", benchmark$choice$end_method),
         subtitle = "benchmark on simulated recordings with known recovery constants") + theme_spark(cfg)
  list(plot = plot, width = 9, height = 4.2, info = list())
}

figure_numerical_checks <- function(results, cfg) {
  checks <- need(results, "supp_checks")
  colours <- figure_palette(cfg, "datasets")
  a <- ggplot(checks$bound, aes(eps, sup_R)) + geom_line(aes(y = bound), linetype = 2, colour = "grey40") + geom_point(size = 1.8) + geom_line() +
    scale_x_log10() + scale_y_log10() +
    labs(x = tex("$\\epsilon = \\sigma/\\sqrt{\\nu}$"), y = tex("Neglected term, $\\sup R(t)$"), title = "Threshold approximation",
         subtitle = "dashed: analytical bound") + theme_spark(cfg)
  b <- ggplot(checks$fitted_regime, aes(Dataset, eps, colour = Dataset)) + geom_jitter(width = 0.15, size = 1.8) + scale_y_log10() +
    scale_colour_manual(values = colours, guide = "none") + labs(x = NULL, y = tex("$\\epsilon$ fitted"), title = "Fitted regime",
                                                                  subtitle = "relative timing jitter per recording") + theme_spark(cfg)
  c <- ggplot(checks$calibration, aes(crit, KS_Stat)) + geom_abline(slope = 1, linetype = 2, colour = "grey40") + geom_point(size = 1.8) +
    labs(x = tex("$1.36/\\sqrt{N}$"), y = tex("KS distance $D$"), title = "Diagnostic calibration",
         subtitle = "well-specified simulations") + theme_spark(cfg)
  list(plot = figure_grid(list(a, b, c), ncol = 3), width = 9, height = 4.4, info = list())
}



# ---- Figure S7: responses by type of maneuver (PRCP) -----------------------------------------------------------------------------
figure_prcp_maneuvers <- function(results, cfg) {
  prcp <- need(results, "prcp", "run analysis/datasets/05_prcp.R")
  if (is.null(prcp$by_kind) || !nrow(prcp$by_kind)) stop("results/prcp.rds has no responses by kind: rerun analysis/datasets/05_prcp.R", call. = FALSE)
  kinds <- intersect(c("slow tilt", "rapid tilt", "stand-up", "tilt"), prcp$by_kind$kind)
  locked <- rbindlist(lapply(kinds, function(k) {
    one <- event_locked(prcp$series_all, "tilt", c(-120, 150), maneuvers = k)
    if (is.null(one) || !nrow(one)) NULL else cbind(maneuver = k, one)
  }))
  variables <- intersect(c("x_p", "x_s", "hr", "map"), unique(locked$variable))
  locked <- locked[variable %in% variables]
  summary <- locked[, .(median = median(value, na.rm = TRUE), q1 = quantile(value, 0.25, na.rm = TRUE), q3 = quantile(value, 0.75, na.rm = TRUE),
                        q10 = quantile(value, 0.10, na.rm = TRUE), q90 = quantile(value, 0.90, na.rm = TRUE), n = .N), by = .(maneuver, variable, u)][n >= 3]
  states <- figure_palette(cfg, "states")
  colour_of <- c(x_p = states[["x_p"]], x_s = states[["x_s"]], hr = figure_palette(cfg, "accent"), map = "#6B6B6B")
  summary[, colour := colour_of[variable]]
  counts <- prcp$by_kind[match(kinds, kind)]
  summary[, `:=`(variable = factor(variable, variables, locked_labels()[variables]),
                 maneuver = factor(maneuver, kinds, sprintf("%s (%d events, %d recordings)", sub("^(.)", "\\U\\1", kinds, perl = TRUE), counts$n_events, counts$n_recordings)))]
  plot <- ggplot(summary, aes(u, median)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey35", linewidth = 0.3) + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_ribbon(aes(ymin = q10, ymax = q90, fill = colour), alpha = 0.14) + geom_ribbon(aes(ymin = q1, ymax = q3, fill = colour), alpha = 0.30) +
    geom_line(aes(colour = colour), linewidth = 0.9) + scale_colour_identity() + scale_fill_identity() +
    facet_grid(variable ~ maneuver, scales = "free_y", switch = "y") +
    labs(x = "Time from the onset of the maneuver (s)", y = NULL, title = "Response by type of maneuver",
         subtitle = "median over events; bands: 25-75% and 10-90%") +
    theme_spark(cfg) + theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_markdown(angle = 90, face = "plain"))
  list(plot = plot, width = 9, height = 8.5, info = list(kinds = kinds))
}


# ---- Registry: PNG name -> builder -----------------------------------------------------------------------------------------------------
FIGURE_REGISTRY <- list(
  fig1_statistical_validation       = figure_fit_diagnostics,
  fig4_benchmarking_comparison      = figure_cwt_comparison,
  fig5_identifiability              = figure_identifiability,
  fig7_prcp                         = figure_prcp_group,
  fig8_fantasia                     = figure_fantasia_groups,
  fig9_exercise                     = figure_exercise,
  fig10_fantasia_physiology         = figure_fantasia_physiology,
  figS1_simulation_recovery         = figure_simulation_recovery,
  figS2_sensitivity                 = figure_sensitivity,
  figS3_recovery_benchmark          = figure_recovery_benchmark,
  figS4_supplement_checks           = figure_numerical_checks,
  figS5_event_detection             = figure_event_detection,
  figS6_landscapes                  = figure_landscapes,
  figS7_maneuvers                   = figure_prcp_maneuvers
)
