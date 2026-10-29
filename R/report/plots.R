# ==============================================================================
# R/report/plots.R -- shared plot building blocks (ggplot2 + patchwork) used by R/report/figures.R.
# Everything that controls the look reads cfg$figures (config.R, section 6).
#
#   style          figure_palette(), theme_spark(), tag_panel(), figure_grid()
#   time series    series_stack()       one panel per variable on a shared time axis
#   event-locked   event_locked(), plot_event_locked()   recordings aligned on the event
#   state plane    plot_state_plane()   paths of (x_p, x_s) in baseline-SD units
#
# Alignment rule: variables that share a time axis are facets of ONE ggplot (facet_grid), so their panel
# areas coincide by construction; separate plots are combined with patchwork (never wrap_elements()).
# ==============================================================================

state_labels <- function() tex_named(c(x_p = "vagal-like $x_p$", x_s = "sympathetic-like $x_s$"))   # function: symbols.R loads after this file


# ---- Style ---------------------------------------------------------------------------------------------
#' A palette group of cfg$figures$palette (groups, datasets, states, phases, models) or the accent colour.
figure_palette <- function(cfg, which = NULL) {
  palette <- cfg$figures$palette
  if (is.null(which)) palette else palette[[which]]
}

theme_spark <- function(cfg) {
  size <- cfg$figures$base_size
  theme_minimal(base_size = size, base_family = figure_font(cfg)) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey91", linewidth = 0.25),
      axis.line = element_line(colour = "grey45", linewidth = 0.3),
      axis.ticks = element_line(colour = "grey45", linewidth = 0.3),
      # Every text element renders HTML (ggtext), so labels written with tex() show italic symbols,
      # sub- and superscripts and Greek letters
      axis.title.x = element_markdown(size = size, margin = margin(t = 4)),
      axis.title.y = element_markdown(size = size, angle = 90, margin = margin(r = 4)),
      axis.text = element_markdown(size = size - 1, colour = "grey20"),
      axis.text.x = element_markdown(size = size - 1, colour = "grey20", margin = margin(t = 2)),
      axis.text.y = element_markdown(size = size - 1, colour = "grey20", margin = margin(r = 2)),
      # The same for every position: a more specific element of the base theme (or of a newer ggplot2) would otherwise win over axis.text.x / axis.text.y
      axis.text.x.bottom = element_markdown(size = size - 1, colour = "grey20", margin = margin(t = 2)),
      axis.text.x.top = element_markdown(size = size - 1, colour = "grey20", margin = margin(b = 2)),
      axis.text.y.left = element_markdown(size = size - 1, colour = "grey20", margin = margin(r = 2)),
      axis.text.y.right = element_markdown(size = size - 1, colour = "grey20", margin = margin(l = 2)),
      plot.title = element_markdown(face = "bold", size = size + 1, hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_markdown(colour = "grey35", size = size - 1, hjust = 0, lineheight = 1.05, margin = margin(b = 5)),
      plot.title.position = "panel",
      strip.text = element_markdown(face = "bold", size = size - 0.5, margin = margin(3, 3, 3, 3)),
      strip.text.x = element_markdown(face = "bold", size = size - 0.5, margin = margin(3, 3, 3, 3)),
      strip.text.y = element_markdown(face = "bold", size = size - 0.5, angle = -90, margin = margin(3, 3, 3, 3)),
      strip.text.x.top = element_markdown(face = "bold", size = size - 0.5, margin = margin(3, 3, 3, 3)),
      strip.text.x.bottom = element_markdown(face = "bold", size = size - 0.5, margin = margin(3, 3, 3, 3)),
      strip.text.y.left = element_markdown(face = "bold", size = size - 0.5, angle = 90, margin = margin(3, 3, 3, 3)),
      strip.text.y.right = element_markdown(face = "bold", size = size - 0.5, angle = -90, margin = margin(3, 3, 3, 3)),
      strip.background = element_rect(fill = "grey95", colour = NA),
      legend.key.size = unit(0.38, "cm"), legend.title = element_markdown(size = size - 1),
      legend.text = element_markdown(size = size - 1), legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(0.15, "cm"),
      plot.tag = element_text(face = "bold", size = size + 4),
      plot.margin = margin(4, 8, 4, 4)
    )
}

#' Panel letter in the top-left corner of a plot.
tag_panel <- function(plot, letter) plot + labs(tag = letter)

#' Several plots in a grid: patchwork keeps the panel areas of one row / column aligned.
#' @param plots list of ggplots (NULL entries are dropped); `letters` tags them A, B, ... in order
figure_grid <- function(plots, ncol, widths = NULL, heights = NULL, letters = TRUE) {
  plots <- Filter(Negate(is.null), plots)
  if (letters) plots <- mapply(tag_panel, plots, LETTERS[seq_along(plots)], SIMPLIFY = FALSE)
  wrap_plots(plots, ncol = ncol, widths = widths, heights = heights)
}

#' Dashed vertical lines at the finite event times.
add_event_lines <- function(plot, times) {
  times <- times[is.finite(times)]
  if (!length(times)) return(plot)
  plot + geom_vline(xintercept = times, linetype = 2, colour = "grey30", linewidth = 0.3)
}


# ---- Time series ---------------------------------------------------------------------------------------
#' Shaded intervals of the perturbation (tilts or the exercise bout) of one recording.
#' @return data.table(xmin, xmax, label, phase)
perturbation_intervals <- function(series, kind) {
  if (kind == "tilt") {
    w <- series$ev$windows
    kinds <- if ("kind" %in% names(w)) w$kind else rep("tilt", nrow(w))
    stand <- kinds == "stand-up"
    label <- character(nrow(w)); label[!stand] <- as.character(seq_len(sum(!stand))); label[stand] <- paste0("S", seq_len(sum(stand)))
    return(data.table(xmin = w$t_onset, xmax = w$t_tilt_end, label = label, phase = ifelse(stand, "stand-up", "tilt")))
  }
  data.table(xmin = series$t_on, xmax = series$t_end, label = "", phase = "exercise")
}

#' Stacked panels on a shared time axis: latent drives (with 95% bands), heart rate and, when present,
#' MAP and tilt angle. The perturbation is shaded; tilts are numbered above the drives.
#' @param kind "tilt" or "exercise"
series_stack <- function(series, kind, cfg, title = NULL, subtitle = NULL) {
  panel_names <- c("Drives (Hz)", "HR (bpm)", "MAP\n(mmHg)", "Tilt angle\n(\u00b0)")
  band <- function(variance) if (is.null(variance)) 0 else 1.96 * sqrt(variance)
  states <- rbind(
    data.table(t = series$t, state = "x_p", v = series$x_p, lo = series$x_p - band(series$vp), hi = series$x_p + band(series$vp)),
    data.table(t = series$t, state = "x_s", v = series$x_s, lo = series$x_s - band(series$vs), hi = series$x_s + band(series$vs)))
  states[, panel := panel_names[1]]
  others <- list(data.table(t = series$t, state = "other", v = 60 / series$rr, lo = NA_real_, hi = NA_real_, panel = panel_names[2]))
  if (!is.null(series$map)) others[[length(others) + 1]] <- data.table(t = series$t, state = "other", v = series$map, lo = NA_real_, hi = NA_real_, panel = panel_names[3])
  if (!is.null(series$angle)) others[[length(others) + 1]] <- data.table(t = series$t, state = "other", v = series$angle, lo = NA_real_, hi = NA_real_, panel = panel_names[4])
  long <- rbind(states, rbindlist(others))
  long[, panel := factor(panel, panel_names)]
  long <- long[panel %in% unique(panel)]
  intervals <- perturbation_intervals(series, kind)
  first_panel <- factor(panel_names[1], levels(long$panel))
  state_colours <- figure_palette(cfg, "states")
  phase_colours <- c(tilt = figure_palette(cfg, "phases")[["tilt"]], `stand-up` = figure_palette(cfg, "phases")[["recovery"]], exercise = figure_palette(cfg, "phases")[["exercise"]])

  plot <- ggplot(long, aes(t, v)) +
    geom_rect(data = intervals, aes(xmin = xmin, xmax = xmax, fill = phase), ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha = 0.13) +
    geom_ribbon(data = long[state != "other"], aes(ymin = lo, ymax = hi, fill = state), alpha = 0.18) +
    geom_line(data = long[state != "other"], aes(colour = state), linewidth = 0.3) +
    geom_line(data = long[state == "other"], colour = figure_palette(cfg, "accent"), linewidth = 0.35) +
    facet_grid(panel ~ ., scales = "free_y", switch = "y") +
    scale_colour_manual(values = state_colours, labels = state_labels(), name = NULL) +
    scale_fill_manual(values = c(phase_colours, state_colours), guide = "none") +
    labs(x = "Time (s)", y = NULL, title = title, subtitle = subtitle) +
    theme_spark(cfg) +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_markdown(angle = 90, face = "plain"),
          legend.position = "top", legend.justification = "left")
  if (kind == "tilt" && nrow(intervals) > 1) {
    labels <- data.table(t = (intervals$xmin + intervals$xmax) / 2, v = Inf, label = intervals$label, panel = first_panel)
    plot <- plot + geom_text(data = labels, mapping = aes(x = t, y = v, label = label), vjust = 1.5, size = 2.9,
                             colour = "grey25", inherit.aes = FALSE)
  }
  plot
}


# ---- Event-locked averages ----------------------------------------------------------------------------------
#' Recordings aligned on their events. For every event the variables are expressed relative to the event's
#' own baseline: the states and (optionally) the CWT band powers in baseline SDs, heart rate and MAP as
#' changes. A recording's curve is the mean over its events; `locked` has one row per recording and time.
#' @param series_list named list of series (analysis/datasets/05_prcp.R or 09_exercise.R)
#' @param kind "tilt": aligned on each tilt onset, shown until the tilt returns;
#'             "exercise": aligned on the end of exercise (recovery), rest as baseline
#' @param window_s time range relative to the alignment point
#' @param with_cwt include the CWT LF / HF power (series carry `cwt`)
#' @param maneuvers kinds to include ("slow tilt", "rapid tilt", "stand-up", "tilt"); NULL: all the tilts
#' @return data.table(id, u, variable, value) with variable in x_p, x_s, hr, map, LF, HF
event_locked <- function(series_list, kind = c("tilt", "exercise"), window_s = c(-120, 150), with_cwt = FALSE, smooth_s = 5, maneuvers = NULL) {
  kind <- match.arg(kind)
  grid <- seq(window_s[1], window_s[2], by = 1)
  one_event <- function(series, t_align, t_base_from, t_base_to, t_limit) {
    series_to_z <- function(t, x, z = TRUE) {
      g <- to_grid(t, x, smooth_s = smooth_s)
      base <- g[t >= t_base_from & t < t_base_to, x]
      if (length(base) < 10) return(NULL)
      value <- approx(g$t, g$x, xout = t_align + grid, rule = 1)$y
      value[t_align + grid > t_limit] <- NA
      centred <- value - mean(base)
      data.table(u = grid, value = if (z) centred / sd(base) else centred)
    }
    pieces <- list(x_p = series_to_z(series$t, series$x_p), x_s = series_to_z(series$t, series$x_s),
                   hr = series_to_z(series$t, 60 / series$rr, z = FALSE),
                   map = if (!is.null(series$map)) series_to_z(series$t, series$map, z = FALSE))
    if (with_cwt && !is.null(series$cwt)) {
      pieces$LF <- series_to_z(series$cwt$t, series$cwt$LF)
      pieces$HF <- series_to_z(series$cwt$t, series$cwt$HF)
    }
    rbindlist(lapply(names(pieces), function(v) if (!is.null(pieces[[v]])) cbind(variable = v, pieces[[v]])))
  }
  events <- lapply(names(series_list), function(id) {
    series <- series_list[[id]]
    rows <- if (kind == "tilt") {
      w <- series$ev$windows
      if ("kind" %in% names(w)) w <- if (is.null(maneuvers)) w[w$kind != "stand-up", ] else w[w$kind %in% maneuvers, ]   # default: the tilts
      lapply(seq_len(nrow(w)), function(i) one_event(series, w$t_onset[i], w$t_base_from[i], w$t_onset[i], w$t_tilt_end[i]))
    } else {
      base_from <- series$t[1] + 20
      list(one_event(series, series$t_end, base_from, series$t_on - 5, Inf))
    }
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) return(NULL)
    cbind(id = id, rbindlist(rows)[, .(value = mean(value, na.rm = TRUE)), by = .(variable, u)][is.finite(value)])
  })
  rbindlist(events)
}

locked_labels <- function(reference = "baseline") tex_named(c(x_p = paste0("$x_p$\n(", reference, " SD)"), x_s = paste0("$x_s$\n(", reference, " SD)"), hr = "$\\Delta$HR\n(bpm)",
                                         map = "$\\Delta$MAP\n(mmHg)", LF = "LF power\n(SD)", HF = "HF power\n(SD)"))

#' Median and interquartile band over recordings of event-locked curves, individual recordings in grey.
#' @param variables which variables to show, in order; @param alignment label of time zero
plot_event_locked <- function(locked, cfg, variables = c("x_p", "x_s", "hr", "map"), alignment = "tilt onset",
                              title = NULL, subtitle = NULL, show_individual = TRUE, individual_max = 30, highlight = NULL, reference = "baseline") {
  data <- locked[variable %in% variables]
  data[, variable := factor(variable, variables, locked_labels(reference)[variables])]
  many <- uniqueN(data$id) > individual_max                    # many recordings: summarize (quantile bands) instead of drawing every line
  summary <- data[, .(median = median(value, na.rm = TRUE), q1 = quantile(value, 0.25, na.rm = TRUE), q3 = quantile(value, 0.75, na.rm = TRUE),
                      q10 = quantile(value, 0.10, na.rm = TRUE), q90 = quantile(value, 0.90, na.rm = TRUE), n = .N), by = .(variable, u)]
  summary <- summary[n >= max(2, 0.5 * uniqueN(data$id))]
  colour_of <- c(x_p = figure_palette(cfg, "states")[["x_p"]], x_s = figure_palette(cfg, "states")[["x_s"]],
                 hr = figure_palette(cfg, "accent"), map = "#6B6B6B", LF = figure_palette(cfg, "states")[["x_s"]],
                 HF = figure_palette(cfg, "states")[["x_p"]])
  summary[, colour := colour_of[match(as.character(variable), locked_labels(reference))] ]
  plot <- ggplot(summary, aes(u, median)) + geom_vline(xintercept = 0, linetype = 2, colour = "grey35", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.25)
  if (many) {
    plot <- plot + geom_ribbon(aes(ymin = q10, ymax = q90, fill = colour), alpha = 0.14) + geom_ribbon(aes(ymin = q1, ymax = q3, fill = colour), alpha = 0.30)
    if (!is.null(highlight) && as.character(highlight) %in% as.character(data$id))
      plot <- plot + geom_line(data = data[as.character(id) == as.character(highlight)], aes(u, value, group = variable), colour = "grey20", linewidth = 0.35)
  } else {
    if (show_individual) plot <- plot + geom_line(data = data, aes(u, value, group = interaction(id, variable)), colour = "grey70", linewidth = 0.25, alpha = 0.6)
    plot <- plot + geom_ribbon(aes(ymin = q1, ymax = q3, fill = colour), alpha = 0.22)
  }
  plot +
    geom_line(aes(colour = colour), linewidth = 0.9) +
    scale_colour_identity() + scale_fill_identity() +
    facet_grid(variable ~ ., scales = "free_y", switch = "y") +
    labs(x = paste0("Time from ", alignment, " (s)"), y = NULL, title = title, subtitle = subtitle) +
    theme_spark(cfg) +
    theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_markdown(angle = 90, face = "plain"))
}

# ---- State plane ---------------------------------------------------------------------------------------
#' Paths in the (x_p, x_s) plane in baseline-SD units, one thin line per recording coloured by phase, with an
#' optional bold median path.
#' @param paths data.table(id, u, v, phase); @param median_path data.table(u, v, phase) or NULL
plot_state_plane <- function(paths, cfg, median_path = NULL, title = NULL, subtitle = NULL, legend = TRUE, individual_max = 30, highlight = NULL, reference = "baseline") {
  colours <- figure_palette(cfg, "phases")
  base <- ggplot() + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.25) + geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.25)
  many <- uniqueN(paths$id) > individual_max
  if (many) {
    # Many recordings: one panel per phase on common axes (the central 99% of the pooled points), so that the panels show where each phase sits in the
    # whole plane: density of the points of the phase (coloured hexagons) over the pooled density of all phases (grey), ellipses that hold 50% and 90% of
    # the points, the median path through the phase (open circle: start) and the highlighted recording
    present <- intersect(names(colours), as.character(unique(paths$phase)))
    paths[, phase := factor(as.character(phase), present)]
    window <- function(x) { q <- quantile(x, c(0.005, 0.995), na.rm = TRUE); q + c(-1, 1) * 0.08 * diff(q) }
    xlim <- window(paths$u); ylim <- window(paths$v)
    inside <- paths[u >= xlim[1] & u <= xlim[2] & v >= ylim[1] & v <= ylim[2]]
    context <- rbindlist(lapply(present, function(p) cbind(inside[, .(u, v)], phase = factor(p, present))))
    bin <- c(diff(xlim) / 28, diff(ylim) / 28)
    paths[, progress := { r <- range(t); if (diff(r) > 0) (t - r[1]) / diff(r) else 0 }, by = .(id, phase)]
    trajectory <- paths[, .(u = median(u, na.rm = TRUE), v = median(v, na.rm = TRUE)), by = .(phase, bin = pmin(floor(progress * 12), 11))][order(phase, bin)]
    phase_title <- function(x) tools::toTitleCase(as.character(x))
    plot <- base + geom_hex(data = context, aes(u, v), binwidth = bin, fill = "grey91") +
      lapply(present, function(p) geom_hex(data = inside[phase == p], aes(u, v, alpha = after_stat(count)), binwidth = bin, fill = colours[[p]])) +
      scale_alpha_continuous(range = c(0.25, 1), trans = "log10", guide = "none") +
      stat_ellipse(data = paths, aes(u, v, colour = phase, group = phase), type = "norm", level = 0.5, linewidth = 0.8) +
      stat_ellipse(data = paths, aes(u, v, colour = phase, group = phase), type = "norm", level = 0.9, linewidth = 0.5, linetype = 2)
    if (!is.null(highlight) && as.character(highlight) %in% as.character(paths$id))
      plot <- plot + geom_path(data = paths[as.character(id) == as.character(highlight)], aes(u, v, group = interaction(id, phase)), colour = "grey25", linewidth = 0.3)
    return(plot +
      geom_path(data = trajectory, aes(u, v, group = phase), colour = "white", linewidth = 2, lineend = "round") +
      geom_path(data = trajectory, aes(u, v, group = phase, colour = phase), linewidth = 1.1, lineend = "round", arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed")) +
      geom_point(data = trajectory[, .SD[1], by = phase], aes(u, v), shape = 21, fill = "white", colour = "grey15", size = 2.4, stroke = 0.8) +
      facet_wrap(~phase, nrow = 1, labeller = labeller(phase = phase_title)) + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
      scale_colour_manual(values = colours, breaks = present, guide = "none") +
      labs(x = tex(paste0("$x_p$ (", reference, " SD)")), y = tex(paste0("$x_s$ (", reference, " SD)")), title = title, subtitle = subtitle) +
      theme_spark(cfg) + theme(legend.position = "none"))
  }
  plot <- base + geom_path(data = paths, aes(u, v, group = interaction(id, phase), colour = phase), linewidth = 0.35, alpha = 0.55)
  if (!is.null(median_path) && nrow(median_path)) {
    ends <- median_path[, .SD[.N], by = phase]
    plot <- plot +
      geom_path(data = median_path, aes(u, v, group = phase), colour = "white", linewidth = 1.7, lineend = "round") +
      geom_path(data = median_path, aes(u, v, group = phase, colour = phase), linewidth = 0.95, lineend = "round") +
      geom_point(data = median_path[1], aes(u, v), shape = 21, fill = "white", colour = figure_palette(cfg, "accent"), size = 2.4, stroke = 0.7) +
      geom_point(data = median_path[.N], aes(u, v), shape = 24, fill = "white", colour = figure_palette(cfg, "accent"), size = 2.4, stroke = 0.7)
  }
  plot +
    scale_colour_manual(values = colours, breaks = intersect(names(colours), as.character(unique(paths$phase))), name = NULL) +
    labs(x = tex(paste0("$x_p$ (", reference, " SD)")), y = tex(paste0("$x_s$ (", reference, " SD)")), title = title, subtitle = subtitle) +
    theme_spark(cfg) + theme(legend.position = if (legend) "bottom" else "none")
}

# ---- Density landscape ---------------------------------------------------------------------------------
#' Density landscape of the (x_p, x_s) plane: filled quasi-potential U = -ln(p / p_max) (dark = basin),
#' contours, basin markers and optional flow arrows. Descriptive map of where the states lie.
#' @param landscape density_landscape() result; @param flow flow_field() result or NULL
plot_landscape <- function(landscape, cfg, flow = NULL, limits = landscape$limits, title = NULL, subtitle = NULL, show_legend = TRUE) {
  style <- cfg$figures$landscape
  plot <- ggplot() +
    geom_raster(data = landscape$grid, aes(u, v, fill = potential), interpolate = TRUE) +
    scale_fill_viridis_c(option = style$palette, direction = 1, limits = c(0, style$potential_max),
                         name = tex("$U = -\\ln(p/p_{max})$"), guide = if (show_legend) "colourbar" else "none")
  if (isTRUE(style$show_contours)) {
    plot <- plot + geom_contour(data = landscape$grid, aes(u, v, z = potential), colour = "white", alpha = 0.35, linewidth = 0.25, bins = 7)
  }
  if (isTRUE(style$show_flow) && !is.null(flow) && nrow(flow)) {
    cell <- c(diff(limits$u), diff(limits$v)) / style$flow_bins
    scale <- 0.85 / max(flow$speed)
    plot <- plot + geom_segment(data = flow, aes(x = u - du * scale * cell[1] / 2, y = v - dv * scale * cell[2] / 2,
                                                 xend = u + du * scale * cell[1] / 2, yend = v + dv * scale * cell[2] / 2),
                                arrow = grid::arrow(length = unit(0.13, "cm"), type = "closed"), colour = "white", alpha = 0.9, linewidth = 0.4)
  }
  if (isTRUE(style$show_basins) && nrow(landscape$basins)) {
    plot <- plot + geom_point(data = landscape$basins, aes(u, v), shape = 21, fill = NA, colour = "white", size = 3.4, stroke = 0.9) +
      geom_point(data = landscape$basins, aes(u, v), shape = 3, colour = "white", size = 1.6, stroke = 0.7)
  }
  plot +
    coord_cartesian(xlim = limits$u, ylim = limits$v, expand = FALSE) +
    labs(x = tex("$x_p$ (baseline SD)"), y = tex("$x_s$ (baseline SD)"), title = title, subtitle = subtitle) +
    theme_spark(cfg) + theme(panel.grid = element_blank())
}
