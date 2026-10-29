# ==============================================================================
# R/signals/landscape.R -- the occupied region of the (x_p, x_s) plane as an energy landscape.
#
# The filtered states of one or several recordings, in units of the baseline (rest) SDs, are smoothed
# into a 2D density p(u, v). The quasi-potential U = -ln(p / p_max) turns high-density regions into
# basins ("attractors": where the states spend their time) and low-density regions into barriers.
# It is a descriptive transform of where the reconstructed states lie, not a fitted potential of the model.
#
#   normalized_path()       one recording -> (u, v, phase) on a 1 Hz grid, in baseline-SD units
#   density_landscape()     weighted 2D kernel density, potential U and the basins (local maxima)
#   flow_field()            mean velocity of the states by cell (arrows point towards attractors)
#   landscape_limits()      plotting limits covering the data
#   phase_modes()           mode of each phase and the distance between phase modes
# ==============================================================================


# ---- Paths ---------------------------------------------------------------------------------------------
#' Phase labels of a recording's beats: baseline/tilt/recovery (perturbation) or rest/exercise/recovery.
#' @param series list with t and either ev (windows: tilt_windows() of every analysed tilt) or t_on / t_end
series_phases <- function(series, kind = c("tilt", "exercise")) {
  kind <- match.arg(kind)
  if (kind == "tilt") return(assign_tilt_phases(series$t, series$ev$windows)$phase)
  events <- list(mode = "exercise", detected = is.finite(series$t_on), t_onset = series$t_on, t_end_peak = series$t_end)
  assign_phase(series$t, events)
}

#' One recording as a path in the (x_p, x_s) plane: 1 Hz grid, `smooth_s` moving average, coordinates
#' in baseline-SD units (u = x_p / sd_p, v = x_s / sd_s) and the phase of each point.
#' @return data.table(t, u, v, phase)
normalized_path <- function(series, kind = c("tilt", "exercise"), smooth_s = 5) {
  kind <- match.arg(kind)
  grid_p <- to_grid(series$t, series$x_p, smooth_s = smooth_s)
  grid_s <- to_grid(series$t, series$x_s, smooth_s = smooth_s)
  n <- min(nrow(grid_p), nrow(grid_s))
  beat_phase <- series_phases(series, kind)
  phase_at <- function(time) beat_phase[findInterval(time, series$t, all.inside = TRUE)]
  data.table(t = grid_p$t[1:n], u = grid_p$x[1:n] / series$sd_p, v = grid_s$x[1:n] / series$sd_s,
             phase = phase_at(grid_p$t[1:n]))
}


# ---- Density, potential, basins ------------------------------------------------------------------------------
#' Plot limits (u and v) covering the central `quantile` mass of the points, padded by 10%.
landscape_limits <- function(u, v, quantile = c(0.005, 0.995)) {
  pad <- function(x) { r <- stats::quantile(x, quantile, na.rm = TRUE); r + c(-1, 1) * 0.1 * diff(r) }
  list(u = pad(u), v = pad(v))
}

#' Weighted Gaussian kernel density of points on a regular grid. Each point carries a weight (use
#' 1 / n_points of its recording so that every recording counts equally). Bandwidth: Silverman's rule
#' per axis on the effective sample size, times `bandwidth_scale`.
#' @return list(u, v, density (n x n matrix, rows = u), bandwidth)
weighted_density <- function(u, v, weights, limits, n_grid = 90, bandwidth_scale = 1) {
  weights <- weights / sum(weights)
  effective_n <- 1 / sum(weights^2)
  weighted_sd <- function(x) sqrt(sum(weights * (x - sum(weights * x))^2))
  bandwidth <- c(u = 1.06 * weighted_sd(u) * effective_n^(-1 / 5), v = 1.06 * weighted_sd(v) * effective_n^(-1 / 5)) *
    bandwidth_scale
  grid_u <- seq(limits$u[1], limits$u[2], length.out = n_grid)
  grid_v <- seq(limits$v[1], limits$v[2], length.out = n_grid)
  kernel_u <- dnorm(outer(grid_u, u, "-") / bandwidth[["u"]]) * rep(weights, each = n_grid)   # grid x points
  kernel_v <- dnorm(outer(grid_v, v, "-") / bandwidth[["v"]])
  list(u = grid_u, v = grid_v, density = kernel_u %*% t(kernel_v) / prod(bandwidth), bandwidth = bandwidth)
}

#' Local maxima of a density matrix above `min_rel` of its peak (8-neighbourhood).
density_basins <- function(landscape, min_rel = 0.2) {
  d <- landscape$density
  n <- nrow(d); m <- ncol(d)
  padded <- matrix(-Inf, n + 2, m + 2)
  padded[2:(n + 1), 2:(m + 1)] <- d
  is_peak <- matrix(TRUE, n, m)
  for (du in -1:1) for (dv in -1:1) {
    if (du == 0 && dv == 0) next
    is_peak <- is_peak & d >= padded[(2:(n + 1)) + du, (2:(m + 1)) + dv]
  }
  is_peak <- is_peak & d >= min_rel * max(d)
  index <- which(is_peak, arr.ind = TRUE)
  data.table(u = landscape$u[index[, 1]], v = landscape$v[index[, 2]], density = d[is_peak])[order(-density)]
}

#' Density landscape of a set of points: density, potential U = -ln(p / p_max) (capped at
#' `potential_max`) as a long table, and its basins.
#' @param points data.table(u, v, weight)
density_landscape <- function(points, limits, n_grid = 90, bandwidth_scale = 1, potential_max = 6, min_rel = 0.2) {
  kde <- weighted_density(points$u, points$v, points$weight, limits, n_grid, bandwidth_scale)
  grid <- CJ(i = seq_len(n_grid), j = seq_len(n_grid))
  grid[, `:=`(u = kde$u[i], v = kde$v[j], density = kde$density[cbind(i, j)])]
  grid[, potential := pmin(-log(density / max(density)), potential_max)]
  list(grid = grid[, .(u, v, density, potential)], basins = density_basins(kde, min_rel),
       bandwidth = kde$bandwidth, limits = limits)
}


# ---- Flow field ----------------------------------------------------------------------------------------------------
#' Mean velocity (baseline SDs per second) of the paths by cell of a `n_bins` x `n_bins` grid over
#' `limits`; cells with fewer than `min_n` steps are dropped. Arrows point in the direction the states
#' move on average from that cell, so they converge on attractors.
#' @param paths data.table(id, t, u, v) on a regular 1 Hz grid
flow_field <- function(paths, limits, n_bins = 14, min_n = 20) {
  steps <- paths[order(id, t), .(u = u[-.N], v = v[-.N], du = diff(u), dv = diff(v), dt = diff(t)), by = id]
  steps <- steps[dt > 0 & dt <= 2]
  steps[, `:=`(du = du / dt, dv = dv / dt)]
  steps[, `:=`(cell_u = findInterval(u, seq(limits$u[1], limits$u[2], length.out = n_bins + 1), all.inside = TRUE),
               cell_v = findInterval(v, seq(limits$v[1], limits$v[2], length.out = n_bins + 1), all.inside = TRUE))]
  edges_u <- seq(limits$u[1], limits$u[2], length.out = n_bins + 1)
  edges_v <- seq(limits$v[1], limits$v[2], length.out = n_bins + 1)
  field <- steps[, .(n = .N, du = mean(du), dv = mean(dv)), by = .(cell_u, cell_v)][n >= min_n]
  field[, `:=`(u = (edges_u[cell_u] + edges_u[cell_u + 1]) / 2, v = (edges_v[cell_v] + edges_v[cell_v + 1]) / 2)]
  field[, speed := sqrt(du^2 + dv^2)][]
}


# ---- Phase modes ---------------------------------------------------------------------------------------------------------
#' Mode (strongest basin) of each phase's landscape and the distance between phase modes, in baseline SDs.
#' @param landscapes named list of density_landscape() results (names = phases)
phase_modes <- function(landscapes) {
  modes <- rbindlist(lapply(names(landscapes), function(phase) {
    basins <- landscapes[[phase]]$basins
    data.table(phase = phase, u = basins$u[1], v = basins$v[1], n_basins = nrow(basins))
  }))
  reference <- modes[1]
  modes[, shift_from_first := sqrt((u - reference$u)^2 + (v - reference$v)^2)][]
}


# ---- Pooled landscapes of several recordings ------------------------------------------------------------------
#' Landscapes by phase pooled over recordings: every recording contributes equally within a phase
#' (weight 1 / number of its points in that phase).
#' @param series_list named list of recording series (see analysis/datasets/05_prcp.R, 09_exercise.R)
#' @param settings cfg$figures$landscape
#' @return list(paths = all points, landscapes = by phase, limits, modes = phase modes and their shifts)
pooled_phase_landscapes <- function(series_list, kind, settings) {
  paths <- rbindlist(lapply(names(series_list), function(id) cbind(id = id, normalized_path(series_list[[id]], kind))))
  paths <- paths[!is.na(phase)]
  paths[, weight := 1 / .N, by = .(id, phase)]
  limits <- landscape_limits(paths$u, paths$v, settings$limits_quantile)
  phases <- levels(droplevels(paths$phase))
  landscapes <- lapply(setNames(phases, phases), function(p) {
    density_landscape(paths[phase == p], limits, settings$n_grid, settings$bandwidth_scale,
                      settings$potential_max, settings$min_rel_density)
  })
  list(paths = paths, landscapes = landscapes, limits = limits, modes = phase_modes(landscapes))
}
