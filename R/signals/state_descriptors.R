# ==============================================================================
# R/signals/state_descriptors.R -- perturbation-response descriptors of the reconstructed states (Castillo-Aguilar et al., STATE preprint, section 3.3),
# with the uncertainty of the filter (supplementary "Perturbation-response descriptors"). Specification fixed before looking at any recording.
#
# State (x_p, x_s); resting region = a disc around the rest mean in coordinates whitened by the stationary SDs of the two drives; distance d_t = Euclidean
# norm of the 5-s moving average. Tolerance delta = sqrt(-2 ln 0.001) = 3.717, the exact 0.999 quantile of the stationary distance of two independent OU drives.
#   Dmax         largest distance inside the perturbation window
#   Trec         first delay after the end of the perturbation from which d stays <= delta for `hold` s (NA: no return within the follow-up)
#   Preturn(H)   Trec <= H (NA when the follow-up is shorter than H and there was no return)
#   Ptransition  the regime of the second-stage HMM (R/signals/regimes.R) in the last 30 s before the horizon differs from the resting regime
# Uncertainty: paths drawn from the filter's marginal posterior N(x_t, v_t) with AR(1) errors (rho = exp(-1 s / tau) of each drive; an approximation, the joint
# filtering law is not stored). The stationary variances that set the scale come from M3 refitted on the resting beats alone: the full-recording fit absorbs the
# perturbation into sigma^2 (validated by simulation; see the supplement).
#
#   states_on_grid()       beat-wise states and variances -> 1-s grid
#   distance_descriptors() Dmax, Trec, Preturn, d_H of one distance series
#   descriptor_draws()     the same over posterior draws (matrix, one row per draw)
#   summarise_draws()      median and 90% interval over draws; share of draws that returned within each horizon
#   regime_transition()    Ptransition with the aware HMM, K = 1..3 by BIC
# ==============================================================================
STATE_DELTA <- sqrt(-2 * log(0.001))

smooth_centered <- function(x, k) if (k <= 1) x else { y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2)); ifelse(is.na(y), x, y) }

#' Beat-wise states and their estimation variances resampled on a 1-s grid (columns t, xp, xs, vp, vs).
states_on_grid <- function(t, xp, xs, vp, vs) {
  g <- seq(ceiling(t[1]), floor(t[length(t)])); f <- function(z) approx(t, z, xout = g, rule = 2)$y
  data.table(t = g, xp = f(xp), xs = f(xs), vp = f(vp), vs = f(vs))
}

first_sustained <- function(inside, hold) { r <- rle(inside); e <- cumsum(r$lengths); s <- e - r$lengths + 1; j <- which(r$values & r$lengths >= hold)[1]; if (is.na(j)) NA_integer_ else s[j] }

#' Descriptors of one distance series d on the grid tg. pert = c(start, end) of the perturbation; H horizon (s); hold dwell time (s); rest_mask: grid points of the rest window.
distance_descriptors <- function(tg, d, pert, H, hold, delta = STATE_DELTA, rest_mask = NULL) {
  t_end <- pert[2]; after <- tg >= t_end; follow <- max(tg) - t_end
  i <- first_sustained(d[after] <= delta, hold); trec <- if (is.na(i)) NA_real_ else tg[after][i] - t_end
  c(Dmax = max(d[tg >= pert[1] & tg <= t_end]), Trec = trec, Preturn = if (!is.na(trec) && trec <= H) 1 else if (follow >= H) 0 else NA_real_,
    d_H = if (follow >= H) mean(d[tg > t_end + H - hold & tg <= t_end + H]) else NA_real_, follow_s = follow,
    rest_exceed = if (is.null(rest_mask)) NA_real_ else mean(d[rest_mask] > delta))
}

#' Whitened, smoothed distance from the rest centre: u, v are the two coordinates already centred and divided by the stationary SDs.
whitened_distance <- function(u, v, smooth_s = 5) sqrt(smooth_centered(u, smooth_s)^2 + smooth_centered(v, smooth_s)^2)

#' Descriptors over posterior draws. g: states_on_grid(); rest, pert: c(start, end) in s; V: stationary variances (Hz^2) of the two drives; tau: their time constants (s).
descriptor_draws <- function(g, rest, pert, V, tau, n_draw = 100, H = 300, hold = 15, smooth_s = 5, seed = 1) {
  r <- g$t >= rest[1] & g$t < rest[2]; if (sum(r) < 30) return(NULL)
  n <- nrow(g); cp <- mean(g$xp[r]); cs <- mean(g$xs[r]); rp <- exp(-1 / tau[1]); rs <- exp(-1 / tau[2])
  ar1 <- function(rho) { z <- rnorm(n); w <- sqrt(1 - rho^2) * z; w[1] <- z[1]; as.numeric(stats::filter(w, rho, method = "recursive")) }   # unit-variance AR(1)
  with_seed(seed, t(vapply(seq_len(n_draw), function(k) {
    u <- (g$xp + sqrt(g$vp) * ar1(rp) - cp) / sqrt(V[1]); v <- (g$xs + sqrt(g$vs) * ar1(rs) - cs) / sqrt(V[2])
    distance_descriptors(g$t, whitened_distance(u, v, smooth_s), pert, H, hold, rest_mask = r) }, numeric(6))))
}

#' Median and 90% interval over the draws; for each horizon in `horizons`, the share of draws that returned within it (NA when the follow-up is shorter and no return).
summarise_draws <- function(m, horizons = c(120, 300, 600)) {
  if (is.null(m)) return(NULL)
  q <- function(x) { x <- x[is.finite(x)]; if (!length(x)) rep(NA_real_, 3) else quantile(x, c(.5, .05, .95), names = FALSE) }
  a <- q(m[, "Dmax"]); b <- q(m[, "Trec"]); follow <- unname(m[1, "follow_s"])
  ret <- vapply(horizons, function(H) { x <- ifelse(is.finite(m[, "Trec"]) & m[, "Trec"] <= H, 1, ifelse(follow >= H, 0, NA_real_)); mean(x, na.rm = TRUE) }, numeric(1))
  c(Dmax = a[1], Dmax_lo = a[2], Dmax_hi = a[3], Trec = b[1], Trec_lo = b[2], Trec_hi = b[3], Preturn = mean(m[, "Preturn"], na.rm = TRUE),
    setNames(ret, paste0("Preturn_", horizons)), follow_s = follow, rest_exceed = mean(m[, "rest_exceed"]))
}

#' Ptransition: K by BIC (1..Kmax) of the aware HMM on the binned states; transition = the regime occupying the last `window_s` before t_end + H differs from the
#' regime of the rest window. Also the mean posterior probability of not being in the resting regime over that window. NA when the follow-up ends before the horizon.
regime_transition <- function(g, rest, pert, V, H = 300, Kmax = 3L, window_s = 30) {
  b <- bin_states(g$t, g$xp, g$xs, g$vp, g$vs); n <- nrow(b); if (n < 30) return(c(K = NA_real_, trans = NA_real_, ptrans = NA_real_))
  Y <- as.matrix(b[, .(y1, y2)]); Vt <- as.matrix(b[, .(v1, v2)])
  fits <- lapply(seq_len(Kmax), function(K) hmm_fit(Y, Vt, V, K)); bic <- vapply(fits, function(f) -2 * f$ll + f$npar * log(n), numeric(1)); K <- which.min(bic); f <- fits[[K]]
  inrest <- b$t >= rest[1] & b$t < rest[2]; post <- b$t > pert[2] + H - window_s & b$t <= pert[2] + H
  if (!any(post) || max(b$t) < pert[2] + H - 5 || sum(inrest) < 3) return(c(K = K, trans = NA_real_, ptrans = NA_real_))
  if (K == 1) return(c(K = 1, trans = 0, ptrans = 0))
  base <- as.integer(names(which.max(table(f$path[inrest]))))
  c(K = K, trans = as.numeric(as.integer(names(which.max(table(f$path[post])))) != base), ptrans = mean(1 - f$gamma[post, base]))
}
