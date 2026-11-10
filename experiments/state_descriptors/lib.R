# ==============================================================================
# experiments/state_descriptors/lib.R -- the perturbation-response descriptors of STATE (Castillo-Aguilar et al. 2026, bioRxiv 10.64898/2026.11.02.763675, section 3.3)
# computed on the latent states of production M3 (nothing in the pipeline changes). Source from the project root:  source("experiments/state_descriptors/lib.R")
#
# State A_t = (x_p, x_s); rest region R_rest = a disc around the rest mean in whitened coordinates (each drive divided by its stationary SD, sqrt(V_p), sqrt(V_s), V_p = lambda_P sigma^2,
# V_s = sigma^2, from the M3 fit); distance d_t = Euclidean norm. Descriptors, all fixed BEFORE looking at a recording (STATE rule: delta, tau, H are part of the definition):
#   Dmax      maximum distance inside the perturbation window
#   Trec      first delay after the end of the perturbation at which d stays <= delta for `hold` seconds (NA = not within the follow-up)
#   Preturn   Trec <= H (NA when the follow-up is shorter than H and there was no return before it)
#   d_H       distance at the horizon (mean over the last `hold` s before t_end + H)
#   Ptransition  the persistent regions other than rest are the regimes of the second-stage HMM of experiments/hmm_states (aware emission, number of regimes by BIC, K = 1..3); the post-perturbation state is the
#             regime occupying the last 30 s before t_end + H; transition = that regime differs from the resting regime (the regime of the rest window); also the posterior probability of leaving the resting regime.
# Two estimators of the same descriptors:
#   plug-in   the filtered mean path; scale = SD of the path in the rest window; delta = 0.999 quantile of d in the rest window (the literal STATE rule). The filtered path is shrunk
#             toward its mean, so this scale is smaller than the stationary one and the descriptors are in "filtered rest SDs".
#   draws     paths drawn from the filter's marginal posterior N(x_t, v_t) with AR(1) errors (rho = exp(-dt / tau) of each drive, an approximation: the joint filtering law is not stored),
#             scale = stationary SD of the model, delta = sqrt(-2 ln 0.001) = 3.717, the exact 0.999 quantile of the stationary distance of the two independent OU drives (chi-square 2).
#             Descriptors are summarised by their median and 90% interval over draws, and Preturn by the share of draws that returned.
# ==============================================================================
suppressMessages(library(data.table)); source("experiments/hmm_states/lib.R")        # second-stage HMM (bin_states, hmm_fit)
DELTA_MODEL <- sqrt(-2 * log(0.001))                                              # 3.717: 0.999 quantile of the stationary distance of two independent OU drives

ma <- function(x, k) if (k <= 1) x else { y <- as.numeric(stats::filter(x, rep(1 / k, k), sides = 2)); ifelse(is.na(y), x, y) }

#' Resample the beat-wise states onto a 1-s grid (columns t, xp, xs, vp, vs).
to_grid1 <- function(t, xp, xs, vp, vs) { g <- seq(ceiling(t[1]), floor(t[length(t)])); f <- function(z) approx(t, z, xout = g, rule = 2)$y; data.table(t = g, xp = f(xp), xs = f(xs), vp = f(vp), vs = f(vs)) }

first_hold <- function(inside, hold) { r <- rle(inside); e <- cumsum(r$lengths); s <- e - r$lengths + 1; j <- which(r$values & r$lengths >= hold)[1]; if (is.na(j)) NA_integer_ else s[j] }

#' Descriptors of one distance series d(t) on the grid `tg`. win_pert = c(start, end) of the perturbation; H horizon (s); hold dwell time (s).
descr_d <- function(tg, d, win_pert, H, hold, delta, rest_mask = NULL) {
  t_end <- win_pert[2]; inpert <- tg >= win_pert[1] & tg <= t_end; after <- tg >= t_end; L <- max(tg) - t_end
  i <- first_hold(d[after] <= delta, hold); trec <- if (is.na(i)) NA_real_ else tg[after][i] - t_end
  ret <- if (!is.na(trec) && trec <= H) 1 else if (L >= H) 0 else NA_real_
  dh <- if (L >= H) mean(d[tg > t_end + H - hold & tg <= t_end + H]) else NA_real_
  c(Dmax = max(d[inpert]), Trec = trec, Preturn = ret, d_H = dh, follow_s = L, rest_exceed = if (is.null(rest_mask)) NA_real_ else mean(d[rest_mask] > delta)) }

#' Whitened, smoothed distance from the rest centre. u, v: coordinates already centred and scaled.
dist_of <- function(u, v, smooth_s = 5) sqrt(ma(u, smooth_s)^2 + ma(v, smooth_s)^2)

#' Plug-in descriptors from the filtered mean path. g: output of to_grid1; rest / pert: c(start, end) windows (s).
descr_plugin <- function(g, rest, pert, H = 300, hold = 15, smooth_s = 5, q = 0.999) {
  r <- g$t >= rest[1] & g$t < rest[2]; if (sum(r) < 30) return(NULL)
  u <- (g$xp - mean(g$xp[r])) / sd(g$xp[r]); v <- (g$xs - mean(g$xs[r])) / sd(g$xs[r]); d <- dist_of(u, v, smooth_s)
  delta <- unname(quantile(d[r], q)); c(descr_d(g$t, d, pert, H, hold, delta), delta = delta) }

#' Rest-window variance of the TRUE state: variance of the filtered mean plus the mean estimation variance (law of total variance), per drive.
rest_V <- function(g, r) c(var(g$xp[r]) + mean(g$vp[r]), var(g$xs[r]) + mean(g$vs[r]))

#' Draws-based descriptors. V = c(V_p, V_s) stationary variances (Hz^2), tau = c(tau_p, tau_s) (s). Returns the per-draw matrix of descriptors.
descr_draws <- function(g, rest, pert, V, tau, n_draw = 100, H = 300, hold = 15, smooth_s = 5, delta = DELTA_MODEL, seed = 1) {
  r <- g$t >= rest[1] & g$t < rest[2]; if (sum(r) < 30) return(NULL); set.seed(seed); n <- nrow(g); cp <- mean(g$xp[r]); cs <- mean(g$xs[r])
  if (is.null(V)) V <- rest_V(g, r)                                              # scale from the rest window alone (law of total variance), independent of the perturbation
  ar1 <- function(rho) { z <- rnorm(n); w <- sqrt(1 - rho^2) * z; w[1] <- z[1]; as.numeric(stats::filter(w, rho, method = "recursive")) }   # unit-variance AR(1)
  rp <- exp(-1 / tau[1]); rs <- exp(-1 / tau[2])
  out <- t(vapply(seq_len(n_draw), function(k) { u <- (g$xp + sqrt(g$vp) * ar1(rp) - cp) / sqrt(V[1]); v <- (g$xs + sqrt(g$vs) * ar1(rs) - cs) / sqrt(V[2])
    descr_d(g$t, dist_of(u, v, smooth_s), pert, H, hold, delta, r) }, numeric(6)))
  out }

#' Median and 90% interval over the draws; Preturn = share of draws that returned (among those with a determinate answer).
summ_draws <- function(m) { if (is.null(m)) return(NULL); f <- function(x) { x <- x[is.finite(x)]; if (!length(x)) c(NA, NA, NA) else quantile(x, c(.5, .05, .95), names = FALSE) }
  a <- f(m[, "Dmax"]); b <- f(m[, "Trec"])
  c(Dmax = a[1], Dmax_lo = a[2], Dmax_hi = a[3], Trec = b[1], Trec_lo = b[2], Trec_hi = b[3], Trec_ret_share = mean(is.finite(m[, "Trec"])), Preturn = mean(m[, "Preturn"], na.rm = TRUE),
    d_H = median(m[, "d_H"], na.rm = TRUE), follow_s = m[1, "follow_s"], rest_exceed = mean(m[, "rest_exceed"])) }

#' Ptransition (STATE section 3.3): does the system sit in a persistent region other than rest at the horizon? V: stationary variances (rest-calibrated), as in descr_draws.
#' Returns K (regimes by BIC), trans (0/1; NA if the recording ends before t_end + H) and ptrans (mean posterior probability of not being in the resting regime over the last 30 s before the horizon).
ptransition <- function(g, rest, pert, V, H = 300, Kmax = 3L, window_s = 30) {
  b <- bin_states(g$t, g$xp, g$xs, g$vp, g$vs); n <- nrow(b); if (n < 30) return(NULL); Y <- as.matrix(b[, .(y1, y2)]); Vt <- as.matrix(b[, .(v1, v2)])
  fits <- lapply(seq_len(Kmax), function(K) hmm_fit(Y, Vt, V, K, "aware")); bic <- vapply(fits, function(f) -2 * f$ll + f$npar * log(n), numeric(1)); K <- which.min(bic); f <- fits[[K]]
  inrest <- b$t >= rest[1] & b$t < rest[2]; post <- b$t > pert[2] + H - window_s & b$t <= pert[2] + H
  if (!any(post) || max(b$t) < pert[2] + H - 5 || sum(inrest) < 3) return(c(K = K, trans = NA_real_, ptrans = NA_real_))
  if (K == 1) return(c(K = 1, trans = 0, ptrans = 0))
  base <- as.integer(names(which.max(table(f$path[inrest])))); c(K = K, trans = as.numeric(as.integer(names(which.max(table(f$path[post])))) != base), ptrans = mean(1 - f$gamma[post, base])) }
