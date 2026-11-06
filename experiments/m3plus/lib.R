# ==============================================================================
# experiments/m3plus/lib.R -- EXPERIMENTAL helpers for the M3 extensions (oscillator, rate-dependent noise, robust weights, smoother, regressors,
# held-out predictive score). Load after analysis/00_setup.R:  source("experiments/m3plus/lib.R")
# ==============================================================================
suppressMessages(Rcpp::sourceCpp("experiments/m3plus/ext_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))

#' The pipeline's artifact multiplier (same formula as the engine and fast_filter.cpp).
gate_mult <- function(dy, ctl) {
  n <- length(dy); jb <- c(0, abs(diff(dy)) / pmax(dy[-n], 1e-6)); jf <- c(abs(diff(dy)) / pmax(dy[-n], 1e-6), 0)
  pmin(1 + (pmax(jb, jf) / ctl$jump_threshold)^ctl$jump_power, 1e10)
}

#' One pass of the extended filter. osc = list(alpha, omega, lo) or NULL; qs = per-beat scale of the drives' process noise; Xextra = extra regressors (N x m).
ext_run <- function(dy, ctl, kp, ks, lp, lR, osc = NULL, qs = NULL, Xextra = NULL, gm = NULL, full = FALSE, smooth = FALSE, beta_fixed = numeric(0), sig2_fixed = -1) {
  n <- length(dy); if (is.null(gm)) gm <- gate_mult(dy, ctl); if (is.null(qs)) qs <- rep(1, n)
  X <- cbind(dy, if (!is.null(Xextra)) Xextra)
  ext_filter(dy, gm, qs, X, kp, ks, lp, lR, if (is.null(osc)) 0 else osc$alpha, if (is.null(osc)) 0 else osc$omega, if (is.null(osc)) 0 else osc$lo,
             ctl$lik_gate_max, full, smooth, beta_fixed, sig2_fixed)
}

#' Dominant frequency (rad/s) of the global wavelet spectrum in the respiratory range 0.12-0.45 Hz (the oscillator's centre).
peak_omega <- function(dy, priors, lo = 0.12, hi = 0.45) {
  g <- global_spectrum_ext(dy, priors); k <- g$f >= lo & g$f <= hi; if (!any(k)) return(2 * pi * 0.25); 2 * pi * g$f[k][which.max(g$p[k])]
}
global_spectrum_ext <- function(dy, priors) {
  time_cum <- cumsum(dy); fs <- priors$cwt_fs; t_grid <- seq(min(time_cum), max(time_cum), by = 1 / fs)
  hr <- splinefun(time_cum, 1 / dy, method = "monoH.FC")(t_grid)
  invisible(capture.output(w <- biwavelet::wt(cbind(t_grid, hr - mean(hr)), do.sig = FALSE)))
  pm <- w$power; for (i in seq_len(ncol(pm))) pm[w$period > w$coi[i], i] <- NA
  p <- rowMeans(pm, na.rm = TRUE); f <- 1 / w$period; o <- order(f); keep <- is.finite(p[o]) & p[o] > 0; list(f = f[o][keep], p = p[o][keep])
}

# ---- Extended simulator -------------------------------------------------------------------------------------------------------------------------
#' Like simulate_rr() plus: a damped oscillator in the vagal pathway (osc = list(alpha, f_hz, lo), V_o = lo * sig2), noise whose scale follows the rate
#' ((rate / nu0)^(-het)), deterministic step responses (steps = list(t0, t1, amp_p, amp_s, tau_p, tau_s)) and ectopic-like artifacts.
sim_ext <- function(duration_s, kp, ks, lp, sig2, lR, nu0 = 1.2, h = 0.02, seed, osc = NULL, het = 0, steps = NULL, artifacts = 0, nu_min = 0.3, m_p_fun = NULL, m_s_fun = NULL) {
  set.seed(seed); n <- ceiling(duration_s / h); time <- (0:n) * h
  a_p <- exp(-kp * h); a_s <- exp(-ks * h); sd_p <- sqrt(lp * sig2 * (1 - a_p^2)); sd_s <- sqrt(sig2 * (1 - a_s^2))
  Vo <- if (is.null(osc)) 0 else osc$lo * sig2; ea <- if (is.null(osc)) 0 else exp(-osc$alpha * h); om <- if (is.null(osc)) 0 else 2 * pi * osc$f_hz; cc <- cos(om * h); ss <- sin(om * h); sd_o <- sqrt(Vo * (1 - ea^2))
  e <- matrix(rnorm(4 * n), n); eph <- rnorm(n + 1)
  xp <- xs <- z1 <- z2 <- numeric(n + 1); xp[1] <- rnorm(1, 0, sqrt(lp * sig2)); xs[1] <- rnorm(1, 0, sqrt(sig2)); if (Vo > 0) { z1[1] <- rnorm(1, 0, sqrt(Vo)); z2[1] <- rnorm(1, 0, sqrt(Vo)) }
  mp <- ms <- numeric(n + 1)
  if (!is.null(steps)) { st <- steps
    prof <- function(amp, tau) { up <- pmin(pmax(time - st$t0, 0), st$t1 - st$t0); m <- amp * (1 - exp(-up / tau)); ifelse(time > st$t1, amp * (1 - exp(-(st$t1 - st$t0) / tau)) * exp(-(time - st$t1) / tau), m) }
    mp <- prof(st$amp_p, st$tau_p); ms <- prof(st$amp_s, st$tau_s) }
  if (!is.null(m_p_fun)) mp <- m_p_fun(time); if (!is.null(m_s_fun)) ms <- m_s_fun(time)
  rate <- numeric(n + 1); rate[1] <- nu0
  for (i in 1:n) {
    g <- if (het != 0) min(max((rate[i] / nu0)^(-het), 0.1), 10) else 1; sg <- sqrt(g)
    xp[i + 1] <- a_p * xp[i] + sd_p * sg * e[i, 1]; xs[i + 1] <- a_s * xs[i] + sd_s * sg * e[i, 2]
    if (Vo > 0) { z1[i + 1] <- ea * (cc * z1[i] - ss * z2[i]) + sd_o * sg * e[i, 3]; z2[i + 1] <- ea * (ss * z1[i] + cc * z2[i]) + sd_o * sg * e[i, 4] }
    rate[i + 1] <- max(nu0 - (xp[i + 1] + mp[i + 1]) + (xs[i + 1] + ms[i + 1]) - z1[i + 1], nu_min)
  }
  increment <- rate * h + sqrt(sig2 * lR * h) * eph; phase <- cumsum(increment)
  beats <- numeric(); index <- 1L; phase_last <- 0; look <- ceiling(6 / h)
  repeat { we <- min(index + look, n + 1L); if (index >= we) break; w <- (index + 1L):we; fc <- which(phase[w] - phase_last >= 1)[1]; if (is.na(fc)) break
    i <- w[fc]; fr <- min(max((1 - (phase[i - 1L] - phase_last)) / increment[i], 0), 1); beats <- c(beats, time[i - 1L] + fr * h); phase_last <- phase[i - 1L] + fr * increment[i]; index <- i }
  at <- function(x) approx(time, x, xout = beats, rule = 2)$y
  out <- data.table(t = beats, rr = diff(c(0, beats)), x_p = at(xp + mp), x_s = at(xs + ms), z1 = at(z1), m_p = at(mp), m_s = at(ms))[-1]
  if (artifacts > 0) { idx <- sample(nrow(out), round(artifacts * nrow(out))); out$rr[idx] <- out$rr[idx] * ifelse(runif(length(idx)) > 0.5, 1.7, 0.45); out$artifact <- seq_len(nrow(out)) %in% idx } else out$artifact <- FALSE
  out
}

# ---- Extended fit ------------------------------------------------------------------------------------------------------------------------------------
#' Causal local rate (mean of 1/RR over the previous `L` beats) relative to the recording mean: the scale variable of the rate-dependent noise.
rate_local <- function(dy, L = 8) { r <- 1 / dy; n <- length(dy); out <- numeric(n); out[1] <- r[1]; cs <- c(0, cumsum(r)); for (k in 2:n) { a <- max(1, k - L); out[k] <- (cs[k] - cs[a]) / (k - a) }; out / mean(r) }

#' MAP fit of M3 and its extensions. spec: osc (damped oscillator; ln alpha, ln lambda_o free, omega at the spectral peak), het (noise scale rate^-gamma, gamma free),
#' robust (Student-t weights instead of the heuristic artifact multiplier), lp (override of the anchor). theta = (ln ks, ln(kp-ks), ln lR, [ln alpha, ln lo], [gamma]).
fit_ext1 <- function(series, spec = list(), omega = 0, th_start = NULL) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; n <- length(dy)
  use_osc <- isTRUE(spec$osc); use_het <- isTRUE(spec$het); lp <- if (is.null(spec$lp)) pri$energy_ratio else spec$lp
  gm0 <- gate_mult(dy, ctl); rho <- rate_local(dy); nu_t <- if (is.null(spec$nu_t)) 4 else spec$nu_t
  decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); lR <- exp(th[3]); i <- 4L; osc <- NULL; qs <- NULL
    if (use_osc) { osc <- list(alpha = exp(th[i]), omega = omega, lo = exp(th[i + 1])); i <- i + 2L }
    if (use_het) { qs <- pmin(pmax(rho^(-th[i]), 0.1), 10); i <- i + 1L }
    list(ks = ks, kp = kp, lR = lR, osc = osc, qs = qs, gamma = if (use_het) th[length(th)] else 0) }
  run <- function(p, gm = gm0, full = FALSE, smooth = FALSE, ...) { if (!is.null(spec$adapt)) { r1 <- ext_run(dy, ctl, p$kp, p$ks, lp, p$lR, osc = p$osc, qs = p$qs, Xextra = if (!is.null(spec$inputs)) input_design(dy, p$kp, p$ks, spec$inputs$events) else NULL, gm = gm, full = TRUE)
      if (is.finite(r1$LL_star)) p$qs <- (if (is.null(p$qs)) 1 else p$qs) * adaptive_scale(r1$std_innov, gm, spec$adapt, ctl$diag_gate_max) }
    ext_run(dy, ctl, p$kp, p$ks, lp, p$lR, osc = p$osc, qs = p$qs, Xextra = if (!is.null(spec$inputs)) input_design(dy, p$kp, p$ks, spec$inputs$events) else NULL, gm = gm, full = full, smooth = smooth, ...) }
  robust_gm <- function(p) { r <- run(p, gm = rep(1, n), full = TRUE); if (!is.finite(r$LL_star)) return(gm0); pmax(1, (nu_t + r$std_innov^2) / (nu_t + 1)) }
  cost <- function(th) {
    p <- decode(th); if (use_osc && (p$osc$alpha > 5 || p$osc$alpha < 1e-3)) return(1e12)
    gm <- if (isTRUE(spec$robust)) robust_gm(p) else gm0; r <- run(p, gm = gm); if (!is.finite(r$LL_star)) return(1e12)
    pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), r, dy, pri)
    ex <- 0; if (use_osc) ex <- ex - 0.5 * ((log(p$osc$alpha) - log(0.2)) / 0.6)^2 - 0.5 * ((log(p$osc$lo) - log(0.3)) / 1.5)^2
    if (use_het) ex <- ex - 0.5 * (p$gamma / 1.5)^2
    v <- -(r$LL_star + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  extras0 <- c(if (use_osc) c(log(0.2), log(0.3)), if (use_het) 0)
  starts <- if (is.null(th_start)) lapply(.m3_starts(pri, ctl), function(st) c(unname(.par_to_theta(st)), extras0)) else list(th_start)
  runs <- lapply(starts, function(th0) tryCatch(optim(th0, cost, method = "BFGS", control = list(maxit = ctl$maxit)), error = function(e) NULL))
  vals <- vapply(runs, function(o) if (is.null(o) || !is.finite(o$value)) Inf else o$value, 0); if (!any(is.finite(vals))) return(list(ok = FALSE))
  k <- which(vals <= min(vals) + 0.01)[1]; th <- runs[[k]]$par; p <- decode(th); gm <- if (isTRUE(spec$robust)) robust_gm(p) else gm0
  list(ok = TRUE, theta = th, p = p, lp = lp, omega = omega, gm = gm, cost = vals[k], run = function(...) run(p, gm = gm, ...), decode = decode, dy = dy, ctl = ctl, n_par = length(th), spec = spec,
       tau_p = 1 / p$kp, tau_s = 1 / p$ks)
}

#' fit_ext1 with the oscillator frequency searched on a grid (0.15-0.40 Hz; profile of the posterior), two damping starts per frequency, from the fit without oscillator.
OSC_GRID_HZ <- c(0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
fit_ext <- function(series, spec = list()) {
  if (!isTRUE(spec$osc)) return(fit_ext1(series, spec))
  base <- fit_ext1(series, modifyList(spec, list(osc = FALSE))); if (!isTRUE(base$ok)) return(base)
  best <- NULL
  for (f in (if (is.null(spec$omega)) OSC_GRID_HZ else spec$omega / (2 * pi))) for (a0 in c(0.1, 0.4)) {
    th_b <- base$theta; th0 <- c(th_b[1:3], log(a0), log(0.3), if (isTRUE(spec$het)) th_b[4])
    fit <- fit_ext1(series, spec, omega = 2 * pi * f, th_start = th0); if (isTRUE(fit$ok) && (is.null(best) || fit$cost < best$cost)) best <- fit }
  if (is.null(best)) base else { best$base_cost <- base$cost; best }
}

# ---- Known inputs (events) as regressors with OU residuals -----------------------------------------------------------------------------------------
#' Response of a first-order drive (rate kappa) to a unit rectangular input on [t0, t1]: g(t); curves at the beat times and its integral over every beat interval.
input_curve <- function(t, kappa, t0, t1) { up <- pmin(pmax(t - t0, 0), t1 - t0); g <- (1 - exp(-kappa * up)) / kappa; ifelse(t > t1, g * exp(-kappa * (t - t1)), g) }
input_design <- function(dy, kp, ks, events, grid = 12) {                  # events: data.frame(t0, t1); two columns (vagal, sympathetic) per event
  t_end <- cumsum(dy); t_start <- c(0, t_end[-length(dy)])
  integ <- function(kappa, t0, t1) { vals <- sapply(seq(0, 1, length.out = grid), function(a) input_curve(t_start + a * dy, kappa, t0, t1)); rowSums((vals[, -1] + vals[, -grid]) / 2) * dy / (grid - 1) }
  do.call(cbind, lapply(seq_len(nrow(events)), function(e) { m <- cbind(-integ(kp, events$t0[e], events$t1[e]), integ(ks, events$t0[e], events$t1[e])); colnames(m) <- paste0(c("vagal", "symp"), e); m }))
}

# ---- Model comparison helpers ---------------------------------------------------------------------------------------------------------------------
#' Estimated states of a fitted model: vagal total (OU + oscillator + event response), sympathetic total (OU + event response), with their variances.
#' smooth = TRUE gives the RTS-smoothed estimates (M3+, M3O+).
state_estimates <- function(f, smooth = FALSE) {
  r <- f$run(full = TRUE, smooth = smooth); sfx <- if (smooth) "_smooth" else ""
  xp <- r[[paste0("X_p", sfx)]]; xs <- r[[paste0("X_s", sfx)]]; xpv <- r[[paste0("X_p_var", sfx)]]; xsv <- r[[paste0("X_s_var", sfx)]]
  if (!is.null(f$p$osc)) { xp <- xp + r[[paste0("X_o", sfx)]]; xpv <- r[[paste0("X_v_var", sfx)]] }
  if (!is.null(f$spec$inputs)) { ev <- f$spec$inputs$events; tt <- cumsum(f$dy); b <- r$beta
    for (e in seq_len(nrow(ev))) { xp <- xp + b[2 * e] * input_curve(tt, f$p$kp, ev$t0[e], ev$t1[e]); xs <- xs + b[2 * e + 1] * input_curve(tt, f$p$ks, ev$t0[e], ev$t1[e]) } }
  list(xp = xp, xpv = xpv, xs = xs, xsv = xsv, r = r)
}
#' Events estimated from the RR series with the pipeline's detector: data.frame(t0, t1) or NULL. exercise: one bout, end by the fitted profile; tilt: all candidate pairs.
events_from_rr <- function(dy, mode, seg_cfg) {
  d <- tryCatch(detect_events(dy, mode, seg_cfg), error = function(e) NULL); if (is.null(d) || !isTRUE(d$detected)) return(NULL)
  if (mode == "exercise") { t1 <- if (is.finite(d$t_end_profile)) d$t_end_profile else if (is.finite(d$t_end_cp)) d$t_end_cp else d$t_recovery   # onset and end of the fitted profile, as the pipeline does
    t0 <- if (is.finite(d$t_onset_profile)) d$t_onset_profile else d$t_onset; ev <- data.frame(t0 = t0, t1 = t1) }
  else ev <- data.frame(t0 = d$events$t_onset, t1 = d$events$t_recovery)
  ev <- ev[is.finite(ev$t0) & is.finite(ev$t1) & ev$t1 > ev$t0, , drop = FALSE]; if (nrow(ev)) ev else NULL
}

# ---- Adaptive noise scale -----------------------------------------------------------------------------------------------------------------------------
#' Local scale of the process noise from the PAST standardized innovations: exponentially weighted mean of z^2 with half-life `hl` beats (clipped to [0.2, 5]);
#' beats flagged as artifacts count as 1. A two-stage approximation of stochastic volatility of the drives (amplitude modulation); the states keep their meaning.
adaptive_scale <- function(z, gm, hl = 25, diag_gate = 1.5) {
  u <- z^2; u[!is.finite(u) | gm >= diag_gate] <- 1; a <- 0.5^(1 / hl); n <- length(u); w <- numeric(n); w[1] <- 1
  for (k in 2:n) w[k] <- (1 - a) * u[k - 1] + a * w[k - 1]
  pmin(pmax(w, 0.2), 5)
}

# ---- Stochastic set-points (ext2 filter) ---------------------------------------------------------------------------------------------------------------
suppressMessages(Rcpp::sourceCpp("experiments/m3plus/ext2_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))
ext2_run <- function(dy, ctl, kp, ks, lp, lR, osc = NULL, qs = NULL, q = c(0, 0), gm = NULL, full = FALSE, smooth = FALSE, beta_fixed = numeric(0), sig2_fixed = -1) {
  n <- length(dy); if (is.null(gm)) gm <- gate_mult(dy, ctl); if (is.null(qs)) qs <- rep(1, n)
  ext2_filter(dy, gm, qs, cbind(dy), kp, ks, lp, lR, if (is.null(osc)) 0 else osc$alpha, if (is.null(osc)) 0 else osc$omega, if (is.null(osc)) 0 else osc$lo, q[1], q[2], ctl$lik_gate_max, full, smooth, beta_fixed, sig2_fixed)
}

#' MAP fit with stochastic set-points. spec: sp = subset of c("p", "s") (which drives get a moving target; ln q free, prior N(ln 0.01, 2.5^2)), osc (damped oscillator, frequency grid),
#' adapt (half-life in beats of the two-stage adaptive noise scale). theta = (ln ks, ln(kp-ks), ln lR, [ln alpha, ln lo], [ln q_p], [ln q_s]).
fit_ext2_1 <- function(series, spec = list(), omega = 0, th_start = NULL) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; n <- length(dy); use_osc <- isTRUE(spec$osc); sp <- spec$sp; lp <- pri$energy_ratio; gm0 <- gate_mult(dy, ctl)
  decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); lR <- exp(th[3]); i <- 4L; osc <- NULL; q <- c(0, 0)
    if (use_osc) { osc <- list(alpha = exp(th[i]), omega = omega, lo = exp(th[i + 1])); i <- i + 2L }
    if ("p" %in% sp) { q[1] <- exp(th[i]); i <- i + 1L }; if ("s" %in% sp) { q[2] <- exp(th[i]); i <- i + 1L }
    list(ks = ks, kp = kp, lR = lR, osc = osc, q = q) }
  run <- function(p, gm = gm0, full = FALSE, smooth = FALSE, ...) { qs <- NULL
    if (!is.null(spec$adapt)) { r1 <- ext2_run(dy, ctl, p$kp, p$ks, lp, p$lR, osc = p$osc, q = p$q, gm = gm, full = TRUE); if (is.finite(r1$LL_star)) qs <- adaptive_scale(r1$std_innov, gm, spec$adapt, ctl$diag_gate_max) }
    ext2_run(dy, ctl, p$kp, p$ks, lp, p$lR, osc = p$osc, qs = qs, q = p$q, gm = gm, full = full, smooth = smooth, ...) }
  cost <- function(th) { p <- decode(th); if (use_osc && (p$osc$alpha > 5 || p$osc$alpha < 1e-3)) return(1e12); if (any(p$q > 1e3)) return(1e12)
    r <- run(p); if (!is.finite(r$LL_star)) return(1e12); pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), r, dy, pri); ex <- 0
    if (use_osc) ex <- ex - 0.5 * ((log(p$osc$alpha) - log(0.2)) / 0.6)^2 - 0.5 * ((log(p$osc$lo) - log(0.3)) / 1.5)^2
    for (qq in p$q[p$q > 0]) ex <- ex - 0.5 * ((log(qq) - log(0.01)) / 2.5)^2
    v <- -(r$LL_star + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  extras0 <- c(if (use_osc) c(log(0.2), log(0.3)), rep(log(0.02), length(sp)))
  starts <- if (is.null(th_start)) lapply(.m3_starts(pri, ctl), function(st) c(unname(.par_to_theta(st)), extras0)) else list(th_start)
  runs <- lapply(starts, function(th0) tryCatch(optim(th0, cost, method = "BFGS", control = list(maxit = ctl$maxit)), error = function(e) NULL))
  vals <- vapply(runs, function(o) if (is.null(o) || !is.finite(o$value)) Inf else o$value, 0); if (!any(is.finite(vals))) return(list(ok = FALSE))
  k <- which(vals <= min(vals) + 0.01)[1]; th <- runs[[k]]$par; p <- decode(th)
  list(ok = TRUE, theta = th, p = p, lp = lp, omega = omega, cost = vals[k], run = function(...) run(p, ...), dy = dy, ctl = ctl, n_par = length(th), spec = spec, tau_p = 1 / p$kp, tau_s = 1 / p$ks)
}
fit_ext2 <- function(series, spec = list()) {
  if (!isTRUE(spec$osc)) return(fit_ext2_1(series, spec))
  base <- fit_ext2_1(series, modifyList(spec, list(osc = FALSE))); if (!isTRUE(base$ok)) return(base); best <- NULL
  for (f in OSC_GRID_HZ) for (a0 in c(0.1, 0.4)) { th_b <- base$theta; th0 <- c(th_b[1:3], log(a0), log(0.3), th_b[-(1:3)])
    fit <- fit_ext2_1(series, spec, omega = 2 * pi * f, th_start = th0); if (isTRUE(fit$ok) && (is.null(best) || fit$cost < best$cost)) best <- fit }
  if (is.null(best)) base else best
}
