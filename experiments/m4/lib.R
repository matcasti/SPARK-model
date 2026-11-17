# ==============================================================================
# experiments/m4/lib.R -- EXPERIMENTAL. The M4 ladder of README.md on top of experiments/m3plus (extended filter), experiments/regimes/v2 (regimes) and
# experiments/features (feature anchor). Source after analysis/00_setup.R.
#   fit_family()   M3, M4a, M4b, M4c (+ M4d-H with eb, + M4d-I with events) for one series; each model = list(kind = "ext" | "reg", fit, lp, ...)
#   weak_taus()    weak-prior fit (for the empirical-Bayes priors of M4d-H)      eb_from()  hyperparameters from the other recordings
#   m_states(), m_diag(), m_heldout(), m_sd()   accessors common to both kinds
# ==============================================================================
source("experiments/regimes/v2/lib2.R"); source("experiments/features/lib.R")
.feat <- readRDS("experiments/features/_out/sim_features.rds")
ANCHOR_MODEL <- with_seed(1L, train_full(.feat, log(.feat$lambda_p), unlist(SETS[c("Sf", "E", "Gf", "F")])))   # target: the OU part (the oscillator is modelled in M4)
feature_anchor <- function(series) { fe <- add_scale_free(as.data.table(t(all_features(series$dy, series$priors)))); min(max(exp(ANCHOR_MODEL(fe)), 0.05), 50) }
ANCHOR_SD <- 0.5; GATE <- 9; REG_ARGS <- list(stay_prior = 7, stay_sd = 0.7, recip = 20)
MS_STARTS <- CJ(tp = c(2, 5, 10, 20), ts = c(60, 200, 600))

better <- function(a, b) isTRUE(b$ok) && (!isTRUE(a$ok) || b$cost < a$cost - 1e-8)
fit_base_ms <- function(series, spec = list()) {
  best <- fit_ext1(series, spec)
  for (j in seq_len(nrow(MS_STARTS))) { th <- c(log(1 / MS_STARTS$ts[j]), log(1 / MS_STARTS$tp[j] - 1 / MS_STARTS$ts[j]), log(0.02))
    f <- tryCatch(fit_ext1(series, spec, th_start = th), error = function(e) list(ok = FALSE)); if (better(best, f)) best <- f }
  best
}
fit_osc_from <- function(series, base, spec = list(), freqs = OSC_GRID_HZ, damp = c(0.1, 0.4)) {
  best <- list(ok = FALSE)
  for (f in freqs) for (a0 in damp) { fit <- tryCatch(fit_ext1(series, modifyList(spec, list(osc = TRUE)), omega = 2 * pi * f, th_start = c(base$theta[1:3], log(a0), log(0.3))), error = function(e) list(ok = FALSE))
    if (better(best, fit)) best <- fit }
  best
}
#' Refit of an oscillator model from a given theta with a changed spec/series (two damping starts and the given time constants).
refit_osc <- function(series, from, spec, tps = NULL) {
  best <- list(ok = FALSE); ths <- list(from$theta, replace(from$theta, 4, log(0.4)))
  for (tp in tps) { ks <- exp(from$theta[1]); if (1 / tp > ks) ths[[length(ths) + 1]] <- replace(from$theta, 2, log(1 / tp - ks)) }
  for (th in ths) { f <- tryCatch(fit_ext1(series, modifyList(spec, list(osc = TRUE)), omega = from$omega, th_start = th), error = function(e) list(ok = FALSE)); if (better(best, f)) best <- f }
  best
}
with_priors <- function(series, eb) { s <- series; if (!is.null(eb)) { s$priors$log_tau_p_mu <- eb$mu_p; s$priors$sigma_log_tau_p <- eb$sd_p; s$priors$log_tau_s_mu <- eb$mu_s; s$priors$sigma_log_tau_s <- eb$sd_s }; s }
with_anchor <- function(series, lp) { s <- series; s$priors$energy_ratio <- lp; s }

#' Weak-prior fit (oscillator, anchor lp): the time constants that the data prefer, for the empirical-Bayes priors.
weak_taus <- function(series, from, lp) {
  w <- with_priors(series, list(mu_p = mean(log(c(0.5, 60))), sd_p = diff(log(c(0.5, 60))) / (2 * 1.96), mu_s = mean(log(c(10, 2000))), sd_s = diff(log(c(10, 2000))) / (2 * 1.96)))
  f <- refit_osc(w, from, list(lp = lp), tps = c(3, 10)); if (isTRUE(f$ok)) c(tau_p = f$tau_p, tau_s = f$tau_s) else c(tau_p = NA_real_, tau_s = NA_real_)
}
eb_from <- function(taus) { lp <- log(taus[, 1]); ls <- log(taus[, 2]); k <- is.finite(lp) & is.finite(ls)
  list(mu_p = mean(lp[k]), sd_p = max(sd(lp[k]), 0.15), mu_s = mean(ls[k]), sd_s = max(sd(ls[k]), 0.15)) }

#' The ladder for one series. eb: hyperparameters for M4d-H (NULL: skip); events: data.frame(t0, t1) on the series' time axis for M4d-I (NULL: skip).
#' stage1 = TRUE also returns the weak-prior time constants (input of eb_from for the other recordings).
fit_family <- function(series, eb = NULL, events = NULL, stage1 = TRUE) {
  out <- list(); lp0 <- series$priors$energy_ratio
  m3 <- fit_ext(series); out$M3 <- list(kind = "ext", fit = m3, lp = lp0)
  base <- fit_base_ms(series); a <- fit_osc_from(series, base); if (!isTRUE(a$ok)) return(out)
  out$M4a <- list(kind = "ext", fit = a, lp = lp0)
  m <- feature_anchor(series); grid <- m * exp(c(-1, -0.5, 0, 0.5, 1))
  cand <- lapply(grid, function(l) refit_osc(series, a, list(lp = l)))
  tot <- vapply(seq_along(grid), function(j) if (isTRUE(cand[[j]]$ok)) cand[[j]]$cost + 0.5 * ((log(grid[j]) - log(m)) / ANCHOR_SD)^2 else Inf, 0)
  if (!any(is.finite(tot))) return(out); j <- which.min(tot); b <- cand[[j]]; lb <- grid[j]
  out$M4b <- list(kind = "ext", fit = b, lp = lb, lam_feature = m)
  regimes <- function(s, lpr = lb, osc = b$p$osc) { sb <- with_anchor(s, lpr); bb <- tryCatch(fit_ext1(sb, list(), th_start = base$theta), error = function(e) list(ok = FALSE)); if (!isTRUE(bb$ok)) return(NULL)
    cp <- cp_search(sb, bb, min_gain = -Inf); gain <- if (is.null(cp)) NA_real_ else cp$gain
    if (is.null(cp) || cp$gain < GATE) return(list(open = FALSE, gain = gain))
    h <- tryCatch(do.call(fit_regime, c(list(sb, 2, cp, bb), REG_ARGS, list(osc = osc))), error = function(e) list(ok = FALSE))
    if (!isTRUE(h$ok)) return(list(open = FALSE, gain = gain)); list(open = TRUE, gain = gain, h = h) }
  rg <- regimes(series)
  out$M4c <- if (isTRUE(rg$open)) list(kind = "reg", fit = rg$h, lp = lb, gate_open = TRUE, gain = rg$gain) else modifyList(out$M4b, list(gate_open = FALSE, gain = if (is.null(rg)) NA_real_ else rg$gain))
  rw <- regimes(series, lp0, a$p$osc)   # post hoc M4c-w: the same regimes with the production anchor and the M4a oscillator (isolates the anchor)
  out$M4cw <- if (isTRUE(rw$open)) list(kind = "reg", fit = rw$h, lp = lp0, gate_open = TRUE, gain = rw$gain) else modifyList(out$M4a, list(gate_open = FALSE, gain = if (is.null(rw)) NA_real_ else rw$gain))
  if (stage1) out$weak <- weak_taus(series, b, lb)
  if (!is.null(eb)) { se <- with_priors(series, eb); fe <- refit_osc(se, b, list(lp = lb), tps = exp(eb$mu_p)); rh <- if (isTRUE(out$M4c$gate_open)) regimes(se) else NULL
    out$M4dH <- if (isTRUE(rh$open)) list(kind = "reg", fit = rh$h, lp = lb, gate_open = TRUE) else if (isTRUE(fe$ok)) list(kind = "ext", fit = fe, lp = lb, gate_open = FALSE) else NULL }
  if (!is.null(events) && nrow(events)) { fi <- refit_osc(series, b, list(lp = lb, inputs = list(events = events)), tps = c(4, 8)); if (isTRUE(fi$ok)) out$M4dI <- list(kind = "ext", fit = fi, lp = lb) }
  out
}

# ---- accessors ------------------------------------------------------------------------------------------------------------------------
m_states <- function(m, smooth = FALSE) {
  if (m$kind == "ext") { e <- state_estimates(m$fit, smooth); return(list(xp = e$xp, xpv = e$xpv, xs = e$xs, xsv = e$xsv)) }
  r <- m$fit$run(smooth = TRUE); o <- function(x) if (is.null(x)) 0 else x
  if (smooth) list(xp = r$Xs_p + o(r$Xs_o), xpv = if (!is.null(r$Xs_v_var)) r$Xs_v_var else r$Xs_p_var, xs = r$Xs_s, xsv = r$Xs_s_var)
  else list(xp = r$X_p + o(r$X_o), xpv = if (!is.null(r$X_v_var)) r$X_v_var else r$X_p_var, xs = r$X_s, xsv = r$X_s_var)
}
z_diag <- function(z, gm, ctl) { k <- which(is.finite(z) & gm < ctl$diag_gate_max); z <- z[k]; n <- length(z)
  c(ks = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(.975) / sqrt(n))) }
m_diag <- function(m, gm, ctl) z_diag(if (m$kind == "ext") m$fit$run(full = TRUE)$innovations_diagnostic else m$fit$run()$z, gm, ctl)
m_tau <- function(m) c(tau_p = m$fit$tau_p, tau_s = m$fit$tau_s)
#' Stationary SDs (Hz): vagal-like total (OU + oscillator) and sympathetic-like.
m_sd <- function(m) { if (m$kind == "ext") { s2 <- m$fit$run()$sig2; lo <- if (is.null(m$fit$p$osc)) 0 else m$fit$p$osc$lo } else { s2 <- m$fit$p$sig2; lo <- if (is.null(m$fit$osc)) 0 else m$fit$osc$lo }
  c(SD_p = sqrt(s2 * (m$lp + lo)), SD_s = sqrt(s2)) }
#' Held-out score: mean one-step log predictive density over beats `idx` of the full series dy, parameters (and GLS amplitudes) from the head fit.
m_heldout <- function(m, dy, gm, ctl, idx) {
  if (m$kind == "reg") return(mean(m$fit$run(dy, gm)$logdens[idx]))
  f <- m$fit; rr <- f$run(); X <- if (!is.null(f$spec$inputs)) input_design(dy, f$p$kp, f$p$ks, f$spec$inputs$events) else NULL
  r <- ext_run(dy, ctl, f$p$kp, f$p$ks, f$lp, f$p$lR, osc = f$p$osc, Xextra = X, gm = gm, full = TRUE, beta_fixed = rr$beta, sig2_fixed = rr$sig2)
  if (!is.finite(r$LL_star)) NA_real_ else mean(r$logdens[idx])
}
MODELS <- c("M3", "M4a", "M4b", "M4c", "M4cw", "M4dH", "M4dI")
