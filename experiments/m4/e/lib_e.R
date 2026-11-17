# ==============================================================================
# experiments/m4/e/lib_e.R -- EXPERIMENTAL. Components of M4e (README.md) for one series, so that the gates and the interval factors can be set afterwards without refitting.
#   fit_components()  M3 (production), B (multi-start, no oscillator), O (B + oscillator), RB / RO (two-regime HMM without / with the O oscillator; fitted when the change-point gain >= 3)
#   pick()            the component that a variant (q_osc, G) selects, from the stored statistics
#   states_c()        states of a component, posterior SDs multiplied by cal (regime components only)
# ==============================================================================
source("experiments/m4/lib.R")
REG_FLOOR <- 3; Q_COMPONENTS <- TRUE; J_COMPONENTS <- TRUE; O_COMPONENTS <- TRUE

fit_components <- function(series, spec_extra = list()) {
  out <- list(); lp0 <- series$priors$energy_ratio
  out$M3 <- list(kind = "ext", fit = fit_ext(series), lp = lp0)
  base <- fit_base_ms(series); if (!isTRUE(base$ok)) return(list(models = out, stats = c(lr_osc = NA, gain = NA)))
  out$B <- list(kind = "ext", fit = base, lp = lp0)
  a <- fit_osc_from(series, base); if (isTRUE(a$ok)) out$O <- list(kind = "ext", fit = a, lp = lp0)
  lr <- if (isTRUE(a$ok)) unname(2 * (base$cost - a$cost)) else NA_real_
  cp <- tryCatch(cp_search(series, base, min_gain = -Inf), error = function(e) NULL); gain <- if (is.null(cp)) NA_real_ else cp$gain
  if (is.finite(gain) && gain >= REG_FLOOR) for (nm in c("RB", "RO")) {
    osc <- if (nm == "RO") (if (isTRUE(a$ok)) a$p$osc else NULL) else NULL; if (nm == "RO" && is.null(osc)) next
    args <- modifyList(c(REG_ARGS, list(osc = osc)), spec_extra)
    h <- tryCatch(do.call(fit_regime, c(list(series, 2, cp, base), args)), error = function(e) list(ok = FALSE))
    if (isTRUE(h$ok)) out[[nm]] <- list(kind = "reg", fit = h, lp = lp0, gate_open = TRUE) }
  if (isTRUE(Q_COMPONENTS) && is.finite(gain) && gain >= REG_FLOOR) for (nm in c("QB", "QO")) {
    osc <- if (nm == "QO") (if (isTRUE(a$ok)) a$p$osc else NULL) else NULL; if (nm == "QO" && is.null(osc)) next
    h <- tryCatch(fit_regime4(series, cp, base, osc = osc), error = function(e) list(ok = FALSE))
    if (isTRUE(h$ok)) out[[nm]] <- list(kind = "reg", fit = h, lp = lp0, gate_open = TRUE) }
  if (isTRUE(J_COMPONENTS) && is.finite(gain) && gain >= REG_FLOOR && isTRUE(a$ok)) {
    h <- tryCatch(fit_regime4(series, cp, base, rates_free = FALSE, omega = a$omega, osc_start = a$p$osc), error = function(e) list(ok = FALSE))
    if (isTRUE(h$ok)) out$JO <- list(kind = "reg", fit = h, lp = lp0, gate_open = TRUE) }
  if (isTRUE(O_COMPONENTS) && is.finite(gain) && gain >= REG_FLOOR) {
    h <- tryCatch(fit_regime4(series, cp, base, rates_free = TRUE, order_on = TRUE), error = function(e) list(ok = FALSE)); if (isTRUE(h$ok)) out$QBo <- list(kind = "reg", fit = h, lp = lp0, gate_open = TRUE)
    if (isTRUE(a$ok)) { h <- tryCatch(fit_regime4(series, cp, base, rates_free = TRUE, order_on = TRUE, omega = a$omega, osc_start = a$p$osc), error = function(e) list(ok = FALSE))
      if (isTRUE(h$ok)) out$QJo <- list(kind = "reg", fit = h, lp = lp0, gate_open = TRUE) } }
  list(models = out, stats = c(lr_osc = lr, gain = gain))
}
#' Component chosen by the variant: oscillator if lr_osc >= q_osc, regimes if gain >= G (falls back to the linear component if the regime fit is missing).

COV_GRID <- seq(1, 2.5, by = 0.1)

# ---- E2: regimes with their own response kinetics (reg4_filter.cpp) ---------------------------------------------------------------------------------------------------------
suppressMessages(Rcpp::sourceCpp("experiments/m4/e/reg4_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))
reg4_run <- function(dy, ctl, p, osc, gm, full = FALSE, smooth = FALSE)
  reg4_filter(dy, gm, p$kp, p$ks, p$lp, p$lR, p$sig2, p$nu0, p$mu, p$rate, p$Pi, c(1, rep(0, nrow(p$mu) - 1)), if (is.null(osc)) 0 else osc$alpha, if (is.null(osc)) 0 else osc$omega, if (is.null(osc)) 0 else osc$lo,
              ctl$lik_gate_max, full, smooth)
RATE_SD <- 1.5; RATE_BOUNDS <- log(c(1 / 600, 2))
#' Two-regime HMM (reference regime 0, active regime 1). theta = (ln ks, ln(kp - ks), ln lR, ln sigma^2, nu0, mu_p, mu_s, stay0, stay1, [ln r_off,p, ln r_off,s, ln r_on,p, ln r_on,s if rates_free],
#' [ln alpha, ln lo if omega is given: the oscillator estimated jointly at angular frequency omega]). rates_free = FALSE: response rates = kappa (the model of fit_regime).
#' Priors as fit_regime (REG_ARGS), ln r ~ N(ln kappa_b, RATE_SD^2), oscillator priors as fit_ext1. osc: a fixed oscillator (used only if omega is NULL).
fit_regime4 <- function(series, cp, base, osc = NULL, rates_free = TRUE, omega = NULL, osc_start = NULL, order_on = FALSE, stay_prior = 7, stay_sd = 0.7, recip = 20, mu_sd = 0.8, mu_max = 6) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; lp <- pri$energy_ratio; gm <- gate_mult(dy, ctl)
  nr <- if (rates_free) 4L else 0L; io <- 9L + nr; jo <- !is.null(omega)
  decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); st <- plogis(th[8:9]); Pi <- rbind(c(st[1], 1 - st[1]), c(1 - st[2], st[2]))
    rate <- if (rates_free) rbind(exp(th[10:11]), exp(c(th[12], if (order_on) th[12] - exp(th[13]) else th[13]))) else rbind(c(kp, ks), c(kp, ks))
    o <- if (jo) list(alpha = exp(th[io + 1]), omega = omega, lo = exp(th[io + 2])) else osc
    list(ks = ks, kp = kp, lR = exp(th[3]), sig2 = exp(th[4]), nu0 = th[5], mu = rbind(c(0, 0), th[6:7]), Pi = Pi, rate = rate, lp = lp, osc = o, th = th) }
  lr_of <- function(th) c(th[10:12], if (order_on) th[12] - exp(th[13]) else th[13])
  cost <- function(th) { if (rates_free && any(lr_of(th) < RATE_BOUNDS[1] | lr_of(th) > RATE_BOUNDS[2])) return(1e12); p <- decode(th)
    if (p$nu0 <= 0.2 || any(abs(p$mu) > mu_max)) return(1e12); if (jo && (p$osc$alpha > 5 || p$osc$alpha < 1e-3 || p$osc$lo > 50)) return(1e12)
    r <- reg4_run(dy, ctl, p, p$osc, gm); if (!is.finite(r$LL)) return(1e12)
    pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), list(sig2 = p$sig2, nu0 = p$nu0), dy, pri)
    ex <- -recip * max(p$mu[2, 1] * p$mu[2, 2], 0) - 0.5 * sum((p$mu / mu_sd)^2) - 0.5 * sum(((th[8:9] - stay_prior) / stay_sd)^2)
    if (rates_free) ex <- ex - 0.5 * sum(((lr_of(th) - log(c(p$kp, p$ks, p$kp, p$ks))) / RATE_SD)^2)
    if (jo) ex <- ex - 0.5 * ((log(p$osc$alpha) - log(0.2)) / 0.6)^2 - 0.5 * ((log(p$osc$lo) - log(0.3)) / 1.5)^2
    v <- -(r$LL + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  m1 <- c(cp$mu_p, cp$mu_s); if (m1[1] * m1[2] > 0) m1 <- c(-0.1, 0.1) * sign(m1[2] - m1[1] + 1e-9); m1 <- pmax(pmin(m1, 3), -3)
  th0 <- c(log(base$p$ks), log(base$p$kp - base$p$ks), log(base$p$lR), log(cp$sig2), cp$nu0, m1, stay_prior, stay_prior)
  rs <- if (rates_free) list(log(c(base$p$kp, base$p$ks, base$p$kp, base$p$ks)), log(c(1 / 30, 1 / 60, 1 / 8, 1 / 20))) else list(numeric(0))
  if (rates_free && order_on) rs <- lapply(rs, function(r) c(r[1:3], log(max(r[3] - r[4], 0.05))))
  os <- if (jo) unique(list(if (!is.null(osc_start)) log(c(osc_start$alpha, max(osc_start$lo, 1e-3))) else log(c(0.2, 0.3)), log(c(0.2, 0.3)))) else list(numeric(0))
  starts <- unlist(lapply(rs, function(r) lapply(os, function(o) c(th0, r, o))), recursive = FALSE)
  ps <- c(1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, rep(1, nr), if (jo) c(1, 1))
  fits <- lapply(starts, function(t0) tryCatch(optim(t0, cost, method = "BFGS", control = list(maxit = 500, parscale = ps)), error = function(e) NULL))
  fits <- Filter(function(f) !is.null(f) && is.finite(f$value) && f$value < 1e11, fits); if (!length(fits)) return(list(ok = FALSE))
  o <- fits[[which.min(vapply(fits, function(f) f$value, 0))]]; p <- decode(o$par)
  list(ok = TRUE, p = p, K = 2, cost = o$value, n_par = length(o$par), dy = dy, ctl = ctl, lp = lp, gm = gm, tau_p = 1 / p$kp, tau_s = 1 / p$ks, osc = p$osc,
       tau_resp = c(on_p = 1 / p$rate[2, 1], off_p = 1 / p$rate[1, 1], on_s = 1 / p$rate[2, 2], off_s = 1 / p$rate[1, 2]),
       run = function(dyy = dy, gmm = gm, full = TRUE, smooth = FALSE) reg4_run(dyy, ctl, p, p$osc, gmm, full, smooth))
}
