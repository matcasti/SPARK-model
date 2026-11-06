# ==============================================================================
# experiments/anchors/lib.R -- EXPERIMENTAL. Alternative spectral anchors for lambda_P = V_p / V_s. The paper's anchor splits the wavelet energy of the
# band by a logistic frequency rule; here lambda_P may depend on the time constants (kp, ks) through the spectrum the model itself implies:
# an OU drive with rate k and stationary variance V has PSD(f) = V * 2k / (k^2 + (2 pi f)^2), so its share of variance inside a band is
# (2/pi) * (atan(2 pi f_hi / k) - atan(2 pi f_lo / k)).
#   hard    the pipeline's anchor (a constant per recording)
#   const   the same value for every recording (1.3): is the data-driven part of the anchor informative?
#   bandfix the pipeline's anchor corrected for the share of each OU variance that lies inside the band: lambda0 * F_s(ks) / F_p(kp)
#   lfhf    LF/HF band-power ratio matched to the model: (V_p F_p,LF + V_s F_s,LF) / (V_p F_p,HF + V_s F_s,HF) = P_LF / P_HF
#   spec    non-negative least squares of the whole global wavelet spectrum on the two OU spectra (relative error)
#   oracle  the true lambda_P (simulations only)
# ==============================================================================
band_fraction <- function(k, lo, hi) (2 / pi) * (atan(2 * pi * hi / k) - atan(2 * pi * lo / k))
clamp_lambda <- function(x) if (!is.finite(x)) 1 else min(max(x, 1e-3), 1e3)

#' Global wavelet spectrum in the anchor band, built exactly as extract_spectral_priors() does (ascending frequencies).
global_spectrum <- function(dy, priors) {
  time_cum <- cumsum(dy); fs <- priors$cwt_fs; t_grid <- seq(min(time_cum), max(time_cum), by = 1 / fs)
  hr <- splinefun(time_cum, 1 / dy, method = "monoH.FC")(t_grid)
  invisible(capture.output(w <- biwavelet::wt(cbind(t_grid, hr - mean(hr)), do.sig = FALSE)))
  pm <- w$power; for (i in seq_len(ncol(pm))) pm[w$period > w$coi[i], i] <- NA
  p <- rowMeans(pm, na.rm = TRUE); f <- 1 / w$period; o <- order(f); f <- f[o]; p <- p[o]
  keep <- f >= priors$lf_lower_hz & f <= priors$hf_upper_hz & is.finite(p) & p > 0
  list(f = f[keep], p = p[keep])
}

#' A rule kp, ks -> lambda_P for the given series context.
make_rule <- function(name, series, lp_true = NA_real_) {
  pri <- series$priors; lo <- pri$lf_lower_hz; mid <- pri$lf_upper_hz; hi <- pri$hf_upper_hz; lam0 <- pri$energy_ratio
  switch(name,
    hard = function(kp, ks) lam0,
    const = function(kp, ks) 1.3,                       # the typical value the pipeline's anchor takes on the simulated series
    oracle = function(kp, ks) lp_true,
    bandfix = function(kp, ks) clamp_lambda(lam0 * band_fraction(ks, lo, hi) / band_fraction(kp, lo, hi)),
    lfhf = { R <- pri$p_lf / pri$p_hf
      function(kp, ks) { fp_lf <- band_fraction(kp, lo, mid); fp_hf <- band_fraction(kp, mid, hi); fs_lf <- band_fraction(ks, lo, mid); fs_hf <- band_fraction(ks, mid, hi)
        den <- R * fp_hf - fp_lf; if (!is.finite(R) || den <= 0) return(1e-3); clamp_lambda((fs_lf - R * fs_hf) / den) } },
    spec = { g <- global_spectrum(series$dy, pri)
      function(kp, ks) { ap <- 2 * kp / (kp^2 + (2 * pi * g$f)^2) / g$p; as <- 2 * ks / (ks^2 + (2 * pi * g$f)^2) / g$p   # columns divided by the data: relative error
        A <- matrix(c(sum(ap * ap), sum(ap * as), sum(ap * as), sum(as * as)), 2); b <- c(sum(ap), sum(as))
        cf <- tryCatch(solve(A, b), error = function(e) c(NA, NA))
        if (!all(is.finite(cf))) return(1)
        if (cf[1] <= 0) cf <- c(0, sum(as) / sum(as * as)) else if (cf[2] <= 0) cf <- c(sum(ap) / sum(ap * ap), 0)           # non-negativity
        clamp_lambda(cf[1] / max(cf[2], 1e-12)) } },
    stop("unknown rule: ", name))
}

#' M3 objective with lambda_P = rule(kp, ks), minimized by the pipeline's three-start scheme. Returns the fit with the same fields as the pipeline's fits.
fit_rule <- function(series, rule) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy
  objective <- function(ks, kp, lR) {
    fl <- filter_eval(dy, kp, ks, rule(kp, ks), lR, ctl); if (!is.finite(fl$LL_star)) return(1e12)
    pen <- .log_prior_terms(c(ks = ks, kp = kp, lR = lR), fl, dy, pri)
    v <- -(fl$LL_star + pen$tau_s + pen$tau_p + pen$topology + pen$lR); if (is.finite(v)) v else 1e12 }
  best <- .optimize_m3(objective, pri, ctl); if (is.null(best)) return(list(ok = FALSE, spec = "M3"))
  par <- list(ks = best$par[["ks"]], kp = best$par[["kp"]], lR = best$par[["lR"]]); par$lp <- rule(par$kp, par$ks)
  f <- filter_eval(dy, par$kp, par$ks, par$lp, par$lR, ctl, full = TRUE)
  list(ok = TRUE, spec = "M3", conv = 0, tau_p = 1 / par$kp, tau_s = 1 / par$ks, par = par, X_p = f$X_p, X_s = f$X_s, X_p_var = f$X_p_var, X_s_var = f$X_s_var,
       share_p = var(f$X_p) / (var(f$X_p) + var(f$X_s)), rmse_ms = 1000 * sqrt(mean((dy - f$implied_rr)^2)), LL_star = f$LL_star, cost = best$value,
       diag = { keep <- which(is.finite(f$innovations_diagnostic) & f$gate_multiplier < ctl$diag_gate_max); z <- f$innovations_diagnostic[keep]; n <- length(z)
         list(ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf_viol = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(0.975) / sqrt(n))) })
}
