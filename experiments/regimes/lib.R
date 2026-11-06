# ==============================================================================
# experiments/regimes/lib.R -- EXPERIMENTAL. Latent Markov regimes for the PNS/SNS targets. Loads the extended filter of experiments/m3plus (for M3, M3E and the simulator).
#   reg_run()     one pass of the regime filter       cp_search()  change-point regimes estimated with the SPARK likelihood (a rectangular regime [t0, t1])
#   fit_regime()  MAP fit of the HMM (K regimes)      sim_regime() simulator with regime-like or exercise-like targets
# ==============================================================================
source("experiments/m3plus/lib.R")
suppressMessages(Rcpp::sourceCpp("experiments/regimes/reg_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))
reg_run <- function(dy, ctl, kp, ks, lp, lR, sig2, nu0, mu, Pi, w0 = c(1, rep(0, nrow(mu) - 1)), gm = NULL, full = FALSE, osc = NULL, smooth = FALSE, eps = 0, kappa = 25) {
  if (is.null(gm)) gm <- gate_mult(dy, ctl); if (is.null(osc) && !smooth && eps == 0) return(reg_filter(dy, gm, kp, ks, lp, lR, sig2, nu0, mu, Pi, w0, ctl$lik_gate_max, full))
  reg2_filter(dy, gm, kp, ks, lp, lR, sig2, nu0, mu, Pi, w0, if (is.null(osc)) 0 else osc$alpha, if (is.null(osc)) 0 else osc$omega, if (is.null(osc)) 0 else osc$lo, ctl$lik_gate_max, full, smooth, eps, kappa) }

# ---- Change-point regime (one rectangular regime [t0, t1]) estimated with the SPARK likelihood, kappas fixed at their M3 values --------------------------------------
#' Grid search over onset and end (profile of the log-likelihood with the amplitudes profiled by GLS), then Nelder-Mead on (t0, t1). Returns NULL if the gain over M3 is below `min_gain`.
cp_search <- function(series, base, min_gain = 10, n_t0 = 28) {
  ctl <- series$control; dy <- series$dy; pri <- series$priors; p <- base$p; tt <- cumsum(dy); T_end <- tt[length(tt)]; gm <- gate_mult(dy, ctl)
  ll <- function(t0, t1) { if (!is.finite(t0) || !is.finite(t1) || t1 - t0 < 20 || t0 < 10 || t1 > T_end + 60) return(-Inf)
    r <- ext_run(dy, ctl, p$kp, p$ks, pri$energy_ratio, p$lR, Xextra = input_design(dy, p$kp, p$ks, data.frame(t0 = t0, t1 = t1)), gm = gm); if (is.finite(r$LL_star)) r$LL_star else -Inf }
  ll0 <- ext_run(dy, ctl, p$kp, p$ks, pri$energy_ratio, p$lR, gm = gm)$LL_star
  grid <- expand.grid(t0 = seq(0.06 * T_end, 0.9 * T_end, length.out = n_t0), d = c(40, 70, 110, 160, 220, 300, 400)); grid$t1 <- pmin(grid$t0 + grid$d, T_end); grid <- unique(grid[, c("t0", "t1")])
  v <- mapply(ll, grid$t0, grid$t1); if (!any(is.finite(v))) return(NULL); b <- which.max(v)
  o <- tryCatch(optim(c(grid$t0[b], grid$t1[b]), function(x) -ll(x[1], x[2]), method = "Nelder-Mead", control = list(maxit = 80, reltol = 1e-8)), error = function(e) NULL)
  t0 <- grid$t0[b]; t1 <- grid$t1[b]; best <- v[b]; if (!is.null(o) && is.finite(o$value) && -o$value > best) { t0 <- o$par[1]; t1 <- o$par[2]; best <- -o$value }
  if (best - ll0 < min_gain) return(NULL)
  r <- ext_run(dy, ctl, p$kp, p$ks, pri$energy_ratio, p$lR, Xextra = input_design(dy, p$kp, p$ks, data.frame(t0 = t0, t1 = t1)), gm = gm)
  list(t0 = t0, t1 = t1, gain = best - ll0, nu0 = r$beta[1], mu_p = r$beta[2] / p$kp, mu_s = r$beta[3] / p$ks, sig2 = r$sig2)
}

# ---- HMM of regimes (GPB1) ---------------------------------------------------------------------------------------------------------------------------------------
#' theta = (ln ks, ln(kp - ks), ln lR, ln sigma^2, nu0, [mu_p, mu_s for regimes 1..K-1], [stay logit for regimes 0..K-1]); regime 0 is the reference (mu = 0).
#' Priors: the pipeline's priors on the time scales and the jitter, mu ~ N(0, 0.8^2) Hz, stay logit ~ N(5, 2^2) (sticky: regimes last ~150 beats a priori).
fit_regime <- function(series, K = 3, cp = NULL, base = NULL, stay_prior = 5, stay_sd = 2, mu_sd = 0.8, mu_max = 6, anchor_sd = NULL, osc = NULL, recip = 0, multi = FALSE) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; lp <- pri$energy_ratio; gm <- gate_mult(dy, ctl); if (is.null(base)) base <- fit_ext(series); rb <- base$run()
  nm <- 5L + 2L * (K - 1L) + K; decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); mu <- rbind(c(0, 0), if (K > 1) matrix(th[6:(5 + 2 * (K - 1))], K - 1, 2, byrow = TRUE))
    stay <- plogis(th[(6 + 2 * (K - 1)):nm]); Pi <- matrix(0, K, K); for (j in 1:K) { Pi[j, ] <- if (K > 1) (1 - stay[j]) / (K - 1) else 0; Pi[j, j] <- stay[j] }
    list(ks = ks, kp = kp, lR = exp(th[3]), sig2 = exp(th[4]), nu0 = th[5], mu = mu, Pi = Pi, th = th) }
  cost <- function(th) { p <- decode(th); if (p$nu0 <= 0.2 || any(abs(p$mu) > mu_max)) return(1e12)
    r <- reg_run(dy, ctl, p$kp, p$ks, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, gm = gm, osc = osc); if (!is.finite(r$LL)) return(1e12)
    pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), list(sig2 = p$sig2, nu0 = p$nu0), dy, pri); ex <- -recip * sum(pmax(p$mu[-1, 1] * p$mu[-1, 2], 0)) - 0.5 * sum((p$mu / mu_sd)^2) - (if (is.null(anchor_sd)) 0 else 0.5 * (((log(p$ks) - log(base$p$ks)) / anchor_sd)^2 + ((log(p$kp) - log(base$p$kp)) / anchor_sd)^2)) - 0.5 * sum(((th[(6 + 2 * (K - 1)):nm] - stay_prior) / stay_sd)^2)
    v <- -(r$LL + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  nu0 <- if (!is.null(cp)) cp$nu0 else rb$nu0; mu_init <- if (K > 1) { m1 <- if (!is.null(cp)) c(cp$mu_p, cp$mu_s) else c(-0.02, 0.02); if (K == 2) m1 else c(m1, 0.4 * m1, rep(0.4 * m1, max(K - 3, 0)) ) } else numeric(0)
  if (recip > 0 && K == 2 && mu_init[1] * mu_init[2] > 0) mu_init <- c(-0.1, 0.1) * sign(mu_init[2] - mu_init[1] + 1e-9) * c(1, 1)
  mu_init <- pmax(pmin(mu_init[seq_len(2 * (K - 1))], min(3, 0.8 * mu_max)), -min(3, 0.8 * mu_max))
  th0 <- c(log(base$p$ks), log(base$p$kp - base$p$ks), log(base$p$lR), log(if (!is.null(cp)) cp$sig2 else rb$sig2), nu0, mu_init, rep(stay_prior, K))
  starts <- list(th0)
  if (multi && K == 2) { kp5 <- 1 / 5; ks0 <- exp(th0[1]); b <- th0; b[2] <- log(max(kp5 - ks0, 1e-3)); starts[[2]] <- replace(b, 6:7, c(-0.2, 0.3)); starts[[3]] <- b }   # tau_p 5 s with typical reciprocal targets / with the change-point targets
  fits <- lapply(starts, function(t0) tryCatch(optim(t0, cost, method = "BFGS", control = list(maxit = 400, parscale = c(1, 1, 1, 1, 0.1, rep(0.1, 2 * (K - 1)), rep(1, K)))), error = function(e) NULL))
  fits <- Filter(function(f) !is.null(f) && is.finite(f$value), fits); o <- if (length(fits)) fits[[which.min(vapply(fits, function(f) f$value, 0))]] else NULL
  if (is.null(o) || !is.finite(o$value) || o$value >= 1e11) return(list(ok = FALSE)); p <- decode(o$par)
  list(ok = TRUE, p = p, K = K, cost = o$value, n_par = nm, dy = dy, ctl = ctl, lp = lp, gm = gm, tau_p = 1 / p$kp, tau_s = 1 / p$ks, cp_used = !is.null(cp),
       osc = osc, run = function(dyy = dy, gmm = gm, full = TRUE, smooth = FALSE) reg_run(dyy, ctl, p$kp, p$ks, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, gm = gmm, full = full, osc = osc, smooth = smooth))
}
#' Times at which the posterior probability of leaving the reference regime crosses 0.5 upwards and stays above for at least `min_len` beats.
regime_onsets <- function(W, tt, min_len = 8) { act <- rowSums(W[, -1, drop = FALSE]) > 0.5; if (!any(act)) return(numeric(0)); r <- rle(act); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1; k <- which(r$values & r$lengths >= min_len); tt[starts[k]] }

# ---- Simulator for the regime experiments --------------------------------------------------------------------------------------------------------------------
#' scenario: "regime" (rectangular targets, rise and decay at the drives' own rates: the model is correct), "exercise" (the pipeline's exercise profile: different on/off constants),
#' "two" (two rectangular regimes, tilt-like and smaller), "rest" (no regime). Truth: states x_p, x_s (with the response) and the event times.
sim_regime <- function(scenario, seed, dur = 420, osc = NULL) {
  set.seed(seed); tau_p <- runif(1, 3.5, 6); tau_s <- runif(1, 120, 250); lp <- exp(runif(1, log(1.5), log(4))); kp <- 1 / tau_p; ks <- 1 / tau_s; ev <- NULL; mp <- ms <- NULL
  if (scenario == "regime") { t_on <- runif(1, 0.18, 0.28) * dur; t_off <- t_on + runif(1, 100, 140); a_p <- runif(1, 0.15, 0.30); a_s <- runif(1, 0.20, 0.40); ev <- data.frame(t0 = t_on, t1 = t_off)
    mp <- prof_exercise(t_on, t_off, -a_p, tau_p, tau_p); ms <- prof_exercise(t_on, t_off, a_s, tau_s, tau_s) }
  if (scenario == "exercise") { pr <- draw_exercise_pars(seed, dur); ev <- data.frame(t0 = pr$t_on, t1 = pr$t_off); mp <- prof_exercise(pr$t_on, pr$t_off, -pr$amp_p, pr$ton_p, pr$toff_p); ms <- prof_exercise(pr$t_on, pr$t_off, pr$amp_s, pr$ton_s, pr$toff_s) }
  if (scenario == "two") { dur <- 640; t1 <- runif(1, 0.12, 0.2) * dur; t2 <- t1 + runif(1, 150, 200) ; t3 <- t2 + runif(1, 90, 140); t4 <- t3 + runif(1, 150, 200); a_p <- runif(1, 0.08, 0.16); a_s <- runif(1, 0.08, 0.16); ev <- data.frame(t0 = c(t1, t3), t1 = c(t2, t4))
    mp <- function(t) prof_exercise(t1, t2, -a_p, tau_p, tau_p)(t) + prof_exercise(t3, t4, -a_p, tau_p, tau_p)(t); ms <- function(t) prof_exercise(t1, t2, a_s, tau_s, tau_s)(t) + prof_exercise(t3, t4, a_s, tau_s, tau_s)(t) }
  s <- sim_ext(dur, kp, ks, lp, 9e-4, 0.02, seed = seed, osc = osc, m_p_fun = mp, m_s_fun = ms); list(s = s, ev = ev, tau_p = tau_p, tau_s = tau_s, lp = lp, off = s$t[1] - s$rr[1], scenario = scenario)
}

#' Onsets with hysteresis: the posterior of leaving the reference regime must exceed `hi` and the episode lasts until it falls below `lo`; episodes shorter than `min_len` beats are ignored.
regime_onsets2 <- function(W, tt, hi = 0.8, lo = 0.2, min_len = 30) { p <- rowSums(W[, -1, drop = FALSE]); state <- FALSE; start <- NA_integer_; out <- numeric(0)
  for (k in seq_along(p)) { if (!state && p[k] > hi) { state <- TRUE; start <- k } else if (state && p[k] < lo) { if (k - start >= min_len) out <- c(out, tt[start]); state <- FALSE } }
  if (state && length(p) - start >= min_len) out <- c(out, tt[start]); out }
