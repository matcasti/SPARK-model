# experiments/regimes/v2/lib2.R -- EXPERIMENTAL (round 2): regime filter with oscillator and smoother (reg2_filter.cpp) on top of regimes/lib.R.
source("experiments/regimes/lib.R"); suppressMessages(Rcpp::sourceCpp("experiments/regimes/v2/reg2_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))

# ---- Graded regime (chain of 3 states 0 - 1 - 2): targets mu_0 = 0, mu_1 = f mu, mu_2 = mu (f in (0, 1)), transitions only between neighbours ---------------------------------
#' A staircase for ramps and slow recoveries: the middle state is the "transition" level. theta = (ln ks, ln(kp - ks), ln lR, ln sigma^2, nu0, mu_p, mu_s, logit f, stay logits of the 3 states).
fit_regime3c <- function(series, cp = NULL, base = NULL, stay_prior = 7, stay_sd = 0.7, mid_stay = 5, osc = NULL, mu_sd = 0.8, mu_max = 6) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; lp <- pri$energy_ratio; gm <- gate_mult(dy, ctl); if (is.null(base)) base <- fit_ext(series); rb <- base$run()
  decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); f <- plogis(th[8]); m <- c(th[6], th[7]); mu <- rbind(c(0, 0), f * m, m); st <- plogis(th[9:11])
    Pi <- matrix(0, 3, 3); Pi[1, 1] <- st[1]; Pi[1, 2] <- 1 - st[1]; Pi[2, 2] <- st[2]; Pi[2, 1] <- Pi[2, 3] <- (1 - st[2]) / 2; Pi[3, 3] <- st[3]; Pi[3, 2] <- 1 - st[3]
    list(ks = ks, kp = kp, lR = exp(th[3]), sig2 = exp(th[4]), nu0 = th[5], mu = mu, Pi = Pi, f = f, th = th) }
  cost <- function(th) { p <- decode(th); if (p$nu0 <= 0.2 || any(abs(p$mu) > mu_max)) return(1e12)
    r <- reg_run(dy, ctl, p$kp, p$ks, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, gm = gm, osc = osc); if (!is.finite(r$LL)) return(1e12)
    pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), list(sig2 = p$sig2, nu0 = p$nu0), dy, pri)
    ex <- -0.5 * sum((c(th[6], th[7]) / mu_sd)^2) - 0.5 * (th[8] / 1)^2 - 0.5 * ((th[9] - stay_prior) / stay_sd)^2 - 0.5 * ((th[10] - mid_stay) / stay_sd)^2 - 0.5 * ((th[11] - stay_prior) / stay_sd)^2
    v <- -(r$LL + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  m1 <- if (!is.null(cp)) c(cp$mu_p, cp$mu_s) else c(-0.02, 0.02); m1 <- pmax(pmin(m1, 3), -3)
  th0 <- c(log(base$p$ks), log(base$p$kp - base$p$ks), log(base$p$lR), log(if (!is.null(cp)) cp$sig2 else rb$sig2), if (!is.null(cp)) cp$nu0 else rb$nu0, m1, 0, stay_prior, mid_stay, stay_prior)
  o <- tryCatch(optim(th0, cost, method = "BFGS", control = list(maxit = 500, parscale = c(1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1))), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$value) || o$value >= 1e11) return(list(ok = FALSE)); p <- decode(o$par)
  list(ok = TRUE, p = p, K = 3, cost = o$value, n_par = 11L, dy = dy, ctl = ctl, lp = lp, gm = gm, tau_p = 1 / p$kp, tau_s = 1 / p$ks, osc = osc,
       run = function(dyy = dy, gmm = gm, full = TRUE, smooth = FALSE) reg_run(dyy, ctl, p$kp, p$ks, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, gm = gmm, full = full, osc = osc, smooth = smooth))
}

# ---- Regime smoothing (Kim's backward recursion on the filtered regime probabilities) and state smoother with soft regime inputs -----------------------------------------------
#' Smoothed regime probabilities: Ws[k, j] = W[k, j] sum_i Pi[j, i] Ws[k + 1, i] / A[k + 1, i], with A the prior (predicted) probabilities of the filter.
regime_smooth <- function(W, A, Pi) { N <- nrow(W); Ws <- W; for (k in (N - 1):1) { m <- Ws[k + 1, ] / pmax(A[k + 1, ], 1e-300); v <- W[k, ] * as.numeric(Pi %*% m); Ws[k, ] <- v / sum(v) }; Ws }
#' Smoothed states of a fitted regime model h (fit_regime / fit_regime3c): regime smoothing, then the soft-input RTS smoother.
regime_states_smooth <- function(h) { r <- h$run(full = TRUE, smooth = TRUE); Ws <- regime_smooth(r$W, r$A, h$p$Pi); p <- h$p; o <- h$osc
  s <- reg2_smooth_states(h$dy, h$gm, p$kp, p$ks, h$lp, p$lR, p$sig2, p$nu0, p$mu, Ws, if (is.null(o)) 0 else o$alpha, if (is.null(o)) 0 else o$omega, if (is.null(o)) 0 else o$lo); s$Ws <- Ws; s$filt <- r; s }

# ---- Regimes with their own drive rates ("rest - active - recovery"): reg3_filter.cpp -----------------------------------------------------------------------------------------
suppressMessages(Rcpp::sourceCpp("experiments/regimes/v2/reg3_filter.cpp", cacheDir = "experiments/m3plus/_out/cache"))
reg3_run <- function(dy, ctl, kp, ks, lp, lR, sig2, nu0, mu, Pi, w0, gm, osc = NULL, full = FALSE)
  reg3_filter(dy, gm, kp, ks, lp, lR, sig2, nu0, mu, Pi, w0, if (is.null(osc)) 0 else osc$alpha, if (is.null(osc)) 0 else osc$omega, if (is.null(osc)) 0 else osc$lo, ctl$lik_gate_max, full)
#' Cycle 0 (rest) -> 1 (active, targets mu) -> 2 (recovery: target 0, drive rates multiplied by c_p, c_s <= 1) -> 0; the active regime may also return directly to rest.
#' theta = (ln ks, ln(kp - ks), ln lR, ln sigma^2, nu0, mu_p, mu_s, ln c_p, ln c_s, stay logits s0, s1, s2, logit q), q = P(1 -> 2 | leaving 1). Priors: c ~ lognormal(ln 0.5, 0.7), stay logits N(stay_prior, stay_sd).
fit_regime_rec <- function(series, cp = NULL, base = NULL, stay_prior = 7, stay_sd = 0.7, osc = NULL, mu_sd = 0.8, mu_max = 6, c_prior = log(0.5), c_sd = 0.7) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; lp <- pri$energy_ratio; gm <- gate_mult(dy, ctl); if (is.null(base)) base <- fit_ext(series); rb <- base$run()
  decode <- function(th) { ks <- exp(th[1]); kp <- ks + exp(th[2]); cp_ <- exp(th[8]); cs_ <- exp(th[9]); m <- c(th[6], th[7]); mu <- rbind(c(0, 0), m, c(0, 0)); st <- plogis(th[10:12]); q <- plogis(th[13])
    Pi <- matrix(0, 3, 3); Pi[1, 1] <- st[1]; Pi[1, 2] <- 1 - st[1]; Pi[2, 2] <- st[2]; Pi[2, 3] <- (1 - st[2]) * q; Pi[2, 1] <- (1 - st[2]) * (1 - q); Pi[3, 3] <- st[3]; Pi[3, 1] <- 1 - st[3]
    list(ks = ks, kp = kp, kpv = c(kp, kp, kp * cp_), ksv = c(ks, ks, ks * cs_), c_p = cp_, c_s = cs_, lR = exp(th[3]), sig2 = exp(th[4]), nu0 = th[5], mu = mu, Pi = Pi, th = th) }
  cost <- function(th) { p <- decode(th); if (p$nu0 <= 0.2 || any(abs(p$mu) > mu_max) || th[8] > 0.7 || th[9] > 0.7 || th[8] < -4 || th[9] < -4) return(1e12)
    r <- reg3_run(dy, ctl, p$kpv, p$ksv, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, c(1, 0, 0), gm, osc); if (!is.finite(r$LL)) return(1e12)
    pen <- .log_prior_terms(c(ks = p$ks, kp = p$kp, lR = p$lR), list(sig2 = p$sig2, nu0 = p$nu0), dy, pri)
    ex <- -0.5 * sum((c(th[6], th[7]) / mu_sd)^2) - 0.5 * sum(((th[8:9] - c_prior) / c_sd)^2) - 0.5 * sum(((th[10:12] - stay_prior) / stay_sd)^2) - 0.5 * (th[13] / 1.5)^2
    v <- -(r$LL + pen$tau_s + pen$tau_p + pen$topology + pen$lR + ex); if (is.finite(v)) v else 1e12 }
  m1 <- if (!is.null(cp)) c(cp$mu_p, cp$mu_s) else c(-0.02, 0.02); m1 <- pmax(pmin(m1, 3), -3)
  th0 <- c(log(base$p$ks), log(base$p$kp - base$p$ks), log(base$p$lR), log(if (!is.null(cp)) cp$sig2 else rb$sig2), if (!is.null(cp)) cp$nu0 else rb$nu0, m1, c_prior, c_prior, stay_prior, stay_prior, stay_prior, 1)
  o <- tryCatch(optim(th0, cost, method = "BFGS", control = list(maxit = 500, parscale = c(1, 1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1))), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$value) || o$value >= 1e11) return(list(ok = FALSE)); p <- decode(o$par)
  list(ok = TRUE, p = p, K = 3, cost = o$value, opt = o, n_par = 13L, dy = dy, ctl = ctl, lp = lp, gm = gm, tau_p = 1 / p$kp, tau_s = 1 / p$ks, osc = osc,
       run = function(dyy = dy, gmm = gm, full = TRUE) reg3_run(dyy, ctl, p$kpv, p$ksv, lp, p$lR, p$sig2, p$nu0, p$mu, p$Pi, c(1, 0, 0), gmm, osc, full))
}

# ---- Hybrid: the regime model dates the transitions, the events model (M3E) estimates the states -----------------------------------------------------------------------------------
#' Episodes where the smoothed probability of the non-reference regime exceeds `thr` for at least `min_len` beats: data.frame(t0, t1) (NULL if none). At most `max_ev` longest episodes are kept.
regime_episodes <- function(Ws, tt, thr = 0.5, min_len = 30, max_ev = 6) { a <- Ws[, 2] > thr; if (!any(a)) return(NULL); r <- rle(a); e <- cumsum(r$lengths); s <- e - r$lengths + 1; k <- which(r$values & r$lengths >= min_len); if (!length(k)) return(NULL)
  k <- k[order(-r$lengths[k])][seq_len(min(max_ev, length(k)))]; k <- sort(k); data.frame(t0 = tt[pmax(s[k], 1)], t1 = tt[e[k]]) }
