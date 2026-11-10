# ==============================================================================
# experiments/ou_logrr/lib.R -- EXPERIMENTAL candidate: two independent OU drives observed in ln(RR), exact per-beat transition,
# scalar double Kalman filter (the level c is profiled), RTS smoother, MAP fit, blockade layer hooks (h). Nothing here touches the
# production engine. Spec: "Hoja de implementacion OU bivariado en log-RR con capa de bloqueo" (2026-11-09).
#
#   y_k = ln RR_k = c + h_p x_p,k + h_s x_s,k + e_k,   h = (1 - b_p, -(1 - b_s)),   e_k ~ N(0, 2 sigma_t^2 / RR_k^2)
#   x_i,k = phi_i,k x_i,k-1 + w_i,k,   phi = exp(-RR_k / tau_i),   w ~ N(0, V_i (1 - phi^2)),   x_0 ~ N(0, diag(V_p, V_s))
# Time unit: seconds. Input: RR in seconds (validated). Source with:  source("experiments/ou_logrr/lib.R")
# ==============================================================================

# ---- Pieces --------------------------------------------------------------------------------------------------------------
#' Exact transition of one OU branch over dt: x' = phi x + N(0, q), parametrised with the stationary variance V.
ou_transition <- function(dt, tau, V) list(phi = exp(-dt / tau), q = V * (-expm1(-2 * dt / tau)))

#' Per-beat measurement variance of ln RR from the R-peak localisation error sigma_t (s): delta method on two time stamps.
obs_noise <- function(rr, sigma_t) 2 * sigma_t^2 / rr^2

#' Doubled scalar Kalman filter. One mean vector is driven by y, the other by the constant 1 (same covariance and gain), so that the level
#' c can be profiled out afterwards (A, B, C, D accumulators). par: tau_p, tau_s, Vp, Vs, sigma_t. skip: beats that only predict.
#' store = TRUE keeps the filtered means, the filtered and predicted covariances and phi for the smoother.
kf_run <- function(y, dt, par, h = c(1, -1), skip = NULL, R = NULL, store = FALSE) {
  n <- length(y)
  if (is.null(R)) R <- obs_noise(dt, par$sigma_t)
  if (is.null(skip)) skip <- rep(FALSE, n)
  my <- m1 <- c(0, 0)
  Ppp <- par$Vp; Pss <- par$Vs; Pps <- 0
  A <- B <- C <- D <- 0; nobs <- 0L
  if (store) { Mf <- Mf1 <- matrix(NA_real_, n, 2); Pf <- Pm <- matrix(NA_real_, n, 3); Phi <- matrix(NA_real_, n, 2) }
  for (k in seq_len(n)) {
    ph <- exp(-dt[k] / par$tau_p); ps <- exp(-dt[k] / par$tau_s)
    my[1] <- ph * my[1]; my[2] <- ps * my[2]; m1[1] <- ph * m1[1]; m1[2] <- ps * m1[2]
    Ppp <- ph^2 * Ppp + par$Vp * (-expm1(-2 * dt[k] / par$tau_p))
    Pss <- ps^2 * Pss + par$Vs * (-expm1(-2 * dt[k] / par$tau_s))
    Pps <- ph * ps * Pps
    if (store) { Pm[k, ] <- c(Ppp, Pss, Pps); Phi[k, ] <- c(ph, ps) }
    if (!skip[k]) {
      vy <- y[k] - (h[1] * my[1] + h[2] * my[2]); v1 <- 1 - (h[1] * m1[1] + h[2] * m1[2])
      a1 <- Ppp * h[1] + Pps * h[2]; a2 <- Pps * h[1] + Pss * h[2]                   # P h'
      S <- h[1] * a1 + h[2] * a2 + R[k]
      K1 <- a1 / S; K2 <- a2 / S
      my <- my + c(K1, K2) * vy; m1 <- m1 + c(K1, K2) * v1
      # Joseph form: P <- (I - K h) P (I - K h)' + R K K'
      la <- 1 - K1 * h[1]; lb <- -K1 * h[2]; lc <- -K2 * h[1]; ld <- 1 - K2 * h[2]
      r1a <- la * Ppp + lb * Pps; r1b <- la * Pps + lb * Pss; r2a <- lc * Ppp + ld * Pps; r2b <- lc * Pps + ld * Pss
      Ppp <- r1a * la + r1b * lb + R[k] * K1^2
      Pps <- r1a * lc + r1b * ld + R[k] * K1 * K2
      Pss <- r2a * lc + r2b * ld + R[k] * K2^2
      A <- A + vy * v1 / S; B <- B + v1^2 / S; C <- C + vy^2 / S; D <- D + log(S); nobs <- nobs + 1L
    }
    if (store) { Mf[k, ] <- my; Mf1[k, ] <- m1; Pf[k, ] <- c(Ppp, Pss, Pps) }
  }
  out <- list(A = A, B = B, C = C, D = D, n = nobs)
  if (store) out <- c(out, list(Mf = Mf, Mf1 = Mf1, Pf = Pf, Pm = Pm, Phi = Phi))
  out
}

#' Profiled level c = A / B and profiled log-likelihood (c integrated out by concentration).
profile_level <- function(acc) list(c = acc$A / acc$B, logLik = -0.5 * (acc$n * log(2 * pi) + acc$D + acc$C - acc$A^2 / acc$B))

#' Rauch-Tung-Striebel smoother on the filtered output of kf_run(store = TRUE) with the profiled level c.
rts_smooth <- function(acc, c_hat) {
  n <- nrow(acc$Mf)
  X <- acc$Mf - c_hat * acc$Mf1                                          # filtered states x_k|k
  X1 <- acc$Mf1; X1s <- X1                                               # the constant-1 filter, smoothed with the same gains (level uncertainty)
  P <- acc$Pf; Xs <- X; Ps <- P
  for (k in (n - 1):1) {
    ph <- acc$Phi[k + 1, ]; Pp <- acc$Pm[k + 1, ]                       # predicted covariance of beat k+1 (Ppp, Pss, Pps)
    det <- Pp[1] * Pp[2] - Pp[3]^2; if (!is.finite(det) || det <= 0) stop("non-positive determinant in the smoother at beat ", k)
    inv <- c(Pp[2], Pp[1], -Pp[3]) / det                               # inverse of the 2x2 (pp, ss, ps)
    # J = P_k F P_pred^-1 with F = diag(ph)
    PF <- c(P[k, 1] * ph[1], P[k, 2] * ph[2], P[k, 3] * ph[2], P[k, 3] * ph[1])  # [PF11, PF22, PF12, PF21]
    J11 <- PF[1] * inv[1] + PF[3] * inv[3]; J12 <- PF[1] * inv[3] + PF[3] * inv[2]
    J21 <- PF[4] * inv[1] + PF[2] * inv[3]; J22 <- PF[4] * inv[3] + PF[2] * inv[2]
    d <- Xs[k + 1, ] - ph * X[k, ]                                      # smoothed minus predicted mean
    Xs[k, ] <- X[k, ] + c(J11 * d[1] + J12 * d[2], J21 * d[1] + J22 * d[2])
    d1 <- X1s[k + 1, ] - ph * X1[k, ]; X1s[k, ] <- X1[k, ] + c(J11 * d1[1] + J12 * d1[2], J21 * d1[1] + J22 * d1[2])
    dP11 <- Ps[k + 1, 1] - Pp[1]; dP22 <- Ps[k + 1, 2] - Pp[2]; dP12 <- Ps[k + 1, 3] - Pp[3]
    # P_s = P + J dP J'
    t11 <- J11 * dP11 + J12 * dP12; t12 <- J11 * dP12 + J12 * dP22; t21 <- J21 * dP11 + J22 * dP12; t22 <- J21 * dP12 + J22 * dP22
    Ps[k, 1] <- P[k, 1] + t11 * J11 + t12 * J12
    Ps[k, 2] <- P[k, 2] + t21 * J21 + t22 * J22
    Ps[k, 3] <- P[k, 3] + t11 * J21 + t12 * J22
  }
  list(x = Xs, P = Ps, x_filt = X, P_filt = P, x1_smooth = X1s, x1_filt = X1)
}

# ---- Priors, objective, fit -----------------------------------------------------------------------------------------------
#' Log-normal-in-tau priors (Normal on ln tau). tau_prior "spark": tau_p 3-6 s, tau_s 60-100 s at +-2 sd; "effector": 0.5-2 s and 10-20 s.
ou_priors <- function(tau_prior = c("spark", "effector"), sigma_t = 0.001, lambda_P = NULL, bounds = NULL) {
  tau_prior <- match.arg(tau_prior)
  b <- if (tau_prior == "spark") list(p = c(3, 6), s = c(60, 100)) else list(p = c(0.5, 2), s = c(10, 20))
  if (!is.null(bounds)) b <- modifyList(b, bounds)                              # e.g. list(p = c(3, 60)) for a wide tau_p prior
  mu <- function(x) mean(log(x)); sd2 <- function(x) diff(log(x)) / 4                # +-2 sd
  list(tau_prior = tau_prior, mu_p = mu(b$p), sd_p = sd2(b$p), mu_s = mu(b$s), sd_s = sd2(b$s), sigma_t = sigma_t, lambda_P = lambda_P)
}

#' theta = (ln tau_p, ln(tau_s - tau_p), ln Vp, ln Vs); with a fixed ratio lambda_P (V_p / V_s) theta has three elements (ln Vs last).
theta_to_par <- function(theta, pri) {
  tp <- exp(theta[1]); ts <- tp + exp(theta[2])
  if (is.null(pri$lambda_P)) { Vp <- exp(theta[3]); Vs <- exp(theta[4]) } else { Vs <- exp(theta[3]); Vp <- pri$lambda_P * Vs }
  list(tau_p = tp, tau_s = ts, Vp = Vp, Vs = Vs, sigma_t = pri$sigma_t)
}

#' Negative log posterior (profiled log-likelihood + log-priors on ln tau_p, ln tau_s). Extreme parameters return a large finite penalty.
neg_log_post <- function(theta, y, dt, pri, h = c(1, -1), skip = NULL, R = NULL, penalty = 1e10) {
  par <- theta_to_par(theta, pri)
  val <- tryCatch({
    acc <- kf_run(y, dt, par, h = h, skip = skip, R = R)
    lp <- profile_level(acc)$logLik + dnorm(log(par$tau_p), pri$mu_p, pri$sd_p, log = TRUE) + dnorm(log(par$tau_s), pri$mu_s, pri$sd_s, log = TRUE)
    -lp
  }, error = function(e) penalty)
  if (!is.finite(val)) penalty else val
}

#' Full fit of one record (rr in seconds). Single start, BFGS, numerical Hessian on the transformed scale. Returns parameters, level,
#' filtered and smoothed states (with 95% half-widths), log-likelihood and the standardised innovations.
fit_candidate <- function(rr, pri = ou_priors(), h = c(1, -1), skip = NULL, R = NULL, start = NULL, smooth = TRUE, hessian = TRUE) {
  stopifnot(all(is.finite(rr)), all(rr > 0))
  if (median(rr) > 20) rr <- rr / 1000                                      # milliseconds
  y <- log(rr); dt <- rr; v <- var(y)
  th0 <- if (!is.null(start)) start else c(log(exp(pri$mu_p)), log(exp(pri$mu_s) - exp(pri$mu_p)),
                                             if (is.null(pri$lambda_P)) c(log(0.5 * v), log(0.5 * v)) else log(v / (1 + pri$lambda_P)))
  opt <- optim(th0, neg_log_post, y = y, dt = dt, pri = pri, h = h, skip = skip, R = R, method = "BFGS", control = list(maxit = 500, reltol = 1e-10))
  par <- theta_to_par(opt$par, pri)
  acc <- kf_run(y, dt, par, h = h, skip = skip, R = R, store = TRUE); prof <- profile_level(acc)
  lev <- function(m1, Pc, B) cbind(Pc[, 1] + m1[, 1]^2 / B, Pc[, 2] + m1[, 2]^2 / B, Pc[, 3] + m1[, 1] * m1[, 2] / B)   # + Var(c_hat) m1 m1' (Var(c_hat) = 1/B)
  out <- list(par = par, theta = opt$par, c = prof$c, logLik = prof$logLik, neg_log_post = opt$value, convergence = opt$convergence, h = h, pri = pri,
              y = y, dt = dt, t = cumsum(dt), x_filt = acc$Mf - prof$c * acc$Mf1, P_filt = acc$Pf, P_filt_c = lev(acc$Mf1, acc$Pf, acc$B), B = acc$B)
  if (hessian) {
    H <- tryCatch(optimHess(opt$par, neg_log_post, y = y, dt = dt, pri = pri, h = h, skip = skip, R = R), error = function(e) NULL)
    out$se <- if (is.null(H)) rep(NA_real_, length(opt$par)) else tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA_real_, length(opt$par)), warning = function(w) rep(NA_real_, length(opt$par)))
  }
  if (smooth) { sm <- rts_smooth(acc, prof$c); out$x_smooth <- sm$x; out$P_smooth <- sm$P; out$P_smooth_c <- lev(sm$x1_smooth, sm$P, acc$B) }
  out
}

#' Standardised one-step innovations with the level fixed: z_k = v_k / sqrt(S_k) (used by the diagnostics).
innovations <- function(fit, skip = NULL) {
  y <- fit$y - fit$c; dt <- fit$dt; par <- fit$par; h <- fit$h; n <- length(y); R <- obs_noise(dt, par$sigma_t)
  m <- c(0, 0); Ppp <- par$Vp; Pss <- par$Vs; Pps <- 0; z <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    ph <- exp(-dt[k] / par$tau_p); ps <- exp(-dt[k] / par$tau_s); m <- c(ph * m[1], ps * m[2])
    Ppp <- ph^2 * Ppp + par$Vp * (-expm1(-2 * dt[k] / par$tau_p)); Pss <- ps^2 * Pss + par$Vs * (-expm1(-2 * dt[k] / par$tau_s)); Pps <- ph * ps * Pps
    if (!is.null(skip) && skip[k]) next
    a1 <- Ppp * h[1] + Pps * h[2]; a2 <- Pps * h[1] + Pss * h[2]; S <- h[1] * a1 + h[2] * a2 + R[k]
    v <- y[k] - (h[1] * m[1] + h[2] * m[2]); z[k] <- v / sqrt(S)
    K <- c(a1, a2) / S; m <- m + K * v
    la <- 1 - K[1] * h[1]; lb <- -K[1] * h[2]; lc <- -K[2] * h[1]; ld <- 1 - K[2] * h[2]
    r1a <- la * Ppp + lb * Pps; r1b <- la * Pps + lb * Pss; r2a <- lc * Ppp + ld * Pps; r2b <- lc * Pps + ld * Pss
    Ppp <- r1a * la + r1b * lb + R[k] * K[1]^2; Pps <- r1a * lc + r1b * ld + R[k] * K[1] * K[2]; Pss <- r2a * lc + r2b * ld + R[k] * K[2]^2
  }
  z
}
