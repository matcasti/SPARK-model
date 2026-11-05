# ==============================================================================
# experiments/parametrization/strategies.R -- EXPERIMENTAL (not part of the pipeline or the manuscript).
# Alternative parametrizations and optimizers for the M3 objective. All of them minimize the SAME
# function (the negative log-posterior of M3 in physical coordinates, see make_objective), so any
# difference between strategies is an optimization effect, not a change of the estimator.
# ==============================================================================

#' The M3 objective in physical coordinates (ks, kp, lR), with an evaluation counter.
make_objective <- function(series) {
  dy <- series$dy; pri <- series$priors; ctl <- series$control
  counter <- new.env(); counter$n <- 0L
  f <- function(ks, kp, lR) {
    ks <- unname(ks); kp <- unname(kp); lR <- unname(lR)
    counter$n <- counter$n + 1L
    if (!is.finite(ks) || !is.finite(kp) || !is.finite(lR) || ks <= 0 || kp <= ks || lR <= 0) return(1e12)
    fl <- filter_eval(dy, kp, ks, pri$energy_ratio, lR, ctl)
    if (!is.finite(fl$LL_star)) return(1e12)
    pen <- .log_prior_terms(c(ks = ks, kp = kp, lR = lR), fl, dy, pri)
    v <- -(fl$LL_star + pen$tau_s + pen$tau_p + pen$topology + pen$lR)
    if (is.finite(v)) v else 1e12
  }
  list(f = f, counter = counter)
}

# ---- Parametrizations: z (unconstrained) <-> physical (ks, kp, lR) ----------------------------------------------
#' Each entry: to(z) -> c(ks, kp, lR); from(ks, kp, lR) -> z. `pri` is the series' priors (for the whitened variant).
parametrizations <- function(pri) {
  mu_p <- pri$log_tau_p_mu; sd_p <- pri$sigma_log_tau_p
  mu_s <- pri$log_tau_s_mu; sd_s <- pri$sigma_log_tau_s
  d0 <- mu_s - mu_p; sd_d <- sqrt(sd_p^2 + sd_s^2) / d0          # ln(tau_s/tau_p) = d, eta = ln d
  list(
    # production: ln ks, ln(kp - ks), ln lR
    base = list(
      to = function(z) { ks <- exp(z[1]); c(ks, ks + exp(z[2]), exp(z[3])) },
      from = function(ks, kp, lR) c(log(ks), log(kp - ks), log(lR))),
    # geometric-mean rate and logit of the separation U = (kp - ks) / (kp + ks) (the coordinate of the Beta prior)
    mean_logit = list(
      to = function(z) { m <- exp(z[1]); U <- plogis(z[2]); c(m * sqrt((1 - U) / (1 + U)), m * sqrt((1 + U) / (1 - U)), exp(z[3])) },
      from = function(ks, kp, lR) c(0.5 * log(kp * ks), qlogis((kp - ks) / (kp + ks)), log(lR))),
    # anchor on the fast (data-informed) rate, separation as a log-ratio: ks = kp * exp(-d), d = exp(eta)
    fast_ratio = list(
      to = function(z) { kp <- exp(z[1]); c(kp * exp(-exp(z[2])), kp, exp(z[3])) },
      from = function(ks, kp, lR) c(log(kp), log(log(kp / ks)), log(lR))),
    # the same, centred and scaled by the prior spread (prior-whitened: the prior is ~N(0, 1) in z1, z2)
    whitened = list(
      to = function(z) {
        ltp <- mu_p + sd_p * z[1]; d <- exp(log(d0) + sd_d * z[2])
        c(exp(-(ltp + d)), exp(-ltp), exp(z[3])) },
      from = function(ks, kp, lR) {
        ltp <- -log(kp); d <- log(kp / ks)
        c((ltp - mu_p) / sd_p, (log(d) - log(d0)) / sd_d, log(lR)) })
  )
}

# ---- Optimizers ---------------------------------------------------------------------------------------------------
#' Run one optimizer in the coordinates of `par`, from the physical start `phys0`. Returns the optimum in physical terms.
run_local <- function(obj, par, phys0, optimizer = c("bfgs", "nlminb", "bfgs_nm", "lbfgsb"), maxit = 1000) {
  optimizer <- match.arg(optimizer)
  fz <- function(z) obj$f(par$to(z)[1], par$to(z)[2], par$to(z)[3])
  z0 <- par$from(phys0[1], phys0[2], phys0[3])
  obj$counter$n <- 0L; t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(switch(optimizer,
    bfgs    = { o <- optim(z0, fz, method = "BFGS", control = list(maxit = maxit)); list(z = o$par, v = o$value, conv = o$convergence) },
    nlminb  = { o <- nlminb(z0, fz, control = list(iter.max = maxit, eval.max = 2 * maxit)); list(z = o$par, v = o$objective, conv = o$convergence) },
    bfgs_nm = { o <- optim(z0, fz, method = "BFGS", control = list(maxit = maxit))
                p <- optim(o$par, fz, method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-12))
                if (p$value < o$value) list(z = p$par, v = p$value, conv = o$convergence) else list(z = o$par, v = o$value, conv = o$convergence) },
    lbfgsb  = { o <- optim(z0, fz, method = "L-BFGS-B", lower = z0 - 12, upper = z0 + 12, control = list(maxit = maxit)); list(z = o$par, v = o$value, conv = o$convergence) }),
    error = function(e) NULL)
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(res) || !is.finite(res$v)) return(list(ok = FALSE, evals = obj$counter$n, secs = secs))
  phys <- par$to(res$z)
  list(ok = TRUE, v = res$v, conv = res$conv, ks = phys[1], kp = phys[2], lR = phys[3], evals = obj$counter$n, secs = secs)
}

#' Coarse grid over (tau_p, tau_s, lR) around the priors; the best node is polished by BFGS in `par`.
run_grid <- function(obj, par, pri, n = 7, maxit = 1000) {
  t0 <- proc.time()[["elapsed"]]; obj$counter$n <- 0L
  tp <- exp(pri$log_tau_p_mu + pri$sigma_log_tau_p * seq(-3, 3, length.out = n))
  ts <- exp(pri$log_tau_s_mu + pri$sigma_log_tau_s * seq(-3, 3, length.out = n))
  lr <- c(1e-4, 1e-3, 1e-2)
  g <- as.matrix(expand.grid(tp = tp, ts = ts, lr = lr)); g <- g[g[, "ts"] > 1.5 * g[, "tp"], , drop = FALSE]
  v <- apply(g, 1, function(r) obj$f(1 / r[2], 1 / r[1], r[3]))
  b <- g[which.min(v), ]
  grid_evals <- obj$counter$n
  out <- run_local(obj, par, c(1 / b[["ts"]], 1 / b[["tp"]], b[["lr"]]), "bfgs", maxit)
  out$evals <- out$evals + grid_evals; out$secs <- proc.time()[["elapsed"]] - t0
  out
}

#' Variable projection: ln lR is optimized out (1-D, inner) so BFGS works in the two kinetic coordinates only.
run_profile <- function(obj, par, phys0, maxit = 1000) {
  t0 <- proc.time()[["elapsed"]]; obj$counter$n <- 0L
  inner <- function(z12, lR_hint) {
    o <- optimize(function(l) { ph <- par$to(c(z12, l)); obj$f(ph[1], ph[2], ph[3]) },
                  interval = c(log(1e-6), log(1)), tol = 1e-6)
    list(v = o$objective, l = o$minimum)
  }
  z0 <- par$from(phys0[1], phys0[2], phys0[3])
  fz <- function(z12) inner(z12)$v
  res <- tryCatch(optim(z0[1:2], fz, method = "BFGS", control = list(maxit = maxit)), error = function(e) NULL)
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(res) || !is.finite(res$value)) return(list(ok = FALSE, evals = obj$counter$n, secs = secs))
  l <- inner(res$par)$l; phys <- par$to(c(res$par, l))
  list(ok = TRUE, v = res$value, conv = res$convergence, ks = phys[1], kp = phys[2], lR = phys[3], evals = obj$counter$n, secs = secs)
}

#' Dispersed physical starts shared by all strategies: (ks, kp, lR). Half near the priors, half wide.
dispersed_starts <- function(pri, n = 6) {
  mk <- function(wide) {
    repeat {
      ltp <- if (wide) runif(1, log(1.5), log(20)) else rnorm(1, pri$log_tau_p_mu, 2 * pri$sigma_log_tau_p)
      lts <- if (wide) runif(1, log(40), log(800)) else rnorm(1, pri$log_tau_s_mu, 2 * pri$sigma_log_tau_s)
      if (exp(lts) > 1.5 * exp(ltp)) break
    }
    c(exp(-lts), exp(-ltp), exp(runif(1, log(1e-4), log(5e-2))))
  }
  lapply(seq_len(n), function(i) mk(wide = i %% 2 == 0))
}
