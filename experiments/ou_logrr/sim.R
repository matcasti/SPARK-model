# ==============================================================================
# experiments/ou_logrr/sim.R -- synthetic data with known truth and parameter / state recovery (step 2 of the roadmap).
#   Rscript experiments/ou_logrr/sim.R [E1|E2] [reps] [duration=<s>] [variant=base|wide|oracle|ratio<m>|wide_oracle|wide_ratio<m>]      writes experiments/ou_logrr/_out/recovery_<scenario>.rds and prints a summary
# Generator G1 (inside the model): exact OU states on a 5 ms grid; a beat closes at the first grid time at which the elapsed time reaches
# exp(c + x_p - x_s); R-peak times get a detection error N(0, sigma_t^2); RR are the differences of the noisy times.
# ==============================================================================
source("experiments/ou_logrr/lib.R")

simulate_beats <- function(par, c0 = log(0.85), duration = 600, offset = NULL, step = 0.005) {
  tt <- seq(0, duration + 6, by = step); n <- length(tt)
  ou <- function(tau, V) { phi <- exp(-step / tau); as.numeric(stats::filter(c(rnorm(1, 0, sqrt(V)), rnorm(n - 1, 0, sqrt(V * (1 - phi^2)))), phi, method = "recursive")) }
  off <- if (is.null(offset)) list(p = rep(0, n), s = rep(0, n)) else offset(tt)
  sp <- ou(par$tau_p, par$Vp) + off$p; ss <- ou(par$tau_s, par$Vs) + off$s
  target <- exp(c0 + sp - ss)
  i <- 1L; idx <- integer(0); win <- 1000L
  repeat {
    j <- (i + 1L):min(n, i + win); hit <- which((j - i) * step >= target[j])[1]
    if (is.na(hit)) break
    i <- j[hit]; if (tt[i] > duration) break; idx <- c(idx, i)
  }
  obs <- c(0, tt[idx]) + rnorm(length(idx) + 1L, 0, par$sigma_t)
  list(rr = diff(obs), t = tt[idx], xp = sp[idx], xs = ss[idx], c = c0)
}

#' E1: stationary, 10 min. E2: E1 plus 300 s rest, 120 s exercise (fast vagal step to -0.25, slow linear sympathetic ramp to +0.30), 300 s recovery.
scenario <- function(name, duration = NULL) {
  par <- list(tau_p = 4.5, tau_s = 75, Vp = 0.00125, Vs = 0.00125, sigma_t = 0.001)
  if (name == "E1") return(list(par = par, duration = if (is.null(duration)) 600 else duration, offset = NULL))
  offset <- function(t) {
    ex <- t >= 300 & t < 420; rec <- t >= 420
    p <- ifelse(ex, -0.25 * (1 - exp(-(t - 300) / 5)), ifelse(rec, -0.25 * (1 - exp(-120 / 5)) * exp(-(t - 420) / 20), 0))
    s <- ifelse(ex, 0.30 * (t - 300) / 120, ifelse(rec, 0.30 * exp(-(t - 420) / 50), 0))
    list(p = p, s = s)
  }
  list(par = par, duration = 720, offset = offset)
}

recover_one <- function(sc, seed, pri = ou_priors()) {
  set.seed(seed); d <- simulate_beats(sc$par, duration = sc$duration, offset = sc$offset)
  secs_fit <- system.time(fit <- fit_candidate(d$rr, pri))[["elapsed"]]
  z <- innovations(fit); lb <- tryCatch(Box.test(z, lag = 20, type = "Ljung-Box")$p.value, error = function(e) NA_real_)
  n <- min(length(d$xp), nrow(fit$x_smooth)); xs <- fit$x_smooth[seq_len(n), ]; sdp <- sqrt(fit$P_smooth[seq_len(n), 1]); sds <- sqrt(fit$P_smooth[seq_len(n), 2])
  tv <- cbind(d$xp[seq_len(n)], d$xs[seq_len(n)])
  chg <- function(x, tb) { rest <- tb >= 150 & tb < 300; endx <- tb >= 405 & tb < 420; mean(x[endx]) - mean(x[rest]) }
  ex <- sc$duration > 700                                                   # E2: change from rest to the end of exercise, estimated and true
  tb <- d$t[seq_len(n)]
  chg_p <- if (ex) c(chg(xs[, 1], tb), chg(tv[, 1], tb)) else c(NA_real_, NA_real_); chg_s <- if (ex) c(chg(xs[, 2], tb), chg(tv[, 2], tb)) else c(NA_real_, NA_real_)
  data.frame(seed = seed, chg_p_est = chg_p[1], chg_p_true = chg_p[2], chg_s_est = chg_s[1], chg_s_true = chg_s[2], n = n, conv = fit$convergence, secs = secs_fit, tau_p = fit$par$tau_p, tau_s = fit$par$tau_s, Vp = fit$par$Vp, Vs = fit$par$Vs, c_hat = fit$c, c_true = d$c,
             rmse_p = sqrt(mean((xs[, 1] - tv[, 1])^2)), rmse_s = sqrt(mean((xs[, 2] - tv[, 2])^2)), cor_p = cor(xs[, 1], tv[, 1]), cor_s = cor(xs[, 2], tv[, 2]),
             cover_p = mean(abs(xs[, 1] - tv[, 1]) <= 1.96 * sdp), cover_s = mean(abs(xs[, 2] - tv[, 2]) <= 1.96 * sds),
             z_mean = mean(z, na.rm = TRUE), z_var = var(z, na.rm = TRUE), lb_p = lb, logLik = fit$logLik)
}

summarise_recovery <- function(r, sc) {
  q <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x, na.rm = TRUE), quantile(x, .1, na.rm = TRUE), quantile(x, .9, na.rm = TRUE))
  tp <- sc$par
  cat(sprintf("\n%d records; converged %d; median fit %.1f s; beats %d\n", nrow(r), sum(r$conv == 0), median(r$secs), round(median(r$n))))
  cat("median [P10, P90] of estimate / truth\n")
  cat(sprintf("  tau_p  %s   (accept 0.67-1.5)\n  tau_s  %s   (accept 0.67-1.5)\n  V_p    %s   (accept 0.75-1.33)\n  V_s    %s   (accept 0.75-1.33)\n",
              q(r$tau_p / tp$tau_p), q(r$tau_s / tp$tau_s), q(r$Vp / tp$Vp), q(r$Vs / tp$Vs)))
  cat(sprintf("  c (estimate - truth)  %s\n", q(r$c_hat - r$c_true)))
  cat("smoothed states vs truth\n")
  cat(sprintf("  RMSE   p %s   s %s\n  cor    p %s   s %s\n  cover  p %s   s %s   (accept 0.90-0.98)\n", q(r$rmse_p), q(r$rmse_s), q(r$cor_p), q(r$cor_s), q(r$cover_p), q(r$cover_s)))
  if (any(is.finite(r$chg_p_true))) cat(sprintf("rest -> end of exercise (true change: p %.3f, s %.3f): estimated / true  p %s   s %s\n", median(r$chg_p_true), median(r$chg_s_true), q(r$chg_p_est / r$chg_p_true), q(r$chg_s_est / r$chg_s_true)))
  cat(sprintf("innovations: |mean| %s (accept < 0.05); variance %s (0.9-1.1); Ljung-Box p > 0.05 in %.0f%%\n", q(abs(r$z_mean)), q(r$z_var), 100 * mean(r$lb_p > 0.05, na.rm = TRUE)))
}

#' Prior variants for diagnosis: base (the spec), wide (tau_p 3-60 s), ratio<m> (V_p / V_s fixed at m times the true ratio; "oracle" = ratio1),
#' and their combinations, e.g. wide_ratio3 or wide_oracle.
make_priors <- function(variant, par) {
  wide <- grepl("wide", variant)
  m <- if (grepl("oracle", variant)) 1 else if (grepl("ratio", variant)) as.numeric(sub(".*ratio", "", variant)) else NA_real_
  ou_priors(bounds = if (wide) list(p = c(3, 60)) else NULL, lambda_P = if (is.finite(m)) m * par$Vp / par$Vs else NULL)
}

if (!interactive() && sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE); name <- if (length(a)) a[1] else "E1"; reps <- if (length(a) > 1) as.integer(a[2]) else 40L
  opt <- function(k, default) { hit <- grep(paste0("^", k, "="), a, value = TRUE); if (length(hit)) sub(".*=", "", hit[1]) else default }
  dur <- opt("duration", NA); dur <- if (is.na(dur)) NULL else as.numeric(dur); variant <- opt("variant", "base")
  sc <- scenario(name, dur); dir.create("experiments/ou_logrr/_out", showWarnings = FALSE, recursive = TRUE)
  pri <- make_priors(variant, sc$par)
  r <- do.call(rbind, lapply(seq_len(reps), function(i) recover_one(sc, seed = 1000 + i, pri = pri)))
  tag <- paste0(name, if (!is.null(dur)) paste0("_", dur) else "", if (variant != "base") paste0("_", variant) else "")
  saveRDS(list(results = r, scenario = name, variant = variant), sprintf("experiments/ou_logrr/_out/recovery_%s.rds", tag))
  cat(sprintf("scenario %s, prior variant %s\n", name, variant)); summarise_recovery(r, sc)
}
