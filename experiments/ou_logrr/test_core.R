# ==============================================================================
# experiments/ou_logrr/test_core.R -- unit tests of the core (step 1 of the roadmap). Run from the project root:
#   Rscript experiments/ou_logrr/test_core.R
# ==============================================================================
source("experiments/ou_logrr/lib.R")
n_fail <- 0L
check <- function(name, value, tol) { ok <- isTRUE(value <= tol); if (!ok) n_fail <<- n_fail + 1L; cat(sprintf("%-6s %-62s %.3g (tol %.0e)\n", if (ok) "ok" else "FAIL", name, value, tol)) }
set.seed(1)
par <- list(tau_p = 4.5, tau_s = 75, Vp = 0.00125, Vs = 0.00125, sigma_t = 0.001)
rr <- pmax(0.4, 0.85 + cumsum(rnorm(200, 0, 0.01)) * 0 + rnorm(200, 0, 0.03)); y <- log(rr)       # arbitrary positive RR series; the tests are algebraic

# 1. Transition: one step of d1 + d2 equals two successive steps (mean factor and variance)
d1 <- 0.7; d2 <- 1.9
for (b in c("p", "s")) {
  tau <- par[[paste0("tau_", b)]]; V <- par[[paste0("V", b)]]
  one <- ou_transition(d1 + d2, tau, V); a <- ou_transition(d1, tau, V); bb <- ou_transition(d2, tau, V)
  check(sprintf("transition %s: phi(d1+d2) = phi(d1) phi(d2)", b), abs(one$phi - a$phi * bb$phi), 1e-10)
  check(sprintf("transition %s: variance of one step = two steps", b), abs(one$q - (bb$phi^2 * a$q + bb$q)), 1e-10)
}

# 2. Likelihood: filter with c known against the dense multivariate normal (N = 200)
tk <- cumsum(rr); R <- obs_noise(rr, par$sigma_t); h <- c(1, -1); cc <- 0.3
Sig <- outer(tk, tk, function(a, b) h[1]^2 * par$Vp * exp(-abs(a - b) / par$tau_p) + h[2]^2 * par$Vs * exp(-abs(a - b) / par$tau_s)) + diag(R)
Lc <- chol(Sig); r <- backsolve(Lc, y - cc, transpose = TRUE)
dense <- -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(Lc))) + sum(r^2))
acc <- kf_run(y, rr, par, h = h)
filt <- -0.5 * (acc$n * log(2 * pi) + acc$D + acc$C - 2 * cc * acc$A + cc^2 * acc$B)
check("likelihood: filter vs dense normal (c fixed)", abs(filt - dense), 1e-8)

# 3. Profile: c = A / B against numerical optimisation over c
f <- function(cv) -0.5 * (acc$n * log(2 * pi) + acc$D + acc$C - 2 * cv * acc$A + cv^2 * acc$B)
num <- optimize(f, c(-2, 2), maximum = TRUE, tol = 1e-12)$maximum
check("profile: c = A/B vs optimize", abs(profile_level(acc)$c - num), 1e-6)

# 4. Covariance symmetric positive definite at every beat (stored filtered and predicted)
st <- kf_run(y, rr, par, h = h, store = TRUE)
pd <- function(M) all(M[, 1] > 0 & M[, 2] > 0 & M[, 1] * M[, 2] - M[, 3]^2 > 0)
check("covariance: filtered and predicted positive definite (0 = yes)", as.numeric(!(pd(st$Pf) && pd(st$Pm))), 0)

# 5. Smoother: variance not above the filtered one; last point identical
sm <- rts_smooth(st, profile_level(st)$c)
check("smoother: max(smoothed var - filtered var) <= 0 (p)", max(sm$P[, 1] - sm$P_filt[, 1]), 1e-12)
check("smoother: max(smoothed var - filtered var) <= 0 (s)", max(sm$P[, 2] - sm$P_filt[, 2]), 1e-12)
n <- length(y)
check("smoother: last point equals the filtered one", max(abs(sm$x[n, ] - sm$x_filt[n, ]), abs(sm$P[n, ] - sm$P_filt[n, ])), 1e-10)

# 6. Artifacts: skipping the update of beat k equals merging its interval with the next (observation k omitted)
k <- 80; skip <- rep(FALSE, length(y)); skip[k] <- TRUE
a1 <- kf_run(y, rr, par, h = h, skip = skip, R = R)
keep <- -k; dt2 <- rr[keep]; dt2[k] <- rr[k] + rr[k + 1]                      # after dropping beat k, the former beat k+1 sits at index k
a2 <- kf_run(y[keep], dt2, par, h = h, R = R[keep])
l1 <- profile_level(a1)$logLik; l2 <- profile_level(a2)$logLik
check("artifacts: skip vs merged interval (logLik)", abs(l1 - l2), 1e-10)
check("artifacts: skip vs merged interval (level c)", abs(profile_level(a1)$c - profile_level(a2)$c), 1e-10)

# 7. Smoother against the dense Gaussian posterior (mean and variance of both states at every beat)
c_hat <- profile_level(st)$c; Si <- chol2inv(chol(Sig)); res <- y - c_hat
err_m <- 0; err_v <- 0
for (i in 1:2) {
  Vi <- if (i == 1) par$Vp else par$Vs; tau <- if (i == 1) par$tau_p else par$tau_s
  Cxy <- sweep(exp(-abs(outer(tk, tk, "-")) / tau) * Vi, 2, h[i], "*")           # Cov(x_i(t_k), y_j)
  post_m <- drop(Cxy %*% Si %*% res); post_v <- Vi - rowSums((Cxy %*% Si) * Cxy)
  err_m <- max(err_m, abs(sm$x[, i] - post_m)); err_v <- max(err_v, abs(sm$P[, i] - post_v))
}
check("smoother: mean vs dense Gaussian posterior", err_m, 1e-8)
check("smoother: variance vs dense Gaussian posterior", err_v, 1e-8)

# 8. Level uncertainty (flat prior on c): B = 1' Sigma^-1 1 and the smoothed constant-1 filter equals Cov(x, y) Sigma^-1 1
one <- rep(1, length(y)); B_dense <- drop(t(one) %*% Si %*% one)
check("level: B vs 1' Sigma^-1 1", abs(st$B - B_dense) / B_dense, 1e-9)
err_g <- 0
for (i in 1:2) { Vi <- if (i == 1) par$Vp else par$Vs; tau <- if (i == 1) par$tau_p else par$tau_s
  Cxy <- sweep(exp(-abs(outer(tk, tk, "-")) / tau) * Vi, 2, h[i], "*"); err_g <- max(err_g, abs(sm$x1_smooth[, i] - drop(Cxy %*% Si %*% one))) }
check("level: smoothed constant-1 filter vs dense Cov(x,y) Sigma^-1 1", err_g, 1e-8)

cat(if (n_fail == 0L) "\nAll core tests passed.\n" else sprintf("\n%d core test(s) FAILED.\n", n_fail))
if (n_fail > 0L) quit(status = 1)
