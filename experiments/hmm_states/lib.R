# ==============================================================================
# experiments/hmm_states/lib.R -- second-stage hidden Markov model on the estimated latent states, with their uncertainty.
# Source from the project root:  source("experiments/hmm_states/lib.R")
#
# Stage 1 is production M3 (nothing changes): filtered states x_p, x_s with their estimation variances. Stage 2 (here): the states, averaged on a 5-s grid, are the observations of a
# K-regime hidden Markov chain. "aware" emission (model-driven): y_t | regime k ~ N(mu_k, V + v_t) per drive, where V is the stationary variance of the drive FROM THE M3 FIT
# (V_p = lambda_P sigma^2, V_s = sigma^2: the fluctuation of the OU drive around a regime mean) and v_t is the estimation variance of the filter at that time: only the regime means
# are free, and a poorly determined state weighs little. "naive" emission: free diagonal variance per regime, uncertainty ignored (a standard Gaussian HMM).
# EM with a sticky Dirichlet prior on the transitions, scaled forward-backward, Viterbi.
# ==============================================================================
suppressMessages(library(data.table))

#' Average the beat-wise states on a regular grid (`by` seconds). The variance of a bin is the mean of the filter variances (the errors within a few seconds are strongly correlated).
bin_states <- function(t, xp, xs, vp, vs, by = 5) { b <- floor((t - min(t)) / by)
  d <- data.table(b, t, xp, xs, vp, vs)[, .(t = mean(t), y1 = mean(xp), y2 = mean(xs), v1 = mean(vp), v2 = mean(vs), n = .N), by = b][n >= 2]; d[, .(t, y1, y2, v1, v2)] }

.emission <- function(Y, Vt, mu, V, Sig, aware) {                                # log N(y_t | k), summed over the two drives; T x K
  K <- nrow(mu); lb <- matrix(0, nrow(Y), K)
  for (k in seq_len(K)) for (d in 1:2) { s2 <- if (aware) V[d] + Vt[, d] else Sig[k, d]; lb[, k] <- lb[, k] + dnorm(Y[, d], mu[k, d], sqrt(s2), log = TRUE) }
  lb }

.fb <- function(lb, A, pi0) {                                                    # scaled forward-backward: loglik, gamma, summed xi
  T <- nrow(lb); K <- ncol(lb); m <- apply(lb, 1, max); B <- exp(lb - m); alpha <- matrix(0, T, K); cs <- numeric(T)
  a <- pi0 * B[1, ]; cs[1] <- sum(a); alpha[1, ] <- a / cs[1]
  for (t in 2:T) { a <- as.numeric(alpha[t - 1, ] %*% A) * B[t, ]; cs[t] <- sum(a); alpha[t, ] <- a / cs[t] }
  beta <- matrix(1, T, K); for (t in (T - 1):1) beta[t, ] <- as.numeric(A %*% (B[t + 1, ] * beta[t + 1, ])) / cs[t + 1]
  gamma <- alpha * beta; gamma <- gamma / rowSums(gamma); xi <- matrix(0, K, K)
  for (t in 1:(T - 1)) { x <- outer(alpha[t, ], B[t + 1, ] * beta[t + 1, ]) * A / cs[t + 1]; xi <- xi + x / sum(x) }
  list(ll = sum(log(cs)) + sum(m), gamma = gamma, xi = xi) }

.viterbi <- function(lb, A, pi0) { T <- nrow(lb); K <- ncol(lb); lA <- log(A); d <- matrix(-Inf, T, K); psi <- matrix(0L, T, K); d[1, ] <- log(pi0) + lb[1, ]
  for (t in 2:T) for (k in 1:K) { v <- d[t - 1, ] + lA[, k]; psi[t, k] <- which.max(v); d[t, k] <- max(v) + lb[t, k] }
  p <- integer(T); p[T] <- which.max(d[T, ]); for (t in (T - 1):1) p[t] <- psi[t + 1, p[t + 1]]; p }

#' Fit a K-regime HMM to the binned states. Y: T x 2 matrix; Vt: T x 2 estimation variances; V: stationary variance of the two drives (aware mode).
#' kappa: pseudo-count of staying (sticky prior); starts: k-means initialisations. Returns loglik, number of parameters, regime means, transitions, posterior (gamma) and Viterbi path.
hmm_fit <- function(Y, Vt, V, K, mode = c("aware", "naive"), kappa = 10, starts = 4, iters = 150, tol = 1e-6) {
  mode <- match.arg(mode); aware <- mode == "aware"; T <- nrow(Y); Z <- cbind(Y[, 1] / sqrt(V[1] + Vt[, 1]), Y[, 2] / sqrt(V[2] + Vt[, 2]))
  if (K == 1) { w <- 1 / (if (aware) cbind(V[1] + Vt[, 1], V[2] + Vt[, 2]) else matrix(1, T, 2)); mu <- matrix(colSums(Y * w) / colSums(w), 1); Sig <- matrix(apply(Y, 2, var), 1)
    lb <- .emission(Y, Vt, mu, V, Sig, aware); return(list(ll = sum(lb), npar = if (aware) 2 else 4, K = 1, mu = mu, A = matrix(1), gamma = matrix(1, T, 1), path = rep(1L, T), mode = mode)) }
  best <- NULL
  for (s in seq_len(starts)) {
    km <- tryCatch(kmeans(Z, K, nstart = 1, iter.max = 30), error = function(e) NULL); if (is.null(km)) next
    cl <- km$cluster; mu <- t(sapply(seq_len(K), function(k) colMeans(Y[cl == k, , drop = FALSE]))); mu[is.na(mu)] <- 0; Sig <- matrix(rep(apply(Y, 2, var), each = K), K)
    A <- matrix((1 - 0.95) / (K - 1), K, K); diag(A) <- 0.95; pi0 <- rep(1 / K, K); ll_old <- -Inf
    for (it in seq_len(iters)) { lb <- .emission(Y, Vt, mu, V, Sig, aware); r <- .fb(lb, A, pi0); g <- r$gamma
      for (k in 1:K) for (d in 1:2) { w <- if (aware) 1 / (V[d] + Vt[, d]) else rep(1, T); den <- sum(g[, k] * w); if (den > 1e-8) mu[k, d] <- sum(g[, k] * w * Y[, d]) / den
        if (!aware) Sig[k, d] <- max(sum(g[, k] * (Y[, d] - mu[k, d])^2) / max(sum(g[, k]), 1e-8), 1e-3 * var(Y[, d])) }
      cnt <- r$xi + kappa * diag(K) + 0.5 * (1 - diag(K)); A <- cnt / rowSums(cnt); pi0 <- (g[1, ] + 1e-3) / sum(g[1, ] + 1e-3)
      if (is.finite(r$ll) && abs(r$ll - ll_old) < tol) break; ll_old <- r$ll }
    lb <- .emission(Y, Vt, mu, V, Sig, aware); r <- .fb(lb, A, pi0)
    if (is.null(best) || r$ll > best$ll) best <- list(ll = r$ll, npar = (if (aware) 2 * K else 4 * K) + K * (K - 1) + (K - 1), K = K, mu = mu, Sig = Sig, A = A, pi0 = pi0, gamma = r$gamma, path = .viterbi(lb, A, pi0), mode = mode) }
  best }

#' Baseline regime = the one occupying most of the first `n0` bins; the first leave (>= `min_run` bins away from it) and the return afterwards, in seconds.
regime_events <- function(fit, t, n0 = 24L, min_run = 6L) { p <- fit$path; if (fit$K < 2) return(list(base = 1L, onset = NA_real_, end = NA_real_, away = 0))
  base <- as.integer(names(which.max(table(p[seq_len(min(n0, length(p)))])))); away <- p != base; r <- rle(away); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  j <- which(r$values & r$lengths >= min_run)[1]; if (is.na(j)) return(list(base = base, onset = NA_real_, end = NA_real_, away = mean(away)))
  list(base = base, onset = t[starts[j]], end = if (ends[j] < length(p)) t[ends[j] + 1] else NA_real_, away = mean(away)) }

#' Gaussian-mixture-free detection statistic: twice the gain in log-likelihood of K = 2 over K = 1 (same emission mode).
detect_lr <- function(b, V, mode = "aware", K = 2) { Y <- as.matrix(b[, .(y1, y2)]); Vt <- as.matrix(b[, .(v1, v2)]); f1 <- hmm_fit(Y, Vt, V, 1, mode); f2 <- hmm_fit(Y, Vt, V, K, mode)
  list(lr = 2 * (f2$ll - f1$ll), fit1 = f1, fit2 = f2) }
