# ==============================================================================
# R/signals/regimes.R -- second-stage hidden Markov model on the filtered states, with their estimation uncertainty (supplementary "Perturbation-response
# descriptors"). The states, averaged on a 5-s grid, are the observations of a K-regime chain; the "aware" emission is y_t | regime k ~ N(mu_k, V + v_t) per drive,
# where V is the stationary variance of the drive (from an M3 fit) and v_t the estimation variance of the filter, so only the regime means are free and a poorly
# determined state weighs little. EM with a sticky Dirichlet prior on the transitions, scaled forward-backward, Viterbi.
#
#   bin_states()   beat-wise states and variances -> 5-s bins (t, y1, y2, v1, v2)
#   hmm_fit()      K-regime fit: log-likelihood, number of parameters, means, transitions, posterior and Viterbi path
# ==============================================================================

#' Average the states on a regular grid (`by` s). The variance of a bin is the mean of the filter variances (errors within a few seconds are strongly correlated).
bin_states <- function(t, xp, xs, vp, vs, by = 5) {
  b <- floor((t - min(t)) / by)
  d <- data.table(b, t, xp, xs, vp, vs)[, .(t = mean(t), y1 = mean(xp), y2 = mean(xs), v1 = mean(vp), v2 = mean(vs), n = .N), by = b][n >= 2]
  d[, .(t, y1, y2, v1, v2)]
}

hmm_emission <- function(Y, Vt, mu, V) {                                       # log N(y_t | k), summed over the two drives (T x K)
  K <- nrow(mu); lb <- matrix(0, nrow(Y), K)
  for (k in seq_len(K)) for (d in 1:2) lb[, k] <- lb[, k] + dnorm(Y[, d], mu[k, d], sqrt(V[d] + Vt[, d]), log = TRUE)
  lb
}

hmm_forward_backward <- function(lb, A, pi0) {                                 # scaled forward-backward: log-likelihood, posterior, summed transition counts
  n <- nrow(lb); K <- ncol(lb); m <- apply(lb, 1, max); B <- exp(lb - m); alpha <- matrix(0, n, K); cs <- numeric(n)
  a <- pi0 * B[1, ]; cs[1] <- sum(a); alpha[1, ] <- a / cs[1]
  for (t in 2:n) { a <- as.numeric(alpha[t - 1, ] %*% A) * B[t, ]; cs[t] <- sum(a); alpha[t, ] <- a / cs[t] }
  beta <- matrix(1, n, K); for (t in (n - 1):1) beta[t, ] <- as.numeric(A %*% (B[t + 1, ] * beta[t + 1, ])) / cs[t + 1]
  gamma <- alpha * beta; gamma <- gamma / rowSums(gamma); xi <- matrix(0, K, K)
  for (t in 1:(n - 1)) { x <- outer(alpha[t, ], B[t + 1, ] * beta[t + 1, ]) * A / cs[t + 1]; xi <- xi + x / sum(x) }
  list(ll = sum(log(cs)) + sum(m), gamma = gamma, xi = xi)
}

hmm_viterbi <- function(lb, A, pi0) {
  n <- nrow(lb); K <- ncol(lb); lA <- log(A); d <- matrix(-Inf, n, K); psi <- matrix(0L, n, K); d[1, ] <- log(pi0) + lb[1, ]
  for (t in 2:n) for (k in 1:K) { v <- d[t - 1, ] + lA[, k]; psi[t, k] <- which.max(v); d[t, k] <- max(v) + lb[t, k] }
  p <- integer(n); p[n] <- which.max(d[n, ]); for (t in (n - 1):1) p[t] <- psi[t + 1, p[t + 1]]; p
}

#' K-regime aware HMM. Y: T x 2 binned states; Vt: T x 2 estimation variances; V: stationary variances of the two drives. kappa: pseudo-count of staying (sticky prior);
#' starts: k-means initialisations. @return list(ll, npar, K, mu, A, gamma, path)
hmm_fit <- function(Y, Vt, V, K, kappa = 10, starts = 4, iters = 150, tol = 1e-6) {
  n <- nrow(Y)
  if (K == 1) { w <- 1 / cbind(V[1] + Vt[, 1], V[2] + Vt[, 2]); mu <- matrix(colSums(Y * w) / colSums(w), 1)
    return(list(ll = sum(hmm_emission(Y, Vt, mu, V)), npar = 2, K = 1, mu = mu, A = matrix(1), gamma = matrix(1, n, 1), path = rep(1L, n))) }
  Z <- cbind(Y[, 1] / sqrt(V[1] + Vt[, 1]), Y[, 2] / sqrt(V[2] + Vt[, 2])); best <- NULL
  for (s in seq_len(starts)) {
    km <- tryCatch(kmeans(Z, K, nstart = 1, iter.max = 30), error = function(e) NULL); if (is.null(km)) next
    mu <- t(sapply(seq_len(K), function(k) colMeans(Y[km$cluster == k, , drop = FALSE]))); mu[is.na(mu)] <- 0
    A <- matrix((1 - 0.95) / (K - 1), K, K); diag(A) <- 0.95; pi0 <- rep(1 / K, K); ll_old <- -Inf
    for (it in seq_len(iters)) {
      r <- hmm_forward_backward(hmm_emission(Y, Vt, mu, V), A, pi0); g <- r$gamma
      for (k in 1:K) for (d in 1:2) { w <- 1 / (V[d] + Vt[, d]); den <- sum(g[, k] * w); if (den > 1e-8) mu[k, d] <- sum(g[, k] * w * Y[, d]) / den }
      cnt <- r$xi + kappa * diag(K) + 0.5 * (1 - diag(K)); A <- cnt / rowSums(cnt); pi0 <- (g[1, ] + 1e-3) / sum(g[1, ] + 1e-3)
      if (is.finite(r$ll) && abs(r$ll - ll_old) < tol) break; ll_old <- r$ll
    }
    lb <- hmm_emission(Y, Vt, mu, V); r <- hmm_forward_backward(lb, A, pi0)
    if (is.null(best) || r$ll > best$ll) best <- list(ll = r$ll, npar = 2 * K + K * (K - 1) + (K - 1), K = K, mu = mu, A = A, gamma = r$gamma, path = hmm_viterbi(lb, A, pi0))
  }
  best
}
