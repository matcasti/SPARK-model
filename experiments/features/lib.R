# ==============================================================================
# experiments/features/lib.R -- EXPERIMENTAL. Spectral, entropy, geometric and fractal features of an RR series, IAAFT surrogates, the simulation
# generator with randomized scale, M3 with a given anchor and the cross-validated elastic-net anchor. Source after analysis/00_setup.R.
#   features_S / _E / _G / _F()   one feature set (named numeric vector); all_features() binds them
#   iaaft()                        surrogate with the same amplitude distribution and periodogram (Schreiber & Schmitz 1996)
#   draw_truth(), simulate_one()   parameters drawn from the design ranges (README) and the simulated recording
#   fit_with_anchor()              M3 with lambda_P fixed at a given value
#   cv_predict()                   out-of-fold elastic-net prediction from a feature set
# ==============================================================================
suppressMessages(library(glmnet))
Rcpp::sourceCpp("experiments/features/features.cpp")

SETS <- list(S = "^S_", E = "^E_", G = "^G_", F = "^F_", Sf = "^Sf_", Gf = "^Gf_")

#' Scale-free versions of S and G (in place): ratios and relative powers only. Entropies (r = 0.2 SD) and DFA exponents are scale-free already.
add_scale_free <- function(D) {
  D[, `:=`(Sf_anchor = S_anchor, Sf_wratio = S_wlf - S_whf, Sf_lfhf = S_lfhf, Sf_lfrel = S_lf - log(exp(S_lf) + exp(S_hf) + exp(S_vlf)),
           Sf_hfrel = S_hf - log(exp(S_lf) + exp(S_hf) + exp(S_vlf)), Sf_vlfrel = S_vlf - log(exp(S_lf) + exp(S_hf) + exp(S_vlf)), Sf_rmssd_sdnn = S_rmssd - S_sdnn,
           Gf_ratio = G_ratio, Gf_csi = G_csi, Gf_guzik = G_guzik, Gf_porta = G_porta)]
  invisible(D)
}

features_S <- function(rr, priors) {
  bp <- function(b) band_power(rr, b)
  lf <- bp(c(0.04, 0.15)); hf <- bp(c(0.15, 0.40)); vlf <- bp(c(0.01, 0.04))
  c(S_anchor = log(priors$energy_ratio), S_wlf = log(priors$p_lf), S_whf = log(priors$p_hf), S_lf = log(lf), S_hf = log(hf), S_vlf = log(vlf),
    S_lfhf = log(lf / hf), S_total = log(var(rr)), S_meanrr = mean(rr), S_sdnn = log(sd(rr)), S_rmssd = log(sqrt(mean(diff(rr)^2))))
}

coarse <- function(x, s) { n <- floor(length(x) / s); if (s == 1) x else colMeans(matrix(x[seq_len(n * s)], nrow = s)) }

features_E <- function(rr, scales = 1:6) {
  r <- 0.2 * sd(rr)
  mse <- vapply(scales, function(s) sampen_cpp(coarse(rr, s), 2L, r), numeric(1))
  ok <- is.finite(mse)
  c(E_sampen = mse[1], E_pe3 = permen_cpp(rr, 3L, 1L), E_pe4 = permen_cpp(rr, 4L, 1L), E_pe3d2 = permen_cpp(rr, 3L, 2L), E_pe4d2 = permen_cpp(rr, 4L, 2L),
    setNames(mse[-1], paste0("E_mse", scales[-1])), E_mse_area = if (all(ok)) sum(mse) else NA_real_,
    E_mse_slope = if (sum(ok) >= 3) unname(coef(lm(mse[ok] ~ scales[ok]))[2]) else NA_real_)
}

features_G <- function(rr) {
  d <- diff(rr); sd1 <- sd(d) / sqrt(2); sd2 <- sqrt(max(2 * var(rr) - sd1^2, 1e-12)); nz <- d != 0
  c(G_sd1 = log(sd1), G_sd2 = log(sd2), G_ratio = log(sd1 / sd2), G_csi = sd2 / sd1, G_cvi = log10(1e6 * sd1 * sd2 * 16),
    G_guzik = sum(d[d > 0]^2) / sum(d[nz]^2), G_porta = mean(d[nz] < 0))
}

dfa_exponent <- function(x, sizes) {
  y <- cumsum(x - mean(x)); N <- length(y)
  Fn <- vapply(sizes, function(n) { nb <- floor(N / n); if (nb < 4) return(NA_real_)
    Y <- matrix(y[seq_len(nb * n)], nrow = n); tc <- seq_len(n) - (n + 1) / 2
    Yc <- sweep(Y, 2, colMeans(Y)); b <- colSums(tc * Yc) / sum(tc^2); sqrt(mean((Yc - outer(tc, b))^2)) }, numeric(1))
  ok <- is.finite(Fn) & Fn > 0; if (sum(ok) < 3) NA_real_ else unname(coef(lm(log(Fn[ok]) ~ log(sizes[ok])))[2])
}

features_F <- function(rr) c(F_alpha1 = dfa_exponent(rr, 4:16), F_alpha2 = dfa_exponent(rr, unique(round(exp(seq(log(16), log(64), length.out = 12))))))

all_features <- function(rr, priors) c(features_S(rr, priors), features_E(rr), features_G(rr), features_F(rr))
nonspectral <- function(rr) c(features_E(rr), features_G(rr), features_F(rr))

#' IAAFT surrogate of x (iterated amplitude-adjusted Fourier transform).
iaaft <- function(x, iters = 100) {
  n <- length(x); sorted <- sort(x); amp <- Mod(fft(x)); s <- sample(x)
  for (i in seq_len(iters)) {
    ph <- fft(s); s <- Re(fft(amp * exp(1i * Arg(ph)), inverse = TRUE)) / n
    s_new <- sorted[rank(s, ties.method = "first")]; if (identical(s_new, s)) break; s <- s_new
  }
  s
}

# ---- simulations ------------------------------------------------------------------------------------------------------------
draw_truth <- function(seed, osc) with_seed(seed, {
  p <- list(tau_p = runif(1, 2.5, 10), tau_s = exp(runif(1, log(60), log(400))), lambda_p = exp(runif(1, log(0.25), log(12))), lR = exp(runif(1, log(0.005), log(0.08))),
            sig2 = exp(runif(1, log(1e-4), log(1e-2))), nu0 = runif(1, 0.8, 1.8), osc = osc)
  if (osc) p$o <- list(alpha = 0.18, f_hz = runif(1, 0.15, 0.40), lo = exp(runif(1, log(0.1), log(1.5))))
  p })

#' Simulated recording. With the oscillator, the true vagal-like state is x_p + z1 (the oscillation enters the rate with the vagal sign) and the
#' anchor M3 needs is the effective ratio lambda_eff = lambda_P + lo; without it lambda_eff = lambda_P.
simulate_one <- function(p, seed, duration_s = 420) {
  if (!p$osc) { s <- simulate_rr(duration_s, nu0 = p$nu0, kp = 1 / p$tau_p, ks = 1 / p$tau_s, sig2 = p$sig2, lp = p$lambda_p, lR = p$lR, h = 0.01, seed = seed); lam <- p$lambda_p }
  else { s <- simulate_rr_osc(duration_s, 1 / p$tau_p, 1 / p$tau_s, p$lambda_p, p$sig2, p$lR, nu0 = p$nu0, h = 0.02, seed = seed, osc = p$o); s$x_p <- s$x_p + s$z1; lam <- p$lambda_p + p$o$lo }
  list(rr = s$rr, x_p = s$x_p, x_s = s$x_s, lambda_eff = lam)
}

fit_with_anchor <- function(series, lp) { s <- series; s$priors$energy_ratio <- min(max(lp, 1e-3), 1e3); fit_variant(s, model_spec("M3")) }

# ---- learner ------------------------------------------------------------------------------------------------------------------
#' Design matrix of the columns matching `patterns`: features and their squares; NA imputed by the training median.
design <- function(D, patterns, train) {
  cols <- unique(unlist(lapply(patterns, function(p) grep(p, names(D), value = TRUE))))
  X <- as.matrix(D[, ..cols])
  for (j in seq_len(ncol(X))) { m <- median(X[train, j], na.rm = TRUE); X[!is.finite(X[, j]), j] <- if (is.finite(m)) m else 0 }
  keep <- apply(X[train, , drop = FALSE], 2, sd) > 0; X <- X[, keep, drop = FALSE]
  cbind(X, setNames(as.data.frame(X^2), paste0(colnames(X), "_sq")) |> as.matrix())
}

#' Out-of-fold predictions of y from the columns matching `patterns` (elastic net, alpha 0.5, lambda by inner 5-fold CV). folds: integer vector.
cv_predict <- function(D, y, patterns, folds, extra = NULL) {
  pred <- rep(NA_real_, nrow(D))
  for (k in sort(unique(folds))) {
    tr <- folds != k; X <- design(D, patterns, tr); if (!is.null(extra)) X <- cbind(X, extra)
    m <- cv.glmnet(X[tr, , drop = FALSE], y[tr], alpha = 0.5, nfolds = 5)
    pred[!tr] <- as.numeric(predict(m, X[!tr, , drop = FALSE], s = "lambda.min"))
  }
  pred
}

#' Model trained on all rows (for applying a sim-trained anchor to real data): returns a function features -> prediction.
train_full <- function(D, y, patterns) {
  tr <- rep(TRUE, nrow(D)); cols <- unique(unlist(lapply(patterns, function(p) grep(p, names(D), value = TRUE))))
  X <- design(D, patterns, tr); med <- vapply(cols, function(c) median(D[[c]], na.rm = TRUE), numeric(1)); kept <- colnames(X)
  m <- cv.glmnet(X, y, alpha = 0.5, nfolds = 5)
  function(newD) { Z <- as.matrix(newD[, ..cols]); for (j in seq_along(cols)) Z[!is.finite(Z[, j]), j] <- med[j]
    Z <- cbind(Z, setNames(as.data.frame(Z^2), paste0(colnames(Z), "_sq")) |> as.matrix())[, kept, drop = FALSE]
    as.numeric(predict(m, Z, s = "lambda.min")) }
}
