# ==============================================================================
# experiments/features/sim.R -- EXPERIMENTAL. Q1 and Q2 of README.md on simulated recordings with known truth.
#   Rscript experiments/features/sim.R            (from the project root; ~10 min on 6 cores)
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/features/lib.R")
n_cores <- max(1L, min(6L, cfg$cores - 1L)); N <- 1000L; out <- "experiments/features/_out"; dir.create(out, showWarnings = FALSE)

# ---- 1. simulate and compute the features --------------------------------------------------------------------------------
jobs <- lapply(seq_len(N), function(i) list(i = i, osc = i > N / 2, seed = id_seed(paste("features", i), cfg$seed)))
one <- function(job) {
  p <- draw_truth(job$seed, job$osc); s <- simulate_one(p, job$seed); ser <- prepare_series(s$rr)
  list(row = data.table(i = job$i, osc = job$osc, lambda_eff = s$lambda_eff, lambda_p = p$lambda_p, tau_p = p$tau_p, tau_s = p$tau_s, sig2 = p$sig2, nu0 = p$nu0,
                        share_true = var(s$x_p) / (var(s$x_p) + var(s$x_s)), n = length(s$rr), t(all_features(s$rr, ser$priors))),
       sim = s, series = ser, truth = p)
}
log_msg("Simulating %d recordings", N)
sims <- drop_failed(map_recordings(jobs, one, n_cores), "simulation")
D <- rbindlist(lapply(sims, `[[`, "row")); setkey(D, i)
add_scale_free(D)   # scale-free versions (lib.R): amplitude features leak lambda through the independent sigma^2 prior
saveRDS(D, file.path(out, "sim_features.rds"))

# ---- 2. Q1: out-of-fold prediction of ln lambda_eff ------------------------------------------------------------------------
folds <- with_seed(cfg$seed, { f <- integer(nrow(D)); for (g in c(FALSE, TRUE)) { k <- which(D$osc == g); f[k] <- sample(rep(1:5, length.out = length(k))) }; f })
y <- log(D$lambda_eff)
COMBOS <- list(S = "S", E = "E", G = "G", F = "F", EGF = c("E", "G", "F"), SE = c("S", "E"), SG = c("S", "G"), SF = c("S", "F"), SEGF = c("S", "E", "G", "F"),
               Sf = "Sf", Gf = "Gf", EGfF = c("E", "Gf", "F"), SfE = c("Sf", "E"), SfGf = c("Sf", "Gf"), SfF = c("Sf", "F"), SfEGfF = c("Sf", "E", "Gf", "F"))
pred <- list(const = vapply(seq_along(y), function(j) median(y[folds != folds[j]]), numeric(1)), wavelet = D$S_anchor)
for (nm in names(COMBOS)) pred[[nm]] <- with_seed(cfg$seed, cv_predict(D, y, unlist(SETS[COMBOS[[nm]]]), folds))
score <- function(p, k) c(med_abs_ln_err = median(abs(p[k] - y[k])), R2 = 1 - sum((p[k] - y[k])^2) / sum((y[k] - mean(y[k]))^2))
q1 <- rbindlist(lapply(names(pred), function(nm) rbindlist(lapply(list(all = TRUE, model = !D$osc, oscillator = D$osc), function(k) as.list(score(pred[[nm]], k))), idcol = "generator")[, set := nm]))
q1 <- dcast(q1, set ~ generator, value.var = c("med_abs_ln_err", "R2"))[order(med_abs_ln_err_all)]
cat("\n== Q1: out-of-fold error of ln lambda_eff (median |ln error|; R2) ==\n"); print(q1[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
saveRDS(list(pred = pred, folds = folds, q1 = q1), file.path(out, "sim_q1.rds"))

# ---- 3. Q2a: M3 refitted with each anchor ----------------------------------------------------------------------------------
ANCHORS <- c("wavelet", "const", "S", "SEGF", "Sf", "SfEGfF", "oracle")
anchor_value <- function(nm, j) switch(nm, oracle = D$lambda_eff[j], exp(pred[[nm]][j]))
fit_job <- function(z) {
  j <- match(z$row$i, D$i)
  rbindlist(lapply(ANCHORS, function(nm) { fit <- fit_with_anchor(z$series, anchor_value(nm, j)); if (!isTRUE(fit$ok)) return(NULL)
    truth <- list(tau_p = z$truth$tau_p, tau_s = z$truth$tau_s, lambda_p = z$sim$lambda_eff)
    recovery_row(fit, truth, z$sim, D$share_true[j])[, anchor := nm] }))[, `:=`(i = z$row$i, osc = z$row$osc)]
}
log_msg("Refitting M3 with %d anchors", length(ANCHORS))
F2 <- rbindlist(drop_failed(map_recordings(sims, fit_job, n_cores), "refit")); F2[, anchor := factor(anchor, ANCHORS)]
saveRDS(F2, file.path(out, "sim_q2a.rds"))
summ <- function(d) d[, .(n = .N, abs_ln_lp = median(abs(ln_err_lp)), abs_ln_tau_p = median(abs(ln_err_tau_p)), abs_ln_tau_s = median(abs(ln_err_tau_s)),
                          cor_xp = median(cor_x_p), cor_xs = median(cor_x_s), rmse_xp = median(rmse_x_p), rmse_xs = median(rmse_x_s),
                          cover_xp = mean(cover_x_p), cover_xs = mean(cover_x_s), abs_share = median(abs(share_err))), by = anchor][order(anchor)]
rnd <- function(d) d[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]
cat("\n== Q2a: M3 with each anchor, all simulations ==\n"); print(rnd(summ(F2[ok == TRUE])), row.names = FALSE)
cat("\n-- model generator --\n"); print(rnd(summ(F2[ok == TRUE & osc == FALSE])), row.names = FALSE)
cat("\n-- oscillator generator --\n"); print(rnd(summ(F2[ok == TRUE & osc == TRUE])), row.names = FALSE)

# ---- 4. Q2b: direct prediction of the time constants and the share, alone and stacked on M3 --------------------------------
w <-F2[anchor == "wavelet"][match(D$i, i)]
M3X <- cbind(m3_ln_tau_p = log(w$tau_p_true) + w$ln_err_tau_p, m3_ln_tau_s = log(w$tau_s_true) + w$ln_err_tau_s, m3_share = w$est_share)
M3X[!is.finite(M3X)] <- 0
targets <- list(ln_tau_p = log(D$tau_p), ln_tau_s = log(D$tau_s), share = D$share_true)
q2b <- rbindlist(lapply(names(targets), function(tg) { yy <- targets[[tg]]; mcol <- switch(tg, ln_tau_p = "m3_ln_tau_p", ln_tau_s = "m3_ln_tau_s", share = "m3_share")
  P <- list(M3_alone = M3X[, mcol],
            S = with_seed(cfg$seed, cv_predict(D, yy, SETS$S, folds)), SEGF = with_seed(cfg$seed, cv_predict(D, yy, unlist(SETS[c("S", "E", "G", "F")]), folds)),
            Sf = with_seed(cfg$seed, cv_predict(D, yy, SETS$Sf, folds)), SfEGfF = with_seed(cfg$seed, cv_predict(D, yy, unlist(SETS[c("Sf", "E", "Gf", "F")]), folds)),
            M3_plus_Sf = with_seed(cfg$seed, cv_predict(D, yy, SETS$Sf, folds, extra = M3X)), M3_plus_SfEGfF = with_seed(cfg$seed, cv_predict(D, yy, unlist(SETS[c("Sf", "E", "Gf", "F")]), folds, extra = M3X)),
            M3_plus_S = with_seed(cfg$seed, cv_predict(D, yy, SETS$S, folds, extra = M3X)), M3_plus_SEGF = with_seed(cfg$seed, cv_predict(D, yy, unlist(SETS[c("S", "E", "G", "F")]), folds, extra = M3X)))
  inp <- D$tau_p >= 3 & D$tau_p <= 6 & D$tau_s >= 100 & D$tau_s <= 300   # truth inside the production priors: no prior-mismatch advantage for the learners
  rbindlist(lapply(names(P), function(nm) data.table(target = tg, predictor = nm, med_abs_err = median(abs(P[[nm]] - yy)), R2 = 1 - sum((P[[nm]] - yy)^2) / sum((yy - mean(yy))^2),
                                                     med_abs_err_in_prior = median(abs(P[[nm]] - yy)[inp]), n_in_prior = sum(inp)))) }))
cat("\n== Q2b: out-of-fold prediction of the time constants and the variance share ==\n"); print(rnd(q2b), row.names = FALSE)
saveRDS(q2b, file.path(out, "sim_q2b.rds"))

# ---- 5. which features carry the signal (full-data elastic net on S+E+G+F, standardized coefficients) -----------------------
X <- design(D, unlist(SETS[c("Sf", "E", "Gf", "F")]), rep(TRUE, nrow(D))); m <- with_seed(cfg$seed, cv.glmnet(scale(X), y, alpha = 0.5)); cf <- as.matrix(coef(m, s = "lambda.1se"))[-1, 1]
cat("\n== Largest standardized coefficients for ln lambda_eff (scale-free Sf+E+Gf+F) ==\n"); print(round(head(cf[order(-abs(cf))], 15), 3))
cat("\nSpearman of ln lambda_eff with each single feature (top 15 by |rho|):\n"); fc <- grep("^(S|E|G|F|Sf|Gf)_", names(D), value = TRUE)
rho <- vapply(fc, function(c) suppressWarnings(cor(D[[c]], y, method = "spearman", use = "complete.obs")), numeric(1)); print(round(head(rho[order(-abs(rho))], 15), 3))
