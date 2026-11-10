# ==============================================================================
# experiments/pooling/run.R -- empirical-Bayes partial pooling of the time constants: does estimating the prior of tau_p and tau_s from the population (instead of fixing it at 3-6 s and
# 100-300 s) make the latent summaries more reliable? Fantasia, 40 subjects, two 10-min windows (A 600-1200 s, B 1200-1800 s; the test-retest design of the paper).
#   Rscript experiments/pooling/run.R
# Settings (M3 refitted in each; everything else as in production): production (fixed priors); weak (tau_p 0.5-60 s, tau_s 10-2000 s); empirical Bayes (population N(mu, sd_between) for
# ln tau_p and ln tau_s, estimated on the OTHER half of the subjects, two folds). Between-subject variance from the covariance of the two windows of the weak fits (test-retest).
# Reliability: ICC(3,1) and Spearman of the log values between windows; known groups: AUC of young vs older.
# ==============================================================================
source("analysis/00_setup.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); w <- cfg$windows$fantasia; dataset <- load_dataset("fantasia", cfg); ids <- names(dataset$rec_full); set.seed(5)
rr_of <- function(id, win) { r <- window_rec(dataset$rec_full[[id]], start_s = win[1], duration_s = win[2] - win[1]); r$rr }
fit_set <- function(rr, bounds_p, bounds_s) { pri <- spark_priors(tau_p_bounds = bounds_p, tau_s_bounds = bounds_s); ctl <- spark_control()
  f <- tryCatch(fit_model_fast(rr, pri, ctl), error = function(e) NULL); if (is.null(f)) return(NULL)
  p <- f$params; data.table(tau_p = 1 / p$kp, tau_s = 1 / p$ks, SD_p = sqrt(p$lp * p$sig2), SD_s = sqrt(p$sig2), lambda_P = p$lp) }
jobs <- CJ(id = ids, win = c("A", "B")); wins <- list(A = w$stability_A, B = w$stability_B)
run_setting <- function(bp, bs, label) { res <- map_recordings(seq_len(nrow(jobs)), function(i) { x <- fit_set(rr_of(jobs$id[i], wins[[jobs$win[i]]]), bp, bs); if (is.null(x)) NULL else cbind(id = jobs$id[i], win = jobs$win[i], x) }, n_cores)
  r <- rbindlist(drop_failed(res, "fit")); r[, setting := label]; r }
prod <- run_setting(c(3, 6), c(100, 300), "production"); weak <- run_setting(c(0.5, 60), c(10, 2000), "weak")
# hyperparameters of ln tau from the weak fits of a set of subjects: mean, and between-subject sd from the test-retest covariance
hyper <- function(sub) { a <- weak[win == "A" & id %in% sub][order(id)]; b <- weak[win == "B" & id %in% sub][order(id)]; ok <- a$id %in% b$id; a <- a[ok]; b <- b[b$id %in% a$id]
  h <- function(v) { la <- log(a[[v]]); lb <- log(b[[v]]); c(mu = mean(c(la, lb)), sd_b = sqrt(max(cov(la, lb), 0.02^2))) }; list(p = h("tau_p"), s = h("tau_s")) }
folds <- split(sample(ids), rep(1:2, length.out = length(ids))); eb <- list()
for (k in 1:2) { H <- hyper(folds[[3 - k]]); bp <- exp(H$p["mu"] + c(-1.96, 1.96) * H$p["sd_b"]); bs <- exp(H$s["mu"] + c(-1.96, 1.96) * H$s["sd_b"])
  cat(sprintf("fold %d (hyperparameters from the other half): tau_p %.1f s [%.1f, %.1f], tau_s %.0f s [%.0f, %.0f]\n", k, exp(H$p["mu"]), bp[1], bp[2], exp(H$s["mu"]), bs[1], bs[2]))
  sub <- jobs[id %in% folds[[k]]]; res <- map_recordings(seq_len(nrow(sub)), function(i) { x <- fit_set(rr_of(sub$id[i], wins[[sub$win[i]]]), bp, bs); if (is.null(x)) NULL else cbind(id = sub$id[i], win = sub$win[i], x) }, n_cores)
  eb[[k]] <- rbindlist(drop_failed(res, "fit")) }
ebf <- rbindlist(eb); ebf[, setting := "empirical Bayes"]; all <- rbind(prod, weak, ebf); saveRDS(all, "experiments/pooling/_out/fits.rds")
icc <- function(a, b) { x <- cbind(a, b); n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1); mse <- sum((x - rowMeans(x) - matrix(colMeans(x) - gm, n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
grp <- fantasia_group(ids, cfg); names(grp) <- ids
auc <- function(x, g) { r <- rank(x); n1 <- sum(g == "young"); n2 <- sum(g == "older"); (sum(r[g == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
out <- rbindlist(lapply(unique(all$setting), function(st) { d <- all[setting == st]; a <- d[win == "A"][order(id)]; b <- d[win == "B"][order(id)]; keep <- intersect(a$id, b$id); a <- a[id %in% keep]; b <- b[id %in% keep]
  rbindlist(lapply(c("tau_p", "tau_s", "SD_p", "SD_s"), function(v) { la <- log(a[[v]]); lb <- log(b[[v]]); data.table(setting = st, variable = v, n = length(keep), icc = icc(la, lb), spearman = cor(la, lb, method = "spearman"),
    auc_young_older = auc(a[[v]], grp[a$id]), sd_log = sd(c(la, lb)), prior_share = if (v == "tau_p") mean(a[[v]] >= 3 & a[[v]] <= 6) else if (v == "tau_s") mean(a[[v]] >= 100 & a[[v]] <= 300) else NA_real_) })) }))
options(width = 200); cat("\nTest-retest (windows A and B) and young-vs-older AUC of M3 refitted under each prior\n"); print(out[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)][order(variable, match(setting, c("production", "weak", "empirical Bayes")))], row.names = FALSE)
