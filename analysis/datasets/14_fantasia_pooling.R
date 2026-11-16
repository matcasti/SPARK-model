# ==============================================================================
# 14_fantasia_pooling.R -- are the time constants the priors? Fantasia windows A and B refitted with the production priors, with weak priors and with empirical-Bayes
# priors estimated on the other half of the subjects (two folds); test-retest reliability (ICC(3,1) of log values) and young-vs-older discrimination of each.
#
# Reads   CSV recordings (data/fantasia)
# Writes  results/pooling.rds   fits (one row per subject x window x setting), reliability table, hyperparameters per fold
# Method  R/steps/extended_validation.R (pooling_fit, pooling_hyperparameters, pooling_reliability); settings in cfg$pooling.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. Load --------------------------------------------------------------------------------------------------------------------
dataset <- load_dataset("fantasia", cfg); ids <- names(dataset$rec_full); w <- cfg$windows$fantasia
windows <- list(A = w$stability_A, B = w$stability_B)
rr_of <- function(id, win) window_rec(dataset$rec_full[[id]], start_s = windows[[win]][1], duration_s = diff(windows[[win]]))$rr
jobs <- CJ(id = ids, win = c("A", "B"))
run_setting <- function(sub, bounds_p, bounds_s, label) {
  res <- map_recordings(seq_len(nrow(sub)), function(i) { x <- pooling_fit(rr_of(sub$id[i], sub$win[i]), bounds_p, bounds_s); if (is.null(x)) NULL else cbind(id = sub$id[i], win = sub$win[i], x) }, cfg$cores)
  out <- rbindlist(drop_failed(res, "fit")); out[, setting := label]; out
}

# ---- 2. Production and weak priors --------------------------------------------------------------------------------------------------
production <- run_setting(jobs, c(3, 6), c(100, 300), "production")
weak <- run_setting(jobs, cfg$pooling$weak_p, cfg$pooling$weak_s, "weak")

# ---- 3. Empirical Bayes: hyperparameters from the other half of the subjects ----------------------------------------------------------
folds <- with_seed(cfg$pooling$seed, split(sample(ids), rep(1:2, length.out = length(ids))))
eb <- list(); hyper <- list()
for (k in 1:2) {
  h <- pooling_hyperparameters(weak, folds[[3 - k]])
  bp <- exp(h$p[["mu"]] + c(-1.96, 1.96) * h$p[["sd_b"]]); bs <- exp(h$s[["mu"]] + c(-1.96, 1.96) * h$s[["sd_b"]])
  hyper[[k]] <- data.table(fold = k, tau_p = exp(h$p[["mu"]]), tau_p_lo = bp[1], tau_p_hi = bp[2], tau_s = exp(h$s[["mu"]]), tau_s_lo = bs[1], tau_s_hi = bs[2])
  eb[[k]] <- run_setting(jobs[id %in% folds[[k]]], bp, bs, "empirical Bayes")
}

# ---- 4. Reliability and save ------------------------------------------------------------------------------------------------------------
fits <- rbind(production, weak, rbindlist(eb))
groups <- setNames(fantasia_group(ids, cfg), ids)
reliability <- pooling_reliability(fits, groups)
write_result(list(fits = fits, reliability = reliability, hyper = rbindlist(hyper), settings = cfg$pooling), "pooling", cfg)
print(reliability[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
