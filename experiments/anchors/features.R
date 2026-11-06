# EXPERIMENTAL. Does the spectrum of a simulated series carry information about the true lambda_P? Features of each simulated recording vs the truth.
# Run from the project root:  Rscript experiments/anchors/features.R
source("analysis/00_setup.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
one <- function(job) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  pri <- extract_spectral_priors(sim$rr, spark_priors(), verbose = FALSE)
  data.table(scenario = job$scenario, sim = job$sim, lp_true = truth$lambda_p, tau_p = truth$tau_p, tau_s = truth$tau_s, anchor = pri$energy_ratio, p_lf = pri$p_lf, p_hf = pri$p_hf,
             sd_rr = sd(sim$rr), mean_rr = mean(sim$rr), n = length(sim$rr), share_true = var(sim$x_p) / (var(sim$x_p) + var(sim$x_s)))
}
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
d <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"))
dir.create("experiments/anchors/_out", showWarnings = FALSE); saveRDS(d, "experiments/anchors/_out/features.rds")
d[, `:=`(y = log(lp_true), lr = log(p_lf / p_hf), la = log(anchor), lt = log(p_lf + p_hf))]
cat("Spearman correlation of the TRUE log lambda_P with:\n"); print(round(sapply(c(anchor = "la", LF_over_HF = "lr", total_band_power = "lt"), function(v) cor(d$y, d[[v]], method = "spearman")), 2))
cat("\nWithin scenario (anchor is ~constant: how much does it vary?):\n"); print(d[, .(anchor_med = round(median(anchor), 2), anchor_sd_log = round(sd(log(anchor)), 2), true_sd_log = round(sd(y), 2), cor_LFHF_true = round(cor(lr, y, method = "spearman"), 2)), by = scenario])
cv <- function(form) { pred <- numeric(nrow(d)); for (sc in unique(d$scenario)) { m <- lm(form, d[scenario != sc]); pred[d$scenario == sc] <- predict(m, d[scenario == sc]) }; c(abs_err_log = median(abs(pred - d$y)), r2 = 1 - sum((pred - d$y)^2) / sum((d$y - mean(d$y))^2)) }
cat("\nLeave-one-scenario-out prediction of ln lambda_P (median absolute error in log units; R2):\n")
print(round(rbind(constant = c(median(abs(d$y - median(d$y))), 0), anchor_as_is = c(median(abs(d$y - d$la)), NA), calib_anchor = cv(y ~ la), calib_LFHF = cv(y ~ lr), calib_LFHF_total = cv(y ~ lr + lt), calib_all = cv(y ~ lr + lt + la)), 2))
