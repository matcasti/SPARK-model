# ==============================================================================
# experiments/blockade/validation/m3_blockade.R -- does the SPARK decomposition (production M3) respond to pharmacological blockade the way the physiology says?
#   Rscript experiments/blockade/validation/m3_blockade.R
# M3 (production fit, own spectral anchor, production priors) is fitted to every 5-min segment of the Lancaster data (7 men; spontaneous and paced 0.25 Hz). Quantities in physical
# units: SD_p = sqrt(lambda_P sigma^2) and SD_s = sqrt(sigma^2) (Hz), tau_p, tau_s, anchor lambda_P, mean heart rate, SD of the filtered states. Ratios to saline of the SAME day
# (day 1 for atropine and double blockade, day 2 for propranolol alone); geometric mean over subjects; bootstrap CI over subjects.
# Predictions fixed before the run: (P1) high-dose atropine lowers SD_p to < 0.3 of saline in at least 5 of 6 subjects (spontaneous); (P2) SD_s does not fall by more than
# 3x under atropine; (P3) propranolol alone does not lower SD_p (ratio >= 0.8 in most subjects); (P4) double blockade leaves both SDs < 0.5 of saline.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/blockade/lib.R")
set.seed(11); seg <- fread("experiments/blockade/_out/protocol.csv"); rows <- list()
for (i in seq_len(nrow(seg))) { r <- seg[i]; tk <- clean_peaks(scan(sprintf("experiments/blockade/_out/rpeaks/rpeaks_Subject%d_day%d.csv", r$subject, r$day), quiet = TRUE))
  k <- tk >= r$t0 + 20 & tk <= r$t1; rr <- diff(tk[k]); if (length(rr) < 150 || mean(rr > 0.25 & rr < 2.5) < 0.97) next; rr <- rr[rr > 0.25 & rr < 2.5]
  s <- prepare_series(rr); f <- fit_variant(s, model_spec("M3")); if (!isTRUE(f$ok)) next
  rows[[length(rows) + 1L]] <- data.table(subject = r$subject, day = r$day, cond = r$cond, breathing = r$breathing, n = length(rr), hr = 60 / mean(rr), tau_p = f$tau_p, tau_s = f$tau_s,
    lambda_P = s$priors$energy_ratio, SD_p = sqrt(s$priors$energy_ratio * f$sig2), SD_s = sqrt(f$sig2), sd_xp = sd(f$X_p), sd_xs = sd(f$X_s), nu0 = f$nu0, jitter_ms = f$rr_jitter_ms) }
m <- rbindlist(rows); saveRDS(m, "experiments/blockade/_out/m3_blockade.rds")
ref <- function(d, cn) d[cond == cn, .(subject, breathing, ref = cn)]
ratios <- rbindlist(lapply(c("atropine_low", "atropine_high", "double", "propranolol"), function(cn) { base <- if (cn == "propranolol") 2L else 1L
  a <- m[cond == cn]; b <- m[cond == "saline" & day == base, .(subject, breathing, SD_p0 = SD_p, SD_s0 = SD_s, sdxp0 = sd_xp, sdxs0 = sd_xs, hr0 = hr, tp0 = tau_p, ts0 = tau_s, lam0 = lambda_P)]
  x <- merge(a, b, by = c("subject", "breathing")); x[, .(subject, breathing, cond, SD_p = SD_p / SD_p0, SD_s = SD_s / SD_s0, sd_xp = sd_xp / sdxp0, sd_xs = sd_xs / sdxs0, d_hr = hr - hr0, tau_p = tau_p, tau_p0 = tp0, tau_s = tau_s, tau_s0 = ts0, lambda_P = lambda_P, lambda0 = lam0)] }))
saveRDS(ratios, "experiments/blockade/_out/m3_blockade_ratios.rds")
gm <- function(x) exp(mean(log(x))); bs <- function(x) { b <- replicate(2000, gm(sample(x, replace = TRUE))); sprintf("%.2f [%.2f, %.2f]", gm(x), quantile(b, .025), quantile(b, .975)) }
cat("\nRatio to saline of the same day (geometric mean over subjects [95% bootstrap CI]); n = subjects\n")
print(ratios[, .(n = .N, SD_p = bs(SD_p), SD_s = bs(SD_s), sd_xp = bs(sd_xp), sd_xs = bs(sd_xs), d_hr_bpm = round(median(d_hr))), by = .(breathing, cond)][order(breathing, match(cond, c("atropine_low", "atropine_high", "double", "propranolol")))], row.names = FALSE)
cat("\nPredictions (spontaneous breathing):\n")
sp <- ratios[breathing == "spontaneous"]
cat(sprintf("  P1 atropine: SD_p < 0.3 of saline in %d of %d subjects\n", sum(sp[cond == "atropine_high", SD_p] < 0.3), nrow(sp[cond == "atropine_high"])))
cat(sprintf("  P2 atropine: SD_s ratio range %.2f-%.2f; fell by more than 3x in %d subjects\n", min(sp[cond == "atropine_high", SD_s]), max(sp[cond == "atropine_high", SD_s]), sum(sp[cond == "atropine_high", SD_s] < 1 / 3)))
cat(sprintf("  P3 propranolol alone: SD_p ratio >= 0.8 in %d of %d subjects\n", sum(sp[cond == "propranolol", SD_p] >= 0.8), nrow(sp[cond == "propranolol"])))
cat(sprintf("  P4 double blockade: both SDs < 0.5 of saline in %d of %d subjects (SD_p < 0.5 in %d, SD_s < 0.5 in %d)\n", sum(sp[cond == "double", SD_p < 0.5 & SD_s < 0.5]), nrow(sp[cond == "double"]), sum(sp[cond == "double", SD_p] < 0.5), sum(sp[cond == "double", SD_s] < 0.5)))
cat("\nSame for paced breathing: SD_p atropine", paste(round(ratios[breathing == "paced_0.25" & cond == "atropine_high", SD_p], 2), collapse = " "), "| SD_s", paste(round(ratios[breathing == "paced_0.25" & cond == "atropine_high", SD_s], 2), collapse = " "), "\n")
cat("\nM3 parameters under saline (spontaneous, day 1) vs high-dose atropine, median: tau_p", round(median(m[cond == "saline" & day == 1 & breathing == "spontaneous", tau_p]), 1), "->", round(median(m[cond == "atropine_high" & breathing == "spontaneous", tau_p]), 1),
    " tau_s", round(median(m[cond == "saline" & day == 1 & breathing == "spontaneous", tau_s])), "->", round(median(m[cond == "atropine_high" & breathing == "spontaneous", tau_s])), " lambda_P anchor", round(median(m[cond == "saline" & day == 1 & breathing == "spontaneous", lambda_P]), 1), "->", round(median(m[cond == "atropine_high" & breathing == "spontaneous", lambda_P]), 2), "\n")
