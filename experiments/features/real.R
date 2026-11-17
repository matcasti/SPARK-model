# ==============================================================================
# experiments/features/real.R -- EXPERIMENTAL. R1 (information beyond the spectrum: IAAFT surrogates) and R2 (blockade test of sim-trained anchors).
#   Rscript experiments/features/real.R           (from the project root, after sim.R; needs results/fit_*.rds and data/blockade/)
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/features/lib.R")
n_cores <- max(1L, min(6L, cfg$cores - 1L)); out <- "experiments/features/_out"; N_SURR <- 19L

# ---- R1: real value against 19 IAAFT surrogates; rejection = real value is the lowest or the highest of the 20 (two-sided, alpha 0.10) --------
surrogate_test <- function(rr, seed) {
  v0 <- nonspectral(rr); S <- with_seed(seed, t(vapply(seq_len(N_SURR), function(k) nonspectral(iaaft(rr)), numeric(length(v0)))))
  data.table(feature = names(v0), real = v0, surr_mean = colMeans(S, na.rm = TRUE), surr_sd = apply(S, 2, sd, na.rm = TRUE),
             reject = vapply(seq_along(v0), function(j) { x <- S[, j]; if (!is.finite(v0[j]) || sum(is.finite(x)) < N_SURR) NA else v0[j] < min(x) || v0[j] > max(x) }, logical(1)))
}
rec <- rbindlist(lapply(c("fantasia", "prcp", "exercise"), function(ds) { st <- readRDS(sprintf("results/fit_%s.rds", ds))$states
  rbindlist(lapply(unique(st$SubjectID), function(id) data.table(dataset = ds, id = id, rr = list(st[SubjectID == id, RR_Interval])))) }))
ctrl <- rbindlist(lapply(1:100, function(i) { idx <- if (i <= 50) i else 450L + i; sd <- id_seed(paste("features", idx), cfg$seed); p <- draw_truth(sd, idx > 500)   # the first 50 of each generator in sim.R
  data.table(dataset = if (i > 50) "sim_oscillator" else "sim_model", id = as.character(i), rr = list(simulate_one(p, sd)$rr)) }))
jobs <- c(lapply(seq_len(nrow(rec)), function(j) rec[j]), lapply(seq_len(nrow(ctrl)), function(j) ctrl[j]))
log_msg("Surrogate test: %d real + %d simulated recordings x %d surrogates", nrow(rec), nrow(ctrl), N_SURR)
R1 <- rbindlist(drop_failed(map_recordings(jobs, function(r) surrogate_test(r$rr[[1]], id_seed(paste("surr", r$dataset, r$id), cfg$seed))[, `:=`(dataset = r$dataset, id = r$id)], n_cores), "surrogate"))
saveRDS(R1, file.path(out, "real_surrogates.rds"))
tab <- dcast(R1[, .(share = round(mean(reject, na.rm = TRUE), 2)), by = .(dataset, feature)], feature ~ dataset, value.var = "share")
cat("\n== R1: share of recordings whose feature lies outside its 19 IAAFT surrogates (expected 0.10 if the series is linear Gaussian) ==\n"); print(tab, row.names = FALSE)
cat("\nAny feature rejected (share of recordings):\n"); print(R1[, .(any = any(reject, na.rm = TRUE)), by = .(dataset, id)][, .(share_any = round(mean(any), 2)), by = dataset], row.names = FALSE)

# ---- R2: blockade segments, anchors trained on all simulations -----------------------------------------------------------------
D <- readRDS(file.path(out, "sim_features.rds")); y <- log(D$lambda_eff)
models <- with_seed(cfg$seed, list(Sf = train_full(D, y, SETS$Sf), SfEGfF = train_full(D, y, unlist(SETS[c("Sf", "E", "Gf", "F")]))))
protocol <- fread(file.path(cfg$paths$blockade, "protocol.csv")); st <- cfg$blockade
seg_job <- function(j) {
  g <- protocol[j]; tk <- blockade_peaks(cfg, g$subject, g$day); k <- tk >= g$t0 + st$trim_s & tk <= g$t1; rr <- diff(tk[k])
  if (length(rr) < st$min_beats || mean(rr > 0.25 & rr < 2.5) < 0.97) return(NULL)
  rr <- rr[rr > 0.25 & rr < 2.5]; ser <- prepare_series(rr); fe <- add_scale_free(as.data.table(t(all_features(rr, ser$priors))))
  anchors <- c(wavelet = ser$priors$energy_ratio, Sf = exp(models$Sf(fe)), SfEGfF = exp(models$SfEGfF(fe)))
  rbindlist(lapply(names(anchors), function(a) { fit <- fit_with_anchor(ser, anchors[[a]]); if (!isTRUE(fit$ok)) return(NULL); lp <- fit$par$lp
    data.table(anchor = a, subject = g$subject, day = g$day, cond = g$cond, breathing = g$breathing, n = length(rr), hr = 60 / mean(rr), sd_rr_ms = 1000 * sd(rr),
               tau_p = fit$tau_p, tau_s = fit$tau_s, lambda_P = lp, SD_p = sqrt(lp * fit$sig2), SD_s = sqrt(fit$sig2)) }))
}
log_msg("Blockade: %d segments x 3 anchors", nrow(protocol))
B <- rbindlist(drop_failed(map_recordings(seq_len(nrow(protocol)), seg_job, n_cores), "segment"))
saveRDS(B, file.path(out, "real_blockade.rds"))
cat("\n== R2: anchor lambda_P by condition (median over segments, spontaneous breathing) ==\n")
print(dcast(B[breathing == "spontaneous", .(lambda = round(median(lambda_P), 2)), by = .(anchor, cond)], cond ~ anchor, value.var = "lambda"), row.names = FALSE)
cat("\n== R2: ratio to saline of the same day (geometric mean [bootstrap 95% CI]); criterion: atropine SD_p < 0.3 and SD_s towards 1; propranolol SD_p not raised ==\n")
for (a in unique(B$anchor)) { s <- blockade_ratio_summary(blockade_ratios(B[anchor == a]), st)
  cat("\n--", a, "--\n"); print(s[, .(breathing, cond, n, SD_p = sprintf("%.2f [%.2f, %.2f]", SD_p, SD_p_lo, SD_p_hi), SD_s = sprintf("%.2f [%.2f, %.2f]", SD_s, SD_s_lo, SD_s_hi))], row.names = FALSE) }
