# Summaries of v2/sim.R and v2/real.R (run from the project root:  Rscript experiments/regimes/v2/summarize.R [sim|real])
what <- commandArgs(TRUE); what <- if (length(what)) what[1] else "both"
if (what %in% c("sim", "both")) local({
suppressMessages(library(data.table)); options(width = 230); r <- readRDS("experiments/regimes/v2/_out/sim.rds"); GATE <- 9
w <- r[model %in% c("M3", "M3O")]; g0 <- r[model %in% c("HMM", "HMMO")]; g0[, fallback := ifelse(model == "HMM", "M3", "M3O")]
g0 <- merge(g0, w[, .(fallback = model, mode, scenario, sim, f_rp = rmse_p, f_rs = rmse_s, f_cp = cover_p, f_cs = cover_s)], by = c("fallback", "mode", "scenario", "sim"), all.x = TRUE)
mk <- function(T) { g <- copy(g0); g[, open := is.finite(gain) & gain >= GATE & is.finite(dll) & dll >= T]; g[, `:=`(rmse_p = ifelse(open, rmse_p, f_rp), rmse_s = ifelse(open, rmse_s, f_rs), cover_p = ifelse(open, cover_p, f_cp), cover_s = ifelse(open, cover_s, f_cs), model = paste0(model, "g", T))]; g }
a <- rbind(r[, .(scenario, sim, model, mode, rmse_p, rmse_s, cover_p, cover_s)], rbindlist(lapply(c(0, 10, 20), function(T) mk(T)[, .(scenario, sim, model, mode, rmse_p, rmse_s, cover_p, cover_s, open)]), fill = TRUE), fill = TRUE)
cat("\nGate opening (%) by scenario (change-point gain >= 9 and, second stage, fitted log-likelihood gain >= T):\n"); print(dcast(rbindlist(lapply(c(0, 10, 20), function(T) mk(T)[model %like% "^HMMg" & mode == "filter", .(T = T, open = round(100 * mean(open))), by = scenario])), scenario ~ T, value.var = "open"))
ref <- a[model == "M3", .(scenario, sim, mode, r0p = rmse_p, r0s = rmse_s)]; a <- merge(a, ref, by = c("scenario", "sim", "mode"))
cat("\nRMSE relative to M3 of the same mode (%, negative better; mean over simulations), median ratio, and 95% coverage\n")
print(a[mode == "filter", .(n = .N, vagal = round(100 * (mean(rmse_p) / mean(r0p) - 1)), symp = round(100 * (mean(rmse_s) / mean(r0s) - 1)), med_v = round(median(rmse_p / r0p), 2), med_s = round(median(rmse_s / r0s), 2), worse_v = round(100 * mean(rmse_p > r0p)), worst_v = round(max(rmse_p / r0p), 1), cover_v = round(mean(cover_p), 2), cover_s = round(mean(cover_s), 2)), by = .(scenario, model)][order(scenario, match(model, c("M3", "M3O", "HMM", "HMMg0", "HMMg10", "HMMg20", "HMMO", "HMMOg0", "HMMOg10", "HMMOg20")))])
cat("\nSmoothed (RTS) estimates, same scale:\n"); print(a[mode == "smooth" & model %in% c("M3", "M3O", "HMMg0", "HMMOg0"), .(vagal = round(100 * (mean(rmse_p) / mean(r0p) - 1)), symp = round(100 * (mean(rmse_s) / mean(r0s) - 1))), by = .(scenario, model)][order(scenario, model)])
})
if (what %in% c("real", "both")) local({
suppressMessages(library(data.table)); options(width = 230); x <- readRDS("experiments/regimes/v2/_out/real.rds"); res <- x$res; held <- x$held; GATE <- x$gate
res[, open := is.finite(gain) & gain >= GATE]
cat("\nGate opening (%):\n"); print(unique(res[, .(ds, id, open)])[, .(n = .N, open_pct = round(100 * mean(open))), by = ds])
# gated models: where the gate is closed the regime models are M3 (HMM*) / M3O (HMMOr)
fb <- res[model %in% c("M3", "M3O"), .(ds, id, fb = model, ks_f = ks_ok, pacf_f = pacf, dLL_f = dLL, tau_p_f = tau_p, t50v_f = t50_vagal, t50s_f = t50_symp, dv_f = d_vagal, dsy_f = d_symp)]
g <- res[model %in% c("HMM0", "HMMr", "HMMOr")]; g[, fb := ifelse(model == "HMMOr", "M3O", "M3")]; g <- merge(g, fb, by = c("ds", "id", "fb"), all.x = TRUE)
g[!(open), `:=`(ks_ok = ks_f, pacf = pacf_f, dLL = dLL_f, tau_p = tau_p_f, t50_vagal = t50v_f, t50_symp = t50s_f, d_vagal = dv_f, d_symp = dsy_f)]; g[, model := paste0(model, "g")]
a <- rbind(res[, .(ds, id, model, ks_ok, pacf, dLL, tau_p, t50_vagal, t50_symp, d_vagal, d_symp, hits, n_annot, false_alarms, n_onsets)], g[, .(ds, id, model, ks_ok, pacf, dLL, tau_p, t50_vagal, t50_symp, d_vagal, d_symp, hits, n_annot, false_alarms, n_onsets)], fill = TRUE)
cat("\nIn-sample (gated regime models; dLL = explicit log-likelihood minus M3's, median):\n")
print(a[, .(n = .N, ks_pass = round(100 * mean(ks_ok)), pacf = round(mean(pacf), 1), dLL = round(median(dLL, na.rm = TRUE), 1), tau_p = round(median(tau_p, na.rm = TRUE), 1)), by = .(ds, model)][order(ds, match(model, c("M3", "M3O", "HMM0", "HMM0g", "HMMr", "HMMrg", "HMMOr", "HMMOrg")))])
cat("\nExercise realism:\n"); print(a[ds == "exercise" & !is.na(t50_vagal) & !is.na(t50_symp), .(n = .N, t50_vagal = round(median(t50_vagal), 1), t50_symp = round(median(t50_symp), 1), vagal_faster = round(100 * mean(t50_vagal < t50_symp)), expected_signs = round(100 * mean(d_vagal < 0 & d_symp > 0))), by = model][order(match(model, c("M3", "M3O", "HMM0", "HMM0g", "HMMr", "HMMrg", "HMMOr", "HMMOrg")))])
cat("\nPRCP (10 recordings): annotated onsets found within 30 s and detections without an annotated change:\n"); print(a[ds == "prcp", .(onsets_found = sum(hits, na.rm = TRUE), of = sum(n_annot, na.rm = TRUE), detections = sum(n_onsets), no_annotation = sum(false_alarms, na.rm = TRUE)), by = model][order(match(model, c("M3", "M3O", "HMM0", "HMM0g", "HMMr", "HMMrg", "HMMOr", "HMMOrg")))])
held[, `:=`(ok_tail = tail_bad < 0.1, open = is.finite(gain_head) & gain_head >= GATE)]; hl <- melt(held, id.vars = c("ds", "id", "m3", "m3o", "gain_head", "tail_bad", "ok_tail", "open"), measure.vars = c("hmm0", "hmmr", "hmmor"), variable.name = "model", value.name = "ll")
hl[, `:=`(vs_m3 = ll - m3, vs_m3o = ll - m3o)]; hl[model == "hmmor", vs_m3 := ll - m3]
cat("\nHeld-out forward score, nats/beat; gate on the head fit (closed gate = same as M3 / M3O). Recordings with a degenerate tail (>10% beats of 2 s or more) excluded. vs M3 (all three) and vs M3O (HMMOr)\n")
hl[, d := ifelse(open, vs_m3, 0)]; hl[, d_o := ifelse(open, vs_m3o, 0)]
print(hl[ok_tail == TRUE, .(n = .N, open_pct = round(100 * mean(open)), mean = round(mean(d, na.rm = TRUE), 3), median_open = round(median(vs_m3[open], na.rm = TRUE), 3), better_of_open_pct = round(100 * mean(vs_m3[open] > 0, na.rm = TRUE)), worst = round(min(d, na.rm = TRUE), 2), mean_vs_M3O = round(mean(d_o, na.rm = TRUE), 3)), by = .(ds, model)][order(ds, model)])
cat("\nUngated (all recordings with a valid tail):\n"); print(hl[ok_tail == TRUE, .(n = .N, mean = round(mean(vs_m3, na.rm = TRUE), 3), median = round(median(vs_m3, na.rm = TRUE), 3), better_pct = round(100 * mean(vs_m3 > 0, na.rm = TRUE)), worst = round(min(vs_m3, na.rm = TRUE), 2)), by = .(ds, model)][order(ds, model)])
cat("\nM3O minus M3 held-out (oscillator alone):\n"); print(held[ok_tail == TRUE, .(n = .N, mean = round(mean(m3o - m3, na.rm = TRUE), 3), median = round(median(m3o - m3, na.rm = TRUE), 3), better_pct = round(100 * mean(m3o > m3, na.rm = TRUE))), by = ds])
})
