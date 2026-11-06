# EXPERIMENTAL. Idea 6 on the PRCP recordings: tilts and stand-ups from the annotation file as known rectangular inputs (onset to return) to the two drives, with one
# amplitude per event and per pathway estimated jointly with the OU residuals (GLS), against M3 without inputs. Run from the project root:  Rscript experiments/m3plus/real_inputs.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
d <- load_dataset("prcp", cfg); ev <- fread("data/prcp_events.csv"); rows <- list()
for (rid in names(d$rr)) {
  id <- rid; rec <- d$rec[[rid]]; dy <- d$rr[[rid]]; off <- rec$time[1] - rec$rr[1]; e <- ev[as.character(ev$id) == rid & is.finite(onset_s)]
  e[, ret := ifelse(is.finite(return_end_s), return_end_s, return_s)]; e <- e[is.finite(ret)]; e[, `:=`(t0 = onset_s - off, t1 = ret - off)]; e <- e[t0 > 0 & t1 < sum(dy)]
  series <- prepare_series(dy); b <- fit_ext(series); f <- fit_ext(series, list(inputs = list(events = as.data.frame(e[, .(t0, t1)]))))
  rb <- b$run(full = TRUE); ri <- f$run(full = TRUE); m <- nrow(e); bp <- ri$beta[seq(2, 2 * m + 1, by = 2)]; bs <- ri$beta[seq(3, 2 * m + 2, by = 2)]
  ok_sign <- !grepl("stand", e$kind)                                                                 # tilts: vagal withdrawal (b_p < 0) and sympathetic rise (b_s > 0)
  kd <- function(r, fit) { k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < fit$ctl$diag_gate_max); z <- r$innovations_diagnostic[k]; c(ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(length(z)), pacf = sum(abs(as.numeric(pacf(z, lag.max = 40, plot = FALSE)$acf)) > 1.96 / sqrt(length(z)))) }
  rows[[id]] <- data.table(id = id, n_events = m, n_tilts = sum(ok_sign), dLL = ri$LL_star - rb$LL_star, extra_par = f$n_par - b$n_par + 2 * m, aic_gain = 2 * (ri$LL_star - rb$LL_star) - 2 * 2 * m,
                           tau_p_base = b$tau_p, tau_p_inputs = f$tau_p, tau_s_base = b$tau_s, tau_s_inputs = f$tau_s, ks_base = kd(rb, b)[["ks_ok"]], ks_inputs = kd(ri, f)[["ks_ok"]], pacf_base = kd(rb, b)[["pacf"]], pacf_inputs = kd(ri, f)[["pacf"]],
                           tilts_expected_signs = sum(bp[ok_sign] < 0 & bs[ok_sign] > 0), tilts_vagal_neg = sum(bp[ok_sign] < 0), tilts_symp_pos = sum(bs[ok_sign] > 0))
}
r <- rbindlist(rows); dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(r, "experiments/m3plus/_out/real_inputs.rds")
print(r[, .(id, n_events, n_tilts, dLL = round(dLL, 1), aic_gain = round(aic_gain, 1), tau_p_base = round(tau_p_base, 1), tau_p_inputs = round(tau_p_inputs, 1), tau_s_base = round(tau_s_base), tau_s_inputs = round(tau_s_inputs), ks_base, ks_inputs, pacf_base, pacf_inputs, tilts_expected_signs)])
cat(sprintf("\nRecordings improved by AIC: %d of %d; median logLik gain %.1f; tilts with both expected signs (b_p<0, b_s>0): %d of %d (vagal<0: %d; sympathetic>0: %d)\n", sum(r$aic_gain > 0), nrow(r), median(r$dLL), sum(r$tilts_expected_signs), sum(r$n_tilts), sum(r$tilts_vagal_neg), sum(r$tilts_symp_pos)))
cat(sprintf("tau_p: median %.1f s without inputs, %.1f s with inputs; tau_s: %.0f vs %.0f s; KS pass %d vs %d of %d; PACF violations %.1f vs %.1f\n", median(r$tau_p_base), median(r$tau_p_inputs), median(r$tau_s_base), median(r$tau_s_inputs), sum(r$ks_base), sum(r$ks_inputs), nrow(r), mean(r$pacf_base), mean(r$pacf_inputs)))
