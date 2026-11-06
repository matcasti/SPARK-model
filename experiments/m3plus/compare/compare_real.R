# EXPERIMENTAL. The same models on the 313 real recordings: fit (log-likelihood, AIC-type penalty, KS, PACF, time constants), agreement of the states with M3 and what the
# exercise and tilt summaries look like under each. M3E uses events estimated from the RR series by the pipeline's detector (exercise: one bout; PRCP: candidate tilts);
# Fantasia (rest) has no events, so M3E = M3 there. Run from the project root:  Rscript experiments/m3plus/compare/compare_real.R [quick]
arg <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else ""; quick <- arg == "quick"; ds_only <- if (arg %in% c("exercise", "prcp", "fantasia")) arg else NULL
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
jobs <- list()
for (ds in if (is.null(ds_only)) cfg$datasets else ds_only) { d <- load_dataset(ds, cfg); ids <- names(d$rr)
  if (quick) { set.seed(5); ids <- sort(sample(ids, min(length(ids), c(prcp = 3, fantasia = 4, exercise = 8)[[ds]]))) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
diag_of <- function(r, ctl) { k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < ctl$diag_gate_max); z <- r$innovations_diagnostic[k]; n <- length(z)
  c(ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf = sum(abs(as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)) > qnorm(0.975) / sqrt(n))) }
exercise_z <- function(e, dy, ev) {                                                         # the pipeline's rest-to-end-of-exercise response, on the states of one model
  tt <- cumsum(dy); st <- data.table(Time = tt, X_p = e$xp, X_s = e$xs, RR_Interval = dy)
  tryCatch(exercise_rest_response(st, "x", ev$t0[1], ev$t1[1], cfg)[, .(z_p = z_d_x_p, z_s = z_d_x_s, rest_sd_p = rest_sd_x_p, rest_sd_s = rest_sd_x_s)], error = function(e) NULL) }
tilt_pattern <- function(e, dy, ev) {                                                       # per candidate tilt: baseline 120 s to 5 s before the onset vs 30 s after the onset to the return
  tt <- cumsum(dy); out <- list()
  for (i in seq_len(nrow(ev))) { b <- tt >= ev$t0[i] - 120 & tt < ev$t0[i] - 5; a <- tt >= ev$t0[i] + 30 & tt <= ev$t1[i]; if (sum(b) < 20 || sum(a) < 5) next
    zp <- (mean(e$xp[a]) - mean(e$xp[b])) / sd(e$xp[b]); zs <- (mean(e$xs[a]) - mean(e$xs[b])) / sd(e$xs[b]); out[[i]] <- data.table(z_p = zp, z_s = zs, expected = zp < 0 & zs > 0) }
  rbindlist(out) }
one <- function(job) {
  dy <- job$rr; series <- prepare_series(dy); dy <- series$dy; mode <- c(exercise = "exercise", prcp = "tilt")[job$ds]; ev <- if (!is.na(mode)) events_from_rr(dy, mode, cfg$segmentation) else NULL
  fits <- list(M3 = fit_ext(series), M3O = fit_ext(series, list(osc = TRUE)))
  fits$M3E <- if (is.null(ev)) fits$M3 else fit_ext(series, list(inputs = list(events = ev)))
  fits$M3OE <- if (is.null(ev) || job$ds == "prcp") NULL else fit_ext(series, list(osc = TRUE, inputs = list(events = ev)))        # PRCP: too slow (12 regressors x oscillator grid), not requested
  cells <- list(M3 = list("M3", FALSE), `M3+` = list("M3", TRUE), M3O = list("M3O", FALSE), `M3O+` = list("M3O", TRUE), M3E = list("M3E", FALSE), `M3OE+` = list("M3OE", TRUE))
  base <- state_estimates(fits$M3, FALSE); rows <- list(); ex <- list(); tp <- list()
  for (nm in names(cells)) { f <- fits[[cells[[nm]][[1]]]]; if (is.null(f) || !isTRUE(f$ok)) next; e <- state_estimates(f, cells[[nm]][[2]]); dg <- diag_of(e$r, f$ctl)
    npar <- f$n_par + if (!is.null(f$spec$inputs)) 2 * nrow(f$spec$inputs$events) else 0
    rows[[nm]] <- data.table(model = nm, LL = e$r$LL_star, n_par = npar, ks_ok = dg[["ks_ok"]], pacf = dg[["pacf"]], tau_p = f$tau_p, tau_s = f$tau_s, cor_p = cor(e$xp, base$xp), cor_s = cor(e$xs, base$xs), sd_ratio_p = sd(e$xp) / sd(base$xp), sd_ratio_s = sd(e$xs) / sd(base$xs), n_events = if (!is.null(f$spec$inputs)) nrow(f$spec$inputs$events) else 0L)
    if (job$ds == "exercise" && !is.null(ev)) { z <- exercise_z(e, dy, ev); if (!is.null(z)) ex[[nm]] <- cbind(model = nm, z) }
    if (job$ds == "prcp" && !is.null(ev)) { z <- tilt_pattern(e, dy, ev); if (nrow(z)) tp[[nm]] <- cbind(model = nm, z) } }
  list(fit = cbind(rbindlist(rows), ds = job$ds, id = job$id, detected = !is.null(ev)), ex = if (length(ex)) cbind(rbindlist(ex), id = job$id), tp = if (length(tp)) cbind(rbindlist(tp), id = job$id))
}
log_msg("Real-data comparison | %d recordings", length(jobs)); out <- drop_failed(map_recordings(jobs, one, n_cores), "recording")
fit <- rbindlist(lapply(out, `[[`, "fit")); ex <- rbindlist(lapply(out, `[[`, "ex")); tp <- rbindlist(lapply(out, `[[`, "tp"))
dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(list(fit = fit, ex = ex, tp = tp), sprintf("experiments/m3plus/_out/compare_real_%s.rds", if (quick) "quick" else if (!is.null(ds_only)) ds_only else "full"))
b <- fit[model == "M3", .(ds, id, LL0 = LL, par0 = n_par)]; f <- merge(fit, b, by = c("ds", "id")); f[, `:=`(dLL = LL - LL0, aic_gain = 2 * (LL - LL0) - 2 * (n_par - par0))]
cat("\n== Fit (all recordings). Smoothed rows share the fit of their filtered twin ==\n")
print(f[model %in% c("M3", "M3O", "M3E", "M3OE+"), .(n = .N, extra_par = median(n_par - par0), dLL_med = round(median(dLL), 1), pct_AIC_better = round(100 * mean(aic_gain > 0), 1), ks_pass = round(100 * mean(ks_ok), 1), pacf_viol = round(mean(pacf), 1), tau_p_med = round(median(tau_p), 1), tau_s_med = round(median(tau_s)), pct_tau_p_in_prior = round(100 * mean(tau_p >= 3 & tau_p <= 6), 1)), by = model][order(model)], row.names = FALSE)
cat("\n== Same, by dataset (AIC-better %, KS pass %) ==\n"); print(dcast(f[model %in% c("M3O", "M3E", "M3OE+")], ds ~ model, value.var = c("aic_gain", "ks_ok"), fun.aggregate = function(x) round(100 * mean(if (is.logical(x)) x else x > 0), 1)), row.names = FALSE)
cat("\n== Agreement of the states with those of M3 (median over recordings) ==\n"); print(f[model != "M3", .(cor_vagal = round(median(cor_p), 3), cor_symp = round(median(cor_s), 3), sd_ratio_vagal = round(median(sd_ratio_p), 2), sd_ratio_symp = round(median(sd_ratio_s), 2)), by = model][order(model)], row.names = FALSE)
if (nrow(ex)) { cat("\n== Exercise (rest to end of exercise, in rest SDs), events from the RR series; median [2.5, 97.5] across recordings ==\n")
  q <- function(x) sprintf("%.1f [%.1f, %.1f]", median(x, na.rm = TRUE), quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
  print(ex[, .(n = .N, z_vagal = q(z_p), z_symp = q(z_s), rest_sd_vagal = signif(median(rest_sd_p), 3), rest_sd_symp = signif(median(rest_sd_s), 3)), by = model][order(model)], row.names = FALSE) }
if (nrow(tp)) { cat("\n== PRCP candidate tilts detected from the RR: share with vagal decrease and sympathetic increase ==\n"); print(tp[, .(events = .N, pct_expected = round(100 * mean(expected), 1), z_vagal_med = round(median(z_p), 1), z_symp_med = round(median(z_s), 1)), by = model][order(model)], row.names = FALSE) }
