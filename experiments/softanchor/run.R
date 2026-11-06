# ==============================================================================
# experiments/softanchor/run.R -- EXPERIMENTAL. "Soft anchor": lambda_P is free with a log-normal prior centred on the spectral anchor
# (width s); s = 0 is M3 (hard anchor), s = Inf is M3 without the anchor. Same objective otherwise; the optimizer is the 3-start scheme.
# Run from the project root:  Rscript experiments/softanchor/run.R [quick|full]
# ==============================================================================
mode <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "quick"
source("analysis/00_setup.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
spec_free <- drop_term(model_spec("M3"), "anchor")                      # lambda_P is the 4th parameter
S_VALUES <- c(0, 0.15, 0.3, 0.6, 1.2, Inf)

diagnostics <- function(f, ctl) {
  keep <- which(is.finite(f$innovations_diagnostic) & f$gate_multiplier < ctl$diag_gate_max); z <- f$innovations_diagnostic[keep]; n <- length(z)
  pac <- as.numeric(pacf(z, lag.max = min(ctl$pacf_lag_max, max(1, n - 5)), plot = FALSE)$acf)
  list(ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n), pacf_viol = sum(abs(pac) > qnorm(0.975) / sqrt(n)))
}
fit_soft <- function(series, s) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; ln_anchor <- log(pri$energy_ratio)
  if (s == 0) {                                                          # M3: hard anchor, three-start scheme of the pipeline
    spec <- model_spec("M3"); o <- optimize_spec(series, spec); th <- c(o$par, ln_anchor); par <- unpack_theta(o$par, spec, pri$energy_ratio); value <- o$value
  } else {
    cost <- function(th) { v <- variant_objective(th, dy, spec_free, pri, ctl); if (v >= 1e11) return(1e12); if (is.finite(s)) v + 0.5 * ((th[4] - ln_anchor) / s)^2 else v }
    runs <- lapply(.m3_starts(pri, ctl), function(st) tryCatch(optim(c(unname(.par_to_theta(st)), ln_anchor), cost, method = "BFGS", control = list(maxit = ctl$maxit)), error = function(e) NULL))
    vals <- vapply(runs, function(o) if (is.null(o) || !is.finite(o$value)) Inf else o$value, 0); k <- which(vals <= min(vals) + 0.01)[1]
    th <- runs[[k]]$par; value <- vals[k]; par <- unpack_theta(th, spec_free, pri$energy_ratio)
  }
  f <- filter_eval(dy, par$kp, par$ks, par$lp, par$lR, ctl, full = TRUE); d <- diagnostics(f, ctl)
  data.table(s = s, cost = value, LL = f$LL_star, ln_ratio = th[4] - ln_anchor, tau_p = 1 / par$kp, tau_s = 1 / par$ks, lR = par$lR, nu0 = f$nu0, sig2 = f$sig2,
             sd_xp = sd(f$X_p), sd_xs = sd(f$X_s), ks_ok = d$ks_ok, pacf_viol = d$pacf_viol)
}
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr)
  if (mode == "quick") { set.seed(3); ids <- sort(sample(ids, min(length(ids), c(prcp = 4, fantasia = 12, exercise = 14)[[ds]]))) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]], group = if (ds == "fantasia") substr(id, 3, 3) else NA_character_) }
log_msg("Soft anchor | %s | %d recordings", mode, length(jobs))
res <- parallel::mclapply(jobs, function(j) tryCatch({ series <- prepare_series(j$rr)
  r <- rbindlist(lapply(S_VALUES, function(s) fit_soft(series, s))); r[, `:=`(ds = j$ds, id = j$id, group = j$group, anchor = series$priors$energy_ratio)] }, error = function(e) { message("FAILED ", j$id, ": ", conditionMessage(e)); NULL }),
  mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
dir.create("experiments/softanchor/_out", showWarnings = FALSE); saveRDS(res, sprintf("experiments/softanchor/_out/%s.rds", mode))
base <- res[s == 0, .(ds, id, LL0 = LL)]; res <- merge(res, base, by = c("ds", "id")); res[, dLL := LL - LL0]
cat("\n== By prior width s (all recordings) ==\n")
print(res[, .(n = .N, dLL_med = round(median(dLL), 2), pct_dLL_gt_1 = round(100 * mean(dLL > 1), 1), pct_dLL_gt_1.92 = round(100 * mean(dLL > 1.92), 1), abs_ln_ratio_med = round(median(abs(ln_ratio)), 2),
              ks_pass = round(100 * mean(ks_ok), 1), pacf_viol = round(mean(pacf_viol), 1), tau_p_med = round(median(tau_p), 1), tau_s_med = round(median(tau_s), 1)), by = s], row.names = FALSE)
cat("\n== Direction of the shift (ln lambda_P hat - ln anchor), s = Inf, by dataset ==\n")
print(res[s == Inf, .(n = .N, median = round(median(ln_ratio), 2), q25 = round(quantile(ln_ratio, .25), 2), q75 = round(quantile(ln_ratio, .75), 2)), by = ds], row.names = FALSE)
