# EXPERIMENTAL. What makes the KS test of the M3 innovations fail? For each recording: KS pass of (a) the diagnostic innovations as they are, (b) after clipping at |z| = 3 (tails),
# (c) after dividing by their local scale (rolling SD over 51 beats: heteroskedasticity), (d) after removing an AR(4) fit (serial dependence), (e) after (c) and (d) together.
# Run from the project root:  Rscript experiments/m3plus/structure/diagnose.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); for (id in names(d$rr)) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
ks_ok <- function(z) { z <- z[is.finite(z)]; n <- length(z); unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n) }
one <- function(j) {
  series <- prepare_series(j$rr); f <- fit_ext(series); r <- f$run(full = TRUE); k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < f$ctl$diag_gate_max); z <- r$innovations_diagnostic[k]; n <- length(z)
  loc <- as.numeric(sqrt(stats::filter(z^2, rep(1 / 51, 51), sides = 2))); loc[is.na(loc)] <- sqrt(mean(z^2)); zl <- z / loc
  ar <- tryCatch(ar(z, order.max = 4, aic = FALSE, method = "yule-walker"), error = function(e) NULL); zw <- if (is.null(ar)) z else as.numeric(na.omit(ar$resid)) / sd(ar$resid, na.rm = TRUE)
  ar2 <- tryCatch(ar(zl, order.max = 4, aic = FALSE, method = "yule-walker"), error = function(e) NULL); zlw <- if (is.null(ar2)) zl else as.numeric(na.omit(ar2$resid)) / sd(ar2$resid, na.rm = TRUE)
  q <- function(x) quantile(x, c(.1, .9)); cv <- sd(loc) / mean(loc)
  data.table(ds = j$ds, id = j$id, n = n, raw = ks_ok(z), clipped = ks_ok(pmax(pmin(z, 3), -3)), local_scale = ks_ok(zl / sd(zl)), ar4 = ks_ok(zw), both = ks_ok(zlw),
             kurt = mean(((z - mean(z)) / sd(z))^4) - 3, skew = mean(((z - mean(z)) / sd(z))^3), scale_cv = cv, scale_ratio = max(loc) / min(loc), lag1 = acf(z, plot = FALSE, lag.max = 1)$acf[2])
}
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording")); dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(res, "experiments/m3plus/_out/diagnose.rds")
cat("\nKS pass rate (%) of the M3 innovations, as they are and after removing one feature at a time:\n")
print(res[, lapply(.SD, function(x) round(100 * mean(x), 1)), .SDcols = c("raw", "clipped", "local_scale", "ar4", "both"), by = ds], row.names = FALSE); print(res[, lapply(.SD, function(x) round(100 * mean(x), 1)), .SDcols = c("raw", "clipped", "local_scale", "ar4", "both")], row.names = FALSE)
cat("\nShape of the innovations (median): excess kurtosis, skewness, CV of the local scale, max/min local scale, lag-1 autocorrelation\n")
print(res[, .(kurt = round(median(kurt), 2), skew = round(median(skew), 2), scale_cv = round(median(scale_cv), 2), scale_ratio = round(median(scale_ratio), 1), lag1 = round(median(lag1), 3)), by = ds], row.names = FALSE)
