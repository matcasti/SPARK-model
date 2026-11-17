# ==============================================================================
# experiments/features/second_order.R -- EXPERIMENTAL. Decisive check for Q1: is the gain of entropy/geometric/DFA features over the band ratios
# information beyond second order, or just a finer description of the autocorrelation? Adds a full scale-free second-order set (A: autocorrelation
# of RR at lags 1-20 and relative log power in 12 log-spaced bins 0.01-0.5 Hz) and compares out-of-fold errors. Same simulations and folds as sim.R.
#   Rscript experiments/features/second_order.R
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/features/lib.R")
n_cores <- max(1L, min(6L, cfg$cores - 1L)); out <- "experiments/features/_out"
D <- readRDS(file.path(out, "sim_features.rds")); q1 <- readRDS(file.path(out, "sim_q1.rds")); folds <- q1$folds
EDGES <- exp(seq(log(0.01), log(0.5), length.out = 13))
second_order <- function(rr) {
  ac <- as.numeric(acf(rr, lag.max = 20, plot = FALSE)$acf)[-1]
  tg <- rr_tachogram(rr); sp <- spec.pgram(ts(tg$x, frequency = 4), taper = 0.1, detrend = TRUE, plot = FALSE, spans = c(3, 3))
  bp <- vapply(seq_len(12), function(b) sum(sp$spec[sp$freq >= EDGES[b] & sp$freq < EDGES[b + 1]]), numeric(1)); bp <- log(pmax(bp, 1e-12) / sum(bp))
  c(setNames(ac, sprintf("A_acf%02d", 1:20)), setNames(bp, sprintf("A_bin%02d", 1:12)))
}
log_msg("Second-order features of %d simulations", nrow(D))
A <- rbindlist(drop_failed(map_recordings(D$i, function(i) { p <- draw_truth(id_seed(paste("features", i), cfg$seed), i > 500); s <- simulate_one(p, id_seed(paste("features", i), cfg$seed))
  data.table(i = i, t(second_order(s$rr))) }, n_cores), "simulation"))
D <- merge(D, A, by = "i"); setkey(D, i); stopifnot(identical(D$i, sort(D$i)))
y <- log(D$lambda_eff); pat <- function(...) unlist(c(SETS[intersect(c(...), names(SETS))], if ("A" %in% c(...)) "^A_"))
SET2 <- list(Sf = "Sf", EGfF = c("E", "Gf", "F"), SfEGfF = c("Sf", "E", "Gf", "F"), A = "A", SfA = c("Sf", "A"), SfA_EGfF = c("Sf", "A", "E", "Gf", "F"))
res <- rbindlist(lapply(names(SET2), function(nm) { p <- with_seed(cfg$seed, cv_predict(D, y, pat(SET2[[nm]]), folds))
  rbindlist(lapply(list(all = TRUE, model = !D$osc, oscillator = D$osc), function(k) data.table(err = median(abs(p[k] - y[k])), R2 = 1 - sum((p[k] - y[k])^2) / sum((y[k] - mean(y[k]))^2))), idcol = "generator")[, set := nm] }))
cat("\n== Out-of-fold median |ln error| of ln lambda_eff (R2), scale-free sets ==\n")
print(dcast(res[, .(set, generator, v = sprintf("%.3f (%.2f)", err, R2))], set ~ generator, value.var = "v")[match(names(SET2), set)], row.names = FALSE)
saveRDS(res, file.path(out, "sim_second_order.rds"))

# ---- Likelihood check: the profiled log-likelihood of M3 over a grid of lambda_P (L) is what the model itself knows about lambda_P. --------------
# If E/G/F add nothing to L, they carry no information about lambda_P beyond the likelihood; if they do, the excess comes from outside the model's
# Gaussian likelihood (the integrate-and-fire observation and the oscillator make the simulated RR slightly non-Gaussian).
GRID <- exp(seq(log(0.25), log(12), length.out = 9))
log_msg("Likelihood profiles: %d simulations x %d values of lambda_P", nrow(D), length(GRID))
L <- rbindlist(drop_failed(map_recordings(D$i, function(i) { sd <- id_seed(paste("features", i), cfg$seed); s <- simulate_one(draw_truth(sd, i > 500), sd); ser <- prepare_series(s$rr)
  ll <- vapply(GRID, function(g) { f <- fit_with_anchor(ser, g); if (isTRUE(f$ok)) f$LL_star else NA_real_ }, numeric(1)); ll <- ll - max(ll, na.rm = TRUE)
  data.table(i = i, t(setNames(ll, sprintf("L_ll%d", seq_along(GRID)))), L_argmax = log(GRID[which.max(ll)])) }, n_cores), "profile"))
D <- merge(D, L, by = "i"); setkey(D, i)
pat2 <- function(...) c(pat(...), if ("L" %in% c(...)) "^L_")
SET3 <- list(L = "L", L_Sf = c("L", "Sf"), L_SfA = c("L", "Sf", "A"), L_EGfF = c("L", "E", "Gf", "F"), L_SfA_EGfF = c("L", "Sf", "A", "E", "Gf", "F"))
res3 <- rbindlist(lapply(names(SET3), function(nm) { p <- with_seed(cfg$seed, cv_predict(D, y, pat2(SET3[[nm]]), folds))
  rbindlist(lapply(list(all = TRUE, model = !D$osc, oscillator = D$osc), function(k) data.table(err = median(abs(p[k] - y[k])), R2 = 1 - sum((p[k] - y[k])^2) / sum((y[k] - mean(y[k]))^2))), idcol = "generator")[, set := nm] }))
cat("\n== With the likelihood profile L (median |ln error| (R2)) ==\n")
print(dcast(res3[, .(set, generator, v = sprintf("%.3f (%.2f)", err, R2))], set ~ generator, value.var = "v")[match(names(SET3), set)], row.names = FALSE)
cat(sprintf("\nProfile flatness: median range of the profiled log-likelihood over lambda_P 0.25-12 = %.2f units; argmax error median |ln| = %.3f\n",
            median(-apply(as.matrix(D[, grep("^L_ll", names(D)), with = FALSE]), 1, min, na.rm = TRUE)), median(abs(D$L_argmax - y))))
saveRDS(list(second_order = res, likelihood = res3, profiles = L), file.path(out, "sim_second_order.rds"))
