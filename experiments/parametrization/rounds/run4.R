# ==============================================================================
# experiments/parametrization/rounds/run4.R -- EXPERIMENTAL, round 4: lR-stratified designed starts; finer "racing" grids and designed starts.
# Racing = run every start for only m iterations, then finish only the best one. Run from the project root.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/parametrization/strategies.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
prev_best <- unique(readRDS("experiments/parametrization/_out/run_full.rds")[, .(ds, id, best)])
r2 <- readRDS("experiments/parametrization/_out/round2_full.rds"); prev_best2 <- unique(r2[, .(ds, id, best_all)])
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); for (id in names(d$rr)) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
all_ids <- vapply(jobs, `[[`, "", "id")

nl <- function(obj, par, z0, iter.max = 1000) {
  fz <- function(z) { p <- par$to(z); obj$f(p[1], p[2], p[3]) }
  o <- tryCatch(nlminb(z0, fz, control = list(iter.max = iter.max, eval.max = 2 * iter.max)), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$objective)) list(v = Inf, z = z0) else list(v = o$objective, z = o$par)
}
one <- function(job) {
  series <- prepare_series(job$rr); pri <- series$priors; ctl <- series$control
  obj <- make_objective(series); W <- parametrizations(pri)$whitened
  idx <- match(job$id, all_ids); set.seed(1000L + idx)
  default <- c(pri$ks_mode, pri$ks_mode + max(pri$kp_mode - pri$ks_mode, ctl$kp_gap_floor), ctl$lR_start)
  disp <- dispersed_starts(pri, 6)
  slow_noisy <- c(1 / 150, 1 / 25, 0.3); low_lR <- c(1 / 160, 1 / 12, 1e-5)
  # starts stratified in the noise-to-signal ratio lR (the direction where the narrow basins were found), tau_p spread alongside
  s_lo <- c(1 / 160, 1 / 12, 1e-6); s_hi <- c(1 / 150, 1 / 25, 0.3); s_vlo <- c(1 / 170, 1 / 14, 1e-7); s_mid <- c(1 / 200, 1 / 6, 1e-2)
  sets <- list(LR3 = list(default, s_lo, s_hi), LR3b = list(default, s_vlo, s_hi), LR4 = list(default, s_lo, s_hi, s_mid), LR4b = list(default, s_vlo, s_lo, s_hi))
  rows <- list()
  for (nm in names(sets)) for (m in c(4, 5, 6, 8)) {
    S <- sets[[nm]]; obj$counter$n <- 0L
    short <- lapply(S, function(s) nl(obj, W, W$from(s[1], s[2], s[3]), iter.max = m))
    b <- which.min(vapply(short, `[[`, 0, "v")); fin <- nl(obj, W, short[[b]]$z)
    rows[[length(rows) + 1]] <- data.table(ds = job$ds, id = job$id, method = sprintf("race_%s_m%d", nm, m), obj = fin$v, evals = obj$counter$n)
  }
  rbindlist(rows)
}
res <- parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) { message("FAILED ", j$id, ": ", conditionMessage(e)); NULL }),
                          mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
res <- merge(res, prev_best2, by = c("ds", "id"))
res[, best_all := pmin(best_all, min(obj[is.finite(obj)])), by = .(ds, id)]
res[, hit := is.finite(obj) & obj - best_all < 0.01]
saveRDS(res, "experiments/parametrization/_out/round4_full.rds")
out <- res[, .(success_pct = round(100 * mean(hit), 1), evals_mean = round(mean(evals)), evals_p95 = round(quantile(evals, .95)), misses = paste(id[!hit], collapse = ",")), by = method]
out[, cost_vs_baseline := round(evals_mean / 273, 2)]
print(out[order(-success_pct, evals_mean)], row.names = FALSE)
