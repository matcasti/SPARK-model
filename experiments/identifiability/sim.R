# ==============================================================================
# experiments/identifiability/sim.R -- where does the vagal / sympathetic decomposition stop being identifiable? Simulated recordings with known drives, every
# specification of the paper (M0 structural constraints only, M1 + time-scale priors, M2 + spectral anchor, M3 = M1 + M2) fitted to each. Axes: the time-scale ratio
# tau_s / tau_p (tau_p = 4 s), the variance ratio lambda_P (V_p / V_s) and the recording length. Metrics: correlation of the filtered states with the truth, share of the variance
# given to the vagal drive against the true share, share of fits that put (almost) all variance in one drive, error of the time constants.
#   Rscript experiments/identifiability/sim.R            -> _out/sim.rds and a summary on screen
# ==============================================================================
source("analysis/00_setup.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); N <- 30L
grid_A <- CJ(axis = "time-scale ratio", ratio = c(1.5, 2, 4, 10, 30, 100), lambda = 2, duration = c(300, 600, 1800))      # tau_s = ratio * 4 s
grid_B <- CJ(axis = "variance ratio", ratio = 30, lambda = c(0.2, 0.5, 2, 6, 20), duration = c(300, 600, 1800))
cells <- rbind(grid_A, grid_B); jobs <- unlist(lapply(seq_len(nrow(cells)), function(i) lapply(seq_len(N), function(j) c(as.list(cells[i]), list(sim = j, seed = id_seed(paste("ident", i, j), cfg$seed))))), recursive = FALSE)
one <- function(job) {
  tau_p <- 4; tau_s <- job$ratio * tau_p
  sim <- simulate_rr(job$duration, kp = 1 / tau_p, ks = 1 / tau_s, lp = job$lambda, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  fits <- fit_specs(prepare_series(sim$rr)); true_share <- var(sim$x_p) / (var(sim$x_p) + var(sim$x_s))
  truth <- list(tau_p = tau_p, tau_s = tau_s, lambda_p = job$lambda)
  rows <- rbindlist(lapply(fits, recovery_row, truth = truth, simulated = sim, true_share = true_share)); rows[, `:=`(axis = job$axis, ratio = job$ratio, lambda = job$lambda, duration = job$duration, sim = job$sim)]
}
log_msg("Identifiability grid: %d cells x %d simulations = %d series, 4 specifications each", nrow(cells), N, length(jobs))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); res <- res[ok == TRUE]
saveRDS(res, "experiments/identifiability/_out/sim.rds")
res[, collapse := est_share < 0.05 | est_share > 0.95]
options(width = 200)
for (ax in c("time-scale ratio", "variance ratio")) { cat(sprintf("\n===== %s: correlation of the filtered sympathetic-like state with the truth, median over simulations (and share of one-process fits) =====\n", ax))
  d <- res[axis == ax]; key <- if (ax == "time-scale ratio") "ratio" else "lambda"
  print(dcast(d[, .(cor_s = sprintf("%.2f (%d%%)", median(cor_x_s), round(100 * mean(collapse)))), by = c(key, "spec", "duration")], as.formula(paste(key, "+ duration ~ spec")), value.var = "cor_s"), row.names = FALSE) }
