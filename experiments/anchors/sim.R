# EXPERIMENTAL. Score the anchor rules of lib.R on the pipeline's simulated scenarios (same seeds) against the truth. Run from the project root:
#   Rscript experiments/anchors/sim.R
source("analysis/00_setup.R"); source("experiments/anchors/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); RULES <- if (length(commandArgs(TRUE))) strsplit(commandArgs(TRUE)[1], ",")[[1]] else c("hard", "bandfix", "lfhf", "spec", "oracle")
one <- function(job) {
  truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]())
  sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  series <- prepare_series(sim$rr); true_share <- var(sim$x_p) / (var(sim$x_p) + var(sim$x_s))
  rbindlist(lapply(RULES, function(r) { fit <- fit_rule(series, make_rule(r, series, truth$lambda_p)); x <- recovery_row(fit, truth, sim, true_share); x[, `:=`(rule = r, ks_ok = fit$diag$ks_ok, pacf_viol = fit$diag$pacf_viol)] }))[, `:=`(scenario = job$scenario, sim = job$sim)]
}
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Anchor rules on %d simulated recordings x %d rules", length(jobs), length(RULES))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); res[, rule := factor(rule, RULES)]
dir.create("experiments/anchors/_out", showWarnings = FALSE); saveRDS(res, if (length(commandArgs(TRUE))) "experiments/anchors/_out/sim_const.rds" else "experiments/anchors/_out/sim.rds")
f <- function(d) d[, .(n = .N, abs_ln_err_lp = round(median(abs(ln_err_lp)), 2), bias_ln_lp = round(median(ln_err_lp), 2), cor_xp = round(mean(cor_x_p), 3), cor_xs = round(mean(cor_x_s), 3), rmse_xp = round(mean(rmse_x_p), 3), rmse_xs = round(mean(rmse_x_s), 3),
                    abs_err_tau_p = round(median(abs(ln_err_tau_p)), 2), abs_err_tau_s = round(median(abs(ln_err_tau_s)), 2), cover_xp = round(mean(cover_x_p), 2), cover_xs = round(mean(cover_x_s), 2), LL = round(median(LL_star), 1)), by = rule][order(rule)]
for (sc in names(RECOVERY_SCENARIOS)) { cat("\n==", sc, "==\n"); print(f(res[scenario == sc]), row.names = FALSE) }
cat("\n== All scenarios ==\n"); print(f(res), row.names = FALSE)
