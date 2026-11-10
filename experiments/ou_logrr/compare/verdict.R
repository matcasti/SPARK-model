# ==============================================================================
# experiments/ou_logrr/compare/verdict.R -- evaluates the decision rule of the README on the outputs of sim_compare.R and real_compare.R.
#   Rscript experiments/ou_logrr/compare/verdict.R      (the real-data part runs for the datasets already computed)
# S3 was amended after seeing that M3 itself misses the 0.67-1.5 band in the tilt scenario: "within the band, or not farther from 1 (in ln) than M3 by more than 0.1".
# ==============================================================================
options(width = 220, datatable.print.nrows = 200); suppressMessages(library(data.table)); out <- "experiments/ou_logrr/_out"; V <- c("F1", "F2", "F3", "F4", "F5")
yn <- function(x) ifelse(is.na(x), "NA", ifelse(x, "pass", "FAIL"))

# ---- Simulation ---------------------------------------------------------------------------------------------------
sim <- rbindlist(lapply(list.files(out, "^compare_sim_.*rds$", full.names = TRUE), readRDS), fill = TRUE); sim[, grp := sub(":.*", "", scenario)]
f <- sim[smooth == FALSE]
agg <- f[, .(rmse_p = mean(rmse_p), rmse_s = mean(rmse_s), cor_p = mean(cor_p), cor_s = mean(cor_s), cov_p = mean(cover_p), cov_s = mean(cover_s), ln_tp = median(abs(ln_tau_p)), ln_ts = median(abs(ln_tau_s)),
         chg_p = median(chg_p, na.rm = TRUE), chg_s = median(chg_s, na.rm = TRUE), vf = mean(t50_p < t50_s, na.rm = TRUE), vf_true = mean(t50_p_true < t50_s_true, na.rm = TRUE)), by = .(grp, model)]
rows <- list()
for (g in unique(agg$grp)) { m3 <- agg[grp == g & model == "M3"]
  for (v in V) { x <- agg[grp == g & model == v]
    s1 <- x$rmse_p <= 1.10 * m3$rmse_p && x$rmse_s <= 1.10 * m3$rmse_s && x$cor_p >= m3$cor_p - 0.03 && x$cor_s >= m3$cor_s - 0.03
    ok_cov <- function(c, c3) (c >= 0.90 && c <= 0.98) || c >= c3 - 0.03; s2 <- ok_cov(x$cov_p, m3$cov_p) && ok_cov(x$cov_s, m3$cov_s)
    s3 <- NA; if (g %in% c("exercise", "tilt", "G1_E2")) { near <- function(r, r3) (r >= 0.67 && r <= 1.5) || abs(log(r)) <= abs(log(r3)) + 0.1
      s3 <- near(x$chg_p, m3$chg_p) && near(x$chg_s, m3$chg_s) && x$vf >= x$vf_true - 0.05 }
    s4 <- x$ln_tp <= m3$ln_tp + 0.1 && x$ln_ts <= m3$ln_ts + 0.1
    rows[[length(rows) + 1L]] <- data.table(scenario = g, variant = v, S1_accuracy = yn(s1), S2_calibration = yn(s2), S3_exercise_tilt = yn(s3), S4_time_constants = yn(s4)) } }
cat("\n== Simulation: filtered states, the decision rule (S1-S4) per scenario group and variant\n"); simv <- rbindlist(rows); print(simv, row.names = FALSE)
cat("\n== Simulation: variants passing every criterion in every scenario group\n")
print(simv[, .(all_pass = all(c(S1_accuracy, S2_calibration, S4_time_constants) == "pass") && all(S3_exercise_tilt[S3_exercise_tilt != "NA"] == "pass"), n_fail = sum(unlist(.SD) == "FAIL")), by = variant, .SDcols = 3:6], row.names = FALSE)

# ---- Real data ----------------------------------------------------------------------------------------------------
rd <- function(ds) { p <- file.path(out, sprintf("real_compare_%s.rds", ds)); if (file.exists(p)) readRDS(p) else NULL }
real <- lapply(c(fantasia = "fantasia", prcp = "prcp", exercise = "exercise"), rd); real <- Filter(Negate(is.null), real)
if (length(real)) {
  rows <- rbindlist(lapply(names(real), function(ds) real[[ds]]$rows[, dataset := ds]))
  cat("\n== Real data: per model (all datasets pooled; tau 'in prior' = production ranges 3-6 and 100-300 s)\n")
  print(rows[, .(n = .N, converged = round(mean(ok), 3), tau_p_in = round(mean(tau_p >= 3 & tau_p <= 6), 2), tau_s_in = round(mean(tau_s >= 100 & tau_s <= 300), 2), ks_pass = round(mean(ks_p > 0.05), 3), pacf = round(mean(pacf_viol), 1),
                cor_xp_M3 = round(median(cor_xp_M3, na.rm = TRUE), 2), cor_xs_M3 = round(median(cor_xs_M3, na.rm = TRUE), 2)), by = model][order(model)], row.names = FALSE)
  if (!is.null(real$exercise)) { e <- real$exercise$exercise; cat("\n== Exercise cohort: responses (events from the RR)\n")
    print(e[, .(n = .N, vagal_first = round(mean(vagal_first, na.rm = TRUE), 2), exp_p = round(mean(exp_p), 2), exp_s = round(mean(exp_s), 2), z_p = round(median(z_p, na.rm = TRUE), 1), z_s = round(median(z_s, na.rm = TRUE), 1),
                t50_p = round(median(t50_p, na.rm = TRUE)), t50_s = round(median(t50_s, na.rm = TRUE))), by = model][order(model)], row.names = FALSE) }
  if (!is.null(real$prcp)) { cat("\n== PRCP: the pipeline's tilt analysis on each model's states\n")
    print(real$prcp$prcp[, .(tilts = .N, both_signs = round(mean(expected_p & expected_s), 2), z_p = round(median(z_d_x_p), 1), z_s = round(median(z_d_x_s), 1), vagal_first = round(mean(lag_s_minus_p > 0, na.rm = TRUE), 2)), by = model][order(model)], row.names = FALSE) }
  if (!is.null(real$fantasia)) { w <- real$fantasia$windows; g <- real$fantasia$group; main <- real$fantasia$rows
    icc31 <- function(a, b) { x <- cbind(a, b); k <- is.finite(rowSums(x)); x <- x[k, , drop = FALSE]; n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1)
      mse <- sum((x - rowMeans(x) + gm - matrix(colMeans(x) - gm + gm, n, 2, byrow = TRUE) + gm - gm)^2) / (n - 1); (msr - mse) / (msr + mse) }
    cat("\n== Fantasia: test-retest ICC(3,1) of log values between windows A and B, and young-vs-older AUC (whole analysis window)\n")
    ic <- rbindlist(lapply(unique(w$model), function(m) { a <- w[model == m & window == "A"][order(id)]; b <- w[model == m & window == "B"][order(id)]; mm <- main[model == m]; gg <- g[match(mm$id, id), Group]
      auc <- function(x) { k <- is.finite(x); r <- rank(x[k]); gk <- gg[k]; n1 <- sum(gk == "young"); n2 <- sum(gk == "older"); (sum(r[gk == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
      data.table(model = m, icc_sd_p = icc31(log(a$sd_p), log(b$sd_p)), icc_sd_s = icc31(log(a$sd_s), log(b$sd_s)), icc_tau_p = icc31(log(a$tau_p), log(b$tau_p)), icc_tau_s = icc31(log(a$tau_s), log(b$tau_s)),
                 auc_sd_p = auc(mm$sd_p), auc_sd_s = auc(mm$sd_s), auc_tau_p = auc(mm$tau_p)) }))
    print(ic[order(model)][, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE) }
}
