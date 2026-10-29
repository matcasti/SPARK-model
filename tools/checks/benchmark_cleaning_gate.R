# Choice of the RR cleaning gate (clean_sample) on simulated series with injected artifacts.
# Artifacts (placed in-place so beat counts are preserved): ectopic pairs (short-long),
# missed-detection-like long intervals (about 2x RR), false-detection-like short intervals
# (about 0.4x RR) and isolated spikes. Reported per (n_mad, n_iter):
#   false_alter  proportion of artifact-free beats changed by > 1 ms (should be small)
#   sens         proportion of artifact beats whose error is at least halved
#   rmssd_ratio  RMSSD of the cleaned series / RMSSD of the artifact-free series (1 = ideal)
Sys.setenv(SPARK_CONFIG = normalizePath("tests/fixtures/config_test.R"))
source("analysis/00_setup.R")
n_sim <- 30; grid <- expand.grid(n_mad = c(2, 3, 3.5, 4, 5), n_iter = c(1, 2, 5))
res <- list()
for (i in seq_len(n_sim)) {
  g <- if (i %% 2) "young" else "older"
  true <- sim_rest(1200 + i, g, 600, m_p = rsa_drive(0.25, if (g == "young") 0.05 else 0.025))$rr
  true <- round(true / 0.004) * 0.004
  set.seed(i); n <- length(true); x <- true; art <- rep(FALSE, n)
  pick <- sort(sample(seq(20, n - 20), round(0.03 * n)))
  for (k in pick) {
    type <- sample(4, 1)
    if (type == 1 && !art[k + 1]) { x[k] <- 0.65 * true[k]; x[k + 1] <- true[k + 1] + 0.35 * true[k]; art[c(k, k + 1)] <- TRUE }
    if (type == 2) { x[k] <- 2 * true[k]; art[k] <- TRUE }
    if (type == 3) { x[k] <- 0.4 * true[k]; art[k] <- TRUE }
    if (type == 4) { x[k] <- true[k] * runif(1, 1.35, 1.7); art[k] <- TRUE }
  }
  for (j in seq_len(nrow(grid))) {
    cl <- clean_sample(x, n_mad = grid$n_mad[j], n_iter = grid$n_iter[j])
    res[[length(res) + 1]] <- data.table(sim = i, n_mad = grid$n_mad[j], n_iter = grid$n_iter[j],
      false_alter = mean(abs(cl - x)[!art] > 0.001), sens = mean((abs(cl - true) <= 0.5 * abs(x - true))[art]),
      rmssd_ratio = sqrt(mean(diff(cl)^2)) / sqrt(mean(diff(true)^2)), sdnn_ratio = sd(cl) / sd(true),
      rmssd_ratio_uncleaned = sqrt(mean(diff(x)^2)) / sqrt(mean(diff(true)^2)))
  }
}
r <- rbindlist(res)[, lapply(.SD, mean), by = .(n_mad, n_iter), .SDcols = c("false_alter", "sens", "rmssd_ratio", "sdnn_ratio", "rmssd_ratio_uncleaned")]
print(r[order(n_iter, n_mad)], digits = 3)
fwrite(r, file.path(Sys.getenv("SPARK_TEST_ROOT", "tests/_run"), "cleaning_gate_benchmark.csv"))
