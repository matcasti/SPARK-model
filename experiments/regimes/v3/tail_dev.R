# EXPERIMENTAL (round 4). Protection of the forward-score tail of the gated regime model, designed on DEVELOPMENT data only: the 154 exercise recordings used before the frozen validation, plus the 10 PRCP recordings.
# Candidates (all forward-pass only, the fit is unchanged): sigma^2 inflation by c (2, 4, 8, 16), sigma^2 floor at rho times M3's (0.25, 0.5, 1), contamination density (eps, kappa) = (0.01, 25), (0.01, 100) applied to both M3 and the regime model.
# Run from the project root:  Rscript experiments/regimes/v3/tail_dev.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); GATE <- 9
d <- load_dataset("exercise", cfg); ids <- names(d$rr); set.seed(5); a <- sample(ids, 100); set.seed(21); b <- sample(ids, 60); hard <- c("104", "099", "211", "208", "207", "138"); set.seed(11); c1 <- unique(c(hard, sample(setdiff(ids, hard), 34))); dev_ids <- sort(union(union(a, b), c1))
stopifnot(!length(intersect(dev_ids, readLines("experiments/regimes/v2/validation/ids.txt")))); jobs <- lapply(dev_ids, function(id) list(ds = "exercise", id = id, rr = d$rr[[id]]))
dp <- load_dataset("prcp", cfg); jobs <- c(jobs, lapply(names(dp$rr), function(id) list(ds = "prcp", id = id, rr = dp$rr[[id]])))
one <- function(job) {
  series <- prepare_series(job$rr); dy <- series$dy; n <- length(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); lp <- sh$priors$energy_ratio
  bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]; hr <- bh$run(); pp <- bh$p
  m3r <- function(eps, kappa) reg_run(dy, ctl, pp$kp, pp$ks, lp, pp$lR, hr$sig2, hr$nu0, matrix(0, 1, 2), matrix(1, 1, 1), w0 = 1, gm = gm, full = TRUE, eps = eps, kappa = kappa)
  r0 <- m3r(0, 25); out <- list(m3 = mean(r0$logdens[idx])); for (k in c(25, 100)) out[[paste0("m3_rob", k)]] <- mean(m3r(0.01, k)$logdens_rob[idx])
  open <- !is.null(cph) && cph$gain >= GATE; ok <- FALSE
  if (open) { h <- fit_regime(sh, 2, cph, bh, 7, 0.7, recip = 20); if (isTRUE(h$ok)) { ok <- TRUE; p <- h$p
      run <- function(sig2 = p$sig2, eps = 0, kappa = 25) reg_run(dy, ctl, p$kp, p$ks, lp, p$lR, sig2, p$nu0, p$mu, p$Pi, gm = gm, full = TRUE, eps = eps, kappa = kappa)
      out$hmm <- mean(run()$logdens[idx]); for (cc in c(2, 4, 8, 16)) out[[paste0("hmm_c", cc)]] <- mean(run(p$sig2 * cc)$logdens[idx]); for (rho in c(0.25, 0.5, 1)) out[[paste0("hmm_floor", rho)]] <- mean(run(max(p$sig2, rho * hr$sig2))$logdens[idx])
      for (k in c(25, 100)) out[[paste0("hmm_rob", k)]] <- mean(run(eps = 0.01, kappa = k)$logdens_rob[idx]) } }
  if (!ok) for (nm in c("hmm", paste0("hmm_c", c(2, 4, 8, 16)), paste0("hmm_floor", c(0.25, 0.5, 1)), paste0("hmm_rob", c(25, 100)))) out[[nm]] <- if (grepl("rob", nm)) out[[sub("hmm", "m3", nm)]] else out$m3
  out$sig2_ratio <- if (ok) h$p$sig2 / hr$sig2 else NA_real_; as.data.table(c(list(ds = job$ds, id = job$id, open = open, tail_bad = mean(dy[idx] > 2)), out)) }
log_msg("Tail protection, development data | %d recordings", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording"), fill = TRUE); saveRDS(res, "experiments/regimes/v3/_out/tail_dev.rds"); cat("Saved", nrow(res), "rows\n")
x <- res[tail_bad < 0.1]; summ <- function(d) list(n = length(d), mean = round(mean(d), 3), median = round(median(d), 3), trimmed = round(mean(d, trim = 0.1), 3), better = round(100 * mean(d > 0)), worst = round(min(d), 2))
cat("\nForward score minus M3 (nats/beat); exercise (valid tails) and PRCP. Robust-density rows are compared with M3 under the same density.\n")
for (dss in c("exercise", "prcp")) { y <- x[ds == dss]; cat("\n==", dss, "\n"); tb <- rbindlist(c(lapply(c("hmm", paste0("hmm_c", c(2, 4, 8, 16)), paste0("hmm_floor", c(0.25, 0.5, 1))), function(nm) c(model = nm, summ(y[[nm]] - y$m3))), lapply(c(25, 100), function(k) c(model = paste0("hmm_rob", k, " vs m3_rob", k), summ(y[[paste0("hmm_rob", k)]] - y[[paste0("m3_rob", k)]])))), fill = TRUE); print(tb) }
