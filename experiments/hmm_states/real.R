# ==============================================================================
# experiments/hmm_states/real.R -- the second-stage HMM on the 313 real recordings (states and estimation variances of results/fit_*.rds, i.e. production M3).
#   Rscript experiments/hmm_states/real.R           (needs _out/null.rds from sim.R for the gates)
# Fantasia (rest, 10 min): how often does the model see a regime where nothing was imposed. PRCP (10 recordings, 4 tilts each): agreement of the regimes with the tilt angle (> 20 deg).
# Exercise (263): onset and the second transition against the events detected from the heart rate by the pipeline's detector.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R"); source("experiments/hmm_states/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); options(width = 200); nul <- readRDS("experiments/hmm_states/_out/null.rds"); gate <- nul[, .(g = quantile(lr, .95)), by = mode]; g <- setNames(gate$g, gate$mode)
trans <- function(path, t, min_run = 6L) { r <- rle(path); ends <- cumsum(r$lengths); keep <- r$lengths >= min_run; idx <- which(keep); if (length(idx) < 2) return(numeric(0)); t[ends[idx[-length(idx)]] + 1] }
one <- function(job) { b <- bin_states(job$t, job$xp, job$xs, job$vp, job$vs); if (nrow(b) < 40) return(NULL); V <- job$V; Y <- as.matrix(b[, .(y1, y2)]); Vt <- as.matrix(b[, .(v1, v2)]); out <- list()
  for (m in c("aware", "naive")) { f1 <- hmm_fit(Y, Vt, V, 1, m); f2 <- hmm_fit(Y, Vt, V, 2, m); f3 <- hmm_fit(Y, Vt, V, 3, m); n <- nrow(b)
    bic <- c(-2 * f1$ll + f1$npar * log(n), -2 * f2$ll + f2$npar * log(n), -2 * f3$ll + f3$npar * log(n)); base <- as.integer(names(which.max(table(f2$path[seq_len(15)]))))
    pb2 <- 1 - f2$gamma[, base]; pb3 <- 1 - f3$gamma[, as.integer(names(which.max(table(f3$path[seq_len(15)]))))]; tr3 <- trans(f3$path, b$t); tr2 <- trans(f2$path, b$t)
    ang <- if (!is.null(job$angle)) approx(job$t, job$angle, xout = b$t, rule = 2)$y else NULL; tilt <- if (!is.null(ang)) ang > 20 else NULL
    auc <- function(score, lab) { if (is.null(lab) || length(unique(lab)) < 2) return(NA_real_); r <- rank(score); n1 <- sum(lab); n0 <- sum(!lab); (sum(r[lab]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
    on <- if (!is.null(ang)) b$t[which(diff(as.integer(tilt)) == 1) + 1] else numeric(0)
    hit <- if (length(on)) mean(sapply(on, function(o) any(abs(tr2 - o) <= 30) || any(abs(tr3 - o) <= 30))) else NA_real_
    out[[m]] <- data.table(ds = job$ds, id = job$id, mode = m, bins = n, lr = 2 * (f2$ll - f1$ll), K_bic = which.min(bic), t1 = tr3[1], t2 = tr3[2], auc2 = auc(pb2, tilt), auc3 = auc(pb3, tilt), onset_hit = hit, n_tilts = length(on), transitions_per_hour = length(tr3) / (diff(range(b$t)) / 3600)) }
  rbindlist(out) }
jobs <- list(); for (ds in c("fantasia", "prcp", "exercise")) { f <- readRDS(sprintf("results/fit_%s.rds", ds)); for (id in names(f$params)) { s <- f$states[SubjectID == id]; p <- f$params[[id]]
  jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, t = s$Time, xp = s$X_p, xs = s$X_s, vp = s$X_p_var, vs = s$X_s_var, V = c(p$lp * p$sig2, p$sig2), angle = if ("angle" %in% names(s)) s$angle else NULL, rr = s$RR_Interval) } }
log_msg("HMM second stage on %d real recordings", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording")); saveRDS(res, "experiments/hmm_states/_out/real.rds")
cat("\n== Rest (Fantasia, 10 min): recordings with a regime where none was imposed (LR above the null 95th percentile of the simulations) ==\n")
print(res[ds == "fantasia", .(n = .N, above_gate = round(mean(lr > g[mode]), 2), lr_median = round(median(lr), 1), K_bic_2plus = round(mean(K_bic >= 2), 2), transitions_per_hour = round(median(transitions_per_hour), 1)), by = mode], row.names = FALSE)
cat("\n== PRCP: agreement of the regimes with the tilt angle (AUC of the probability of leaving the baseline regime), onsets matched by a transition within 30 s ==\n")
print(res[ds == "prcp", .(n = .N, above_gate = round(mean(lr > g[mode]), 2), K_bic_median = median(K_bic), auc_K2 = round(median(auc2, na.rm = TRUE), 2), auc_K3 = round(median(auc3, na.rm = TRUE), 2), onset_hit = round(mean(onset_hit, na.rm = TRUE), 2), transitions_per_hour = round(median(transitions_per_hour), 1)), by = mode], row.names = FALSE)
# exercise: against the events of the heart-rate detector
ex <- res[ds == "exercise"]; ev <- rbindlist(lapply(jobs[sapply(jobs, `[[`, "ds") == "exercise"], function(j) { e <- events_from_rr(j$rr, "exercise", cfg$segmentation); if (is.null(e)) NULL else data.table(id = j$id, e_on = e$t0[1], e_off = e$t1[1]) }))
m <- merge(ex, ev, by = "id"); m[, `:=`(on_err = t1 - e_on, off_err = t2 - e_off)]
cat("\n== Exercise cohort: HMM transitions against the events of the pipeline's detector (heart-rate based), K = 3 solution ==\n")
print(m[, .(n = .N, above_gate = round(mean(lr > g[mode]), 2), K_bic_3 = round(mean(K_bic == 3), 2), onset_hit_30s = round(mean(abs(on_err) <= 30, na.rm = TRUE), 2), onset_err_median_s = round(median(on_err, na.rm = TRUE), 1), second_transition_hit_60s = round(mean(abs(off_err) <= 60, na.rm = TRUE), 2), second_err_median_s = round(median(off_err, na.rm = TRUE), 1)), by = mode], row.names = FALSE)
