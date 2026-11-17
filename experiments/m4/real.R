# ==============================================================================
# experiments/m4/real.R -- EXPERIMENTAL. The M4 ladder on real recordings (README.md, "Real recordings"): 60 exercise (the sample of regimes/v2), 40 Fantasia, 10 PRCP.
#   Rscript experiments/m4/real.R
# Parts per recording: "full" (in-sample diagnostics, time constants, exercise realism), "head" (first 70% of the beats; scored on the last 30%),
# and for Fantasia the stability windows "A" and "B" (test-retest). Stage 2 = M4d-H with empirical-Bayes priors from the other recordings of the same dataset and part.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/m4/lib.R"); n_cores <- max(1L, min(6L, cfg$cores - 1L)); out <- "experiments/m4/_out"
ann <- fread("data/prcp_events.csv"); jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- names(d$rr); if (ds == "exercise") { set.seed(21); ids <- sort(sample(ids, 60)) }
  for (id in ids) { rr <- d$rr[[id]]; ev <- NULL
    if (ds == "prcp") { off <- d$rec[[id]]$time[1] - d$rec[[id]]$rr[1]; a <- ann[as.character(ann$id) == id & is.finite(ann$onset_s) & is.finite(ann$return_s)]; if (nrow(a)) ev <- data.frame(t0 = a$onset_s - off, t1 = a$return_s - off) }
    for (part in c("full", "head")) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, part = part, rr = rr, ev = ev) }
  if (ds == "fantasia") { w <- cfg$windows$fantasia; for (id in ids) for (part in c("A", "B")) { win <- if (part == "A") w$stability_A else w$stability_B
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, part = part, rr = window_rec(d$rec_full[[id]], start_s = win[1], duration_s = win[2] - win[1])$rr, ev = NULL) } } }

t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2
  k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
run_job <- function(job, eb = NULL) {
  dy_full <- prepare_series(job$rr)$dy; n <- length(dy_full); n_head <- floor(0.7 * n)
  series <- prepare_series(if (job$part == "head") dy_full[seq_len(n_head)] else dy_full); dy <- series$dy; ctl <- series$control; gm <- gate_mult(dy, ctl); tt <- cumsum(dy)
  fam <- fit_family(series, eb = eb, events = job$ev, stage1 = is.null(eb)); keep <- if (is.null(eb)) setdiff(intersect(MODELS, names(fam)), "M4dH") else intersect("M4dH", names(fam))
  ex_ev <- if (job$ds == "exercise" && job$part == "full") events_from_rr(dy, "exercise", cfg$segmentation) else NULL
  gmf <- gate_mult(dy_full, ctl); idx <- (n_head + 1):n; idx <- idx[gmf[idx] < ctl$lik_gate_max]
  rows <- rbindlist(lapply(keep, function(nm) { m <- fam[[nm]]; dg <- m_diag(m, gm, ctl); sd2 <- m_sd(m)
    row <- data.table(ds = job$ds, id = job$id, part = job$part, model = nm, ks_ok = dg[["ks"]], pacf = dg[["pacf"]], tau_p = m$fit$tau_p, tau_s = m$fit$tau_s, SD_p = sd2[["SD_p"]], SD_s = sd2[["SD_s"]],
                      lp = m$lp, gate_open = isTRUE(m$gate_open))
    if (job$part == "head") row[, `:=`(heldout = tryCatch(m_heldout(m, dy_full, gmf, ctl, idx), error = function(e) NA_real_), tail_bad = mean(dy_full[idx] > 2))]
    if (!is.null(ex_ev)) { e <- m_states(m); t_on <- ex_ev$t0[1]; t_end <- ex_ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) row[, `:=`(t50_vagal = t50(e$xp, tt, t_on, mean(e$xp[rest]), mean(e$xp[endw])), t50_symp = t50(e$xs, tt, t_on, mean(e$xs[rest]), mean(e$xs[endw])),
                                                     d_vagal = mean(e$xp[endw]) - mean(e$xp[rest]), d_symp = mean(e$xs[endw]) - mean(e$xs[rest]))] }
    row }), fill = TRUE)
  list(rows = rows, weak = if (is.null(eb)) data.table(ds = job$ds, id = job$id, part = job$part, t(fam$weak)) else NULL)
}
log_msg("M4 real, stage 1: %d jobs", length(jobs)); s1 <- drop_failed(map_recordings(jobs, run_job, n_cores), "job")
R1 <- rbindlist(lapply(s1, `[[`, "rows"), fill = TRUE); W <- rbindlist(lapply(s1, `[[`, "weak"))
log_msg("Stage 2 (empirical-Bayes priors)")
R2 <- rbindlist(lapply(drop_failed(map_recordings(jobs, function(job) { w <- W[ds == job$ds & part == job$part & id != job$id]; run_job(job, eb = eb_from(as.matrix(w[, .(tau_p, tau_s)]))) }, n_cores), "job"), `[[`, "rows"), fill = TRUE)
R <- rbind(R1, R2, fill = TRUE); R[, model := factor(model, MODELS)]; saveRDS(list(res = R, weak = W), file.path(out, "real.rds"))

f <- R[part == "full"]
cat("\n== In-sample (full records) ==\n")
print(f[, .(n = .N, ks_pass = round(100 * mean(ks_ok), 1), pacf = round(mean(pacf), 1), tau_p_med = round(median(tau_p), 1), tau_p_in = round(100 * mean(tau_p >= 3 & tau_p <= 6)), tau_s_med = round(median(tau_s)),
           lp_med = round(median(lp), 2), gate = round(100 * mean(gate_open))), by = .(ds, model)][order(ds, model)], row.names = FALSE)
cat("\n== Exercise realism (rest -> end of exercise; filtered states) ==\n")
print(f[ds == "exercise" & is.finite(d_vagal), .(n = .N, signs = round(100 * mean(d_vagal < 0 & d_symp > 0)), vagal_first = round(100 * mean(t50_vagal < t50_symp, na.rm = TRUE)),
                                                t50_vagal = round(median(t50_vagal, na.rm = TRUE)), t50_symp = round(median(t50_symp, na.rm = TRUE))), by = model][order(model)], row.names = FALSE)
h <- R[part == "head" & tail_bad < 0.1]; b <- h[model == "M3", .(ds, id, m3 = heldout)]; h <- merge(h, b, by = c("ds", "id")); h[, d := heldout - m3]
boot_lo <- function(x) { x <- x[is.finite(x)]; if (length(x) < 5) return(NA_real_); quantile(with_seed(7L, replicate(2000, mean(sample(x, replace = TRUE)))), 0.025) }
cat("\n== Held-out score minus M3 (nats per beat; last 30%; degenerate tails excluded) ==\n")
print(h[model != "M3", .(n = sum(is.finite(d)), mean = round(mean(d, na.rm = TRUE), 3), lo95 = round(boot_lo(d), 3), median = round(median(d, na.rm = TRUE), 3), better = round(100 * mean(d > 0, na.rm = TRUE))), by = model][order(model)], row.names = FALSE)
print(h[model != "M3", .(mean = round(mean(d, na.rm = TRUE), 3), better = round(100 * mean(d > 0, na.rm = TRUE))), by = .(ds, model)][order(ds, model)], row.names = FALSE)
icc <- function(a, b) { x <- cbind(a, b); x <- x[stats::complete.cases(x), , drop = FALSE]; n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1)
  mse <- sum((x - rowMeans(x) - matrix(colMeans(x) - gm, n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
tw <- dcast(R[part %in% c("A", "B")], id + model ~ part, value.var = c("tau_p", "tau_s", "SD_p", "SD_s"))
cat("\n== Fantasia test-retest ICC(3,1) of the log values (windows A and B) ==\n")
print(tw[, .(n = .N, tau_p = round(icc(log(tau_p_A), log(tau_p_B)), 2), tau_s = round(icc(log(tau_s_A), log(tau_s_B)), 2), SD_p = round(icc(log(SD_p_A), log(SD_p_B)), 2), SD_s = round(icc(log(SD_s_A), log(SD_s_B)), 2)), by = model][order(model)], row.names = FALSE)
