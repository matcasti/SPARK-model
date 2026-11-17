# ==============================================================================
# experiments/m4/e/real_e.R -- EXPERIMENTAL. M4e components on real recordings (README.md, data split).
#   Rscript experiments/m4/e/real_e.R dev|conf [tag]
# dev: the 60 exercise recordings of ../real.R, Fantasia full/head/A/B, PRCP full/head. conf: 60 other exercise recordings, Fantasia full/head/C/D.
# Per component: KS pass, PACF violations, time constants, stationary SDs, held-out score of the last 30% (head fits), exercise realism (full exercise records).
# ==============================================================================
args <- commandArgs(TRUE); phase <- args[1]; tag <- if (length(args) > 1) args[2] else ""; stopifnot(phase %in% c("dev", "conf"))
source("analysis/00_setup.R"); source("experiments/m4/e/lib_e.R"); n_cores <- max(1L, min(7L, cfg$cores)); out <- "experiments/m4/e/_out"; dir.create(out, showWarnings = FALSE)
SPEC_EXTRA <- if (exists("SPEC_EXTRA_E")) SPEC_EXTRA_E else list()
ann <- fread("data/prcp_events.csv"); jobs <- list()
for (ds in cfg$datasets) { if (phase == "conf" && ds == "prcp") next
  d <- load_dataset(ds, cfg); ids <- names(d$rr)
  if (ds == "exercise") { set.seed(21); dev_ids <- sort(sample(ids, 60)); if (phase == "dev") ids <- dev_ids else { set.seed(22); ids <- sort(sample(setdiff(ids, dev_ids), 60)) } }
  for (id in ids) for (part in c("full", "head")) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, part = part, rr = d$rr[[id]])
  if (ds == "fantasia") { wins <- if (phase == "dev") list(A = c(600, 1200), B = c(1200, 1800)) else list(C = c(3600, 4200), D = c(4200, 4800))
    for (id in ids) for (w in names(wins)) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, part = w, rr = window_rec(d$rec_full[[id]], start_s = wins[[w]][1], duration_s = diff(wins[[w]]))$rr) } }

t50 <- function(x, t, t_on, rest_mean, end_mean) { x <- as.numeric(stats::filter(x, rep(1 / 5, 5), sides = 2)); x[is.na(x)] <- 0; half <- (rest_mean + end_mean) / 2
  k <- which(t >= t_on & sign(end_mean - rest_mean) * (x - half) >= 0)[1]; if (is.na(k)) NA_real_ else t[k] - t_on }
run_job <- function(job) {
  dy_full <- prepare_series(job$rr)$dy; n <- length(dy_full); n_head <- floor(0.7 * n)
  series <- prepare_series(if (job$part == "head") dy_full[seq_len(n_head)] else dy_full); dy <- series$dy; ctl <- series$control; gm <- gate_mult(dy, ctl); tt <- cumsum(dy)
  fc <- fit_components(series, SPEC_EXTRA)
  ex_ev <- if (job$ds == "exercise" && job$part == "full") events_from_rr(dy, "exercise", cfg$segmentation) else NULL
  gmf <- gate_mult(dy_full, ctl); idx <- (n_head + 1):n; idx <- idx[gmf[idx] < ctl$lik_gate_max]
  rows <- rbindlist(lapply(names(fc$models), function(nm) { m <- fc$models[[nm]]; dg <- tryCatch(m_diag(m, gm, ctl), error = function(e) c(ks = NA, pacf = NA)); sd2 <- m_sd(m)
    row <- data.table(ds = job$ds, id = job$id, part = job$part, comp = nm, ks_ok = dg[["ks"]], pacf = dg[["pacf"]], tau_p = m$fit$tau_p, tau_s = m$fit$tau_s, SD_p = sd2[["SD_p"]], SD_s = sd2[["SD_s"]], cost = unname(m$fit$cost))
    if (job$part == "head") row[, `:=`(heldout = tryCatch(m_heldout(m, dy_full, gmf, ctl, idx), error = function(e) NA_real_), tail_bad = mean(dy_full[idx] > 2))]
    if (!is.null(ex_ev)) { e <- m_states(m); t_on <- ex_ev$t0[1]; t_end <- ex_ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end
      if (sum(rest) > 20 && sum(endw) > 5) row[, `:=`(t50_vagal = t50(e$xp, tt, t_on, mean(e$xp[rest]), mean(e$xp[endw])), t50_symp = t50(e$xs, tt, t_on, mean(e$xs[rest]), mean(e$xs[endw])),
                                                     d_vagal = mean(e$xp[endw]) - mean(e$xp[rest]), d_symp = mean(e$xs[endw]) - mean(e$xs[rest]))] }
    row }), fill = TRUE)
  list(rows = rows, stats = data.table(ds = job$ds, id = job$id, part = job$part, lr_osc = fc$stats[["lr_osc"]], gain = fc$stats[["gain"]]))
}
log_msg("M4e real (%s%s): %d jobs", phase, tag, length(jobs)); s1 <- drop_failed(map_recordings(jobs, run_job, n_cores), "job")
res <- list(rows = rbindlist(lapply(s1, `[[`, "rows"), fill = TRUE), stats = rbindlist(lapply(s1, `[[`, "stats")))
saveRDS(res, file.path(out, sprintf("real_%s%s.rds", phase, tag))); log_msg("saved real_%s%s.rds (%d rows)", phase, tag, nrow(res$rows))
