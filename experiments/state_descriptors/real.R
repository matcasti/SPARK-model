# ==============================================================================
# experiments/state_descriptors/real.R -- the STATE descriptors on the real recordings (states and estimation variances of results/fit_*.rds, i.e. production M3).
#   Rscript experiments/state_descriptors/real.R         ->  _out/real.rds and a printed summary
# Exercise (263 rest - exercise - recovery recordings; perturbation = exercise as detected by the pipeline; H = 300 s, dwell 15 s) and PRCP (10 people x up to 6 tilts / stand-ups; H = 100 s,
# the recovery window of the pipeline; dwell 15 s). Specification fixed before looking at the data: smooth 5 s, delta 0.999 (plug-in: quantile of the rest window; draws: chi-square 2).
# Questions: (1) how do the descriptors compare with the recovery times of the pipeline; (2) sensitivity to delta, dwell, horizon; (3) H2-like: do they carry what static baseline HRV does not;
# (4) H3-like: are they reproducible within person across repeated tilts (PRCP: 2 repeats of each kind).
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/state_descriptors/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); options(width = 200)
H_EX <- 300; H_TILT <- 100; HOLD <- 15
ex <- readRDS("results/exercise.rds"); fe <- readRDS("results/fit_exercise.rds"); pr <- readRDS("results/prcp.rds"); fp <- readRDS("results/fit_prcp.rds")
static_rest <- function(t, rr, rest) { i <- t >= rest[1] & t < rest[2]; x <- rr[i]; c(hr = 60 / mean(x), sdnn = 1000 * sd(x), rmssd = 1000 * sqrt(mean(diff(x)^2))) }
one <- function(job) { s <- job$s; g <- to_grid1(s$t, s$x_p, s$x_s, s$vp, s$vs); if (!is.na(job$t_max)) g <- g[t <= job$t_max]; V_full <- c(job$par$lp * job$par$sig2, job$par$sig2)
  # scale and time constants from M3 refitted on the resting beats alone (validated in sim.R, variant C): the full-recording fit absorbs the perturbation into sigma^2 (variance inflated 6-7x in exercise simulations)
  k <- s$t >= job$rest[1] & s$t < job$rest[2]; if (sum(k) < 60) return(NULL); sr <- prepare_series(s$rr[k]); fr <- fit_variant(sr, model_spec("M3")); if (!isTRUE(fr$ok)) return(NULL)
  V <- c(sr$priors$energy_ratio * fr$sig2, fr$sig2); tau <- c(fr$tau_p, fr$tau_s)
  pl <- descr_plugin(g, job$rest, job$pert, job$H, HOLD); dr <- summ_draws(descr_draws(g, job$rest, job$pert, V, tau, n_draw = 100, H = job$H, hold = HOLD, seed = job$seed)); if (is.null(pl) || is.null(dr)) return(NULL)
  g_pt <- if (job$kind == "exercise") g else g[t >= job$rest[1]]                  # tilts: only the event window (baseline, tilt, recovery), not the earlier events of the recording
  pt <- ptransition(g_pt, job$rest, job$pert, V, job$H); if (is.null(pt)) pt <- c(K = NA_real_, trans = NA_real_, ptrans = NA_real_)
  grid <- NULL
  if (job$grid) { grid <- rbindlist(lapply(c(0.99, 0.999), function(q) rbindlist(lapply(c(5, 15, 30), function(h) rbindlist(lapply(c(120, 300, 600), function(H) { a <- descr_plugin(g, job$rest, job$pert, H, h, q = q)
      data.table(q = q, hold = h, H = H, Dmax = a["Dmax"], Trec = a["Trec"], Preturn = a["Preturn"]) }))))))
    grid[, id := job$id] }
  st <- static_rest(s$t, s$rr, job$rest)
  list(row = data.table(id = job$id, kind = job$kind, H = job$H, rest_s = diff(job$rest), t(setNames(pl, paste0(names(pl), "_pl"))), t(setNames(dr, paste0(names(dr), "_dr"))), t(pt), t(st), sd_p = s$sd_p, sd_s = s$sd_s, V_p = V[1], V_s = V[2], V_p_full = V_full[1], V_s_full = V_full[2], tau_p = tau[1], tau_s = tau[2]), grid = grid) }
rows <- function(res, what) rbindlist(lapply(Filter(Negate(is.null), res), `[[`, what), fill = TRUE)
# ---- exercise ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
jobs <- lapply(names(ex$series_all), function(id) { s <- ex$series_all[[id]]; list(id = id, kind = "exercise", s = s, par = fe$params[[id]], rest = c(cfg$exercise$rest_skip_s, s$t_on - 5), pert = c(s$t_on, s$t_end), H = H_EX,
  t_max = NA_real_, seed = id_seed(paste("stex", id), cfg$seed), grid = TRUE) })
log_msg("Exercise: %d recordings", length(jobs)); rex <- map_recordings(jobs, one, n_cores); ex_tab <- rows(drop_failed(rex, "exercise"), "row"); ex_grid <- rows(drop_failed(rex, "exercise"), "grid")
# ---- PRCP tilts ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tl <- pr$tilts[usable == TRUE]; jobs <- lapply(seq_len(nrow(tl)), function(i) { e <- tl[i]; id <- e$SubjectID; list(id = paste(id, e$event, sep = "_"), subject = id, kind = e$kind, s = pr$series_all[[id]], par = fp$params[[id]], rest = c(e$t_base_from, e$t_onset - 5),
  pert = c(e$t_onset, e$t_tilt_end), H = H_TILT, t_max = e$t_rec_to, seed = id_seed(paste("sttl", id, e$event), cfg$seed), grid = FALSE) })
log_msg("PRCP tilts: %d events", length(jobs)); rtl <- map_recordings(jobs, one, n_cores); tl_tab <- rows(drop_failed(rtl, "tilt"), "row"); tl_tab[, subject := sub("_.*", "", id)]
saveRDS(list(exercise = ex_tab, grid = ex_grid, tilt = tl_tab, H_EX = H_EX, H_TILT = H_TILT), "experiments/state_descriptors/_out/real.rds")

# ---- summaries ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Real data have no known truth: the descriptors are summarised, compared with the pipeline's recovery times, and checked for redundancy with static baseline HRV and for within-person reproducibility.
med <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x, na.rm = TRUE), quantile(x, .1, na.rm = TRUE), quantile(x, .9, na.rm = TRUE))
sp <- function(a, b) { k <- is.finite(a) & is.finite(b); if (sum(k) < 5) NA_real_ else round(cor(a[k], b[k], method = "spearman"), 2) }
icc1 <- function(x, g) { d <- data.table(x, g)[is.finite(x)]; d <- d[g %in% d[, .N, by = g][N >= 2, g]]; if (uniqueN(d$g) < 3) return(NA_real_); a <- anova(lm(x ~ factor(g), d)); k <- mean(d[, .N, by = g]$N); (a$`Mean Sq`[1] - a$`Mean Sq`[2]) / (a$`Mean Sq`[1] + (k - 1) * a$`Mean Sq`[2]) }
cat("\n== Exercise ==\n"); cat(sprintf("n = %d; rest window %.0f s (median); follow-up after the end of exercise %.0f s (median); follow-up >= H = %d s: %d recordings\n", nrow(ex_tab), median(ex_tab$rest_s), median(ex_tab$follow_s_pl), H_EX, sum(ex_tab$follow_s_pl >= H_EX)))
cat("plug-in (filtered rest SDs): Dmax", med(ex_tab$Dmax_pl), "| delta", med(ex_tab$delta_pl), "| returned within H:", sprintf("%.2f", mean(ex_tab$Preturn_pl, na.rm = TRUE)), "of", sum(is.finite(ex_tab$Preturn_pl)), "\n")
cat("draws   (stationary SDs):   Dmax", med(ex_tab$Dmax_dr), "| Trec (s)", med(ex_tab$Trec_dr), "| share of draws returned within H", med(ex_tab$Preturn_dr), "| rest exceedance of delta", signif(mean(ex_tab$rest_exceed_dr), 2), "\n")
rc <- ex$recovery[definition == ex$summary$recovery_definition & status == "recovered", .(id = SubjectID, state, time)]; w <- dcast(rc, id ~ state, value.var = "time"); m <- merge(ex_tab, w, by = "id")
cat(sprintf("Spearman, Trec (draws, 2D) with the pipeline's recovery time (%s) of x_p: %s, x_s: %s (n = %d)\n", ex$summary$recovery_definition, sp(m$Trec_dr, m$p), sp(m$Trec_dr, m$s), nrow(m)))
cat("\nSensitivity (plug-in; median Trec [s], share returned) by delta quantile, dwell and horizon:\n"); print(ex_grid[, .(Trec = round(median(Trec, na.rm = TRUE)), returned = round(mean(Preturn, na.rm = TRUE), 2), n = sum(is.finite(Preturn))), by = .(q, hold, H)][order(q, hold, H)], row.names = FALSE)
cat("\nRedundancy with static baseline HRV (Spearman with the rest-window values) and R^2 of a linear model on the three static indices (ranks):\n")
st <- c("hr", "sdnn", "rmssd", "sd_p", "sd_s"); for (v in c("Dmax_dr", "Trec_dr", "Dmax_pl", "Trec_pl")) cat(sprintf("%-8s", v), paste(sprintf("%s %.2f", st, sapply(st, function(s) sp(ex_tab[[v]], ex_tab[[s]]))), collapse = " | "), sprintf("| R2 %.2f\n", { d <- ex_tab[is.finite(get(v))]; if (nrow(d) < 20) NA else summary(lm(rank(d[[v]]) ~ rank(d$hr) + rank(d$sdnn) + rank(d$rmssd)))$r.squared }))
cat("\n== PRCP tilts and stand-ups ==\n"); cat(sprintf("%d events, %d people; follow-up %.0f s (median); returned within H = %d s: plug-in %.2f, draws (median share) %.2f\n", nrow(tl_tab), uniqueN(tl_tab$subject), median(tl_tab$follow_s_pl), H_TILT, mean(tl_tab$Preturn_pl, na.rm = TRUE), median(tl_tab$Preturn_dr, na.rm = TRUE)))
print(tl_tab[, .(n = .N, Dmax_draws = round(median(Dmax_dr), 1), Trec_draws_s = round(median(Trec_dr, na.rm = TRUE)), Dmax_plug = round(median(Dmax_pl), 1), Trec_plug_s = round(median(Trec_pl, na.rm = TRUE)), returned_plug = round(mean(Preturn_pl, na.rm = TRUE), 2)), by = kind], row.names = FALSE)
cat("Reproducibility (ICC(1) across the two repeats of each kind within person; 10 people):\n"); print(rbindlist(lapply(c("Dmax_dr", "Dmax_pl", "Trec_dr", "d_H_dr"), function(v) tl_tab[, .(descriptor = v, ICC = round(icc1(get(v), subject), 2), n_people = uniqueN(subject[is.finite(get(v))])), by = kind])), row.names = FALSE)
cat("Same ICC across all events of a person, kind removed by centring:\n"); print(sapply(c("Dmax_dr", "Dmax_pl", "Trec_dr", "d_H_dr"), function(v) { d <- copy(tl_tab)[is.finite(get(v))]; d[, y := get(v) - median(get(v)), by = kind]; round(icc1(d$y, d$subject), 2) }))

cat("\n== Ptransition (HMM regimes as persistent regions; regime at the horizon different from the resting regime) ==\n")
cat(sprintf("exercise: regimes by BIC K = 1 / 2 / 3: %s; Ptransition at H = %d s: %.2f of %d recordings with follow-up (mean posterior probability %.2f)\n", paste(table(factor(ex_tab$K, 1:3)), collapse = " / "), H_EX, mean(ex_tab$trans, na.rm = TRUE), sum(is.finite(ex_tab$trans)), mean(ex_tab$ptrans, na.rm = TRUE)))
cat("agreement with the distance-based return (draws: share of draws returned within H > 0.5):\n"); print(table(transition = ex_tab$trans, returned = ex_tab$Preturn_dr > 0.5, useNA = "no"))
cat("PRCP:\n"); print(tl_tab[, .(n = .N, K2plus = round(mean(K > 1, na.rm = TRUE), 2), Ptransition = round(mean(trans, na.rm = TRUE), 2), n_with_follow = sum(is.finite(trans))), by = kind], row.names = FALSE)
