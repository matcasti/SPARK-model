# ==============================================================================
# experiments/m3plus/compare/go_nogo.R -- does M3 with the respiratory oscillator (M3O) become the default? The decision rule of the comparison with the log-RR candidate, applied to
# M3O against production M3 (filtered states, what production reports; M3O+ against M3+ for the smoothed ones).   Rscript experiments/m3plus/compare/go_nogo.R
# Sources: the simulations and real fits already stored in _out (compare_sim_*.rds, real_full.rds, compare_real_full.rds); new here: S3 (exercise and tilt change ratios) and
# the test-retest and known-group checks of Fantasia (windows A and B).
# Rule: S1 mean standardised RMSE not above M3's by more than 10% and correlation not lower by more than 0.03 (both states); S2 95% coverage in 0.90-0.98 or not lower than M3's by more
# than 0.03; S3 (exercise, tilt) change of each state rest -> end of event, estimated / true, within 0.67-1.5 or not farther from 1 (in ln) than M3's by more than 0.1;
# S4 median |ln error| of tau_p and tau_s not above M3's by more than 0.1. R1 KS pass not lower by more than 0.05, PACF violations not above +1, tau_p inside the production prior
# not lower by more than 0.10; R2 forward score (held-out nats per beat) difference with 95% bootstrap CI above -0.01; R3 PRCP tilts with both expected signs not lower by more
# than 0.05; R4 test-retest ICC of SD_p, SD_s and tau_p not lower by more than 0.05 and young-vs-older AUC of SD_p not lower by more than 0.03.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R"); options(width = 200)
n_cores <- max(1L, min(5L, cfg$cores - 2L)); out <- "experiments/m3plus/_out"
# ---- S1, S2, S4 from the stored simulations --------------------------------------------------------------------------------
sim <- rbindlist(lapply(c("stationary", "rsa", "exercise", "tilt"), function(s) readRDS(sprintf("%s/compare_sim_%s.rds", out, s))), fill = TRUE); sim[, group := sub(":.*", "", scenario)]
agg <- sim[model %in% c("M3", "M3O", "M3+", "M3O+"), .(rmse = mean(rmse), cor = mean(cor), cover = mean(cover)), by = .(group, model, state)]
tau <- sim[model %in% c("M3", "M3O", "M3+", "M3O+") & state == "vagal", .(ln_tp = median(abs(tau_p_err)), ln_ts = median(abs(tau_s_err))), by = .(group, model)]
ok_cov <- function(c, c3) (c >= 0.90 && c <= 0.98) || c >= c3 - 0.03
res <- rbindlist(lapply(unique(agg$group), function(g) rbindlist(lapply(list(c("M3", "M3O"), c("M3+", "M3O+")), function(pr) { a <- agg[group == g & model == pr[1]]; b <- agg[group == g & model == pr[2]]
  s1 <- all(sapply(c("vagal", "sympathetic"), function(st) b[state == st, rmse] <= 1.10 * a[state == st, rmse] && b[state == st, cor] >= a[state == st, cor] - 0.03))
  s2 <- all(sapply(c("vagal", "sympathetic"), function(st) ok_cov(b[state == st, cover], a[state == st, cover])))
  ta <- tau[group == g & model == pr[1]]; tb <- tau[group == g & model == pr[2]]; s4 <- tb$ln_tp <= ta$ln_tp + 0.1 && tb$ln_ts <= ta$ln_ts + 0.1
  data.table(group = g, comparison = paste(pr[2], "vs", pr[1]), S1 = s1, S2 = s2, S4 = s4, rmse_vagal = sprintf("%.2f -> %.2f", a[state == "vagal", rmse], b[state == "vagal", rmse]), rmse_symp = sprintf("%.2f -> %.2f", a[state == "sympathetic", rmse], b[state == "sympathetic", rmse]),
             cover = sprintf("%.2f/%.2f -> %.2f/%.2f", a[state == "vagal", cover], a[state == "sympathetic", cover], b[state == "vagal", cover], b[state == "sympathetic", cover])) }))))
cat("== S1, S2, S4 on the stored simulations (80 per scenario; stationary: the pipeline's four scenarios) ==\n"); print(res[order(comparison, group)], row.names = FALSE)

# ---- S3: exercise and tilt, change of each state from rest to the end of the event ----------------------------------------
DUR <- 420; osc_draw <- function() list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1))))
build <- function(sc, seed) { set.seed(seed); d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))), osc = osc_draw())
  mk <- function(mp, ms) sim_ext(DUR, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = d$osc, m_p_fun = mp, m_s_fun = ms)
  if (sc == "exercise") { p <- draw_exercise_pars(seed, DUR); return(c(d, list(true = c(p$t_on, p$t_off), s = mk(prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s))))) }
  t_on <- runif(1, 0.28, 0.45) * DUR; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
  c(d, list(true = c(t_on, t_off), s = mk(prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))), prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100))))))) }
chg <- function(x, t, t0, t1) { pre <- t >= t0 - 120 & t < t0 - 10; en <- t >= t1 - 15 & t < t1; if (sum(pre) < 5 || sum(en) < 5) NA_real_ else mean(x[en]) - mean(x[pre]) }
t50 <- function(x, t, t0, t1) { pre <- t >= t0 - 120 & t < t0 - 10; en <- t >= t1 - 15 & t < t1; if (sum(pre) < 5 || sum(en) < 5) return(NA_real_); a <- mean(x[pre]); b <- mean(x[en]); k <- which(t >= t0 & (x - a) / (b - a) >= 0.5)[1]; if (is.na(k)) NA_real_ else t[k] - t0 }
one <- function(job) { z <- build(job$scenario, job$seed); s <- z$s; tp <- s$x_p + s$z1; ts <- s$x_s; series <- prepare_series(s$rr); tt <- s$t; rows <- list()
  for (m in c("M3", "M3O")) { f <- fit_ext(series, if (m == "M3O") list(osc = TRUE) else list()); if (!isTRUE(f$ok)) next; e <- state_estimates(f, FALSE)
    rows[[m]] <- data.table(scenario = job$scenario, sim = job$sim, model = m, chg_p = chg(e$xp, tt, z$true[1], z$true[2]) / chg(tp, tt, z$true[1], z$true[2]), chg_s = chg(e$xs, tt, z$true[1], z$true[2]) / chg(ts, tt, z$true[1], z$true[2]),
      vagal_first = t50(e$xp, tt, z$true[1], z$true[2]) < t50(e$xs, tt, z$true[1], z$true[2]), true_vagal_first = t50(tp, tt, z$true[1], z$true[2]) < t50(ts, tt, z$true[1], z$true[2])) }
  rbindlist(rows) }
jobs <- unlist(lapply(c("exercise", "tilt"), function(sc) lapply(1:60, function(i) list(scenario = sc, sim = i, seed = id_seed(paste("gonogo", sc, i), cfg$seed)))), recursive = FALSE)
s3 <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(s3, sprintf("%s/go_nogo_s3.rds", out))
cat("\n== S3: change from rest to the end of the event, estimated / true (median); vagal reaches half of its change first (share; truth in brackets) ==\n")
print(s3[, .(chg_vagal = round(median(chg_p, na.rm = TRUE), 2), chg_symp = round(median(chg_s, na.rm = TRUE), 2), vagal_first = sprintf("%.2f (%.2f)", mean(vagal_first, na.rm = TRUE), mean(true_vagal_first, na.rm = TRUE))), by = .(scenario, model)][order(scenario, model)], row.names = FALSE)

# ---- Real data: stored fits ------------------------------------------------------------------------------------------------
rf <- readRDS(sprintf("%s/real_full.rds", out)); b <- rf[variant == "base"]; o <- rf[variant == "osc"]
cat("\n== R1 diagnostics and tau_p, 313 recordings ==\n"); print(rf[variant %in% c("base", "osc"), .(KS_pass = round(mean(ks_ok), 3), PACF_violations = round(mean(pacf_viol), 1), tau_p_median_s = round(median(tau_p), 1), tau_p_in_prior = round(mean(tau_p >= 3 & tau_p <= 6), 2)), by = variant], row.names = FALSE)
m <- merge(b[, .(ds, id, t0 = tail_logdens)], o[, .(ds, id, t1 = tail_logdens)], by = c("ds", "id")); m <- m[is.finite(t0) & is.finite(t1)]; set.seed(1); d <- m$t1 - m$t0
bsd <- replicate(2000, mean(sample(d, replace = TRUE))); cat(sprintf("\n== R2 forward score (held-out nats per beat), M3O minus M3: %.4f [%.4f, %.4f]; better in %.0f%% of %d recordings\n", mean(d), quantile(bsd, .025), quantile(bsd, .975), 100 * mean(d > 0), length(d)))
print(m[, .(difference = round(mean(t1 - t0), 4)), by = ds], row.names = FALSE)
cr <- readRDS(sprintf("%s/compare_real_full.rds", out)); cat("\n== R3 PRCP tilts with both expected signs, and the exercise responses (rest SDs, median) ==\n")
print(cr$tp[model %in% c("M3", "M3O"), .(tilts = .N, expected = round(mean(expected), 3), z_p = round(median(z_p), 1), z_s = round(median(z_s), 1)), by = model], row.names = FALSE)
print(cr$ex[model %in% c("M3", "M3O"), .(recordings = .N, z_p = round(median(z_p, na.rm = TRUE), 1), z_s = round(median(z_s, na.rm = TRUE), 1)), by = model], row.names = FALSE)

# ---- R4: test-retest and known groups, Fantasia -------------------------------------------------------------------------------
w <- cfg$windows$fantasia; dataset <- load_dataset("fantasia", cfg); ids <- names(dataset$rec_full); grp <- fantasia_group(ids, cfg); names(grp) <- ids
fitw <- function(i) { win <- if (jobs_w$win[i] == "A") w$stability_A else w$stability_B; rr <- window_rec(dataset$rec_full[[jobs_w$id[i]]], start_s = win[1], duration_s = win[2] - win[1])$rr; series <- prepare_series(rr)
  rbindlist(lapply(c("M3", "M3O"), function(m) { f <- fit_ext(series, if (m == "M3O") list(osc = TRUE) else list()); if (!isTRUE(f$ok)) return(NULL); r <- f$run(full = TRUE)
    data.table(id = jobs_w$id[i], win = jobs_w$win[i], model = m, tau_p = f$tau_p, tau_s = f$tau_s, SD_p = sqrt(f$lp * r$sig2), SD_s = sqrt(r$sig2)) })) }
jobs_w <- CJ(id = ids, win = c("A", "B")); fw <- rbindlist(drop_failed(map_recordings(seq_len(nrow(jobs_w)), fitw, n_cores), "window")); saveRDS(fw, sprintf("%s/go_nogo_fantasia.rds", out))
icc <- function(a, b) { x <- cbind(a, b); n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1); mse <- sum((x - rowMeans(x) - matrix(colMeans(x) - gm, n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
auc <- function(x, g) { r <- rank(x); n1 <- sum(g == "young"); n2 <- sum(g == "older"); (sum(r[g == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
r4 <- rbindlist(lapply(c("M3", "M3O"), function(m) { a <- fw[model == m & win == "A"][order(id)]; bb <- fw[model == m & win == "B"][order(id)]; k <- intersect(a$id, bb$id); a <- a[id %in% k]; bb <- bb[id %in% k]
  rbindlist(lapply(c("tau_p", "tau_s", "SD_p", "SD_s"), function(v) data.table(model = m, variable = v, icc = icc(log(a[[v]]), log(bb[[v]])), auc_young_older = auc(a[[v]], grp[a$id])))) }))
cat("\n== R4 Fantasia test-retest ICC(3,1) of the log values (windows A and B) and young-vs-older AUC ==\n"); print(dcast(r4, variable ~ model, value.var = c("icc", "auc_young_older"))[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE)
