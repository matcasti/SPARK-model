# EXPERIMENTAL. Persistence of the regimes: the stay-logit prior (mean, sd) of the 2-regime HMM, on PRCP (annotated onsets; long recordings, where the HMM chattered) and on simulations.
# Run from the project root:  Rscript experiments/regimes/run_tune2.R
source("analysis/00_setup.R"); source("experiments/regimes/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); ev_ann <- fread("data/prcp_events.csv"); GATE <- 9
VARS <- list(`7_sd2` = c(7, 2), `7_sd0.7` = c(7, 0.7), `9_sd0.7` = c(9, 0.7), `11_sd0.5` = c(11, 0.5))
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
prcp_job <- function(j) { series <- prepare_series(j$rr); dy <- series$dy; n <- length(dy); tt <- cumsum(dy); ctl <- series$control; gm <- gate_mult(dy, ctl); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf)
  an <- ev_ann[as.character(id) == j$id & is.finite(onset_s), onset_s] - j$off; ll0 <- base$run()$LL_star - sum(gm < ctl$lik_gate_max) / 2 * (log(2 * pi) + 1)
  n_head <- floor(0.7 * n); sh <- prepare_series(dy[seq_len(n_head)]); bh <- fit_ext(sh); cph <- cp_search(sh, bh, min_gain = -Inf); idx <- (n_head + 1):n; idx <- idx[gm[idx] < ctl$lik_gate_max]; hr <- bh$run()
  m3t <- mean(ext_run(dy, ctl, bh$p$kp, bh$p$ks, sh$priors$energy_ratio, bh$p$lR, gm = gm, full = TRUE, beta_fixed = hr$nu0, sig2_fixed = hr$sig2)$logdens[idx])
  rbindlist(lapply(names(VARS), function(v) { h <- fit_regime(series, 2, cp, base, VARS[[v]][1], VARS[[v]][2]); if (!isTRUE(h$ok)) return(NULL); r <- h$run(); on <- regime_onsets2(r$W, tt); hh <- fit_regime(sh, 2, cph, bh, VARS[[v]][1], VARS[[v]][2])
    data.table(variant = v, id = j$id, dLL = r$LL - ll0, n_annot = length(an), hits = if (length(on)) sum(vapply(an, function(a) min(abs(on - a)) <= 30, TRUE)) else 0L, detections = length(on), false_alarms = if (length(on)) sum(vapply(on, function(o) all(abs(an - o) > 30), TRUE)) else 0L,
               held = if (isTRUE(hh$ok)) mean(hh$run(dy, gm)$logdens[idx]) - m3t else NA_real_) })) }
sim_job <- function(job) { z <- sim_regime(job$scenario, job$seed); s <- z$s; series <- prepare_series(s$rr); tt <- cumsum(series$dy); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); e <- state_estimates(base, FALSE)
  rows <- list(data.table(variant = "M3", rp = score(e$xp, e$xpv, s$x_p)[["rmse"]], rs = score(e$xs, e$xsv, s$x_s)[["rmse"]], on = 0))
  if (!is.null(cp) && cp$gain >= GATE) for (v in names(VARS)) { h <- fit_regime(series, 2, cp, base, VARS[[v]][1], VARS[[v]][2]); if (!isTRUE(h$ok)) next; r <- h$run(); rows[[length(rows) + 1]] <- data.table(variant = v, rp = score(r$X_p, r$X_p_var, s$x_p)[["rmse"]], rs = score(r$X_s, r$X_s_var, s$x_s)[["rmse"]], on = length(regime_onsets2(r$W, tt))) }
  rbindlist(rows)[, `:=`(scenario = job$scenario, sim = job$sim, n_true = if (is.null(z$ev)) 0L else nrow(z$ev))] }
d <- load_dataset("prcp", cfg); pj <- lapply(names(d$rr), function(id) list(id = id, rr = d$rr[[id]], off = d$rec[[id]]$time[1] - d$rec[[id]]$rr[1]))
pr <- rbindlist(drop_failed(map_recordings(pj, prcp_job, n_cores), "recording")); saveRDS(pr, "experiments/regimes/_out/tune2_prcp.rds")
cat("\n== PRCP (10 recordings, 62 annotated onsets): persistence prior ==\n"); print(pr[, .(fits = .N, dLL = round(median(dLL), 1), hits = sum(hits), hit_pct = round(100 * sum(hits) / sum(n_annot)), detections = sum(detections), false_alarms = sum(false_alarms), held_gain = round(mean(held, na.rm = TRUE), 3), held_better = sum(held > 0, na.rm = TRUE)), by = variant][order(match(variant, names(VARS)))], row.names = FALSE)
sj <- unlist(lapply(c("regime", "exercise", "two", "rest"), function(sc) lapply(seq_len(40), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("regimes", sc, i), cfg$seed)))), recursive = FALSE)
sr <- rbindlist(drop_failed(map_recordings(sj, sim_job, n_cores), "simulation")); saveRDS(sr, "experiments/regimes/_out/tune2_sim.rds")
m3 <- sr[variant == "M3", .(scenario, sim, rp0 = rp, rs0 = rs)]; x <- merge(sr[variant != "M3"], m3, by = c("scenario", "sim"))
cat("\n== Simulations (gate 9): RMSE vs M3 (%) and detected onsets (true: regime 1, exercise 1, two 2) ==\n"); print(x[, .(n = .N, vagal = round(100 * (mean(rp) / mean(rp0) - 1)), symp = round(100 * (mean(rs) / mean(rs0) - 1)), onsets = round(mean(on), 2), true = round(mean(n_true), 1)), by = .(scenario, variant)][order(scenario, match(variant, names(VARS)))], row.names = FALSE)
