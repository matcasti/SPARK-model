# ==============================================================================
# experiments/m4/sim.R -- EXPERIMENTAL. The M4 ladder on simulated recordings with known truth (README.md, "Simulation with truth").
#   Rscript experiments/m4/sim.R [N per scenario = 40]
# Stage 1 fits M3, M4a, M4b, M4c, M4d-I and the weak-prior time constants; stage 2 refits with the empirical-Bayes priors of the other simulations
# of the same scenario (leave-one-out) for M4d-H.
# ==============================================================================
args <- commandArgs(TRUE); N <- if (length(args)) as.integer(args[1]) else 40L
source("analysis/00_setup.R"); source("experiments/m4/lib.R"); n_cores <- max(1L, min(6L, cfg$cores - 1L)); out <- "experiments/m4/_out"; dir.create(out, showWarnings = FALSE)
OSC <- list(alpha = 0.2, f_hz = 0.25, lo = 0.3)
SCEN <- list(rest = list(gen = "rest", osc = NULL), rsa = list(gen = "rest", osc = OSC), regime = list(gen = "regime", osc = OSC), exercise = list(gen = "exercise", osc = OSC), two = list(gen = "two", osc = OSC))
jobs <- unlist(lapply(names(SCEN), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("m4", sc, i), cfg$seed)))), recursive = FALSE)
make <- function(job) { z <- sim_regime(SCEN[[job$scenario]]$gen, job$seed, osc = SCEN[[job$scenario]]$osc); s <- z$s
  list(z = z, series = prepare_series(s$rr), truth_p = s$x_p + s$z1, truth_s = s$x_s, ev = if (!is.null(z$ev)) data.frame(t0 = z$ev$t0 - z$off, t1 = z$ev$t1 - z$off) else NULL) }
score <- function(m, d, model) { gm <- gate_mult(d$series$dy, d$series$control); dg <- m_diag(m, gm, d$series$control)
  rbindlist(lapply(c(FALSE, TRUE), function(sm) { e <- m_states(m, sm)
    data.table(model = model, mode = if (sm) "smooth" else "filter", rmse_p = sqrt(mean((e$xp - d$truth_p)^2)) / sd(d$truth_p), rmse_s = sqrt(mean((e$xs - d$truth_s)^2)) / sd(d$truth_s),
               cor_p = cor(e$xp, d$truth_p), cor_s = cor(e$xs, d$truth_s), cover_p = mean(abs(e$xp - d$truth_p) <= 1.96 * sqrt(pmax(e$xpv, 0))), cover_s = mean(abs(e$xs - d$truth_s) <= 1.96 * sqrt(pmax(e$xsv, 0))),
               err_tau_p = abs(log(m$fit$tau_p / d$z$tau_p)), err_tau_s = abs(log(m$fit$tau_s / d$z$tau_s)), err_lp = abs(log(m$lp / d$z$lp)), gate_open = isTRUE(m$gate_open), ks_ok = dg[["ks"]], pacf = dg[["pacf"]]) })) }
stage1 <- function(job) { d <- make(job); fam <- fit_family(d$series, events = d$ev)
  rows <- rbindlist(lapply(intersect(MODELS, names(fam)), function(nm) score(fam[[nm]], d, nm)))
  list(rows = rows[, `:=`(scenario = job$scenario, sim = job$sim)], weak = data.table(scenario = job$scenario, sim = job$sim, t(fam$weak))) }
log_msg("M4 simulations, stage 1: %d jobs", length(jobs)); s1 <- drop_failed(map_recordings(jobs, stage1, n_cores), "simulation")
R1 <- rbindlist(lapply(s1, `[[`, "rows")); W <- rbindlist(lapply(s1, `[[`, "weak"))
stage2 <- function(job) { w <- W[scenario == job$scenario & sim != job$sim]; eb <- eb_from(as.matrix(w[, .(tau_p, tau_s)])); d <- make(job)
  fam <- fit_family(d$series, eb = eb, stage1 = FALSE); if (is.null(fam$M4dH)) return(NULL); score(fam$M4dH, d, "M4dH")[, `:=`(scenario = job$scenario, sim = job$sim)] }
log_msg("Stage 2 (empirical-Bayes priors)"); R2 <- rbindlist(drop_failed(map_recordings(jobs, stage2, n_cores), "simulation"))
R <- rbind(R1, R2, fill = TRUE); R[, model := factor(model, MODELS)]; saveRDS(list(res = R, weak = W), file.path(out, "sim.rds"))

rnd <- function(d) d[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]
summ <- function(d) d[, .(n = .N, rmse_p = mean(rmse_p), rmse_s = mean(rmse_s), cor_p = median(cor_p), cor_s = median(cor_s), cover_p = mean(cover_p), cover_s = mean(cover_s),
                          err_tau_p = median(err_tau_p), err_tau_s = median(err_tau_s), err_lp = median(err_lp), gate = mean(gate_open), ks = mean(ks_ok), pacf = mean(pacf)), by = .(scenario, model)][order(scenario, model)]
cat("\n== Filtered states ==\n"); print(rnd(summ(R[mode == "filter"])), row.names = FALSE)
cat("\n== Smoothed states ==\n"); print(rnd(summ(R[mode == "smooth"])[, .(scenario, model, rmse_p, rmse_s, cover_p, cover_s)]), row.names = FALSE)
cat("\n== RMSE relative to M3 (mean of per-simulation ratios - 1, %, filtered) and share of simulations worse than M3 (vagal) ==\n")
b <- R[mode == "filter" & model == "M3", .(scenario, sim, b_p = rmse_p, b_s = rmse_s)]; x <- merge(R[mode == "filter"], b, by = c("scenario", "sim"))
print(x[model != "M3", .(vagal_pct = round(100 * (mean(rmse_p / b_p) - 1)), symp_pct = round(100 * (mean(rmse_s / b_s) - 1)), worse_vagal = round(100 * mean(rmse_p > 1.05 * b_p)), worse_symp = round(100 * mean(rmse_s > 1.05 * b_s))), by = .(scenario, model)][order(scenario, model)], row.names = FALSE)
