# EXPERIMENTAL. M3 vs M3+ (RTS) vs M3O (oscillator) vs M3O+ (oscillator + RTS) vs M3E (events estimated from the RR series) on simulated recordings with known truth.
# Extras flagged as such: M3E_true (true event times, an upper bound for M3E) and M3OE+ (oscillator + events + RTS). Scenarios: stationary (the pipeline's four), rsa, exercise, tilt.
# Truth for the vagal state = OU + oscillator + event response; for the sympathetic state = OU + event response. Run from the project root:  Rscript experiments/m3plus/compare/compare_sim.R [scenario]
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
only <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else NULL
n_cores <- max(1L, min(5L, cfg$cores - 2L)); N <- 80L; DUR <- 420
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cor = cor(est, truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
osc_draw <- function() list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1))))
build <- function(sc, seed) {
  set.seed(seed); d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))), osc = osc_draw())
  mk <- function(osc, mp = NULL, ms = NULL) sim_ext(DUR, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = osc, m_p_fun = mp, m_s_fun = ms)
  if (sc == "rsa") return(c(d, list(s = mk(d$osc), mode = NULL)))
  if (sc == "exercise") { p <- draw_exercise_pars(seed, DUR); return(c(d, list(mode = "exercise", true = c(p$t_on, p$t_off), s = mk(d$osc, prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s))))) }
  if (sc == "tilt") { t_on <- runif(1, 0.28, 0.45) * DUR; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
    return(c(d, list(mode = "tilt", true = c(t_on, t_off), s = mk(d$osc, prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))), prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100)))))))) }
  stop("scenario")
}
one <- function(job) {
  stationary <- job$scenario == "stationary"; z <- if (stationary) list(mode = NULL) else build(job$scenario, job$seed); s <- z$s
  if (!stationary) { series <- prepare_series(s$rr); dy <- series$dy; off <- s$t[1] - s$rr[1] }
  if (stationary) { truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$sub]]()); sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
    series <- prepare_series(sim$rr); dy <- series$dy; tp <- sim$x_p; ts <- sim$x_s; z$tau_p <- truth$tau_p; z$tau_s <- truth$tau_s; z$mode <- NULL
  } else { tp <- s$x_p + s$z1; ts <- s$x_s }
  ev_hat <- if (!is.null(z$mode)) events_from_rr(dy, z$mode, cfg$segmentation) else NULL
  ev_true <- if (!is.null(z$mode)) data.frame(t0 = z$true[1] - off, t1 = z$true[2] - off) else NULL
  if (!is.null(ev_hat) && z$mode == "tilt") ev_hat <- ev_hat[1, , drop = FALSE]                      # the first candidate, as for a single protocol tilt
  fits <- list(M3 = fit_ext(series), M3O = fit_ext(series, list(osc = TRUE)))
  if (!is.null(z$mode)) { fits$M3E <- if (is.null(ev_hat)) fits$M3 else fit_ext(series, list(inputs = list(events = ev_hat))); fits$M3E_true <- fit_ext(series, list(inputs = list(events = ev_true)))
                          fits$M3OE <- fit_ext(series, list(osc = TRUE, inputs = list(events = if (is.null(ev_hat)) ev_true else ev_hat))) }
  cells <- list(M3 = list("M3", FALSE), `M3+` = list("M3", TRUE), M3O = list("M3O", FALSE), `M3O+` = list("M3O", TRUE), M3E = list("M3E", FALSE), M3E_true = list("M3E_true", FALSE), `M3OE+` = list("M3OE", TRUE))
  rows <- list()
  for (nm in names(cells)) { f <- fits[[cells[[nm]][[1]]]]; if (is.null(f) || !isTRUE(f$ok)) next
    e <- state_estimates(f, cells[[nm]][[2]]); k <- which(is.finite(e$xp) & is.finite(e$xs))
    rows[[nm]] <- data.table(model = nm, state = c("vagal", "sympathetic"), rbind(score(e$xp[k], e$xpv[k], tp[k]), score(e$xs[k], e$xsv[k], ts[k])), tau_p_err = log(f$tau_p / z$tau_p), tau_s_err = log(f$tau_s / z$tau_s),
                             LL = e$r$LL_star, n_par = f$n_par + if (!is.null(f$spec$inputs)) 2 * nrow(f$spec$inputs$events) else 0) }
  out <- rbindlist(rows); out[, `:=`(scenario = if (stationary) paste0("stationary:", job$sub) else job$scenario, sim = job$sim,
        detected = !is.null(ev_hat), onset_err = if (!is.null(ev_hat)) ev_hat$t0[1] - ev_true$t0 else NA_real_, end_err = if (!is.null(ev_hat)) ev_hat$t1[1] - ev_true$t1 else NA_real_)]
  out
}
dir.create("experiments/m3plus/_out", showWarnings = FALSE)
scen <- if (is.null(only)) c("stationary", "rsa", "exercise", "tilt") else only; all <- list()
for (sc in scen) {
  jobs <- if (sc == "stationary") unlist(lapply(names(RECOVERY_SCENARIOS), function(sub) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = "stationary", sub = sub, sim = i, seed = id_seed(paste(sub, i), cfg$seed)))), recursive = FALSE)
         else lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("compare", sc, i), cfg$seed)))
  log_msg("Scenario %s: %d simulated recordings", sc, length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); all[[sc]] <- res
  saveRDS(res, sprintf("experiments/m3plus/_out/compare_sim_%s.rds", sc))
}
res <- rbindlist(all, fill = TRUE); res[, group := sub(":.*", "", scenario)]
for (g in unique(res$group)) { cat(sprintf("\n===== %s =====\n", g))
  d <- res[group == g]; print(d[, .(n = .N, rmse = round(mean(rmse), 3), cor = round(mean(cor), 3), cover95 = round(mean(cover), 3)), by = .(state, model)][order(state, model)], row.names = FALSE)
  print(d[state == "vagal", .(abs_ln_err_tau_p = round(median(abs(tau_p_err)), 2), abs_ln_err_tau_s = round(median(abs(tau_s_err)), 2), LL_med = round(median(LL), 1), n_par = round(median(n_par), 1)), by = model][order(model)], row.names = FALSE) }
cat("\nEvent detection from the RR series (exercise / tilt): detected, median error of onset and end (s)\n"); print(res[group %in% c("exercise", "tilt") & model == "M3" & state == "vagal", .(detected_pct = round(100 * mean(detected), 1), onset_err_med = round(median(onset_err, na.rm = TRUE), 1), onset_err_iqr = paste(round(quantile(onset_err, c(.25, .75), na.rm = TRUE), 1), collapse = " to "), end_err_med = round(median(end_err, na.rm = TRUE), 1)), by = group], row.names = FALSE)
