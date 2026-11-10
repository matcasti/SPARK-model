# ==============================================================================
# experiments/state_descriptors/sim.R -- can the STATE descriptors be recovered from RR? Simulated exercise-like and tilt-like recordings with known true states (600 s, the generator of
# experiments/hmm_states/sim.R): descriptors of the TRUE path against descriptors of the M3 reconstruction (plug-in and posterior draws), same specification.
#   Rscript experiments/state_descriptors/sim.R test      unit checks of the descriptor code
#   Rscript experiments/state_descriptors/sim.R           60 recordings per scenario (several minutes)  ->  _out/sim.rds
# ==============================================================================
source("experiments/state_descriptors/lib.R")
args <- commandArgs(TRUE)
if (length(args) && args[1] == "test") {
  tg <- 0:599; d <- ifelse(tg < 100, 0.5, ifelse(tg < 220, 8, ifelse(tg < 300, 4, 0.5))); r <- descr_d(tg, d, c(100, 220), H = 200, hold = 15, delta = 3.7)
  stopifnot(r["Dmax"] == 8, abs(r["Trec"] - 80) < 1e-9 || r["Trec"] == 80, r["Preturn"] == 1, r["follow_s"] == 379)
  r2 <- descr_d(tg, d, c(100, 220), H = 60, hold = 15, delta = 3.7); stopifnot(r2["Preturn"] == 0, is.na(r2["Trec"]) || r2["Trec"] == 80)
  r3 <- descr_d(tg[1:330], d[1:330], c(100, 220), H = 200, hold = 15, delta = 3.7); stopifnot(r3["Trec"] == 80)                    # (follow-up 110 s: return at 80 s is seen)
  r4 <- descr_d(tg[1:250], d[1:250], c(100, 220), H = 200, hold = 15, delta = 3.7); stopifnot(is.na(r4["Trec"]), is.na(r4["Preturn"]))   # not returned and follow-up shorter than H: unknown
  # AR(1) draws have unit marginal variance; stationary distance of two independent OU drives exceeds DELTA_MODEL 0.1% of the time
  set.seed(2); n <- 2e5; z <- function(rho) { w <- sqrt(1 - rho^2) * rnorm(n); w[1] <- rnorm(1); as.numeric(stats::filter(w, rho, method = "recursive")) }
  cat(sprintf("AR(1) variance %.3f %.3f; share of stationary distance above delta %.4f (0.001)\n", var(z(.8)), var(z(.99)), mean(sqrt(z(.8)^2 + z(.99)^2) > DELTA_MODEL)))
  cat("descriptor checks passed\n"); quit(save = "no") }

source("analysis/00_setup.R"); source("experiments/m3plus/lib.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); DUR <- 600; N_PER <- 60L
H_SC <- c(exercise = 300, tilt = 120, shift = 300); HOLD <- 15
osc_draw <- function() list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1))))
build <- function(sc, seed) { set.seed(seed); d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))), osc = osc_draw())
  mk <- function(mp, ms) sim_ext(DUR, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = d$osc, m_p_fun = mp, m_s_fun = ms)
  if (sc %in% c("exercise", "shift")) { p <- draw_exercise_pars(seed, DUR); if (sc == "shift") p$toff_p <- p$toff_s <- 1e6      # "shift": the displacement persists after the perturbation (no recovery): a true regime shift
    return(c(d, list(true = c(p$t_on, p$t_off), s = mk(prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s))))) }
  t_on <- runif(1, 0.28, 0.45) * DUR; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
  c(d, list(true = c(t_on, t_off), s = mk(prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))), prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100))))))) }

job_fun <- function(job) { z <- build(job$scenario, job$seed); s <- z$s; H <- H_SC[[job$scenario]]; rest <- c(20, z$true[1] - 5); pert <- z$true
  series <- prepare_series(s$rr); f <- fit_variant(series, model_spec("M3")); if (!isTRUE(f$ok)) return(NULL)
  # truth: whitened by the TRUE stationary SDs, rest centre from the true path itself
  gt <- to_grid1(s$t, s$x_p + s$z1, s$x_s, 0 * s$t, 0 * s$t); r <- gt$t >= rest[1] & gt$t < rest[2]; Vt <- c(z$lp * 9e-4, 9e-4)
  dt <- dist_of((gt$xp - mean(gt$xp[r])) / sqrt(Vt[1]), (gt$xs - mean(gt$xs[r])) / sqrt(Vt[2])); tr <- descr_d(gt$t, dt, pert, H, HOLD, DELTA_MODEL, r)
  ge <- to_grid1(s$t, f$X_p, f$X_s, f$X_p_var, f$X_s_var); Ve <- c(series$priors$energy_ratio * f$sig2, f$sig2)
  # variant C: scale from M3 refitted on the resting beats alone (STATE: calibrate from the resting regime, independently of the perturbation)
  rr_rest <- s$rr[s$t < z$true[1] - 5]; sr <- prepare_series(rr_rest); fr <- fit_variant(sr, model_spec("M3")); if (!isTRUE(fr$ok)) return(NULL); VeC <- c(sr$priors$energy_ratio * fr$sig2, fr$sig2)
  drC <- summ_draws(descr_draws(ge, rest, pert, VeC, c(fr$tau_p, fr$tau_s), n_draw = 100, H = H, hold = HOLD, seed = job$seed))
  pt <- ptransition(ge, rest, pert, VeC, H); if (is.null(pt)) return(NULL)
  pl <- descr_plugin(ge, rest, pert, H, HOLD); dr <- summ_draws(descr_draws(ge, rest, pert, Ve, c(f$tau_p, f$tau_s), n_draw = 100, H = H, hold = HOLD, seed = job$seed))
  if (is.null(pl) || is.null(dr)) return(NULL)
  data.table(scenario = job$scenario, sim = job$sim, H = H, t_true = tr["Trec"], D_true = tr["Dmax"], P_true = tr["Preturn"], dH_true = tr["d_H"], ex_true = tr["rest_exceed"],
             D_plug = pl["Dmax"], T_plug = pl["Trec"], P_plug = pl["Preturn"], delta_plug = pl["delta"], ex_plug = pl["rest_exceed"],
             D_med = dr["Dmax"], D_lo = dr["Dmax_lo"], D_hi = dr["Dmax_hi"], T_med = dr["Trec"], T_lo = dr["Trec_lo"], T_hi = dr["Trec_hi"], T_share = dr["Trec_ret_share"], P_draw = dr["Preturn"], dH_draw = dr["d_H"], ex_draw = dr["rest_exceed"],
             D_medC = drC["Dmax"], D_loC = drC["Dmax_lo"], D_hiC = drC["Dmax_hi"], T_medC = drC["Trec"], T_loC = drC["Trec_lo"], T_hiC = drC["Trec_hi"], P_drawC = drC["Preturn"], ex_drawC = drC["rest_exceed"],
             K_hmm = pt["K"], trans = pt["trans"], ptrans = pt["ptrans"], V_ratio_p = Ve[1] / Vt[1], V_ratio_s = Ve[2] / Vt[2], V_ratioC_p = VeC[1] / Vt[1], V_ratioC_s = VeC[2] / Vt[2]) }

jobs <- unlist(lapply(c("exercise", "tilt", "shift"), function(sc) lapply(seq_len(N_PER), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("stdesc", sc, i), cfg$seed)))), recursive = FALSE)
log_msg("State descriptors, simulations: %d recordings", length(jobs)); res <- rbindlist(drop_failed(map_recordings(jobs, job_fun, n_cores), "state descriptors")); saveRDS(res, "experiments/state_descriptors/_out/sim.rds")
