# ==============================================================================
# experiments/ou_logrr/compare/sim_compare.R -- production M3 against the candidate variants (F1-F5) on simulated recordings with known truth.
#   Rscript experiments/ou_logrr/compare/sim_compare.R [scenario ...]      scenarios: stationary rsa exercise tilt G1_E1 G1_E2 (default: all)
# Truth in the production generators (Hz) for stationary / rsa / exercise / tilt (the scenarios of experiments/m3plus/compare), and in the candidate's own
# generator G1 (ln RR) for G1_E1 / G1_E2. Filtered states = what production reports; smoothed = M3+RTS against the candidate's RTS.
# Output: experiments/ou_logrr/_out/compare_sim_<scenario>.rds
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/ou_logrr/compare/lib.R"); source("experiments/ou_logrr/sim.R")
ONLY <- Sys.getenv("ONLY", ""); SUFFIX <- if (nzchar(ONLY)) paste0("_", ONLY) else ""           # ONLY=F6: run just that variant (M3 and the others are in the earlier files of the same seeds)
args <- commandArgs(TRUE); scen <- if (length(args)) args else c("stationary", "rsa", "exercise", "tilt", "G1_E1", "G1_E2")
n_cores <- max(1L, min(5L, cfg$cores - 2L)); DUR <- 420
N <- c(stationary = 30L, rsa = 60L, exercise = 60L, tilt = 60L, G1_E1 = 40L, G1_E2 = 40L)

osc_draw <- function() list(alpha = 0.18, f_hz = runif(1, 0.2, 0.35), lo = exp(runif(1, log(0.2), log(1))))
build <- function(sc, seed) {                                              # as experiments/m3plus/compare/compare_sim.R
  set.seed(seed); d <- list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4))), osc = osc_draw())
  mk <- function(osc, mp = NULL, ms = NULL) sim_ext(DUR, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = osc, m_p_fun = mp, m_s_fun = ms)
  if (sc == "rsa") return(c(d, list(s = mk(d$osc), mode = NULL)))
  if (sc == "exercise") { p <- draw_exercise_pars(seed, DUR); return(c(d, list(mode = "exercise", true = c(p$t_on, p$t_off), s = mk(d$osc, prof_exercise(p$t_on, p$t_off, -p$amp_p, p$ton_p, p$toff_p), prof_exercise(p$t_on, p$t_off, p$amp_s, p$ton_s, p$toff_s))))) }
  t_on <- runif(1, 0.28, 0.45) * DUR; t_off <- t_on + runif(1, 100, 130); ap <- runif(1, 0.08, 0.14); as <- runif(1, 0.06, 0.12)
  c(d, list(mode = "tilt", true = c(t_on, t_off), s = mk(d$osc, prof_exercise(t_on, t_off, -ap, runif(1, 4, 9), exp(runif(1, log(10), log(40)))), prof_exercise(t_on, t_off, as, runif(1, 18, 40), exp(runif(1, log(25), log(100)))))))
}

#' Time (s after onset) at which an estimated state first reaches half of its change between the pre-event mean and the end-of-event mean.
t50 <- function(x, t, t0, t1) { pre <- t >= t0 - 120 & t < t0 - 10; en <- t >= t1 - 15 & t < t1; if (sum(pre) < 5 || sum(en) < 5) return(NA_real_)
  a <- mean(x[pre]); b <- mean(x[en]); k <- which(t >= t0 & (x - a) / (b - a) >= 0.5)[1]; if (is.na(k)) NA_real_ else t[k] - t0 }

one <- function(job) {
  sc <- job$scenario; units <- "hz"; ev <- NULL; extra <- list()
  if (sc == "stationary") {
    truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$sub]]()); sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
    rr <- sim$rr; tp <- sim$x_p; ts <- sim$x_s; tt <- cumsum(rr); tau <- c(truth$tau_p, truth$tau_s); lab <- paste0("stationary:", job$sub)
  } else if (grepl("^G1", sc)) {
    units <- "ln"; s <- scenario(sub("G1_", "", sc)); set.seed(job$seed); d <- simulate_beats(s$par, duration = s$duration, offset = s$offset)
    rr <- d$rr[seq_along(d$xp)]; tp <- d$xp; ts <- d$xs; tt <- d$t; tau <- c(s$par$tau_p, s$par$tau_s); lab <- sc
    if (sc == "G1_E2") ev <- c(300, 420)
  } else {
    z <- build(sc, job$seed); s <- z$s; rr <- s$rr; tp <- s$x_p + s$z1; ts <- s$x_s; tt <- s$t; tau <- c(z$tau_p, z$tau_s); lab <- sc
    if (!is.null(z$mode)) ev <- z$true
  }
  series <- prepare_series(rr); n <- length(series$dy); stopifnot(length(tp) == n)
  rows <- list()
  add <- function(model, st, smooth, fit_tau = c(NA_real_, NA_real_), n_par = NA_real_) {
    if (is.null(st)) return(NULL)
    sp <- score(st$xp, st$xpv, tp); ss <- score(st$xs, st$xsv, ts); ex <- list()
    if (!is.null(ev)) { ex <- list(chg_p = change_ratio(st$xp, tp, tt, ev[1], ev[2]), chg_s = change_ratio(st$xs, ts, tt, ev[1], ev[2]),
                                   t50_p = t50(st$xp, tt, ev[1], ev[2]), t50_s = t50(st$xs, tt, ev[1], ev[2])) }
    rows[[length(rows) + 1L]] <<- data.table(scenario = lab, sim = job$sim, model = model, smooth = smooth, rmse_p = sp[["rmse"]], cor_p = sp[["cor"]], cover_p = sp[["cover"]],
        rmse_s = ss[["rmse"]], cor_s = ss[["cor"]], cover_s = ss[["cover"]], ln_tau_p = log(fit_tau[1] / tau[1]), ln_tau_s = log(fit_tau[2] / tau[2]),
        chg_p = ex$chg_p %||% NA_real_, chg_s = ex$chg_s %||% NA_real_, t50_p = ex$t50_p %||% NA_real_, t50_s = ex$t50_s %||% NA_real_,
        t50_p_true = if (!is.null(ev)) t50(tp, tt, ev[1], ev[2]) else NA_real_, t50_s_true = if (!is.null(ev)) t50(ts, tt, ev[1], ev[2]) else NA_real_)
  }
  if (!nzchar(ONLY)) for (sm in c(FALSE, TRUE)) { m <- tryCatch(m3_states(series, smooth = sm, units = units), error = function(e) NULL)
    if (!is.null(m)) add("M3", m, sm, c(m$fit$tau_p, m$fit$tau_s)) }
  for (v in if (nzchar(ONLY)) ONLY else names(CAND_VARIANTS)) { f <- fit_cand_variant(series, CAND_VARIANTS[[v]]); if (is.null(f)) next
    for (sm in c(FALSE, TRUE)) add(v, cand_states(f, smooth = sm, units = units), sm, c(f$par$tau_p, f$par$tau_s)) }
  rbindlist(rows)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

dir.create("experiments/ou_logrr/_out", showWarnings = FALSE, recursive = TRUE)
for (sc in scen) {
  jobs <- if (sc == "stationary") unlist(lapply(names(RECOVERY_SCENARIOS), function(sub) lapply(seq_len(N[["stationary"]]), function(i) list(scenario = "stationary", sub = sub, sim = i, seed = id_seed(paste("ou", sub, i), cfg$seed)))), recursive = FALSE)
         else lapply(seq_len(N[[sc]]), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("ou", sc, i), cfg$seed)))
  log_msg("Scenario %s: %d simulated recordings", sc, length(jobs))
  res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"), fill = TRUE)
  saveRDS(res, sprintf("experiments/ou_logrr/_out/compare_sim_%s%s.rds", sc, SUFFIX))
}
