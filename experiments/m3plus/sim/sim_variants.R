# EXPERIMENTAL. Ideas 1 (oscillator), 4 (rate-dependent noise), 5 (robust weights) and 6 (known inputs) on simulations whose truth contains the feature,
# scored against the true states (smoothed estimates). Run from the project root:  Rscript experiments/m3plus/sim_variants.R [scenario]
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
only <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else NULL
n_cores <- max(1L, min(5L, cfg$cores - 2L)); N <- 80L
draw <- function(seed) { set.seed(seed); list(tau_p = runif(1, 3.5, 6), tau_s = runif(1, 120, 250), lp = exp(runif(1, log(1.5), log(4)))) }
SCEN <- list(
  rsa = list(variants = list(base = list(), osc = list(osc = TRUE)),
    sim = function(seed) { d <- draw(seed); set.seed(seed + 1); o <- list(alpha = 0.18, f_hz = runif(1, 0.18, 0.35), lo = exp(runif(1, log(0.2), log(1)))); c(d, list(o = o, s = sim_ext(600, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, osc = o))) }),
  het = list(variants = list(base = list(), het = list(het = TRUE)),
    sim = function(seed) { d <- draw(seed); st <- list(t0 = 150, t1 = 450, amp_p = -0.3, amp_s = 0.3, tau_p = 6, tau_s = 40); c(d, list(s = sim_ext(700, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, het = 1.5, steps = st))) }),
  artifacts = list(variants = list(base = list(), robust = list(robust = TRUE)),
    sim = function(seed) { d <- draw(seed); c(d, list(s = sim_ext(600, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, artifacts = 0.02))) }),
  inputs = list(variants = list(base = list(), inputs = list(inputs = TRUE)),
    sim = function(seed) { d <- draw(seed); st <- list(t0 = 150, t1 = 450, amp_p = -0.3, amp_s = 0.3, tau_p = 6, tau_s = 40); c(d, list(st = st, s = sim_ext(700, 1 / d$tau_p, 1 / d$tau_s, d$lp, 9e-4, 0.02, seed = seed, steps = st))) }))
score <- function(est, var, truth) c(rmse = sqrt(mean((est - truth)^2)) / sd(truth), cor = cor(est, truth), cover = mean(abs(est - truth) <= 1.96 * sqrt(pmax(var, 0))))
one <- function(job) {
  sc <- SCEN[[job$scenario]]; z <- sc$sim(job$seed); s <- z$s; series <- prepare_series(s$rr); off <- s$t[1] - s$rr[1]; rows <- list()
  for (v in names(sc$variants)) {
    spec <- sc$variants[[v]]; if (isTRUE(spec$inputs)) spec$inputs <- list(events = data.frame(t0 = z$st$t0 - off, t1 = z$st$t1 - off))
    f <- tryCatch(fit_ext(series, spec), error = function(e) list(ok = FALSE)); if (!isTRUE(f$ok)) next
    r <- f$run(full = TRUE, smooth = TRUE); xp <- r$X_p_smooth; xpv <- r$X_p_var_smooth; xs <- r$X_s_smooth; xsv <- r$X_s_var_smooth
    truth_p <- s$x_p; if (!is.null(z$o)) { xp <- xp + r$X_o_smooth; xpv <- xpv + r$X_o_var_smooth; truth_p <- s$x_p + s$z1 }                       # vagal total = OU + oscillator
    if (!is.null(spec$inputs)) { tt <- cumsum(series$dy); b <- r$beta; xp <- xp + b[2] * input_curve(tt, f$p$kp, spec$inputs$events$t0, spec$inputs$events$t1); xs <- xs + b[3] * input_curve(tt, f$p$ks, spec$inputs$events$t0, spec$inputs$events$t1) }
    ok <- if (isTRUE(spec$robust) || job$scenario == "artifacts") !s$artifact else rep(TRUE, length(xp))
    sp <- score(xp[ok], xpv[ok], truth_p[ok]); ss <- score(xs[ok], xsv[ok], s$x_s[ok])
    extra <- list(osc_hz_err = if (!is.null(f$p$osc)) f$omega / 2 / pi - z$o$f_hz else NA_real_, osc_alpha = if (!is.null(f$p$osc)) f$p$osc$alpha else NA_real_, gamma = f$p$gamma,
                  b_p = if (!is.null(spec$inputs)) r$beta[2] else NA_real_, b_s = if (!is.null(spec$inputs)) r$beta[3] else NA_real_)
    rows[[v]] <- data.table(scenario = job$scenario, sim = job$sim, variant = v, state = c("vagal", "sympathetic"), rbind(sp, ss), tau_p_err = log(f$tau_p / z$tau_p), tau_s_err = log(f$tau_s / z$tau_s),
                            LL = r$LL_star, n_par = f$n_par, ks_ok = { k <- which(is.finite(r$innovations_diagnostic) & r$gate_multiplier < f$ctl$diag_gate_max); unname(ks.test(pnorm(r$innovations_diagnostic[k]), "punif")$statistic) < 1.36 / sqrt(length(k)) }, as.data.table(extra))
  }
  rbindlist(rows)
}
scen <- if (is.null(only)) names(SCEN) else only
for (sc in scen) {
  jobs <- lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste("m3plus", sc, i), cfg$seed)))
  log_msg("Scenario %s: %d simulated recordings", sc, N); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation"))
  dir.create("experiments/m3plus/_out", showWarnings = FALSE); saveRDS(res, sprintf("experiments/m3plus/_out/sim_%s.rds", sc))
  cat(sprintf("\n===== %s =====\n", sc))
  print(res[, .(n = .N, rmse = round(mean(rmse), 3), cor = round(mean(cor), 3), cover95 = round(mean(cover), 3)), by = .(state, variant)][order(state, variant)], row.names = FALSE)
  print(res[, .(abs_ln_err_tau_p = round(median(abs(tau_p_err)), 2), abs_ln_err_tau_s = round(median(abs(tau_s_err)), 2), ks_pass = round(100 * mean(ks_ok), 1), extra_par = median(n_par - 3), LL_med = round(median(LL), 1)), by = variant], row.names = FALSE)
  if (sc == "rsa") print(res[variant == "osc" & state == "vagal", .(freq_err_hz_med = round(median(osc_hz_err), 3), freq_err_abs_med = round(median(abs(osc_hz_err)), 3), alpha_med = round(median(osc_alpha), 2))])
  if (sc == "het") print(res[variant == "het" & state == "vagal", .(gamma_med = round(median(gamma), 2), gamma_q25 = round(quantile(gamma, .25), 2), gamma_q75 = round(quantile(gamma, .75), 2))])
}
