# ==============================================================================
# experiments/softanchor/blockade_prior.R -- EXPERIMENTAL. A soft anchor informed by the blockade data: lambda_P free with a log-normal prior whose CENTRE comes from the structural
# estimate of experiments/blockade/ou_fit.R (about 3.5 for tau_s 100-300 s) instead of from the recording's own spectrum, with the between-subject spread as width.
#   Rscript experiments/softanchor/blockade_prior.R
# Variants (centre, width in ln): hard = production anchor; const_s1.4 = ln 3.5, 1.4 (blockade value, anchor ignored); const_s0.7 = ln 3.5, 0.7; anchor2x_s0.6 = twice the own anchor, 0.6;
# mix_s1.0 = geometric mean of the own anchor and 3.5, 1.0; free = no anchor (the earlier experiments). Simulations: the four scenarios of the pipeline (truth known); real: the 313
# recordings and Fantasia test-retest of two windows (ICC of the log values, young-vs-older AUC).
# ==============================================================================
source("analysis/00_setup.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); options(width = 200)
spec_free <- drop_term(model_spec("M3"), "anchor"); LAM_B <- 3.5
VARS <- list(hard = list(), const_s1.4 = list(c = "const", s = 1.4), const_s0.7 = list(c = "const", s = 0.7), anchor2x_s0.6 = list(c = "anchor2", s = 0.6), mix_s1.0 = list(c = "mix", s = 1.0), free = list(c = "anchor", s = Inf))
fit_v <- function(series, v, lp_true = NULL) {
  pri <- series$priors; ctl <- series$control; dy <- series$dy; la <- log(pri$energy_ratio)
  if (!length(v)) { spec <- model_spec("M3"); o <- optimize_spec(series, spec); par <- unpack_theta(o$par, spec, pri$energy_ratio); lam <- pri$energy_ratio
  } else { ctr <- switch(v$c, const = log(LAM_B), anchor2 = la + log(2), mix = 0.5 * la + 0.5 * log(LAM_B), anchor = la)
    cost <- function(th) { x <- variant_objective(th, dy, spec_free, pri, ctl); if (x >= 1e11) return(1e12); if (is.finite(v$s)) x + 0.5 * ((th[4] - ctr) / v$s)^2 else x }
    runs <- lapply(.m3_starts(pri, ctl), function(st) tryCatch(optim(c(unname(.par_to_theta(st)), ctr), cost, method = "BFGS", control = list(maxit = ctl$maxit)), error = function(e) NULL))
    vals <- vapply(runs, function(o) if (is.null(o) || !is.finite(o$value)) Inf else o$value, 0); k <- which(vals <= min(vals) + 0.01)[1]; par <- unpack_theta(runs[[k]]$par, spec_free, pri$energy_ratio); lam <- par$lp }
  f <- filter_eval(dy, par$kp, par$ks, par$lp, par$lR, ctl, full = TRUE); z <- f$innovations_diagnostic[which(is.finite(f$innovations_diagnostic) & f$gate_multiplier < ctl$diag_gate_max)]; n <- length(z)
  list(ok = TRUE, spec = "M3", conv = 0, tau_p = 1 / par$kp, tau_s = 1 / par$ks, par = par, X_p = f$X_p, X_s = f$X_s, X_p_var = f$X_p_var, X_s_var = f$X_s_var, share_p = var(f$X_p) / (var(f$X_p) + var(f$X_s)), rmse_ms = 1000 * sqrt(mean((dy - f$implied_rr)^2)),
       LL_star = f$LL_star, lambda = lam, sig2 = f$sig2, ks_ok = unname(ks.test(pnorm(z), "punif")$statistic) < 1.36 / sqrt(n)) }
# ---- simulations with known truth --------------------------------------------------------------------------------------------
one <- function(job) { truth <- with_seed(job$seed, RECOVERY_SCENARIOS[[job$scenario]]()); sim <- simulate_rr(cfg$simulation$duration_s, kp = 1 / truth$tau_p, ks = 1 / truth$tau_s, lp = truth$lambda_p, sig2 = 9e-4, lR = 0.02, h = cfg$simulation$h, seed = job$seed)
  series <- prepare_series(sim$rr); ts_ <- var(sim$x_p) / (var(sim$x_p) + var(sim$x_s))
  rbindlist(lapply(names(VARS), function(v) { r <- recovery_row(fit_v(series, VARS[[v]]), truth, sim, ts_); r[, `:=`(variant = v, scenario = job$scenario, anchor_err = log(series$priors$energy_ratio / truth$lambda_p))] })) }
jobs <- unlist(lapply(names(RECOVERY_SCENARIOS), function(sc) lapply(seq_len(cfg$simulation$n_per_scenario), function(i) list(scenario = sc, sim = i, seed = id_seed(paste(sc, i), cfg$seed)))), recursive = FALSE)
log_msg("Blockade-informed soft anchor: %d simulated recordings x %d variants", length(jobs), length(VARS)); sm <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); sm[, variant := factor(variant, names(VARS))]
dir.create("experiments/softanchor/_out", showWarnings = FALSE); saveRDS(sm, "experiments/softanchor/_out/blockade_prior_sim.rds")
cat("\n== Simulations (truth known), median over recordings ==\n"); print(sm[, .(abs_ln_err_lambda = round(median(abs(ln_err_lp)), 2), cor_xp = round(median(cor_x_p), 3), cor_xs = round(median(cor_x_s), 3), rmse_xp = round(median(rmse_x_p), 3), rmse_xs = round(median(rmse_x_s), 3),
                                abs_ln_tau_p = round(median(abs(ln_err_tau_p)), 2), abs_ln_tau_s = round(median(abs(ln_err_tau_s)), 2), cover_xp = round(mean(cover_x_p), 2), cover_xs = round(mean(cover_x_s), 2)), by = variant][order(variant)], row.names = FALSE)
cat("\nBy scenario, |ln error| of lambda_P and RMSE of x_s:\n"); print(dcast(sm[, .(e = sprintf("%.2f / %.2f", median(abs(ln_err_lp)), median(rmse_x_s))), by = .(scenario, variant)], variant ~ scenario, value.var = "e"), row.names = FALSE)
# ---- real data -----------------------------------------------------------------------------------------------------------------
rj <- list(); for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); for (id in names(d$rr)) rj[[length(rj) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
real <- rbindlist(drop_failed(map_recordings(rj, function(j) { s <- prepare_series(j$rr); rbindlist(lapply(names(VARS), function(v) { f <- fit_v(s, VARS[[v]]); data.table(ds = j$ds, id = j$id, variant = v, lambda = f$lambda, anchor = s$priors$energy_ratio, tau_p = f$tau_p, tau_s = f$tau_s, ks_ok = f$ks_ok, LL = f$LL_star, SD_p = sqrt(f$lambda * f$sig2), SD_s = sqrt(f$sig2)) })) }, n_cores), "recording"))
real[, variant := factor(variant, names(VARS))]; saveRDS(real, "experiments/softanchor/_out/blockade_prior_real.rds")
cat("\n== 313 real recordings ==\n"); print(real[, .(lambda_median = round(median(lambda), 2), lambda_p10 = round(quantile(lambda, .1), 2), lambda_p90 = round(quantile(lambda, .9), 1), KS_pass = round(mean(ks_ok), 3), tau_p_median = round(median(tau_p), 1), tau_p_in_prior = round(mean(tau_p >= 3 & tau_p <= 6), 2),
                                dLL_vs_hard_median = NA_real_), by = variant][order(variant)][, !"dLL_vs_hard_median"], row.names = FALSE)
h <- real[variant == "hard", .(ds, id, LL0 = LL)]; print(merge(real, h, by = c("ds", "id"))[variant != "hard", .(dLL_median = round(median(LL - LL0), 1), dLL_p10 = round(quantile(LL - LL0, .1), 1), share_better = round(mean(LL - LL0 > 2), 2)), by = variant][order(variant)], row.names = FALSE)
# ---- Fantasia test-retest ----------------------------------------------------------------------------------------------------
w <- cfg$windows$fantasia; dataset <- load_dataset("fantasia", cfg); ids <- names(dataset$rec_full); grp <- fantasia_group(ids, cfg); names(grp) <- ids; jw <- CJ(id = ids, win = c("A", "B"))
fw <- rbindlist(drop_failed(map_recordings(seq_len(nrow(jw)), function(i) { win <- if (jw$win[i] == "A") w$stability_A else w$stability_B; rr <- window_rec(dataset$rec_full[[jw$id[i]]], start_s = win[1], duration_s = win[2] - win[1])$rr; s <- prepare_series(rr)
  rbindlist(lapply(names(VARS), function(v) { f <- fit_v(s, VARS[[v]]); data.table(id = jw$id[i], win = jw$win[i], variant = v, tau_p = f$tau_p, tau_s = f$tau_s, lambda = f$lambda, SD_p = sqrt(f$lambda * f$sig2), SD_s = sqrt(f$sig2)) })) }, n_cores), "window"))
icc <- function(a, b) { x <- cbind(a, b); n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1); mse <- sum((x - rowMeans(x) - matrix(colMeans(x) - gm, n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
auc <- function(x, g) { r <- rank(x); n1 <- sum(g == "young"); n2 <- sum(g == "older"); (sum(r[g == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2) }
tr <- rbindlist(lapply(names(VARS), function(v) { a <- fw[variant == v & win == "A"][order(id)]; b <- fw[variant == v & win == "B"][order(id)]; k <- intersect(a$id, b$id); a <- a[id %in% k]; b <- b[id %in% k]
  rbindlist(lapply(c("tau_p", "tau_s", "SD_p", "SD_s", "lambda"), function(x) data.table(variant = v, variable = x, icc = icc(log(a[[x]]), log(b[[x]])), auc = auc(a[[x]], grp[a$id])))) }))
cat("\n== Fantasia test-retest ICC (log values) and young-vs-older AUC ==\n"); print(dcast(tr, variable ~ variant, value.var = "icc")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE); print(dcast(tr, variable ~ variant, value.var = "auc")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE)
