# ==============================================================================
# experiments/ou_logrr/compare/real_compare.R -- production M3 against the candidate variants (F1-F5) on the 313 real recordings.
#   Rscript experiments/ou_logrr/compare/real_compare.R [fantasia] [prcp] [exercise]       (default: all three)
# Same windows and RR series as the pipeline (the RR of results/fit_<dataset>.rds). Per recording and model: parameters in physical units, innovation KS and PACF,
# fit error, agreement with M3; exercise: ordering and size of the responses (events from the RR as the pipeline does); PRCP: the pipeline's tilt analysis run on
# each model's states; Fantasia: known groups (young vs older) and test-retest of two windows.   Output: experiments/ou_logrr/_out/real_compare_<dataset>.rds
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/ou_logrr/compare/lib.R")
args <- commandArgs(TRUE); datasets <- if (length(args)) args else c("fantasia", "prcp", "exercise"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
MODELS <- c("M3", names(CAND_VARIANTS))
innov_diag <- function(z) { z <- z[is.finite(z)]; n <- length(z)
  c(ks_p = suppressWarnings(ks.test(pnorm(z), "punif")$p.value), pacf_viol = sum(abs(pacf(z, lag.max = 20, plot = FALSE)$acf) > 1.96 / sqrt(n))) }

#' All models on one RR series. Returns a named list: per model parameters, filtered states in Hz (+ variances), innovations diagnostics, implied-RR error.
fit_models <- function(rr) {
  series <- prepare_series(rr); out <- list()
  m3 <- tryCatch(fit_model_fast(rr, series$priors, series$control), error = function(e) NULL)
  if (!is.null(m3)) { st <- m3$states; p <- m3$params; d <- innov_diag(m3$innovations)
    out$M3 <- list(tau_p = 1 / p$kp, tau_s = 1 / p$ks, sd_p = sqrt(p$lp * p$sig2), sd_s = sqrt(p$sig2), nu0 = p$nu0, lp = p$lp, xp = st[, "X_p"], xs = st[, "X_s"], xpv = st[, "X_p_var"], xsv = st[, "X_s_var"],
                   ks_p = d[["ks_p"]], pacf_viol = d[["pacf_viol"]], rmse_ms = 1000 * sqrt(mean(m3$time_residuals^2, na.rm = TRUE)), ok = TRUE) }
  for (v in names(CAND_VARIANTS)) { f <- fit_cand_variant(series, CAND_VARIANTS[[v]]); if (is.null(f)) next
    s <- cand_states(f); nu <- exp(-f$c); z <- innovations(f); d <- innov_diag(z)
    out[[v]] <- list(tau_p = f$par$tau_p, tau_s = f$par$tau_s, sd_p = nu * sqrt(f$par$Vp), sd_s = nu * sqrt(f$par$Vs), nu0 = nu, lp = f$par$Vp / f$par$Vs, xp = s$xp, xs = s$xs, xpv = s$xpv, xsv = s$xsv,
                     ks_p = d[["ks_p"]], pacf_viol = d[["pacf_viol"]], rmse_ms = 1000 * sqrt(mean((exp(f$c + f$x_filt[, 1] - f$x_filt[, 2]) - rr)^2)), ok = f$convergence == 0) }
  out
}

#' Exercise: onset / end from the RR (as the pipeline), time to half of the change and standardised change for each state.
exercise_metrics <- function(models, dy) {
  ev <- events_from_rr(dy, "exercise", cfg$segmentation); if (is.null(ev)) return(NULL); tt <- cumsum(dy); t0 <- ev$t0[1]; t1 <- ev$t1[1]
  rows <- lapply(names(models), function(m) { x <- models[[m]]; rest <- tt >= max(0, t0 - 120) & tt < t0 - 10; en <- tt >= t1 - 15 & tt < t1; if (sum(rest) < 20 || sum(en) < 5) return(NULL)
    f <- function(v) { a <- mean(v[rest]); b <- mean(v[en]); k <- which(tt >= t0 & (v - a) / (b - a) >= 0.5)[1]
                       c(chg = b - a, z = (b - a) / sd(v[rest]), t50 = if (is.na(k)) NA_real_ else tt[k] - t0) }
    p <- f(x$xp); s <- f(x$xs)
    data.table(model = m, chg_p = p[["chg"]], chg_s = s[["chg"]], z_p = p[["z"]], z_s = s[["z"]], t50_p = p[["t50"]], t50_s = s[["t50"]], vagal_first = p[["t50"]] < s[["t50"]], exp_p = p[["chg"]] < 0, exp_s = s[["chg"]] > 0) })
  rbindlist(rows)
}

per_recording <- function(job) {
  m <- fit_models(job$rr); if (!length(m)) return(NULL)
  rows <- rbindlist(lapply(names(m), function(nm) { x <- m[[nm]]; r3 <- m$M3
    data.table(model = nm, id = job$id, n = length(job$rr), tau_p = x$tau_p, tau_s = x$tau_s, sd_p = x$sd_p, sd_s = x$sd_s, nu0 = x$nu0, lp = x$lp, ks_p = x$ks_p, pacf_viol = x$pacf_viol, rmse_ms = x$rmse_ms, ok = x$ok,
               cor_xp_M3 = if (!is.null(r3)) cor(x$xp, r3$xp) else NA_real_, cor_xs_M3 = if (!is.null(r3)) cor(x$xs, r3$xs) else NA_real_) }))
  list(rows = rows, exercise = if (job$dataset == "exercise") exercise_metrics(m, job$rr) else NULL, states = if (job$dataset == "prcp") lapply(m, function(x) x[c("xp", "xs", "xpv", "xsv")]) else NULL)
}

icc31 <- function(a, b) { x <- cbind(a, b); k <- is.finite(rowSums(x)); x <- x[k, , drop = FALSE]; n <- nrow(x); if (n < 5) return(NA_real_)
  gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1); mse <- sum((x - rowMeans(x) + mean(x) - matrix(colMeans(x), n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
auc <- function(x, g) { k <- is.finite(x); x <- x[k]; g <- g[k]; r <- rank(x); n1 <- sum(g == "young"); n2 <- sum(g == "older"); (sum(r[g == "young"]) - n1 * (n1 + 1) / 2) / (n1 * n2) }

dir.create("experiments/ou_logrr/_out", showWarnings = FALSE, recursive = TRUE)
for (ds in datasets) {
  fits <- read_result(paste0("fit_", ds), cfg); ids <- fits$index$SubjectID
  jobs <- lapply(ids, function(id) list(id = id, dataset = ds, rr = fits$states[SubjectID == id, RR_Interval]))
  log_msg("%s: %d recordings x %d models", ds, length(jobs), length(MODELS))
  res <- drop_failed(map_recordings(jobs, per_recording, n_cores), "recording")
  out <- list(rows = rbindlist(lapply(res, `[[`, "rows")), exercise = rbindlist(lapply(res, `[[`, "exercise")))
  if (ds == "prcp") {                                                                    # the pipeline's tilt analysis on each model's states
    ev <- read_prcp_events(cfg); out$prcp <- rbindlist(lapply(MODELS, function(m) { st <- copy(fits$states)
      for (i in seq_along(ids)) { k <- which(st$SubjectID == ids[i]); s <- res[[i]]$states[[m]]; if (is.null(s)) next
        st[k, `:=`(X_p = s$xp, X_s = s$xs, X_p_var = s$xpv, X_s_var = s$xsv)] }
      r <- drop_failed(map_recordings(ids, prcp_analyse_recording, n_cores, states = st, cfg = cfg, maneuvers = ev), "prcp")
      resp <- bind_results(r, "response"); if (!is.null(resp)) resp <- resp[kind != "stand-up"]; if (is.null(resp) || !nrow(resp)) return(NULL); resp[, model := m]; resp }), fill = TRUE) }
  if (ds == "fantasia") {                                                                # known groups and test-retest of two windows
    dataset <- load_dataset("fantasia", cfg); w <- cfg$windows$fantasia
    win <- function(a) lapply(ids, function(id) { r <- window_rec(dataset$rec_full[[id]], start_s = a[1], duration_s = a[2] - a[1]); list(id = id, dataset = "x", rr = r$rr) })
    wa <- drop_failed(map_recordings(win(w$stability_A), per_recording, n_cores), "window"); wb <- drop_failed(map_recordings(win(w$stability_B), per_recording, n_cores), "window")
    out$windows <- rbind(rbindlist(lapply(wa, `[[`, "rows"))[, window := "A"], rbindlist(lapply(wb, `[[`, "rows"))[, window := "B"])
    out$group <- fits$index[, .(id = SubjectID, Group)] }
  saveRDS(out, sprintf("experiments/ou_logrr/_out/real_compare_%s.rds", ds))
}
