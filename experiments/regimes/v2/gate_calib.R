# EXPERIMENTAL (round 2b). Calibration of the gate (change-point log-likelihood gain >= 9). Null distributions of the gain: (a) parametric bootstrap from each Fantasia recording's M3 fit,
# (b) parametric bootstrap from its M3O fit (respiratory oscillator present, M3 and the search ignore it), (c) the observed gain in the 40 Fantasia recordings (real rest, no regime expected).
# Run from the project root:  Rscript experiments/regimes/v2/gate_calib.R [recordings] [B]
args <- commandArgs(TRUE); NR <- if (length(args) >= 1) as.integer(args[1]) else 20L; B <- if (length(args) >= 2) as.integer(args[2]) else 6L
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L)); d <- load_dataset("fantasia", cfg); ids <- names(d$rr); set.seed(31); ids <- sort(sample(ids, NR))
jobs <- lapply(ids, function(id) list(id = id, rr = d$rr[[id]]))
gain_of <- function(rr) { s <- prepare_series(rr); b <- fit_ext(s); cp <- cp_search(s, b, min_gain = -Inf); if (is.null(cp)) NA_real_ else cp$gain }
one <- function(job) {
  series <- prepare_series(job$rr); dy <- series$dy; dur <- sum(dy); base <- fit_ext(series); rb <- base$run(); fo <- fit_ext(series, list(osc = TRUE)); lp <- series$priors$energy_ratio; out <- list()
  out[[1]] <- data.table(id = job$id, kind = "observed", b = 0L, gain = gain_of(job$rr), n = length(dy))
  for (i in seq_len(B)) { seed <- id_seed(paste("gate", job$id, i), cfg$seed)
    s <- sim_ext(dur, base$p$kp, base$p$ks, lp, rb$sig2, base$p$lR, nu0 = rb$nu0, seed = seed); out[[length(out) + 1]] <- data.table(id = job$id, kind = "boot_M3", b = i, gain = tryCatch(gain_of(s$rr), error = function(e) NA_real_), n = nrow(s))
    if (isTRUE(fo$ok)) { ro <- fo$run(); s <- sim_ext(dur, fo$p$kp, fo$p$ks, lp, ro$sig2, fo$p$lR, nu0 = ro$nu0, seed = seed + 1L, osc = list(alpha = fo$p$osc$alpha, f_hz = fo$p$osc$omega / (2 * pi), lo = fo$p$osc$lo))
      out[[length(out) + 1]] <- data.table(id = job$id, kind = "boot_M3O", b = i, gain = tryCatch(gain_of(s$rr), error = function(e) NA_real_), n = nrow(s)) } }
  rbindlist(out) }
log_msg("Gate calibration | %d Fantasia recordings x %d bootstraps x 2 nulls", length(jobs), B); res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "recording")); saveRDS(res, "experiments/regimes/v2/_out/gate_calib.rds")
q <- function(x) as.list(round(quantile(x, c(.5, .9, .95, .99, 1), na.rm = TRUE), 1)); print(res[, c(list(n_rec = uniqueN(id), n = .N), q(gain)), by = kind]); cat("\nExceedance of the current gate (9) and of the bootstrap 95th percentile:\n")
r1 <- readRDS("experiments/regimes/_out/real.rds")$res[model == "M3" & ds == "fantasia"]; print(res[, .(share_ge_9 = round(100 * mean(gain >= 9, na.rm = TRUE))), by = kind]); cat("Round-1 Fantasia (40): share >= 9 =", round(100 * mean(r1$gain >= 9, na.rm = TRUE)), "%; quantiles", unlist(q(r1$gain)), "\n")
