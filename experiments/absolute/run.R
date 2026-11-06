# EXPERIMENTAL. Absolute autonomic levels from fitted states and parameters: what is identified, what needs an external anchor, and whether it is stable and coherent.
# Part 1 (Fantasia, 40 resting recordings): M3 fitted on the whole record and on 4 consecutive windows; parameters recorded: mean rate, kappa_p, kappa_s, sigma^2, lambda_P (spectral anchor).
# Part 2 (exercise, 100 recordings): M3E (events) fit; rest-window variability, rest rate and vagal withdrawal, for the coherence bound (withdrawal <= implied resting vagal tone).
# Run from the project root:  Rscript experiments/absolute/run.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R"); n_cores <- max(1L, min(5L, cfg$cores - 2L))
pars <- function(f, series, dy) { r <- f$run(); data.table(rate = 1 / mean(dy), nu0 = r$nu0, kp = f$p$kp, ks = f$p$ks, lp = series$priors$energy_ratio, sig2 = r$sig2, lR = f$p$lR, n = length(dy)) }
d <- load_dataset("fantasia", cfg); jobsF <- lapply(names(d$rr), function(id) list(id = id, rr = d$rr[[id]]))
oneF <- function(job) { series <- prepare_series(job$rr); dy <- series$dy; tt <- cumsum(dy); w <- cut(tt, 4, labels = FALSE); out <- list()
  f <- fit_ext(series); out[[1]] <- cbind(data.table(id = job$id, window = 0L), pars(f, series, dy))
  for (k in 1:4) { s <- prepare_series(dy[w == k]); fk <- fit_ext(s); out[[k + 1]] <- cbind(data.table(id = job$id, window = k), pars(fk, s, s$dy)) }; rbindlist(out) }
log_msg("Absolute levels, Fantasia windows | %d recordings", length(jobsF)); resF <- rbindlist(drop_failed(map_recordings(jobsF, oneF, n_cores), "recording")); resF[, group := ifelse(grepl("^f.y", id), "young", "older")]
set.seed(5); de <- load_dataset("exercise", cfg); ids <- sort(sample(names(de$rr), 100)); jobsE <- lapply(ids, function(id) list(id = id, rr = de$rr[[id]]))
oneE <- function(job) { series <- prepare_series(job$rr); dy <- series$dy; tt <- cumsum(dy); ev <- events_from_rr(dy, "exercise", cfg$segmentation); if (is.null(ev)) return(NULL); f <- fit_ext(series, list(inputs = list(events = ev))); if (!isTRUE(f$ok)) return(NULL)
  e <- state_estimates(f, FALSE); t_on <- ev$t0[1]; t_end <- ev$t1[1]; rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end; if (sum(rest) < 20 || sum(endw) < 5) return(NULL)
  cbind(data.table(id = job$id, rest_rate = 1 / mean(dy[rest]), end_rate = 1 / mean(dy[endw]), d_vagal = mean(e$xp[endw]) - mean(e$xp[rest]), d_symp = mean(e$xs[endw]) - mean(e$xs[rest])), pars(f, series, dy)) }
log_msg("Absolute levels, exercise | %d recordings", length(jobsE)); resE <- rbindlist(drop_failed(map_recordings(jobsE, oneE, n_cores), "recording")); saveRDS(list(fantasia = resF, exercise = resE), "experiments/absolute/_out/run.rds"); cat("Saved", nrow(resF), "Fantasia rows and", nrow(resE), "exercise rows\n")
