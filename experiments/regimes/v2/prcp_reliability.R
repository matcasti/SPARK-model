# EXPERIMENTAL (round 2d). Repeated manoeuvres as a test-retest: each PRCP recording has 6-7 annotated tilts / stand-ups. Index per block = mean state in [onset + 30 s, return] minus mean state in the
# supine window [onset - 90 s, onset - 10 s], for the vagal-like and sympathetic-like states (and heart rate as a reference). Per model: share of blocks with the expected sign (vagal down, sympathetic up),
# signal-to-noise of the block effect (mean / SD across the 62 blocks), and the one-way ICC across recordings (how much of the block-to-block variance is between people).
# Models: M3, M3E_det (detector), HMMr (smoothed), HYB. States: filtered and smoothed. Run from the project root:  Rscript experiments/regimes/v2/prcp_reliability.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); ev_ann <- fread("data/prcp_events.csv"); d <- load_dataset("prcp", cfg); n_cores <- max(1L, min(5L, cfg$cores - 2L))
icc1 <- function(y, g) { g <- factor(g); k <- table(g); k0 <- (sum(k) - sum(k^2) / sum(k)) / (nlevels(g) - 1); m <- tapply(y, g, mean); msb <- sum(k * (m - mean(y))^2) / (nlevels(g) - 1); msw <- sum((y - m[as.character(g)])^2) / (length(y) - nlevels(g)); (msb - msw) / (msb + (k0 - 1) * msw) }
one <- function(rid) {
  series <- prepare_series(d$rr[[rid]]); dy <- series$dy; tt <- cumsum(dy); off <- d$rec[[rid]]$time[1] - d$rec[[rid]]$rr[1]; a1 <- ev_ann[as.character(ev_ann$id) == rid & is.finite(onset_s) & is.finite(return_s)]; if (!nrow(a1)) return(NULL)
  blocks <- data.table(on = a1$onset_s - off, ret = a1$return_s - off, kind = a1$kind); base <- fit_ext(series); cp <- cp_search(series, base, min_gain = -Inf); st <- list(); hr <- 60 / dy
  for (sm in c(FALSE, TRUE)) { e <- state_estimates(base, sm); st[[paste0("M3|", sm)]] <- list(xp = e$xp, xs = e$xs) }
  ev_det <- events_from_rr(dy, "tilt", cfg$segmentation); if (!is.null(ev_det)) { f <- fit_ext(series, list(inputs = list(events = ev_det))); if (isTRUE(f$ok)) for (sm in c(FALSE, TRUE)) { e <- state_estimates(f, sm); st[[paste0("M3E_det|", sm)]] <- list(xp = e$xp, xs = e$xs) } }
  if (!is.null(cp) && cp$gain >= 9) { h <- fit_regime(series, 2, cp, base, 7, 0.7, recip = 20); if (isTRUE(h$ok)) { r <- h$run(smooth = TRUE); st[["HMMr|FALSE"]] <- list(xp = r$X_p, xs = r$X_s); st[["HMMr|TRUE"]] <- list(xp = r$Xs_p, xs = r$Xs_s)
      Ws <- regime_smooth(r$W, r$A, h$p$Pi); ev_h <- regime_episodes(Ws, tt, max_ev = 8); if (!is.null(ev_h)) { f <- tryCatch(fit_ext(series, list(inputs = list(events = ev_h))), error = function(e) NULL); if (!is.null(f) && isTRUE(f$ok)) for (sm in c(FALSE, TRUE)) { e <- state_estimates(f, sm); st[[paste0("HYB|", sm)]] <- list(xp = e$xp, xs = e$xs) } } } }
  if (is.null(st[["HMMr|TRUE"]])) { st[["HMMr|FALSE"]] <- st[["M3|FALSE"]]; st[["HMMr|TRUE"]] <- st[["M3|TRUE"]] }; if (is.null(st[["HYB|TRUE"]])) { st[["HYB|FALSE"]] <- st[["M3|FALSE"]]; st[["HYB|TRUE"]] <- st[["M3|TRUE"]] }
  rbindlist(lapply(seq_len(nrow(blocks)), function(i) { b <- blocks[i]; w1 <- tt >= b$on + 30 & tt <= b$ret; w0 <- tt >= b$on - 90 & tt <= b$on - 10; if (sum(w1) < 20 || sum(w0) < 20) return(NULL)
    rbindlist(c(list(data.table(id = rid, block = i, kind = b$kind, model = "HR", smooth = NA, dvagal = NA_real_, dsymp = NA_real_, dhr = mean(hr[w1]) - mean(hr[w0]))),
      lapply(names(st), function(nm) { p <- strsplit(nm, "|", fixed = TRUE)[[1]]; data.table(id = rid, block = i, kind = b$kind, model = p[1], smooth = p[2] == "TRUE", dvagal = mean(st[[nm]]$xp[w1]) - mean(st[[nm]]$xp[w0]), dsymp = mean(st[[nm]]$xs[w1]) - mean(st[[nm]]$xs[w0]), dhr = NA_real_) }))) })) }
log_msg("PRCP block reliability"); res <- rbindlist(drop_failed(map_recordings(as.list(names(d$rr)), one, n_cores), "recording")); saveRDS(res, "experiments/regimes/v2/_out/prcp_reliability.rds")
f <- function(x) list(n = sum(is.finite(x)), mean = signif(mean(x, na.rm = TRUE), 3), snr = round(mean(x, na.rm = TRUE) / sd(x, na.rm = TRUE), 2))
print(res[model == "HR", c(list(), f(dhr), icc = round(icc1(dhr, id), 2))])
out <- res[model != "HR", { c(list(n = .N), vagal_snr = f(-dvagal)$snr, vagal_exp_sign = round(100 * mean(dvagal < 0)), vagal_icc = round(icc1(dvagal, id), 2), symp_snr = f(dsymp)$snr, symp_exp_sign = round(100 * mean(dsymp > 0)), symp_icc = round(icc1(dsymp, id), 2)) }, by = .(model, smooth)]
print(out[order(smooth, match(model, c("M3", "M3E_det", "HMMr", "HYB")))])
