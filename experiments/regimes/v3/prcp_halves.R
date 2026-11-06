# EXPERIMENTAL (round 4). Split-half reliability of person-level indices in PRCP (10 people, 6-7 manoeuvres each, three per half). Each half is fitted independently (M3 and the gated regime model).
# Indices: (a) mean block effect on the vagal and sympathetic states (as in round 2d) from M3 and from the regime model's states; (b) regime parameters (targets mu_p, mu_s and the rate change -mu_p + mu_s) from the
# regime model, when its gate opens in both halves. Reliability = correlation between halves across people (Pearson and Spearman, n = number of people with both halves). Run: Rscript experiments/regimes/v3/prcp_halves.R
source("analysis/00_setup.R"); source("experiments/regimes/v2/lib2.R"); ev_ann <- fread("data/prcp_events.csv"); d <- load_dataset("prcp", cfg)
half_idx <- function(rid, tt, off) { a <- ev_ann[as.character(ev_ann$id) == rid & is.finite(onset_s) & is.finite(return_s)][order(onset_s)]; m <- floor(nrow(a) / 2); Ts <- (a$return_s[m] + a$onset_s[m + 1]) / 2 - off; list(Ts = Ts, a = a, off = off, k = which(tt >= Ts)[1]) }
block_idx <- function(st, tt, blocks) { r <- lapply(seq_len(nrow(blocks)), function(i) { w1 <- tt >= blocks$on[i] + 30 & tt <= blocks$ret[i]; w0 <- tt >= blocks$on[i] - 90 & tt <= blocks$on[i] - 10; if (sum(w1) < 20 || sum(w0) < 20) return(c(NA, NA)); c(mean(st$xp[w1]) - mean(st$xp[w0]), mean(st$xs[w1]) - mean(st$xs[w0])) }); m <- do.call(rbind, r); colMeans(m, na.rm = TRUE) }
rows <- list()
for (rid in names(d$rr)) { series <- prepare_series(d$rr[[rid]]); dy <- series$dy; tt <- cumsum(dy); off <- d$rec[[rid]]$time[1] - d$rec[[rid]]$rr[1]; hs <- half_idx(rid, tt, off); a <- hs$a
  for (half in 1:2) { seg <- if (half == 1) dy[seq_len(hs$k)] else dy[(hs$k + 1):length(dy)]; shift <- if (half == 1) 0 else tt[hs$k]; sel <- if (half == 1) seq_len(floor(nrow(a) / 2)) else (floor(nrow(a) / 2) + 1):nrow(a)
    blocks <- data.table(on = a$onset_s[sel] - off - shift, ret = a$return_s[sel] - off - shift); s <- prepare_series(seg); ds <- s$dy; ts <- cumsum(ds); base <- fit_ext(s); e <- state_estimates(base, FALSE); cp <- cp_search(s, base, min_gain = -Inf)
    m3 <- block_idx(list(xp = e$xp, xs = e$xs), ts, blocks); row <- data.table(id = rid, half = half, m3_vagal = m3[1], m3_symp = m3[2], gain = if (is.null(cp)) NA_real_ else cp$gain)
    if (!is.null(cp) && cp$gain >= 9) { h <- fit_regime(s, 2, cp, base, 7, 0.7, recip = 20); if (isTRUE(h$ok)) { r <- h$run(smooth = TRUE); hm <- block_idx(list(xp = r$Xs_p, xs = r$Xs_s), ts, blocks)
        row[, `:=`(hmm_vagal = hm[1], hmm_symp = hm[2], mu_p = h$p$mu[2, 1], mu_s = h$p$mu[2, 2], rate_change = -h$p$mu[2, 1] + h$p$mu[2, 2], frac_regime = mean(r$W[, 2]))] } }
    rows[[length(rows) + 1]] <- row } }
res <- rbindlist(rows, fill = TRUE); saveRDS(res, "experiments/regimes/v3/_out/prcp_halves.rds"); w <- dcast(res, id ~ half, value.var = setdiff(names(res), c("id", "half")))
rel <- function(v) { a <- w[[paste0(v, "_1")]]; b <- w[[paste0(v, "_2")]]; ok <- is.finite(a) & is.finite(b); if (sum(ok) < 4) return(data.table(index = v, n = sum(ok), pearson = NA_real_, spearman = NA_real_)); data.table(index = v, n = sum(ok), pearson = round(cor(a[ok], b[ok]), 2), spearman = round(cor(a[ok], b[ok], method = "spearman"), 2)) }
print(rbindlist(lapply(c("m3_vagal", "m3_symp", "hmm_vagal", "hmm_symp", "mu_p", "mu_s", "rate_change", "frac_regime"), rel)))
cat("\nGate open in both halves:", sum(is.finite(w$mu_p_1) & is.finite(w$mu_p_2)), "of", nrow(w), "people\n")
