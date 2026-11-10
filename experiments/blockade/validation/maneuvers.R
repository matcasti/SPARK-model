# ==============================================================================
# experiments/blockade/validation/maneuvers.R -- do the SPARK states move with the sympathetic nerve signal (MSNA) and with the maneuvers that drive it?
#   Rscript experiments/blockade/validation/maneuvers.R
# Per session (7 men x 2 days), saline period only: M3 fitted to the RR of the whole period (baseline, paced breathing, hyperventilation, breath-hold; the same fit the paper uses);
# filtered x_p and x_s (Hz) against the 1-s means of the integrated nerve signal. Expectations fixed before the run: (H1) MSNA is higher in the second half of the breath-hold than at
# baseline (>= 10 of 14 sessions); (H2) x_s is higher too (>= 10 of 14); (H3) across sessions the change of MSNA correlates more with the change of x_s than with that of x_p (Spearman).
# Changes are z-scores against the spontaneous baseline segment of the same session. Also the lagged correlation of x_s and x_p with MSNA over the whole saline period.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/blockade/lib.R")
mv <- fread("experiments/blockade/_out/maneuvers.csv"); pr <- fread("experiments/blockade/_out/protocol.csv"); rows <- list(); lag_rows <- list(); LAGS <- seq(-40, 40, by = 5)
for (i in seq_len(nrow(mv))) { m <- mv[i]; s <- m$subject; d <- m$day
  base <- pr[subject == s & day == d & cond == "saline" & breathing == "spontaneous"][1]; if (!nrow(base) || any(is.na(c(m$hv0, m$bh0, m$bh1, m$drug0)))) next
  tk <- clean_peaks(scan(sprintf("experiments/blockade/_out/rpeaks/rpeaks_Subject%d_day%d.csv", s, d), quiet = TRUE)); k <- tk >= base$t0 & tk < m$drug0; tb <- tk[k][-1]; rr <- diff(tk[k]); ok <- rr > 0.25 & rr < 2.5; tb <- tb[ok]; rr <- rr[ok]
  f <- fit_variant(prepare_series(rr), model_spec("M3")); if (!isTRUE(f$ok)) next
  ns <- fread(sprintf("experiments/blockade/_out/msna/msna_Subject%d_day%d.csv", s, d)); grid <- seq(ceiling(min(tb)), floor(max(tb)))
  ms <- approx(ns$t, ns$mean, xout = grid)$y; ms10 <- as.numeric(stats::filter(ms, rep(1 / 10, 10), sides = 2))
  xs <- approx(tb, f$X_s, xout = grid)$y; xp <- approx(tb, f$X_p, xout = grid)$y; hr <- approx(tb, 60 / rr, xout = grid)$y
  zb <- function(x) { b <- grid >= base$t0 + 20 & grid <= base$t1; (x - mean(x[b], na.rm = TRUE)) / sd(x[b], na.rm = TRUE) }
  Z <- list(MSNA = zb(ms10), x_s = zb(xs), x_p = zb(xp), HR = zb(hr)); win <- list(hyperventilation = c(m$hv0, m$bh0), breath_hold_late = c(m$bh0 + (m$bh1 - m$bh0) / 2, m$bh1))
  for (w in names(win)) { a <- grid >= win[[w]][1] & grid <= win[[w]][2]; rows[[length(rows) + 1L]] <- data.table(subject = s, day = d, window = w, signal = names(Z), change = sapply(Z, function(z) mean(z[a], na.rm = TRUE))) }
  for (sg in c("x_s", "x_p", "HR")) { cc <- sapply(LAGS, function(L) { j <- seq_along(grid); jj <- j + L; v <- jj >= 1 & jj <= length(grid); cor(Z[[sg]][j[v]], Z$MSNA[jj[v]], use = "complete.obs") })
    lag_rows[[length(lag_rows) + 1L]] <- data.table(subject = s, day = d, signal = sg, lag = LAGS, r = cc) } }
ch <- rbindlist(rows); lg <- rbindlist(lag_rows); saveRDS(list(change = ch, lagged = lg), "experiments/blockade/_out/maneuvers_results.rds")
w <- dcast(ch, subject + day + window ~ signal, value.var = "change"); options(width = 200)
cat(sprintf("\n%d sessions analysed\n", uniqueN(ch[, paste(subject, day)])))
for (wn in c("breath_hold_late", "hyperventilation")) { x <- w[window == wn]; cat(sprintf("\n== %s: change from baseline, z-scores (median [min, max]); sessions with an increase\n", wn))
  for (sg in c("MSNA", "x_s", "x_p", "HR")) cat(sprintf("  %-5s %5.2f [%5.2f, %5.2f]  up in %d of %d\n", sg, median(x[[sg]]), min(x[[sg]]), max(x[[sg]]), sum(x[[sg]] > 0), nrow(x))) }
x <- w[window == "breath_hold_late"]; cat(sprintf("\nH3, breath-hold: Spearman across sessions of the change of MSNA with x_s %.2f, with x_p %.2f, with HR %.2f\n", cor(x$MSNA, x$x_s, method = "spearman"), cor(x$MSNA, x$x_p, method = "spearman"), cor(x$MSNA, x$HR, method = "spearman")))
cat("\nLagged correlation with MSNA over the whole saline period (positive lag: MSNA later than the state); median over sessions, lag of the maximum |r|\n")
print(lg[, .(r_lag0 = round(median(r[lag == 0], na.rm = TRUE), 2), r_max = round(median(sapply(split(r, paste(subject, day)), function(v) v[which.max(abs(v))])), 2), best_lag_median = median(sapply(split(seq_along(r), paste(subject, day)), function(ix) lag[ix][which.max(abs(r[ix]))]))), by = signal], row.names = FALSE)
