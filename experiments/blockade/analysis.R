# ==============================================================================
# experiments/blockade/analysis.R -- empirical vagal / sympathetic shares per frequency bin from the Lancaster blockade data, against the production partition.
#   Rscript experiments/blockade/analysis.R [rate|lnrr]        needs experiments/blockade/_out/segments.rds (see lib.R: build_segments)
# Reference for atropine and double blockade: saline of day 1; for propranolol alone: saline of day 2 (different day). Ratios are of summed power over subjects
# (power weighted) with a bootstrap over subjects. Output: _out/shares_<units>_<breathing>.rds
# ==============================================================================
source("experiments/blockade/lib.R")
args <- commandArgs(TRUE); u_sel <- if (length(args)) args[1] else "rate"; B <- 2000; set.seed(1)
seg <- readRDS("experiments/blockade/_out/segments.rds")[units == u_sel]
P <- function(cnd, dy, br) seg[cond == cnd & day == dy & breathing == br, .(subject, bin, power)]
slope <- -log(0.05 / 0.95) / (0.15 - 0.10); L <- function(f) 1 / (1 + exp(-slope * (f - 0.10)))
prod_w <- function(f) 0.25 + 0.75 * L(f)                                         # production vagal weight; the sympathetic weight is 1 - w

shares <- function(br) {
  d <- merge(merge(merge(merge(merge(P("saline", 1, br)[, .(subject, bin, sal1 = power)], P("saline", 2, br)[, .(subject, bin, sal2 = power)], all = TRUE),
          P("atropine_high", 1, br)[, .(subject, bin, atrH = power)], all = TRUE), P("atropine_low", 1, br)[, .(subject, bin, atrL = power)], all = TRUE),
          P("double", 1, br)[, .(subject, bin, dbl = power)], all = TRUE), P("propranolol", 2, br)[, .(subject, bin, prop = power)], all = TRUE)
  ratio <- function(num, den, dd) { k <- dd[is.finite(get(num)) & is.finite(get(den))]; ids <- unique(k$subject)
    est <- function(s) { kk <- k[subject %in% s]; kk[, sum(get(num)), by = bin][order(bin)]$V1 / kk[, sum(get(den)), by = bin][order(bin)]$V1 }
    bs <- replicate(B, est(sample(ids, replace = TRUE))); list(est = est(ids), lo = apply(bs, 1, quantile, .025, na.rm = TRUE), hi = apply(bs, 1, quantile, .975, na.rm = TRUE), n = length(ids)) }
  out <- rbindlist(lapply(list(list("atrH", "sal1"), list("atrL", "sal1"), list("dbl", "sal1"), list("prop", "sal2"), list("sal2", "sal1")), function(z) { r <- ratio(z[[1]], z[[2]], d)
    data.table(ratio = paste0(z[[1]], "/", z[[2]]), n = r$n, bin = BINS$bin, f = BINS$f, est = r$est, lo = r$lo, hi = r$hi) }))
  # sympathetic share left under vagal blockade: (atrH - dbl) / sal1, same subjects
  k <- d[is.finite(atrH) & is.finite(dbl) & is.finite(sal1)]; ids <- unique(k$subject)
  est <- function(s) { kk <- k[subject %in% s]; (kk[, sum(atrH), by = bin][order(bin)]$V1 - kk[, sum(dbl), by = bin][order(bin)]$V1) / kk[, sum(sal1), by = bin][order(bin)]$V1 }
  bs <- replicate(B, est(sample(ids, replace = TRUE)))
  out <- rbind(out, data.table(ratio = "(atrH-dbl)/sal1", n = length(ids), bin = BINS$bin, f = BINS$f, est = est(ids), lo = apply(bs, 1, quantile, .025), hi = apply(bs, 1, quantile, .975)))
  out[, breathing := br]; list(table = out, data = d)
}
res <- lapply(c(spontaneous = "spontaneous", paced = "paced_0.25"), shares)
for (br in names(res)) { cat(sprintf("\n===== %s breathing (%s): ratio of power per bin, pooled over subjects; [95%% bootstrap CI] =====\n", br, u_sel)); tb <- res[[br]]$table
  w <- dcast(tb, f ~ ratio, value.var = c("est", "lo", "hi")); setorder(w, f)
  f3 <- function(e, l, h) sprintf("%.3f [%.3f, %.3f]", e, l, h)
  print(data.table(f = round(w$f, 3), atrH_over_sal = f3(w$`est_atrH/sal1`, w$`lo_atrH/sal1`, w$`hi_atrH/sal1`), dbl_over_sal = f3(w$`est_dbl/sal1`, w$`lo_dbl/sal1`, w$`hi_dbl/sal1`),
                   prop_over_sal2 = sprintf("%.2f", w$`est_prop/sal2`), sal2_over_sal1 = sprintf("%.2f", w$`est_sal2/sal1`), symp_residual = sprintf("%.4f", w$`est_(atrH-dbl)/sal1`),
                   prod_w_symp = sprintf("%.3f", 1 - prod_w(w$f))), row.names = FALSE) }
saveRDS(res, sprintf("experiments/blockade/_out/shares_%s.rds", u_sel))
