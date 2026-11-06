# Analysis of the frozen validation (see protocol.md): paired percentile bootstrap over recordings, Holm adjustment per model. Run: Rscript experiments/regimes/v2/validation/analyze.R
suppressMessages(library(data.table)); options(width = 200); set.seed(1); B <- 10000L
x <- readRDS("experiments/regimes/v2/_out/validation.rds"); res <- x$res; held <- x$held[tail_bad < 0.1]
w <- function(v) dcast(res, id ~ model, value.var = v)
pac <- w("pacf"); ks <- w("ks_ok"); tv <- w("t50_vagal"); ts_ <- w("t50_symp"); dv <- w("d_vagal"); dsy <- w("d_symp"); mods <- c("HMMg", "HYB")
boot <- function(d, f = mean) { d <- d[is.finite(d)]; n <- length(d); b <- replicate(B, f(d[sample.int(n, n, TRUE)])); list(est = f(d), lo = quantile(b, .025), hi = quantile(b, .975), b = b, n = n) }
ep <- function(m) { out <- list()
  if (m == "HMMg") { d <- held$hmmg - held$m3; b <- boot(d); out$E1 <- c(b, p = mean(b$b <= -0.01), margin = -0.01, kind = "NI") }
  b <- boot(pac[[m]] - pac$M3); out$E2 <- c(b, p = mean(b$b >= 0), margin = 0, kind = "SUP")
  b <- boot(ks[[m]] - ks$M3); out$E3 <- c(b, p = mean(b$b <= 0), margin = 0, kind = "SUP")
  sgn <- function(m_) (dv[[m_]] < 0 & dsy[[m_]] > 0); fast <- function(m_) (tv[[m_]] < ts_[[m_]]); ok <- is.finite(tv[[m]]) & is.finite(ts_[[m]]) & is.finite(tv$M3) & is.finite(ts_$M3)
  b <- boot(ifelse(ok, as.numeric(sgn(m)) - as.numeric(sgn("M3")), NA)); out$E4 <- c(b, p = mean(b$b <= -0.05), margin = -0.05, kind = "NI")
  b <- boot(ifelse(ok, as.numeric(fast(m)) - as.numeric(fast("M3")), NA)); out$E5 <- c(b, p = mean(b$b <= -0.05), margin = -0.05, kind = "NI")
  tab <- rbindlist(lapply(names(out), function(e) { o <- out[[e]]; data.table(model = m, endpoint = e, rule = o$kind, n = o$n, estimate = round(as.numeric(o$est), 3), lo = round(as.numeric(o$lo), 3), hi = round(as.numeric(o$hi), 3), margin = as.numeric(o$margin), p = as.numeric(o$p)) }))
  tab[, p_holm := p.adjust(p, "holm")]; tab[, met := p_holm < 0.05]; tab }
tabs <- rbindlist(lapply(mods, ep)); verdict <- tabs[, { ni <- all(met[rule == "NI"]); sup <- any(met[rule == "SUP"]); .(verdict = if (sup && ni) "superior" else if (ni) "non-inferior" else "not established") }, by = model]
cat("\nPrimary endpoints (paired differences vs M3; p one-sided bootstrap, Holm within model; NI = non-inferiority against the margin, SUP = superiority against 0):\n"); print(tabs); cat("\nVerdict:\n"); print(verdict)
cat("\nSecondary:\n"); print(res[, .(n = .N, fallback_pct = round(100 * mean(fallback)), ks_pass = round(100 * mean(ks_ok)), pacf = round(mean(pacf), 1), dLL = round(median(dLL), 1), tau_p = round(median(tau_p), 1), t50_vagal = round(median(t50_vagal, na.rm = TRUE), 1), t50_symp = round(median(t50_symp, na.rm = TRUE), 1), secs = round(median(secs), 0)), by = model][order(match(model, c("M3", "M3E_det", "HMMg", "HYB")))])
d <- held$hmmg - held$m3; cat("\nForward difference HMMg - M3: n =", length(d), " median", round(median(d), 3), " trimmed mean", round(mean(d, trim = 0.1), 3), " better in", round(100 * mean(d > 0)), "% worst", round(min(d), 2), "; gate open on the head in", round(100 * mean(held$open_head)), "%\n")
cat("Gate open (full record):", round(100 * mean(unique(res[, .(id, open)])$open)), "% of", uniqueN(res$id), "recordings; failed recordings:", x$n_failed, "\n")
fwrite(tabs, "experiments/regimes/v2/_out/validation_endpoints.csv")
