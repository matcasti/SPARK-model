# ==============================================================================
# experiments/blockade/extras/figure_ou.R -- spectra of ln RR under saline, atropine and double blockade (geometric mean over subjects) against the fitted model spectrum
# and its components (tau_s = 100 s, kappa = 1). Rscript experiments/blockade/extras/figure_ou.R  -> _out/ou_fit_spectra.png
# ==============================================================================
source("experiments/blockade/lib.R"); suppressMessages({ library(ggplot2); library(patchwork) })
r <- readRDS("experiments/blockade/_out/ou_fit.rds")[kappa_mode == "fixed 1" & tau_s == 100]; seg <- fread("experiments/blockade/_out/protocol.csv"); edges <- exp(seq(log(0.011), log(0.5), length.out = 22))
one <- function(br) { out <- list()
  for (s in r[breathing == br, subject]) { p <- r[breathing == br & subject == s]
    g <- function(cn) { x <- seg[subject == s & cond == cn & day == 1 & breathing == br]; tk <- clean_peaks(scan(sprintf("experiments/blockade/_out/rpeaks/rpeaks_Subject%d_day1.csv", s), quiet = TRUE)); pgram(tk, x$t0[1], x$t1[1]) }
    d <- list(saline = g("saline"), atropine = g("atropine_high"), double = g("double")); f <- d$saline$f
    mod <- list(saline = lor(f, p$V_p, p$tau_p) + lor(f, p$V_s, 100) + peak(f, p$A, p$f_r, p$w) + p$floor_sal, atropine = lor(d$atropine$f, p$V_s, 100) + p$floor_atr, double = rep(p$floor_dbl, length(d$double$f)))
    comp <- list(`vagal OU` = lor(f, p$V_p, p$tau_p), `respiratory peak` = peak(f, p$A, p$f_r, p$w), `slow OU (sympathetic-like)` = lor(f, p$V_s, 100))
    for (cn in names(d)) out[[length(out) + 1L]] <- data.table(subject = s, cond = cn, f = d[[cn]]$f, obs = d[[cn]]$I, fit = mod[[cn]])
    for (cm in names(comp)) out[[length(out) + 1L]] <- data.table(subject = s, cond = cm, f = f, obs = NA_real_, fit = comp[[cm]]) }
  x <- rbindlist(out); x[, b := cut(f, edges, labels = FALSE)]; x[, .(f = exp(mean(log(f))), obs = exp(mean(log(obs))), fit = exp(mean(log(fit)))), by = .(cond, b)][!is.na(b)][, breathing := br] }
d <- rbind(one("spontaneous"), one("paced_0.25")); main <- d[cond %in% c("saline", "atropine", "double")]; comp <- d[!cond %in% c("saline", "atropine", "double")]
p1 <- ggplot(main, aes(f, colour = cond)) + geom_point(aes(y = obs)) + geom_line(aes(y = fit)) + facet_wrap(~ breathing) + scale_x_log10() + scale_y_log10() + labs(x = "Frequency (Hz)", y = "Spectrum of ln RR (1/Hz)", colour = NULL, title = "Observed (points) and fitted model spectrum (lines)") + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p2 <- ggplot(rbind(comp[, .(cond, f, fit, breathing)], main[cond == "saline", .(cond = "total (saline)", f, fit, breathing)]), aes(f, fit, colour = cond)) + geom_line(linewidth = 0.9) + facet_wrap(~ breathing) + scale_x_log10() + scale_y_log10(limits = c(1e-6, NA)) +
  labs(x = "Frequency (Hz)", y = NULL, colour = NULL, title = "Fitted components (saline)") + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("experiments/blockade/_out/ou_fit_spectra.png", p1 / p2, width = 10, height = 9, dpi = 150)
