# ==============================================================================
# experiments/identifiability/summary.R -- table and figure of the identifiability grid (_out/sim.rds from sim.R).   Rscript experiments/identifiability/summary.R
# ==============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) }); options(width = 200)
r <- readRDS("experiments/identifiability/_out/sim.rds"); r[, collapse := est_share < 0.05 | est_share > 0.95]
m3 <- r[spec == "M3", .(cor_xp = median(cor_x_p), cor_xs = median(cor_x_s), cover_xs = mean(cover_x_s), abs_ln_tau_p = median(abs(ln_err_tau_p)), abs_ln_tau_s = median(abs(ln_err_tau_s)), share_err = median(abs(share_err)), n = .N), by = .(axis, ratio, lambda, duration)]
cat("M3 (the model of the paper), median over 30 simulations per cell\n"); print(m3[axis == "time-scale ratio"][order(ratio, duration)][, !c("axis", "lambda", "n")][, lapply(.SD, function(x) round(x, 2))], row.names = FALSE)
print(m3[axis == "variance ratio"][order(lambda, duration)][, !c("axis", "ratio", "n")][, lapply(.SD, function(x) round(x, 2))], row.names = FALSE)
cat("\nShare of fits that put (almost) all variance in one drive, by specification (all cells with tau_s/tau_p >= 10 and 10-30 min vs < 10 or 5 min)\n")
r[, regime := ifelse(axis == "time-scale ratio" & ratio >= 10 & duration >= 600, "ratio >= 10, >= 10 min", ifelse(axis == "time-scale ratio" & ratio < 10, "ratio < 10", NA_character_))]
print(r[!is.na(regime), .(collapse = sprintf("%d%%", round(100 * mean(collapse))), cor_xs = round(median(cor_x_s), 2)), by = .(regime, spec)][order(regime, spec)], row.names = FALSE)
spec_col <- c(M0 = "#bdbdbd", M1 = "#fdae61", M2 = "#74add1", M3 = "#d73027")
a <- r[axis == "time-scale ratio", .(cor = median(cor_x_s)), by = .(ratio, duration, spec)]; b <- r[axis == "time-scale ratio", .(col = mean(collapse)), by = .(ratio, duration, spec)]
p1 <- ggplot(a, aes(ratio, cor, colour = spec)) + geom_line() + geom_point() + facet_wrap(~ factor(paste(duration / 60, "min"), levels = c("5 min", "10 min", "30 min"))) + scale_x_log10(breaks = c(1.5, 4, 10, 30, 100)) + scale_colour_manual(values = spec_col) + ylim(0, 1) +
  labs(x = expression(tau[s] / tau[p]), y = "Correlation of the sympathetic-like state with the truth", colour = NULL, title = "A. Time-scale ratio (lambda_P = 2)") + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
p2 <- ggplot(b, aes(ratio, col, colour = spec)) + geom_line() + geom_point() + facet_wrap(~ factor(paste(duration / 60, "min"), levels = c("5 min", "10 min", "30 min"))) + scale_x_log10(breaks = c(1.5, 4, 10, 30, 100)) + scale_colour_manual(values = spec_col) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = expression(tau[s] / tau[p]), y = "Fits with one drive taking all the variance", colour = NULL, title = "B. Collapse to one process") + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("experiments/identifiability/_out/identifiability.png", p1 / p2, width = 9, height = 8, dpi = 150)
