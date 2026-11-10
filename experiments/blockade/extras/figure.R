# ==============================================================================
# experiments/blockade/extras/figure.R -- figure of the comparison: non-vagal share R(f) of the blockade data against the production partition and the fitted forms (A),
# and the anchor lambda_P implied by each weighting on the 310 production recordings (B). Rscript experiments/blockade/extras/figure.R  -> _out/blockade_shares.png
# ==============================================================================
suppressMessages({ library(data.table); library(ggplot2); library(patchwork) })
mk <- function(tag, lab) { o <- readRDS(sprintf("experiments/blockade/_out/forms_%s.rds", tag)); o$m[, setting := lab]; o }
a <- mk("lnrr_spontaneous", "ln RR, spontaneous (6 subjects)"); b <- mk("lnrr_paced_0.25", "ln RR, paced 0.25 Hz (7 subjects)")
emp <- rbind(a$m, b$m); fl <- melt(a$m[, .(f, `production logistic` = prod, `free logistic` = best_logistic, `Lorentzian + floor` = lorentz_floor)], id.vars = "f", variable.name = "form", value.name = "R")
pa <- ggplot() + geom_ribbon(data = emp, aes(f, ymin = lo, ymax = hi, fill = setting), alpha = 0.15) + geom_point(data = emp, aes(f, R_emp, colour = setting), size = 2) +
  geom_line(data = fl, aes(f, R, linetype = form), linewidth = 0.8) + scale_x_log10() + scale_y_log10(limits = c(1e-4, 1)) +
  scale_colour_manual(values = c("#1b7837", "#762a83"), name = NULL) + scale_fill_manual(values = c("#1b7837", "#762a83"), name = NULL) + scale_linetype_manual(values = c(2, 1, 3), name = NULL) +
  geom_vline(xintercept = c(0.04, 0.40), colour = "grey50", linewidth = 0.3) + annotate("text", x = 0.045, y = 3e-4, label = "production range", hjust = 0, size = 3, colour = "grey40") +
  labs(x = "Frequency (Hz)", y = "Non-vagal share R(f) = P(atropine) / P(saline)", title = "A. What atropine leaves of the power", subtitle = "points: geometric mean over subjects; band: P10-P90 across subjects") + theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.box = "vertical")
r <- readRDS("experiments/blockade/_out/anchor.rds")$recordings; lam <- melt(r[, .(dataset, `production (stored)` = stored, `blockade: free logistic` = logistic, `blockade: Lorentzian + floor` = lorentz)], id.vars = "dataset", variable.name = "weights", value.name = "lambda")
pb <- ggplot(lam, aes(weights, lambda, fill = weights)) + geom_boxplot(outlier.size = 0.6, show.legend = FALSE) + scale_y_log10() + facet_wrap(~ dataset) + scale_fill_manual(values = c("#bdbdbd", "#fdae61", "#d73027")) +
  labs(x = NULL, y = expression(lambda[P] == V[p] / V[s] ~ "(log scale)"), title = "B. Anchor implied on the 310 recordings") + theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave("experiments/blockade/_out/blockade_shares.png", pa + pb + plot_layout(widths = c(1.15, 1)), width = 12, height = 5.6, dpi = 160)
