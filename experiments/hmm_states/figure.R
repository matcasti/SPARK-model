# ==============================================================================
# experiments/hmm_states/figure.R -- two examples of the second-stage HMM (uncertainty-aware): one PRCP recording (K = 2, tilt angle on top) and one exercise recording (K = 3).
#   Rscript experiments/hmm_states/figure.R  -> _out/hmm_examples.png
# ==============================================================================
source("experiments/hmm_states/lib.R"); suppressMessages({ library(ggplot2); library(patchwork) })
ex <- function(ds, id, K, title) { f <- readRDS(sprintf("results/fit_%s.rds", ds)); s <- f$states[SubjectID == id]; p <- f$params[[id]]; b <- bin_states(s$Time, s$X_p, s$X_s, s$X_p_var, s$X_s_var); V <- c(p$lp * p$sig2, p$sig2)
  fit <- hmm_fit(as.matrix(b[, .(y1, y2)]), as.matrix(b[, .(v1, v2)]), V, K, "aware"); b[, regime := factor(fit$path)]; r <- rle(as.integer(fit$path)); e <- cumsum(r$lengths); st <- e - r$lengths + 1
  seg <- data.table(x0 = b$t[st], x1 = b$t[e], regime = factor(r$values)); long <- rbind(b[, .(t, state = "vagal-like x_p (Hz)", y = y1, lo = y1 - 1.96 * sqrt(v1), hi = y1 + 1.96 * sqrt(v1))], b[, .(t, state = "sympathetic-like x_s (Hz)", y = y2, lo = y2 - 1.96 * sqrt(v2), hi = y2 + 1.96 * sqrt(v2))])
  pl <- ggplot() + geom_rect(data = seg, aes(xmin = x0, xmax = x1, ymin = -Inf, ymax = Inf, fill = regime), alpha = 0.18) + geom_ribbon(data = long, aes(t, ymin = lo, ymax = hi), alpha = 0.25) + geom_line(data = long, aes(t, y)) +
    facet_wrap(~ state, ncol = 1, scales = "free_y") + scale_fill_brewer(palette = "Set1") + labs(x = "Time (s)", y = NULL, title = title, fill = "regime") + theme_minimal(base_size = 10) + theme(legend.position = "bottom")
  if ("angle" %in% names(s) && any(is.finite(s$angle))) pl <- pl + geom_vline(xintercept = s$Time[which(diff(as.integer(s$angle > 20)) == 1) + 1], linetype = 3, colour = "grey30"); pl }
ggsave("experiments/hmm_states/_out/hmm_examples.png", ex("prcp", "12734", 2, "PRCP 12734, two regimes (dotted: tilt onsets from the angle)") | ex("exercise", "196", 3, "Exercise 196, three regimes"), width = 12, height = 6.5, dpi = 150)
