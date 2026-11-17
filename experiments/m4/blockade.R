# ==============================================================================
# experiments/m4/blockade.R -- EXPERIMENTAL. The M4 ladder on the Lancaster blockade segments (README.md, criterion B). SD_p = vagal-like total (OU + oscillator).
#   Rscript experiments/m4/blockade.R
# M4d-H: empirical-Bayes priors from the weak-prior fits of the segments of the OTHER subjects.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/m4/lib.R"); n_cores <- max(1L, min(6L, cfg$cores - 1L)); out <- "experiments/m4/_out"
protocol <- fread(file.path(cfg$paths$blockade, "protocol.csv")); st <- cfg$blockade
seg_rr <- function(g) { tk <- blockade_peaks(cfg, g$subject, g$day); k <- tk >= g$t0 + st$trim_s & tk <= g$t1; rr <- diff(tk[k])
  if (length(rr) < st$min_beats || mean(rr > 0.25 & rr < 2.5) < 0.97) NULL else rr[rr > 0.25 & rr < 2.5] }
run_seg <- function(j, eb = NULL) { g <- protocol[j]; rr <- seg_rr(g); if (is.null(rr)) return(NULL); series <- prepare_series(rr)
  fam <- fit_family(series, eb = eb, stage1 = is.null(eb)); keep <- if (is.null(eb)) setdiff(intersect(MODELS, names(fam)), c("M4dH", "M4dI")) else intersect("M4dH", names(fam))
  rows <- rbindlist(lapply(keep, function(nm) { m <- fam[[nm]]; s <- m_sd(m)
    data.table(model = nm, subject = g$subject, day = g$day, cond = g$cond, breathing = g$breathing, n = length(rr), hr = 60 / mean(rr), sd_rr_ms = 1000 * sd(rr),
               tau_p = m$fit$tau_p, tau_s = m$fit$tau_s, lambda_P = m$lp, SD_p = s[["SD_p"]], SD_s = s[["SD_s"]], gate_open = isTRUE(m$gate_open)) }))
  list(rows = rows, weak = if (is.null(eb)) data.table(j = j, subject = g$subject, t(fam$weak)) else NULL) }
log_msg("Blockade: %d segments, stage 1", nrow(protocol)); s1 <- drop_failed(map_recordings(seq_len(nrow(protocol)), run_seg, n_cores), "segment")
R1 <- rbindlist(lapply(s1, `[[`, "rows")); W <- rbindlist(lapply(s1, `[[`, "weak"))
R2 <- rbindlist(lapply(drop_failed(map_recordings(W$j, function(jj) { sj <- protocol$subject[jj]; w <- W[subject != sj]; run_seg(jj, eb_from(as.matrix(w[, .(tau_p, tau_s)]))) }, n_cores), "segment"), `[[`, "rows"))
B <- rbind(R1, R2); saveRDS(B, file.path(out, "blockade.rds"))
cat("\n== Ratio to saline of the same day (geometric mean [bootstrap 95% CI]) ==\n")
S <- rbindlist(lapply(intersect(MODELS, unique(B$model)), function(mm) blockade_ratio_summary(blockade_ratios(B[model == mm]), st)[, model := mm]))
print(S[, .(model, breathing, cond, n, SD_p = sprintf("%.2f [%.2f, %.2f]", SD_p, SD_p_lo, SD_p_hi), SD_s = sprintf("%.2f [%.2f, %.2f]", SD_s, SD_s_lo, SD_s_hi))][order(breathing, cond, match(model, MODELS))], row.names = FALSE)
cat("\nGate open (share of segments):", B[model == "M4c", round(mean(gate_open), 2)], "\n")
saveRDS(S, file.path(out, "blockade_summary.rds"))
