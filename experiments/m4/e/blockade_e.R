# ==============================================================================
# experiments/m4/e/blockade_e.R -- EXPERIMENTAL. M4e components on the 81 Lancaster blockade segments (confirmation only; README.md).
#   Rscript experiments/m4/e/blockade_e.R [tag]
# ==============================================================================
args <- commandArgs(TRUE); tag <- if (length(args)) args[1] else ""
source("analysis/00_setup.R"); source("experiments/m4/e/lib_e.R"); n_cores <- max(1L, min(7L, cfg$cores)); out <- "experiments/m4/e/_out"
SPEC_EXTRA <- if (exists("SPEC_EXTRA_E")) SPEC_EXTRA_E else list()
protocol <- fread(file.path(cfg$paths$blockade, "protocol.csv")); st <- cfg$blockade
seg_rr <- function(g) { tk <- blockade_peaks(cfg, g$subject, g$day); k <- tk >= g$t0 + st$trim_s & tk <= g$t1; rr <- diff(tk[k])
  if (length(rr) < st$min_beats || mean(rr > 0.25 & rr < 2.5) < 0.97) NULL else rr[rr > 0.25 & rr < 2.5] }
run_seg <- function(jj) { g <- protocol[jj]; rr <- seg_rr(g); if (is.null(rr)) return(NULL); fc <- fit_components(prepare_series(rr), SPEC_EXTRA)
  rbindlist(lapply(names(fc$models), function(nm) { m <- fc$models[[nm]]; s <- m_sd(m)
    data.table(comp = nm, seg = jj, subject = g$subject, day = g$day, cond = g$cond, breathing = g$breathing, n = length(rr), hr = 60 / mean(rr), sd_rr_ms = 1000 * sd(rr),
               tau_p = m$fit$tau_p, tau_s = m$fit$tau_s, lambda_P = m$lp, SD_p = s[["SD_p"]], SD_s = s[["SD_s"]], cost = unname(m$fit$cost), lr_osc = fc$stats[["lr_osc"]], gain = fc$stats[["gain"]]) })) }
log_msg("M4e blockade%s: %d segments", tag, nrow(protocol)); B <- rbindlist(drop_failed(map_recordings(seq_len(nrow(protocol)), run_seg, n_cores), "segment"))
saveRDS(B, file.path(out, sprintf("blockade%s.rds", tag))); log_msg("saved blockade%s.rds (%d rows)", tag, nrow(B))
