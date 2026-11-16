# ==============================================================================
# 15_blockade.R -- does the reconstruction respond to pharmacological autonomic blockade, and does the sympathetic-like state follow muscle sympathetic nerve
# activity? Public Lancaster "Blocker data" (7 healthy volunteers; saline, atropine 0.04 mg/kg, propranolol 0.2 mg/kg; 5-min segments of spontaneous and
# 0.25-Hz paced breathing; hyperventilation and breath-hold with the integrated nerve signal in the saline period).
#
# Reads   data/blockade/ (tools/prep/blockade.R): protocol.csv, maneuvers.csv, rpeaks/, msna/  -- optional: without them the result records that they are absent
# Writes  results/blockade.rds   available, fits (one row per segment), ratios to saline, summary, maneuver changes and lagged correlations
# Method  R/steps/blockade.R; settings in cfg$blockade. Predictions fixed before the first run (supplementary material).
# ==============================================================================
source("analysis/00_setup.R")

if (!blockade_available(cfg)) {
  log_msg("Blockade data not found in %s: run tools/prep/blockade.R after downloading the Lancaster 'Blocker data' (see README). Recording that they are absent.", cfg$paths$blockade)
  write_result(list(available = FALSE), "blockade", cfg)
} else {
  # ---- 1. Segments: M3 on every 5-min segment -----------------------------------------------------------------------------------
  protocol <- fread(file.path(cfg$paths$blockade, "protocol.csv"))
  peaks <- lapply(split(protocol, by = c("subject", "day")), function(d) blockade_peaks(cfg, d$subject[1], d$day[1]))
  segments <- map_recordings(seq_len(nrow(protocol)), function(i) { r <- protocol[i]; blockade_segment_fit(r, peaks[[paste(r$subject, r$day, sep = ".")]], cfg$blockade) }, cfg$cores)
  fits <- rbindlist(drop_failed(segments, "segment"))
  ratios <- blockade_ratios(fits)
  summary <- blockade_ratio_summary(ratios, cfg$blockade)

  # ---- 2. Maneuvers: states against the nerve signal ----------------------------------------------------------------------------------
  maneuvers <- fread(file.path(cfg$paths$blockade, "maneuvers.csv"))
  sessions <- drop_failed(map_recordings(seq_len(nrow(maneuvers)), function(i) blockade_maneuver_session(maneuvers[i], protocol, cfg), cfg$cores), "session")
  sessions <- Filter(Negate(is.null), sessions)
  change <- rbindlist(lapply(sessions, `[[`, "change")); lagged <- rbindlist(lapply(sessions, `[[`, "lagged"))

  write_result(list(available = TRUE, n_subjects = uniqueN(fits$subject), n_segments = nrow(fits), n_segments_total = nrow(protocol),
                    fits = fits, ratios = ratios, summary = summary, change = change, lagged = lagged), "blockade", cfg)
  print(summary[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
  log_msg("Maneuver sessions analysed: %d", uniqueN(change[, paste(subject, day)]))
}
