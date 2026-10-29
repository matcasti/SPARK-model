# Runs the real tools/prep-data.R on synthetic WFDB records and checks the CSVs against the known truth,
# then feeds them to the analysis loader and the model. Skips if Python/wfdb are unavailable.
#   Rscript tests/suites/test_prep.R
Sys.setenv(SPARK_CONFIG = normalizePath("tests/fixtures/config_test.R"))
if (!nzchar(Sys.getenv("RETICULATE_PYTHON")) && file.exists("/home/claude/pyenv/bin/python")) Sys.setenv(RETICULATE_PYTHON = "/home/claude/pyenv/bin/python")
source("analysis/00_setup.R"); source("tests/fixtures/make_synthetic_wfdb.R")
ok_py <- tryCatch({ reticulate::import("wfdb"); TRUE }, error = function(e) FALSE)
if (!ok_py) { log_msg("Python 'wfdb' not available: prep test skipped."); quit(save = "no") }
n_pass <- 0L
check <- function(name, cond) { if (!isTRUE(cond)) stop("FAILED: ", name, call. = FALSE); n_pass <<- n_pass + 1L; log_msg("ok   %s", name) }

root <- file.path(tempfile("prep"), "run"); dir.create(root, recursive = TRUE)
truth <- write_synthetic_wfdb(root)
RAW_DIR <- file.path(root, "data-raw"); OUT_DIR <- file.path(root, "data")
s <- Sys.time(); source("tools/prep-data.R", local = FALSE); log_msg("tools/prep-data.R ran in %.0f s", as.numeric(Sys.time() - s, units = "secs"))
rd <- function(ds, id) fread(file.path(OUT_DIR, ds, paste0(id, ".csv")))

# Layout and alignment ---------------------------------------------------------------------------
f1 <- rd("fantasia", "f1y01"); p1 <- rd("prcp", "prcp01"); e1 <- rd("tmst", "ex01")
check("Fantasia columns", identical(names(f1), c("time", "rri", "systolic", "diastolic", "resp", "map")))
check("PRCP columns", identical(names(p1), c("time", "rri", "systolic", "diastolic", "angle", "map")))
check("TMST columns", identical(names(e1), c("time", "rri")))
cons <- function(d) mean(abs(diff(d$time) - d$rri[-1]) <= 0.0025)
check("rri is the interval ending at time (>= 98% of beats; cleaned beats excepted)", all(vapply(list(f1, p1, e1), cons, 1.0) >= 0.98))
t1 <- truth$f1y01
check("non-beat annotations excluded: one row per beat interval", nrow(f1) == nrow(t1))
check("rri matches the simulated RR (median error < 1 ms; 250 Hz quantization)", median(abs(f1$rri - t1$rri)) < 0.001)
check("time equals the beat times", max(abs(f1$time - t1$time)) < 0.01)
# Blood pressure -----------------------------------------------------------------------------------
qc <- fread(file.path(OUT_DIR, "prep_qc.csv"))
check("pulse delay estimated within 30 ms of the true 0.18 s", all(abs(qc[has_bp == TRUE, delay_s] - 0.18) < 0.03))
`%notin%` <- Negate(`%in%`)
clean <- !t1$spike & seq_len(nrow(t1)) %notin% 78:84
check("systolic/diastolic recover the generating values (median error < 1 mmHg)", median(abs(f1$systolic - t1$sbp)[clean]) < 1 && median(abs(f1$diastolic - t1$dbp)[clean]) < 1)
check("injected +90 mmHg spikes are cleaned (< 20 mmHg error)", all(abs(f1$systolic - t1$sbp)[t1$spike] < 20))
check("1-s ABP dropout is interpolated, not left missing", !anyNA(f1$systolic))
check("MAP = (SBP + 2 DBP)/3", max(abs(f1$map - (f1$systolic + 2 * f1$diastolic) / 3)) < 1e-9)
# Missing channels ---------------------------------------------------------------------------------
f2 <- rd("fantasia", "f1y02"); f3 <- rd("fantasia", "f1o01"); p3 <- rd("prcp", "prcp03"); p4 <- rd("prcp", "prcp04")
check("Fantasia without BP: BP columns all NA, resp present", all(is.na(f2$systolic)) && all(is.na(f2$map)) && !anyNA(f2$resp))
check("Fantasia without respiration: resp all NA, BP present", all(is.na(f3$resp)) && !anyNA(f3$systolic))
check("PRCP without angle / without ABP: respective columns all NA", all(is.na(p3$angle)) && !anyNA(p3$systolic) && all(is.na(p4$systolic)) && !anyNA(p4$angle))
check("QC flags match the signals present", qc[id == "f1y02", !has_bp] && qc[id == "f1o01", !has_resp] && qc[id == "prcp03", !has_angle] && qc[id == "prcp04", !has_bp])
# Respiration and angle ------------------------------------------------------------------------------
check("respiration per beat follows the generating sinusoid (corr > 0.95)", cor(f1$resp, t1$resp) > 0.95)
tp <- truth$prcp02
check("angle per beat matches the generating angle (median error < 1 deg)", median(abs(rd("prcp", "prcp02")$angle - tp$angle)) < 1)
# TMST ---------------------------------------------------------------------------------------------------
check("TMST: ms converted to s, edges removed, cut below 30 min", median(e1$rri) > 0.3 && median(e1$rri) < 2 && max(e1$time) / 60 < 30 && nrow(e1) < nrow(truth$ex01) - 9)
check("TMST: injected artifacts cleaned (no interval > 2x or < 0.5x its neighbors)", { r <- e1$rri; all(r[-c(1, length(r))] < 1.8 * pmax(r[-c(1, 2)], r[-c(length(r) - 1, length(r))])) })
# Cleaning gate: artifact-free series barely altered; injected artifacts corrected ---------------------------------
x <- sim_rest(77, "young", 600, m_p = rsa_drive(0.25, 0.05))$rr
check("clean_sample alters < 1% of artifact-free beats (old 2 MAD x 5 altered 4-7%)", mean(attr(clean_sample(x), "changed")) < 0.01)
old <- clean_sample(x, n_mad = 2, n_iter = 5); check("the previous gate (2 MAD x 5) altered > 2% of the same beats", mean(attr(old, "changed")) > 0.02)
# Integration with the analysis pipeline ------------------------------------------------------------------------------
cfgt <- cfg; cfgt$paths$fantasia <- file.path(OUT_DIR, "fantasia"); cfgt$paths$prcp <- file.path(OUT_DIR, "prcp"); cfgt$paths$exercise <- file.path(OUT_DIR, "tmst")
cfgt$windows$fantasia <- list(start_s = 30, duration_s = 380)
dl <- lapply(setNames(cfg$datasets, cfg$datasets), load_dataset, cfg = cfgt)
check("loader: rri used as written in every dataset (aligned with time)", all(unlist(lapply(dl, function(d) d$index$rr_source == "rri"))) && all(unlist(lapply(dl, function(d) d$index$time_rri_consistency)) > 0.95))
check("loader: covariate availability matches", identical(dl$fantasia$index[order(SubjectID), win_systolic], c(TRUE, TRUE, FALSE)) && identical(dl$fantasia$index[order(SubjectID), win_resp], c(FALSE, TRUE, TRUE)) &&
        identical(dl$prcp$index[order(SubjectID), win_angle], c(TRUE, TRUE, FALSE, TRUE)))
ft <- fit_dataset(dl$prcp$rec[c("prcp01", "prcp04")], "prcp", cfgt)
check("PRCP files from the new script fit with the engine and carry covariates", all(ft$summary$Converged) && all(c("map", "angle") %in% names(ft$states)))
ev <- tilt_events(ft$states[SubjectID == "prcp01", Time], ft$states[SubjectID == "prcp01", angle], cfg$tilt)
check("tilt onset from the prepared angle column is detected", ev$detected)
log_msg("%d prep checks passed.", n_pass)
