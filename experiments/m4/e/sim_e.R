# ==============================================================================
# experiments/m4/e/sim_e.R -- EXPERIMENTAL. M4e components on simulations with known truth (README.md).
#   Rscript experiments/m4/e/sim_e.R dev|conf [N per scenario = 40] [tag]
# Same five scenarios as ../sim.R; seeds "m4e-dev" or "m4e-conf". Per component and mode: normalized RMSE, coverage on the factor grid COV_GRID, |ln error| of the time constants.
# ==============================================================================
args <- commandArgs(TRUE); phase <- args[1]; N <- if (length(args) > 1) as.integer(args[2]) else 40L; tag <- if (length(args) > 2) args[3] else ""
stopifnot(phase %in% c("dev", "conf"))
source("analysis/00_setup.R"); source("experiments/m4/e/lib_e.R"); n_cores <- max(1L, min(7L, cfg$cores)); out <- "experiments/m4/e/_out"; dir.create(out, showWarnings = FALSE)
SPEC_EXTRA <- if (exists("SPEC_EXTRA_E")) SPEC_EXTRA_E else list()
OSC <- list(alpha = 0.2, f_hz = 0.25, lo = 0.3)
SCEN <- list(rest = list(gen = "rest", osc = NULL), rsa = list(gen = "rest", osc = OSC), regime = list(gen = "regime", osc = OSC), exercise = list(gen = "exercise", osc = OSC), two = list(gen = "two", osc = OSC))
jobs <- unlist(lapply(names(SCEN), function(sc) lapply(seq_len(N), function(i) list(scenario = sc, sim = i, seed = id_seed(paste0("m4e-", phase, " ", sc, " ", i), cfg$seed)))), recursive = FALSE)
make <- function(job) { z <- sim_regime(SCEN[[job$scenario]]$gen, job$seed, osc = SCEN[[job$scenario]]$osc); s <- z$s
  list(z = z, series = prepare_series(s$rr), truth_p = s$x_p + s$z1, truth_s = s$x_s) }
score <- function(m, d, comp) rbindlist(lapply(c(FALSE, TRUE), function(sm) { e <- m_states(m, sm); ep <- abs(e$xp - d$truth_p) / sqrt(pmax(e$xpv, 1e-12)); es <- abs(e$xs - d$truth_s) / sqrt(pmax(e$xsv, 1e-12))
  r <- data.table(comp = comp, mode = if (sm) "smooth" else "filter", rmse_p = sqrt(mean((e$xp - d$truth_p)^2)) / sd(d$truth_p), rmse_s = sqrt(mean((e$xs - d$truth_s)^2)) / sd(d$truth_s),
                  err_tau_p = abs(log(m$fit$tau_p / d$z$tau_p)), err_tau_s = abs(log(m$fit$tau_s / d$z$tau_s)), cost = unname(m$fit$cost), err_lp = abs(log(m$lp / d$z$lp)))
  for (c in COV_GRID) r[, (sprintf("cp%02d", round(10 * c))) := mean(ep <= 1.96 * c)][, (sprintf("cs%02d", round(10 * c))) := mean(es <= 1.96 * c)]; r }))
run1 <- function(job) { d <- make(job); fc <- fit_components(d$series, SPEC_EXTRA)
  rows <- rbindlist(lapply(names(fc$models), function(nm) score(fc$models[[nm]], d, nm)))
  list(rows = rows[, `:=`(scenario = job$scenario, sim = job$sim)], stats = data.table(scenario = job$scenario, sim = job$sim, lr_osc = fc$stats[["lr_osc"]], gain = fc$stats[["gain"]], have = paste(names(fc$models), collapse = ","))) }
log_msg("M4e simulations (%s%s): %d jobs", phase, tag, length(jobs)); s1 <- drop_failed(map_recordings(jobs, run1, n_cores), "simulation")
res <- list(rows = rbindlist(lapply(s1, `[[`, "rows")), stats = rbindlist(lapply(s1, `[[`, "stats")))
saveRDS(res, file.path(out, sprintf("sim_%s%s.rds", phase, tag))); log_msg("saved sim_%s%s.rds (%d rows)", phase, tag, nrow(res$rows))
