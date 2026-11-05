# ==============================================================================
# experiments/parametrization/run.R -- EXPERIMENTAL. Compares parametrizations/optimizers of the M3 fit.
# Run from the project root:  Rscript experiments/parametrization/run.R [quick|full] [out_name]
# Output: experiments/parametrization/_out/<out_name>.rds (git-ignored) and a summary on screen.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) args[1] else "quick"
out_name <- if (length(args) >= 2) args[2] else paste0("run_", mode)
source("analysis/00_setup.R")
source("experiments/parametrization/strategies.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))

# ---- Recordings ----------------------------------------------------------------------------------------------------
jobs <- list()
for (ds in cfg$datasets) {
  d <- load_dataset(ds, cfg)
  ids <- names(d$rr)
  if (mode == "quick") { set.seed(1); ids <- sort(sample(ids, min(length(ids), c(prcp = 5, fantasia = 12, exercise = 18)[[ds]]))) }
  for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]])
}
log_msg("Parametrization experiment | %s | %d recordings | %d cores", mode, length(jobs), n_cores)

# ---- Strategies: parametrization | optimizer ---------------------------------------------------------------------
LOCAL <- list(c("base", "bfgs"), c("mean_logit", "bfgs"), c("fast_ratio", "bfgs"), c("whitened", "bfgs"),
              c("base", "nlminb"), c("whitened", "nlminb"), c("base", "bfgs_nm"), c("base", "lbfgsb"))
GLOBAL <- list(c("base", "grid"), c("whitened", "grid"), c("base", "profile"), c("whitened", "profile"))

one_recording <- function(job) {
  series <- prepare_series(job$rr)
  pri <- series$priors; ctl <- series$control
  obj <- make_objective(series); P <- parametrizations(pri)
  idx <- match(job$id, vapply(jobs, `[[`, "", "id"))
  set.seed(1000L + idx)
  default <- c(pri$ks_mode, pri$ks_mode + max(pri$kp_mode - pri$ks_mode, ctl$kp_gap_floor), ctl$lR_start)
  starts <- c(list(default), dispersed_starts(pri, 6))
  rows <- list()
  add <- function(strategy, start, r) {
    rows[[length(rows) + 1]] <<- data.table(ds = job$ds, id = job$id, n_beats = length(series$dy), strategy = strategy, start = start,
      ok = isTRUE(r$ok), obj = if (isTRUE(r$ok)) r$v else NA_real_, conv = if (isTRUE(r$ok)) r$conv else NA_integer_,
      evals = r$evals, secs = r$secs, tau_p = if (isTRUE(r$ok)) 1 / r$kp else NA_real_, tau_s = if (isTRUE(r$ok)) 1 / r$ks else NA_real_,
      lR = if (isTRUE(r$ok)) r$lR else NA_real_)
  }
  for (s in LOCAL) for (j in seq_along(starts))
    add(paste(s, collapse = "|"), j - 1L, run_local(obj, P[[s[1]]], starts[[j]], s[2]))
  for (s in GLOBAL) {
    r <- if (s[2] == "grid") run_grid(obj, P[[s[1]]], pri) else run_profile(obj, P[[s[1]]], default)
    add(paste(s, collapse = "|"), 0L, r)
  }
  # best-known reference: BFGS + Nelder-Mead polish from 24 further dispersed starts
  set.seed(5000L + idx)
  for (j in seq_len(24)) add("ref", j, run_local(obj, P$base, dispersed_starts(pri, 1)[[1]], "bfgs_nm"))
  rbindlist(rows)
}

res <- parallel::mclapply(jobs, function(j) tryCatch(one_recording(j), error = function(e) { message("FAILED ", j$id, ": ", conditionMessage(e)); NULL }),
                          mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
res[, best := min(obj[ok]), by = .(ds, id)]
res[, gap := obj - best]
dir.create("experiments/parametrization/_out", showWarnings = FALSE)
saveRDS(res, file.path("experiments/parametrization/_out", paste0(out_name, ".rds")))
log_msg("Saved %d rows | %d recordings", nrow(res), uniqueN(res[, paste(ds, id)]))
source("experiments/parametrization/summary.R")
print_summary(res)
