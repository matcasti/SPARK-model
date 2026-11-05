# ==============================================================================
# experiments/parametrization/rounds/run2.R -- EXPERIMENTAL, round 2: cheaper ways to get the multi-start result.
# Baseline to beat: whitened|nlminb from the default start + 3 dispersed starts (100% of best-known, ~263 evaluations).
# Run from the project root:  Rscript experiments/parametrization/rounds/run2.R [full|quick]
# ==============================================================================
mode <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "full"
source("analysis/00_setup.R"); source("experiments/parametrization/strategies.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
prev <- readRDS("experiments/parametrization/_out/run_full.rds")
prev_best <- unique(prev[, .(ds, id, best)])

jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); for (id in names(d$rr)) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
all_ids <- vapply(jobs, `[[`, "", "id")
if (mode == "quick") { set.seed(2); jobs <- jobs[sort(sample(length(jobs), 40))] }
log_msg("Round 2 | %s | %d recordings | %d cores", mode, length(jobs), n_cores)

# nlminb in the prior-whitened coordinates; returns the point so it can be continued
nl <- function(obj, par, z0, iter.max = 1000, rel.tol = 1e-10) {
  fz <- function(z) { p <- par$to(z); obj$f(p[1], p[2], p[3]) }
  o <- tryCatch(nlminb(z0, fz, control = list(iter.max = iter.max, eval.max = 2 * iter.max, rel.tol = rel.tol)), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$objective)) return(list(ok = FALSE, v = Inf, z = z0))
  list(ok = TRUE, v = o$objective, z = o$par)
}
phys_of <- function(par, z) par$to(z)

# data-informed start: AR(1) time constant of the detrended RR series for tau_p, defaults elsewhere
data_start <- function(series, default) {
  dy <- series$dy; k <- 21
  tr <- as.numeric(stats::filter(dy, rep(1 / k, k), sides = 2)); x <- (dy - tr)[is.finite(tr)]
  r1 <- suppressWarnings(acf(x, lag.max = 1, plot = FALSE)$acf[2]); if (!is.finite(r1)) r1 <- 0.5
  kp <- -log(min(max(r1, 0.05), 0.95)) / mean(dy)
  c(default[1], max(kp, 1.6 * default[1]), default[3])
}

one <- function(job) {
  series <- prepare_series(job$rr); pri <- series$priors; ctl <- series$control
  obj <- make_objective(series); W <- parametrizations(pri)$whitened
  idx <- match(job$id, all_ids); set.seed(1000L + idx)
  default <- c(pri$ks_mode, pri$ks_mode + max(pri$kp_mode - pri$ks_mode, ctl$kp_gap_floor), ctl$lR_start)
  starts <- c(list(default), dispersed_starts(pri, 6))           # identical to round 1
  rows <- list()
  record <- function(method, res_v, evals, secs) rows[[length(rows) + 1]] <<- data.table(ds = job$ds, id = job$id, method = method, obj = res_v, evals = evals, secs = secs)
  timed <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- force(expr); list(v = v, secs = proc.time()[["elapsed"]] - t0) }
  fz0 <- function(s) W$from(s[1], s[2], s[3])

  # ---- full runs from each start (also used post hoc for the adaptive stopping rule) ----
  full <- lapply(starts, function(s) { obj$counter$n <- 0L; r <- nl(obj, W, fz0(s)); r$evals <- obj$counter$n; r })
  for (j in seq_along(full)) record(sprintf("full_start%d", j - 1L), full[[j]]$v, full[[j]]$evals, NA)

  # ---- racing: K starts x m iterations, continue only the best ----
  for (cfgr in list(c(4, 4), c(4, 8), c(6, 4), c(6, 8), c(4, 15))) {
    K <- cfgr[1]; m <- cfgr[2]; obj$counter$n <- 0L; t0 <- proc.time()[["elapsed"]]
    short <- lapply(starts[seq_len(K)], function(s) nl(obj, W, fz0(s), iter.max = m))
    b <- which.min(vapply(short, `[[`, 0, "v")); fin <- nl(obj, W, short[[b]]$z)
    record(sprintf("race_K%d_m%d", K, m), fin$v, obj$counter$n, proc.time()[["elapsed"]] - t0)
  }
  # ---- loose tolerance on 4 starts, polish the best ----
  for (tol in c(1e-4, 1e-6)) {
    obj$counter$n <- 0L; t0 <- proc.time()[["elapsed"]]
    rough <- lapply(starts[1:4], function(s) nl(obj, W, fz0(s), rel.tol = tol))
    b <- which.min(vapply(rough, `[[`, 0, "v")); fin <- nl(obj, W, rough[[b]]$z)
    record(sprintf("loose4_tol%g", tol), fin$v, obj$counter$n, proc.time()[["elapsed"]] - t0)
  }
  # ---- designed (deterministic) second/third starts for the two failure modes seen in round 1 ----
  des <- list(default, c(1 / 150, 1 / 25, 0.3), c(1 / 160, 1 / 12, 1e-5))
  obj$counter$n <- 0L; t0 <- proc.time()[["elapsed"]]
  rs <- lapply(des, function(s) nl(obj, W, fz0(s))); record("design3", min(vapply(rs, `[[`, 0, "v")), obj$counter$n, proc.time()[["elapsed"]] - t0)
  obj$counter$n <- 0L; rs2 <- lapply(des[1:2], function(s) nl(obj, W, fz0(s))); record("design2_slow_noisy", min(vapply(rs2, `[[`, 0, "v")), obj$counter$n, NA)
  # ---- data-informed start ----
  ds0 <- data_start(series, default)
  obj$counter$n <- 0L; r <- nl(obj, W, fz0(ds0)); record("data_start_only", r$v, obj$counter$n, NA)
  obj$counter$n <- 0L; r2 <- lapply(list(default, ds0), function(s) nl(obj, W, fz0(s))); record("default_plus_data_start", min(vapply(r2, `[[`, 0, "v")), obj$counter$n, NA)
  rbindlist(rows)
}

res <- parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) { message("FAILED ", j$id, ": ", conditionMessage(e)); NULL }),
                          mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
res <- merge(res, prev_best, by = c("ds", "id"))
res[, best_all := pmin(best, min(obj[is.finite(obj)])), by = .(ds, id)]
res[, gap := obj - best_all]
saveRDS(res, sprintf("experiments/parametrization/_out/round2_%s.rds", mode))
log_msg("Saved round 2 | %d rows", nrow(res))
source("experiments/parametrization/rounds/summary2.R"); print_round2(res)
