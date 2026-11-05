# ==============================================================================
# experiments/parametrization/rounds/run5.R -- EXPERIMENTAL, round 5: out-of-sample check of the lR-stratified racing design.
# New series = first 70% and last 70% of every recording (different likelihood surfaces; the designed starts never saw them).
# Reference "best known" = best of all methods plus 24 dispersed whitened|nlminb starts. Run from the project root.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/parametrization/strategies.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg)
  for (id in names(d$rr)) { rr <- d$rr[[id]]; n <- length(rr)
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = paste0(id, "_head"), rr = rr[seq_len(floor(0.7 * n))])
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = paste0(id, "_tail"), rr = rr[(n - floor(0.7 * n) + 1):n]) } }
nl <- function(obj, par, z0, iter.max = 1000) {
  fz <- function(z) { p <- par$to(z); obj$f(p[1], p[2], p[3]) }
  o <- tryCatch(nlminb(z0, fz, control = list(iter.max = iter.max, eval.max = 2 * iter.max)), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$objective)) list(v = Inf, z = z0) else list(v = o$objective, z = o$par)
}
one <- function(job, idx) {
  series <- prepare_series(job$rr); pri <- series$priors; ctl <- series$control
  obj <- make_objective(series); P <- parametrizations(pri); W <- P$whitened
  set.seed(7000L + idx)
  default <- c(pri$ks_mode, pri$ks_mode + max(pri$kp_mode - pri$ks_mode, ctl$kp_gap_floor), ctl$lR_start)
  disp <- dispersed_starts(pri, 24)
  s_lo <- c(1 / 160, 1 / 12, 1e-6); s_vlo <- c(1 / 170, 1 / 14, 1e-7); s_hi <- c(1 / 150, 1 / 25, 0.3)
  rows <- list(); add <- function(m, v, e) rows[[length(rows) + 1]] <<- data.table(ds = job$ds, id = job$id, method = m, obj = v, evals = e)
  zf <- function(s) W$from(s[1], s[2], s[3])
  r <- run_local(obj, P$base, default, "bfgs");              add("production", if (r$ok) r$v else Inf, r$evals)
  obj$counter$n <- 0L; r <- nl(obj, W, zf(default));           add("whitened_nlminb_1", r$v, obj$counter$n)
  obj$counter$n <- 0L; v <- min(vapply(c(list(default), disp[1:3]), function(s) nl(obj, W, zf(s))$v, 0)); add("fixed_4runs", v, obj$counter$n)
  race <- function(S, m) { obj$counter$n <- 0L; sh <- lapply(S, function(s) nl(obj, W, zf(s), iter.max = m)); b <- which.min(vapply(sh, `[[`, 0, "v")); f <- nl(obj, W, sh[[b]]$z); list(v = f$v, e = obj$counter$n) }
  for (cfgr in list(list("race_LR3_m6", list(default, s_lo, s_hi), 6), list("race_LR3b_m5", list(default, s_vlo, s_hi), 5), list("race_LR3b_m6", list(default, s_vlo, s_hi), 6),
                    list("race_D4_m8", c(list(default), disp[1:3]), 8))) { x <- race(cfgr[[2]], cfgr[[3]]); add(cfgr[[1]], x$v, x$e) }
  v24 <- vapply(disp, function(s) nl(obj, W, zf(s))$v, 0); add("ref24", min(v24), NA_integer_)
  rbindlist(rows)
}
res <- parallel::mclapply(seq_along(jobs), function(i) tryCatch(one(jobs[[i]], i), error = function(e) { message("FAILED ", jobs[[i]]$id, ": ", conditionMessage(e)); NULL }),
                          mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
res[, best_all := min(obj[is.finite(obj)]), by = .(ds, id)]
res[, hit := is.finite(obj) & obj - best_all < 0.01]
saveRDS(res, "experiments/parametrization/_out/round5_full.rds")
out <- res[method != "ref24", .(success_pct = round(100 * mean(hit), 1), evals_mean = round(mean(evals)), evals_p95 = round(quantile(evals, .95)), n_miss = sum(!hit),
                                 worst_gap = signif(max((obj - best_all)[is.finite(obj)]), 2)), by = method]
cat("\nOut-of-sample series:", uniqueN(res$id), "\n"); print(out[order(-success_pct, evals_mean)], row.names = FALSE)
