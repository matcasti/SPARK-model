# ==============================================================================
# experiments/parametrization/rounds/run6.R -- EXPERIMENTAL, round 6: the racing scheme with R's BFGS (what the web app also implements),
# in the prior-whitened coordinates. 313 recordings + 626 head/tail-70% series. Run from the project root.
# ==============================================================================
source("analysis/00_setup.R"); source("experiments/parametrization/strategies.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
best_orig <- unique(readRDS("experiments/parametrization/_out/round4_full.rds")[, .(ds, id, best_all)])
best_ht <- unique(readRDS("experiments/parametrization/_out/round5_full.rds")[, .(ds, id, best_all)])
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg)
  for (id in names(d$rr)) { rr <- d$rr[[id]]; n <- length(rr); k <- floor(0.7 * n)
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = rr)
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = paste0(id, "_head"), rr = rr[seq_len(k)])
    jobs[[length(jobs) + 1]] <- list(ds = ds, id = paste0(id, "_tail"), rr = rr[(n - k + 1):n]) } }
bf <- function(obj, W, z0, maxit) {
  fz <- function(z) { p <- W$to(z); obj$f(p[1], p[2], p[3]) }
  o <- tryCatch(optim(z0, fz, method = "BFGS", control = list(maxit = maxit)), error = function(e) NULL)
  if (is.null(o) || !is.finite(o$value)) list(v = Inf, z = z0) else list(v = o$value, z = o$par)
}
one <- function(job) {
  series <- prepare_series(job$rr); pri <- series$priors; ctl <- series$control
  obj <- make_objective(series); P <- parametrizations(pri); W <- P$whitened
  default <- c(pri$ks_mode, pri$ks_mode + max(pri$kp_mode - pri$ks_mode, ctl$kp_gap_floor), ctl$lR_start)
  s_lo <- c(1 / 160, 1 / 12, 1e-6); s_vlo <- c(1 / 170, 1 / 14, 1e-7); s_hi <- c(1 / 150, 1 / 25, 0.3)
  zf <- function(s) W$from(s[1], s[2], s[3])
  rows <- list(); add <- function(m, v, e) rows[[length(rows) + 1]] <<- data.table(ds = job$ds, id = job$id, method = m, obj = v, evals = e)
  r <- run_local(obj, P$base, default, "bfgs"); add("production", if (r$ok) r$v else Inf, r$evals)
  obj$counter$n <- 0L; r <- bf(obj, W, zf(default), 1000); add("whitened_bfgs_1", r$v, obj$counter$n)
  for (set in list(list("LR3", list(default, s_lo, s_hi)), list("LR3b", list(default, s_vlo, s_hi)))) {
    S <- set[[2]]
    obj$counter$n <- 0L; v <- min(vapply(S, function(s) bf(obj, W, zf(s), 1000)$v, 0)); add(paste0("full3_", set[[1]]), v, obj$counter$n)
    for (m in c(4, 6, 8, 10, 15, 20)) {
      obj$counter$n <- 0L; sh <- lapply(S, function(s) bf(obj, W, zf(s), m)); b <- which.min(vapply(sh, `[[`, 0, "v")); f <- bf(obj, W, sh[[b]]$z, 1000)
      add(sprintf("race_%s_bfgs_m%d", set[[1]], m), f$v, obj$counter$n)
    }
  }
  rbindlist(rows)
}
res <- parallel::mclapply(jobs, function(j) tryCatch(one(j), error = function(e) { message("FAILED ", j$id, ": ", conditionMessage(e)); NULL }),
                          mc.cores = n_cores, mc.preschedule = FALSE)
res <- rbindlist(res[!vapply(res, is.null, TRUE)])
res <- merge(res, rbind(best_orig, best_ht), by = c("ds", "id"), all.x = TRUE)
res[, best_all := pmin(best_all, min(obj[is.finite(obj)])), by = .(ds, id)]
res[, set := ifelse(grepl("_(head|tail)$", id), "out_of_sample", "original")]
res[, hit := is.finite(obj) & obj - best_all < 0.01]
saveRDS(res, "experiments/parametrization/_out/round6_full.rds")
out <- res[, .(success_pct = round(100 * mean(hit), 1), evals_mean = round(mean(evals)), evals_p95 = round(quantile(evals, .95)), n_miss = sum(!hit),
               worst_gap = signif(max((obj - best_all)[is.finite(obj)]), 2)), by = .(set, method)]
for (s in c("original", "out_of_sample")) { cat("\n==", s, "==\n"); print(out[set == s][order(-success_pct, evals_mean)][, !"set"], row.names = FALSE) }
