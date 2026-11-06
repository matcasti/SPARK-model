# EXPERIMENTAL. Idea 7 (automatic differentiation): is the numerical Hessian good enough? Standard errors of ln tau_p and ln tau_s from the Hessian of the M3 objective
# computed with three finite-difference settings (Richardson default, a finer and a coarser step) on 60 recordings, plus the cost of one fit. Run from the project root.
source("analysis/00_setup.R")
set.seed(9); jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); ids <- sample(names(d$rr), min(length(d$rr), c(prcp = 5, fantasia = 15, exercise = 40)[[ds]])); for (id in ids) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]]) }
one <- function(j) {
  series <- prepare_series(j$rr); spec <- model_spec("M3"); t0 <- proc.time()[["elapsed"]]; fit <- fit_variant(series, spec, full = FALSE); secs <- proc.time()[["elapsed"]] - t0
  f <- function(th) variant_objective(th, series$dy, spec, series$priors, series$control)
  se <- function(args) { H <- tryCatch(numDeriv::hessian(f, fit$theta, method.args = args), error = function(e) NULL); if (is.null(H)) return(c(NA, NA)); C <- tryCatch(MASS::ginv(H), error = function(e) NULL); if (is.null(C)) return(c(NA, NA))
    ks <- exp(fit$theta[1]); kp <- ks + exp(fit$theta[2]); J <- rbind(c(1, 0, 0), c(ks / kp, exp(fit$theta[2]) / kp, 0)); sqrt(pmax(diag(J %*% C[1:3, 1:3] %*% t(J)), 0)) }   # SD of ln ks and ln kp
  a <- se(list()); b <- se(list(eps = 1e-5, d = 1e-3, r = 6)); c <- se(list(eps = 1e-3, d = 0.05, r = 3))
  data.table(ds = j$ds, id = j$id, secs = secs, rel_fine = abs(b / a - 1), rel_coarse = abs(c / a - 1), which = c("ln_ks", "ln_kp"))
}
res <- rbindlist(drop_failed(map_recordings(jobs, one, max(1L, min(5L, cfg$cores - 2L))), "recording"))
cat("\nTime of one M3 fit (3 starts, compiled filter): median", round(median(res[which == "ln_ks", secs]), 3), "s; max", round(max(res$secs), 2), "s\n")
cat("Relative difference of the standard error between finite-difference settings (median, 95th percentile, max):\n")
print(res[, .(fine_med = signif(median(rel_fine, na.rm = TRUE), 2), fine_q95 = signif(quantile(rel_fine, .95, na.rm = TRUE), 2), fine_max = signif(max(rel_fine, na.rm = TRUE), 2), coarse_med = signif(median(rel_coarse, na.rm = TRUE), 2), coarse_q95 = signif(quantile(rel_coarse, .95, na.rm = TRUE), 2), coarse_max = signif(max(rel_coarse, na.rm = TRUE), 2)), by = which], row.names = FALSE)
