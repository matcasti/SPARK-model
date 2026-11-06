# EXPERIMENTAL. Does the recording's own spectral anchor fit its own data better than a constant or than the anchor of ANOTHER recording?
# Same number of parameters in every variant (lambda_P is fixed), so the profiled log-likelihood is directly comparable.
# Run from the project root:  Rscript experiments/anchors/real.R
source("analysis/00_setup.R"); source("experiments/anchors/lib.R")
n_cores <- max(1L, min(5L, cfg$cores - 2L))
jobs <- list()
for (ds in cfg$datasets) { d <- load_dataset(ds, cfg); for (id in names(d$rr)) jobs[[length(jobs) + 1]] <- list(ds = ds, id = id, rr = d$rr[[id]], group = if (ds == "fantasia") substr(id, 3, 3) else NA_character_) }
anchors <- unlist(drop_failed(map_recordings(jobs, function(j) extract_spectral_priors(j$rr, spark_priors(), verbose = FALSE)$energy_ratio, n_cores), "anchor"))
for (i in seq_along(jobs)) jobs[[i]]$anchor <- anchors[[i]]
set.seed(11)
for (ds in cfg$datasets) { idx <- which(vapply(jobs, `[[`, "", "ds") == ds); other <- sample(idx)
  for (k in seq_along(idx)) { j <- idx[k]; o <- other[k]; if (o == j) o <- other[(k %% length(idx)) + 1]; jobs[[j]]$perm <- jobs[[o]]$anchor; jobs[[j]]$const <- median(vapply(jobs[idx], `[[`, 0, "anchor")) } }
log_msg("Own vs constant vs permuted anchor | %d recordings", length(jobs))
fit_value <- function(series, value) { r <- fit_rule(series, function(kp, ks) value); data.table(LL = r$LL_star, tau_p = r$tau_p, tau_s = r$tau_s, ks_ok = r$diag$ks_ok, pacf_viol = r$diag$pacf_viol) }
res <- rbindlist(drop_failed(map_recordings(jobs, function(j) { series <- prepare_series(j$rr)
  rbindlist(list(own = fit_value(series, j$anchor), const = fit_value(series, j$const), permuted = fit_value(series, j$perm)), idcol = "variant")[, `:=`(ds = j$ds, id = j$id, group = j$group, anchor = j$anchor)] }, n_cores), "recording"))
dir.create("experiments/anchors/_out", showWarnings = FALSE); saveRDS(res, "experiments/anchors/_out/real.rds")
w <- dcast(res, ds + id + group + anchor ~ variant, value.var = "LL"); w[, `:=`(d_const = own - const, d_perm = own - permuted)]
cat("\n== Log-likelihood of the recording's OWN anchor minus that of the alternative (positive: own anchor fits better) ==\n")
print(w[, .(n = .N, vs_constant_median = round(median(d_const), 2), pct_own_better = round(100 * mean(d_const > 0), 1), pct_own_better_by_2 = round(100 * mean(d_const > 2), 1), pct_const_better_by_2 = round(100 * mean(d_const < -2), 1)), by = ds], row.names = FALSE)
print(w[, .(n = .N, vs_permuted_median = round(median(d_perm), 2), pct_own_better = round(100 * mean(d_perm > 0), 1), pct_own_better_by_2 = round(100 * mean(d_perm > 2), 1), pct_perm_better_by_2 = round(100 * mean(d_perm < -2), 1)), by = ds], row.names = FALSE)
cat("\n== Diagnostics by variant ==\n"); print(res[, .(ks_pass = round(100 * mean(ks_ok), 1), pacf_viol = round(mean(pacf_viol), 1), tau_p_med = round(median(tau_p), 1), pct_tau_p_in_prior = round(100 * mean(tau_p >= 3 & tau_p <= 6), 1)), by = variant], row.names = FALSE)
cat("\nAnchor values (median, range):\n"); print(w[, .(med = round(median(anchor), 2), lo = round(min(anchor), 2), hi = round(max(anchor), 2)), by = ds])
