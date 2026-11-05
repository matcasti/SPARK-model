# Summary of round 2: success (within tol of the best known optimum), evaluations per recording, misses.
print_round2 <- function(res, tol = 0.01) {
  r <- copy(res); r[, hit := is.finite(gap) & gap < tol]
  fulls <- r[grepl("^full_start", method)]
  fulls[, j := as.integer(sub("full_start", "", method))]
  pool <- function(label, js) { d <- fulls[j %in% js, .(hit = any(hit), evals = sum(evals)), by = .(ds, id)]; d[, method := label]; d }
  # adaptive stopping: run default, then dispersed starts in order; stop when two runs agree on the best (within tol) or all used
  adaptive <- function(label, max_runs, agree = 2L) {
    out <- fulls[order(ds, id, j), {
      v <- obj; e <- evals; stop_at <- length(v)
      for (k in 2:min(length(v), max_runs)) { b <- min(v[1:k]); if (sum(v[1:k] - b < tol) >= agree) { stop_at <- k; break }; stop_at <- k }
      list(hit = min(v[1:stop_at]) - min(best_all[1]) < tol, evals = sum(e[1:stop_at]), runs = stop_at)
    }, by = .(ds, id)]
    out[, method := label]; out }
  ad <- list(adaptive("adaptive_agree2_max4", 4L), adaptive("adaptive_agree2_max5", 5L), adaptive("adaptive_agree2_max7", 7L))
  tab <- rbindlist(c(list(pool("fixed_1run (default)", 0), pool("fixed_2runs", 0:1), pool("fixed_3runs", 0:2), pool("fixed_4runs  [baseline]", 0:3)),
                     lapply(ad, function(d) d[, .(ds, id, hit, evals, method)]),
                     list(r[!grepl("^full_start", method), .(ds, id, hit, evals, method)])), fill = TRUE)
  out <- tab[, .(success_pct = round(100 * mean(hit), 1), evals_mean = round(mean(evals)), evals_median = median(evals), evals_p95 = quantile(evals, .95),
                 misses = paste(id[!hit], collapse = ",")), by = method]
  base_e <- out[method == "fixed_4runs  [baseline]", evals_mean]
  out[, cost_vs_baseline := round(evals_mean / base_e, 2)]
  print(out[order(-success_pct, evals_mean), .(method, success_pct, evals_mean, evals_p95, cost_vs_baseline, misses)], row.names = FALSE, digits = 3)
  if ("runs" %in% names(tab)) { cat("\nAdaptive rule: share of recordings by number of runs used\n"); print(tab[grepl("adaptive_agree2_max5", method), .N, by = runs][order(runs)]) }
}
