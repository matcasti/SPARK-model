# Summary of a parametrization experiment: success = within `tol` of the best known optimum of the recording.
print_summary <- function(res, tol = 0.01) {
  r <- copy(res[strategy != "ref"])
  r[, success := ok & gap < tol]
  r[, converged := ok & conv == 0]
  local <- r[start > 0 | grepl("grid|profile", strategy) == FALSE]
  tab <- function(d) d[, .(runs = .N, success = round(100 * mean(success), 1), converged = round(100 * mean(converged), 1),
                           fail = round(100 * mean(!ok), 1), evals_med = median(evals), secs_med = round(median(secs), 3),
                           gap_max = signif(max(gap, na.rm = TRUE), 2)), by = strategy]
  cat("\n== From the PRODUCTION start (start 0) / global strategies ==\n"); print(tab(r[start == 0])[order(-success, evals_med)], row.names = FALSE, digits = 3)
  cat("\n== From 6 DISPERSED starts (robustness) ==\n"); print(tab(r[start > 0])[order(-success, evals_med)], row.names = FALSE)
  # all dispersed starts of a strategy reach the same optimum?
  agree <- r[start > 0, .(all_agree = all(success)), by = .(ds, id, strategy)][, .(recordings_all_6_succeed_pct = round(100 * mean(all_agree), 1)), by = strategy]
  cat("\n== Recordings where ALL 6 dispersed starts reach the best optimum (%) ==\n"); print(agree[order(-recordings_all_6_succeed_pct)], row.names = FALSE)
  ref <- res[strategy == "ref", .(ref_hits = round(100 * mean(gap < tol, na.rm = TRUE), 1)), by = .(ds)]
  cat("\n== Reference runs (24 starts, BFGS+NM) reaching the best optimum (%), by dataset ==\n"); print(ref, row.names = FALSE)
}
