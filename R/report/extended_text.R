# ==============================================================================
# R/report/extended_text.R -- text of the analyses that probe the identification structure from outside the main datasets (pharmacological blockade and
# muscle sympathetic nerve activity, the identifiability boundary, partial pooling of the time constants, perturbation-response descriptors).
# Sourced by manuscript/manuscript.qmd and manuscript/supplementary.qmd after R/report/manuscript_text.R (it uses its helpers: guarded(), num(), pct(),
# count_of(), join_and()). Every function returns markdown and degrades to "[PENDING: ...]" when its result is missing.
#
#   Blockade     abstract_blockade(), narrate_blockade(), blockade_discussion_sentence(), blockade_limitation_sentence(), blockade_table()
#   Boundary     narrate_boundary()
#   Pooling      narrate_pooling(), pooling_table()
#   Descriptors  narrate_descriptors(), descriptors_supplement()
# ==============================================================================

# ---- Helpers -------------------------------------------------------------------------------------------------------------------------
ci_text <- function(est, lo, hi, digits = 2) sprintf("%s (95%% CI %s–%s)", num(est, digits), num(lo, digits), num(hi, digits))
range_text <- function(x, digits = 2) { x <- x[is.finite(x)]; if (!length(x)) "[PENDING]" else if (diff(range(x)) < 0.5 * 10^(-digits)) num(x[1], digits) else paste0(num(min(x), digits), "–", num(max(x), digits)) }
blockade_row <- function(blockade, condition, mode = "spontaneous") { s <- blockade$summary; s[s$cond == condition & s$breathing == mode] }
spearman <- function(a, b) { k <- is.finite(a) & is.finite(b); if (sum(k) < 5) NA_real_ else cor(a[k], b[k], method = "spearman") }


# ---- Pharmacological blockade and muscle sympathetic nerve activity -------------------------------------------------------------------------
abstract_blockade <- function(blockade) guarded("results/blockade.rds", {
  if (!isTRUE(blockade$available)) return(NULL)
  at <- blockade_row(blockade, "atropine_high"); pr <- blockade_row(blockade, "propranolol")
  sprintf("in a public pharmacological-blockade study (%d volunteers), atropine reduced the variability of the vagal-like drive to %s of saline whereas propranolol alone did not (%s), but the sympathetic-like drive also fell under atropine and did not follow muscle sympathetic nerve activity",
          blockade$n_subjects, num(at$SD_p), num(pr$SD_p))
})

narrate_blockade <- function(blockade, figure = "@fig-blockade") guarded("results/blockade.rds", {
  if (!isTRUE(blockade$available)) return(PENDING("blockade data (tools/prep/blockade.R)"))
  at <- blockade_row(blockade, "atropine_high"); db <- blockade_row(blockade, "double"); pr <- blockade_row(blockade, "propranolol")
  atp <- blockade_row(blockade, "atropine_high", "paced_0.25")
  w <- dcast(blockade$change, subject + day + window ~ signal, value.var = "change"); bh <- w[window == "breath_hold_late"]; hv <- w[window == "hyperventilation"]
  lag0 <- blockade$lagged[lag == 0, .(r = median(r, na.rm = TRUE)), by = signal]
  paste0(
    sprintf("We fitted $M_3$, unchanged, to %d five-minute segments of %d volunteers who received saline, atropine (0.04 mg/kg) and propranolol (0.2 mg/kg). ", blockade$n_segments, blockade$n_subjects),
    sprintf("Relative to saline of the same session (geometric mean over volunteers), atropine reduced the standard deviation (SD) of the vagal-like drive to %s during spontaneous breathing (below 0.3 in %s volunteers) and to %s with paced breathing, while heart rate rose by a median of %s bpm; with propranolol added it fell to %s. ",
            ci_text(at$SD_p, at$SD_p_lo, at$SD_p_hi), count_of(at$n_SD_p_below_0.3, at$n), ci_text(atp$SD_p, atp$SD_p_lo, atp$SD_p_hi), num(at$d_hr, 0), ci_text(db$SD_p, db$SD_p_lo, db$SD_p_hi)),
    sprintf("Propranolol alone did not reduce it (ratio %s, 95%% CI %s\u2013%s; at least 0.8 of saline in %s volunteers). ", num(pr$SD_p), num(pr$SD_p_lo), num(pr$SD_p_hi), count_of(pr$n_SD_p_at_least_0.8, pr$n)),
    sprintf("The SD of the sympathetic-like drive, however, also fell under atropine, to %s (95%% CI %s\u2013%s; by more than a factor of three in %s volunteers), and to %s under the double blockade: once the vagal part of the spectrum is removed, the spectral anchor re-assigns the variance, so the slow drive is not independent of the fast one in this reconstruction. ",
            num(at$SD_s), num(at$SD_s_lo), num(at$SD_s_hi), count_of(at$n_SD_s_below_third, at$n), ci_text(db$SD_s, db$SD_s_lo, db$SD_s_hi)),
    sprintf("In the same volunteers the integrated muscle sympathetic nerve signal (MSNA) was higher in the second half of a breath-hold than at baseline in %s sessions, but the sympathetic-like state rose in only %s (median change %s baseline SDs), and across sessions the change of MSNA correlated with that of $x_s$ at Spearman %s ($x_p$ %s, heart rate %s); ",
            count_of(sum(bh$MSNA > 0), nrow(bh)), count_of(sum(bh$x_s > 0), nrow(bh)), num(median(bh$x_s), 1), num(spearman(bh$MSNA, bh$x_s)), num(spearman(bh$MSNA, bh$x_p)), num(spearman(bh$MSNA, bh$HR))),
    sprintf("during hyperventilation $x_s$ rose in %s sessions, $x_p$ fell in %s and heart rate rose in %s, and the median lag-0 correlation of $x_s$ with MSNA over the whole saline period was %s. ",
            count_of(sum(hv$x_s > 0), nrow(hv)), count_of(sum(hv$x_p < 0), nrow(hv)), count_of(sum(hv$HR > 0), nrow(hv)), num(lag0[signal == "x_s", r])),
    "The blockade therefore supports the identity of the fast drive, and the maneuvers show that the slow drive behaves as an acceleratory component of the heart rate, as it is defined, and not as a measure of sympathetic outflow (", figure, ").")
})

blockade_discussion_sentence <- function(blockade) guarded("results/blockade.rds", {
  if (!isTRUE(blockade$available)) return("")
  at <- blockade_row(blockade, "atropine_high"); pr <- blockade_row(blockade, "propranolol")
  sprintf("The pharmacological data add an external test that the main datasets lack: the vagal-like drive lost most of its variability under atropine (%s of saline) and not under propranolol alone (%s), as a vagal component should, whereas the sympathetic-like drive lost variability under atropine too and did not follow muscle sympathetic nerve activity during breath-holds; the label of the fast drive is therefore supported by blockade, and that of the slow drive remains a statement about its time scale.",
          num(at$SD_p), num(pr$SD_p))
})

blockade_anchor_sentence <- function(blockade) guarded("results/blockade.rds", {
  if (!isTRUE(blockade$available)) return("The pharmacological evidence on which the spectral weights are based was drawn from the literature and was not generated here.")
  at <- blockade_row(blockade, "atropine_high")
  sprintf("The pharmacological evidence on which the spectral weights are based was drawn from the literature. In the one blockade dataset analyzed here, atropine left %s of the RR variance of saline (geometric mean), consistent with reports that at rest vagal blockade removes most of the power at low as well as high frequencies [@clemson2022beyond]; the weights instead assign most low-frequency variance to the sympathetic branch, which helps explain why the sympathetic-like amplitude fell under vagal blockade and makes the allocation of low-frequency variance the weakest part of the anchor.",
          pct(at$rr_var, 1))
})

blockade_limitation_sentence <- function(blockade) guarded("results/blockade.rds", {
  if (!isTRUE(blockade$available)) return("The latent states have not been tested under pharmacological blockade or compared with sympathetic nerve activity.")
  sprintf("The comparison with pharmacological blockade and with muscle sympathetic nerve activity rests on %d young volunteers studied supine and at rest in a single public dataset; blockade is not additive and may be incomplete, and it supported only the fast drive: the slow drive lost variability under vagal blockade and did not follow the nerve signal, so it is a slow component of the heart rate and not a measurement of sympathetic outflow.",
          blockade$n_subjects)
})

#' Supplementary table: ratios to saline by breathing mode and condition.
blockade_table <- function(blockade) {
  if (is.null(blockade) || !isTRUE(blockade$available)) return(NULL)
  s <- blockade$summary[order(match(breathing, c("spontaneous", "paced_0.25")), match(cond, c("atropine_high", "double", "propranolol")))]
  data.table(Breathing = ifelse(s$breathing == "spontaneous", "spontaneous", "paced 0.25 Hz"),
             Condition = c(atropine_high = "atropine", double = "atropine + propranolol", propranolol = "propranolol")[s$cond], n = s$n,
             `SD of $x_p$ (ratio)` = sprintf("%s [%s, %s]", num(s$SD_p), num(s$SD_p_lo), num(s$SD_p_hi)),
             `SD of $x_s$ (ratio)` = sprintf("%s [%s, %s]", num(s$SD_s), num(s$SD_s_lo), num(s$SD_s_hi)),
             `Heart rate change (bpm)` = num(s$d_hr, 0))
}


# ---- Identifiability boundary -----------------------------------------------------------------------------------------------------------------
narrate_boundary <- function(ib) guarded("results/ident_boundary.rds", {
  c3 <- ib$by_cell[axis == "time-scale ratio" & spec == "M3"]; v3 <- ib$by_cell[axis == "variance ratio" & spec == "M3"]; g <- ib$by_regime
  col <- function(reg, sp) g[regime == reg & spec == sp, collapse]
  paste0(
    sprintf("Simulations from the model itself locate where the decomposition can be recovered (%d recordings per condition; Supplementary Figure S11): with $M_3$, the correlation between the filtered sympathetic-like state and the truth was %s when $\\tau_s/\\tau_p$ was 4 or less, %s at 10 and %s at 30–100 with recordings of 10–30 min. ",
            ib$settings$n_per_cell, range_text(c3[ratio <= 4, cor_x_s]), range_text(c3[ratio == 10, cor_x_s]), range_text(c3[ratio >= 30 & duration >= 600, cor_x_s])),
    sprintf("Without physiological information ($M_0$), one drive absorbed all the variance in %s of the fits when the ratio was below 10 and in %s when it was 10 or more (at least 10 min); with $M_3$, in %s and %s, so the priors prevent the collapse but cannot create identifiability where the time scales overlap. ",
            pct(col("ratio < 10", "M0")), pct(col("ratio >= 10, >= 10 min", "M0")), pct(col("ratio < 10", "M3")), pct(col("ratio >= 10, >= 10 min", "M3"))),
    sprintf("The weaker drive is the one that is lost: with $\\lambda_P$ = %s the correlation of the sympathetic-like state was %s, and with $\\lambda_P$ = %s that of the vagal-like state was %s.",
            num(max(v3$lambda), 0), range_text(v3[lambda == max(lambda), cor_x_s]), num(min(v3$lambda), 1), range_text(v3[lambda == min(lambda), cor_x_p])))
})


# ---- Partial pooling of the time constants ---------------------------------------------------------------------------------------------------------------
narrate_pooling <- function(pooling) guarded("results/pooling.rds", {
  r <- pooling$reliability; v <- function(st, var, col = "icc") r[setting == st & variable == var][[col]]
  sprintf("When the priors of the time constants were replaced by weak ones or by empirical-Bayes priors estimated on the other half of the subjects, the test-retest ICC of the amplitudes stayed at %s, that of $\\tau_p$ was %s, %s and %s and that of $\\tau_s$ %s, %s and %s (production, weak, empirical Bayes; supplementary material), so the slow time constant carries recording-specific information that the fixed prior suppresses but remains the least reliable summary.",
          range_text(r[variable %in% c("SD_p", "SD_s"), icc]), num(v("production", "tau_p")), num(v("weak", "tau_p")), num(v("empirical Bayes", "tau_p")),
          num(v("production", "tau_s")), num(v("weak", "tau_s")), num(v("empirical Bayes", "tau_s")))
})

pooling_table <- function(pooling) {
  if (is.null(pooling)) return(NULL)
  r <- pooling$reliability[order(match(variable, c("SD_p", "SD_s", "tau_p", "tau_s")), match(setting, c("production", "weak", "empirical Bayes")))]
  lab <- c(SD_p = "SD of the vagal-like drive", SD_s = "SD of the sympathetic-like drive", tau_p = "$\\tau_p$", tau_s = "$\\tau_s$")
  data.table(Quantity = lab[r$variable], Prior = r$setting, `ICC(3,1)` = num(r$icc), Spearman = num(r$spearman), `AUC young vs older` = num(r$auc_young_older),
             `SD of log values` = num(r$sd_log), `Inside production prior` = ifelse(is.finite(r$in_production_prior), pct(r$in_production_prior), ""))
}


# ---- Perturbation-response descriptors -------------------------------------------------------------------------------------------------------------------
narrate_descriptors <- function(sd) guarded("results/state_descriptors.rds", {
  ex <- sd$exercise; m <- merge(ex, sd$recovery_main, by.x = "id", by.y = "id", all.x = TRUE)
  q <- quantile(ex$Trec, c(0.5, 0.1, 0.9), na.rm = TRUE)
  sprintf("Descriptors of the perturbation response computed from the states with their uncertainty, as proposed for latent autonomic trajectories [@castillo2026state], gave a median time to a sustained return to the resting region of %s s (10th–90th percentile %s–%s) after exercise, which correlated with the recovery times above at Spearman %s ($x_p$) and %s ($x_s$); the probability of return depended mainly on the horizon (%s at 120, 300 and 600 s), and in simulations the recovery time was recovered in rank but not exactly (Supplementary Figure S12).",
          num(q[[1]], 0), num(q[[2]], 0), num(q[[3]], 0), num(spearman(m$Trec, m$p)), num(spearman(m$Trec, m$s)),
          join_and(num(vapply(sd$settings$horizons, function(H) mean(ex[[paste0("Preturn_", H)]], na.rm = TRUE), numeric(1)))))
})

descriptors_supplement <- function(sd) guarded("results/state_descriptors.rds", {
  s <- sd$sim; ex <- sd$exercise; tl <- sd$tilt
  by <- function(sc, f) f(s[scenario == sc])
  rho <- function(d) spearman(d$Trec_true, d$Trec); dmax_ratio <- function(d) median(d$Dmax / d$Dmax_true, na.rm = TRUE)
  cover <- function(d) { k <- is.finite(d$Trec_true) & is.finite(d$Trec_lo); mean(d$Trec_true[k] >= d$Trec_lo[k] & d$Trec_true[k] <= d$Trec_hi[k]) }
  bias <- function(d) median(d$Trec - d$Trec_true, na.rm = TRUE)
  r2 <- { d <- ex[is.finite(Trec)]; summary(lm(rank(Trec) ~ rank(rest_hr) + rank(rest_sdnn_ms) + rank(rest_rmssd_ms), d))$r.squared }
  paste0(
    sprintf("*Simulations.* In %d simulated recordings per scenario with known states (a respiratory oscillator that the model does not represent was added), the stationary variance of the full-recording fit was %s (vagal-like) and %s (sympathetic-like) times the true one in exercise-like recordings and %s and %s in tilt-like ones: the fit absorbs the perturbation. Refitting $M_3$ on the resting beats alone gave %s and %s (exercise), and this rest-only scale is the one used. ",
            nrow(s[scenario == "exercise"]), num(by("exercise", function(d) median(d$V_ratio_full_p)), 1), num(by("exercise", function(d) median(d$V_ratio_full_s)), 1),
            num(by("tilt", function(d) median(d$V_ratio_full_p)), 1), num(by("tilt", function(d) median(d$V_ratio_full_s)), 1),
            num(by("exercise", function(d) median(d$V_ratio_rest_p)), 1), num(by("exercise", function(d) median(d$V_ratio_rest_s)), 1)),
    sprintf("With it, the maximum displacement was %s (exercise) and %s (tilt) of the true value (median ratio), the recovery time correlated with the truth at Spearman %s and %s with a median bias of %s s and %s s, and its 90%% interval contained the true value in %s and %s of the recordings. ",
            num(by("exercise", dmax_ratio)), num(by("tilt", dmax_ratio)), num(by("exercise", rho)), num(by("tilt", rho)), num(by("exercise", bias), 0), num(by("tilt", bias), 0),
            pct(by("exercise", cover)), pct(by("tilt", cover))),
    sprintf("A regime transition at the horizon was flagged in %s of the simulated persistent shifts, %s of the exercise recoveries and %s of the tilt recoveries. ",
            pct(by("shift", function(d) mean(d$transition, na.rm = TRUE))), pct(by("exercise", function(d) mean(d$transition, na.rm = TRUE))), pct(by("tilt", function(d) mean(d$transition, na.rm = TRUE)))),
    sprintf("*Recordings.* In the %d exercise recordings the recovery time was nearly unrelated to resting heart rate, SDNN and RMSSD ($R^2$ of the ranks %s), and a transition to a regime other than rest was present at 300 s in %s; after the %d PRCP maneuvers, transitions were flagged in %s, close to the simulated false-alarm rate. These descriptors are conditional on the model, and the tolerance, dwell time and horizon are part of their definition.",
            nrow(ex), num(r2), pct(mean(ex$transition, na.rm = TRUE)), nrow(tl), pct(mean(tl$transition, na.rm = TRUE))))
})
