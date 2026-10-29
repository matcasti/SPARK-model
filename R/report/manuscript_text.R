# ==============================================================================
# R/report/manuscript_text.R -- text of the manuscript generated from the computed results.
# Sourced by manuscript/manuscript.qmd and manuscript/supplementary.qmd; every function returns a
# markdown string and degrades to "[PENDING: ...]" when its results are missing.
#
# The text is DESCRIPTIVE: it reports numbers and states, with explicit and visible rules, whether
# they agree with an expected pattern or differ between conditions; the wording changes with the
# outcome (it does not assume a result). It makes no causal claims. The thresholds are in
# NARRATIVE below. Read the generated sentences and edit the functions if you want other wording.
#
#   helpers          txt(), pct(), median_range_text(), compare_phrase(), join_and(), ...
#   Results          narrate_identifiability(), narrate_constraints(), narrate_prcp(),
#                    narrate_fantasia(), narrate_exercise(), landscape_sentence()
#   Abstract         abstract_identifiability() abstract_prcp() abstract_aging() abstract_exercise()
#   Discussion       findings_summary() (number-free; the numbers are in the Abstract and the Results), sensitivity_sentence(), ...
# ==============================================================================

# Rules used to choose wording (shown in the sentences where they matter)
NARRATIVE <- list(
  agree_share = 0.75,       # share of recordings above which a pattern is called consistent
  disagree_share = 0.25,    # below which it is called not consistent
  similar_rel = 0.20,       # two values within this relative difference are called similar
  effect_size = 0.33,       # |Cliff's delta| at or above which a group difference is called medium or larger
  sync_s = 5                # onsets within this many seconds are called simultaneous
)


# ---- Helpers ---------------------------------------------------------------------------------------------------
PENDING <- function(what) sprintf("[PENDING: %s]", what)

#' Evaluate `e`; "[PENDING]" when it fails, is empty or is not finite.
txt <- function(e) {
  tryCatch({
    x <- as.character(e)
    if (!length(x) || all(is.na(x)) || any(grepl("(^|[^A-Za-z])-?(Inf|NaN|NA)([^A-Za-z]|$)", x))) "[PENDING]" else x
  }, error = function(err) "[PENDING]")
}

#' Run a narrative builder; a missing input or any error becomes "[PENDING: <what>]".
guarded <- function(what, expr) {
  tryCatch({
    out <- expr
    if (is.null(out) || !length(out) || is.na(out[1])) PENDING(what) else out
  }, error = function(e) {
    if (isTRUE(getOption("spark.narrative.debug"))) stop(e)      # options(spark.narrative.debug = TRUE) shows the error
    PENDING(what)
  })
}

# Vectorized: a non-finite value becomes "[PENDING]"
pct <- function(x, digits = 0) if (!length(x)) "[PENDING]" else ifelse(is.finite(x), sprintf(paste0("%.", digits, "f%%"), 100 * x), "[PENDING]")
num <- function(x, digits = 2) {
  if (!length(x)) return("[PENDING]")
  x <- ifelse(is.finite(x) & abs(x) < 0.5 * 10^(-digits), 0, x)        # no "-0"
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), x), "[PENDING]")
}
count_of <- function(k, n) sprintf("%d of %d", as.integer(k), as.integer(n))

#' "median [2.5th, 97.5th percentile]" of a numeric vector.
median_range_text <- function(x, digits = 2, unit = "") {
  x <- x[is.finite(x)]
  if (!length(x)) return("[PENDING]")
  q <- quantile(x, c(0.5, 0.025, 0.975), names = FALSE)
  sprintf(paste0("%.", digits, "f%s [%.", digits, "f, %.", digits, "f]"), q[1], unit, q[2], q[3])
}

#' "a, b and c"
join_and <- function(items) {
  items <- items[nzchar(items)]
  if (length(items) <= 1) return(paste(items, collapse = ""))
  paste0(paste(head(items, -1), collapse = ", "), " and ", tail(items, 1))
}

#' "lower than", "similar to" or "higher than": a relative to b (within NARRATIVE$similar_rel = similar).
compare_phrase <- function(a, b, rel = NARRATIVE$similar_rel) {
  if (!is.finite(a) || !is.finite(b)) return("[PENDING]")
  scale <- max(abs(a), abs(b))
  if (scale == 0 || abs(a - b) <= rel * scale) "similar to" else if (a < b) "lower than" else "higher than"
}

#' Fixed-notation numbers for tables (never scientific notation): `digits` significant digits; p-values below 0.001 as "<0.001".
minus_sign <- function(x) gsub("-", "\u2212", x, fixed = TRUE)        # a true minus sign does not break across lines in tables
fmt_fixed <- function(x, digits = 3) minus_sign(vapply(x, function(v) if (!is.finite(v)) "" else if (v == 0) "0" else sprintf(paste0("%.", max(0, digits - 1 - floor(log10(abs(signif(v, digits))))), "f"), v), character(1)))
fmt_p <- function(p) ifelse(is.finite(p), ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)), "")


#' Display names of the summaries used in the group contrasts.
VARIABLE_LABELS <- c(
  Nu_0 = "$\\nu_0$", Tau_P = "$\\tau_p$", Tau_S = "$\\tau_s$", lambda_P = "$\\lambda_P$", Sigma_sys2 = "$\\hat\\sigma^2_{sys}$",
  rr_jitter_ms = "the timing jitter", sd_x_p = "the SD of $x_p$", sd_x_s = "the SD of $x_s$",
  sd_ratio = "the SD ratio $x_p/x_s$", log_kappa_ratio = "$\\ln(\\kappa_p/\\kappa_s)$", share_p = "the share of $x_p$ variance",
  mean_hr = "mean heart rate", SDNN = "SDNN", RMSSD = "RMSSD", LF_HF = "LF/HF", f_resp = "the respiratory frequency",
  resp_sd_x_p = "the respiratory-band SD of $x_p$", alpha_LF = "the baroreflex alpha index (LF)", sbp = "systolic pressure",
  map = "MAP", pulse_pressure = "pulse pressure"
)
label_of <- function(v) ifelse(v %in% names(VARIABLE_LABELS), VARIABLE_LABELS[v], v)


#' Readable name of a sensitivity setting ("tau_p_bounds=narrow") and of a prior term ("tau").
SETTING_LABELS <- c(tau_p_bounds = "the $\\tau_p$ prior interval", tau_s_bounds = "the $\\tau_s$ prior interval",
                    topo_target_ratio = "the separation target", topo_prior_weight = "the separation weight",
                    jitter_bounds_ms = "the timing-jitter interval", crossover_hz = "the spectral crossover",
                    sns_leakage_at_boundary = "the sympathetic leakage", lower_lf_vagal_fraction = "the vagal floor")
TERM_LABELS <- c(tau = "the time-scale priors", topo = "the separation prior", jitter = "the timing-jitter prior", anchor = "the spectral anchor")
setting_label <- function(setting) {
  parts <- strsplit(setting, "=", fixed = TRUE)[[1]]
  base <- if (parts[1] %in% names(SETTING_LABELS)) SETTING_LABELS[[parts[1]]] else parts[1]
  if (length(parts) > 1) paste0(base, " (", parts[2], ")") else base
}
term_label <- function(setting) {
  term <- sub("M3 without ", "", setting)
  if (term %in% names(TERM_LABELS)) TERM_LABELS[[term]] else term
}

spec_median <- function(metrics, spec, column) {
  stats::median(metrics[[column]][metrics$spec == spec], na.rm = TRUE)
}


# ---- Result 1: is the decomposition identified by RR alone? ----------------------------------------------------------
narrate_identifiability <- function(ident) guarded("results/identifiability.rds", {
  m <- as.data.frame(ident$metrics)
  value <- function(spec, column) spec_median(m, spec, column)
  collapse <- function(spec) mean(m$collapse_frac[m$spec == spec], na.rm = TRUE)
  n_recordings <- length(unique(m$SubjectID))
  per_spec <- function(spec) sprintf(
    "a median of %s distinct solutions per recording and one process absorbed essentially all variability in %s of the fits",
    num(value(spec, "n_solutions"), 1), pct(collapse(spec)))

  more_dispersed <- value("M0", "sd_ln_tau_s_equiv") > value("M3", "sd_ln_tau_s_equiv") * (1 + NARRATIVE$similar_rel) ||
    value("M0", "n_solutions") > value("M3", "n_solutions") * (1 + NARRATIVE$similar_rel) ||
    collapse("M0") > collapse("M3") + NARRATIVE$similar_rel
  datasets <- sort(unique(m$Dataset))
  larger_in <- datasets[vapply(datasets, function(d) {
    d_m <- m[m$Dataset == d, ]
    spec_median(d_m, "M0", "sd_ln_tau_s_equiv") > spec_median(d_m, "M3", "sd_ln_tau_s_equiv") * (1 + NARRATIVE$similar_rel)
  }, logical(1))]
  sd_both <- suppressWarnings(max(value("M0", "sd_ln_tau_s_equiv"), value("M3", "sd_ln_tau_s_equiv")))
  small_both <- is.finite(sd_both) && sd_both < 0.1
  ratio <- value("M3", "best_rmse_ms") / value("M0", "best_rmse_ms")

  consequence <- if (more_dispersed) {
    paste0("In these recordings, the RR series alone therefore left the latent decomposition less well determined than the ",
           "constrained specification, consistent with the premise that physiological information is needed to restrict the ",
           "admissible solutions. This shows that the constraints reduce the ambiguity; it does not show that the retained solution is the correct decomposition.")
  } else {
    paste0("In these recordings, the two specifications did not differ materially in the dispersion of solutions, so this analysis ",
           "by itself does not show that physiological information is required to restrict the latent solution space; that ",
           "premise then rests on the remaining analyses.")
  }
  paste0(
    "Without physiological information ($M_0$; ", n_recordings, " recordings), the multi-start fits of a recording reached ", per_spec("M0"),
    ". Under the full specification ($M_3$), the number of distinct solutions per recording was ", compare_phrase(value("M3", "n_solutions"), value("M0", "n_solutions")), " that under $M_0$ (median ",
    num(value("M3", "n_solutions"), 1), " versus ", num(value("M0", "n_solutions"), 1), "), and one process absorbed essentially all variability in ", pct(collapse("M3")), " versus ", pct(collapse("M0")), " of the fits (@tbl-ident-metrics). ",
    "The dispersion of $\\ln\\tau_s$ among near-equivalent fits under $M_3$ was ", compare_phrase(value("M3", "sd_ln_tau_s_equiv"), value("M0", "sd_ln_tau_s_equiv")),
    " that under $M_0$", if (length(larger_in)) paste0(" (larger under $M_0$ in ", count_of(length(larger_in), length(datasets)), " datasets: ", join_and(larger_in), ")") else "",
    if (small_both) "; it was small under both specifications, so the contrast rests mainly on the number of distinct solutions and on one-process solutions" else "",
    ". The best RR error under $M_3$ was ", compare_phrase(value("M3", "best_rmse_ms"), value("M0", "best_rmse_ms")), " that under $M_0$ (median ",
    num(value("M3", "best_rmse_ms")), " versus ", num(value("M0", "best_rmse_ms")), " ms",
    if (is.finite(ratio) && ratio > 1 + NARRATIVE$similar_rel) paste0("; a factor of ", num(ratio, 1), ", small in absolute terms but not in relative terms") else "", "). ", consequence)
})

# ---- Result 2: effect of the constraints, sensitivity, simulations -------------------------------------------------------
#' part = "main": the effect of the constraints, the prior-posterior contrast and the simulations; part = "supplement": the one-at-a-time sensitivity and the prior-posterior contrast by dataset.
narrate_constraints <- function(ident, sensitivity = NULL, sim_recovery = NULL, part = c("main", "supplement")) guarded("results/identifiability.rds", {
  part <- match.arg(part)
  extra <- ""
  m <- as.data.frame(ident$metrics)
  specs <- c("M0", "M1", "M2", "M3")
  collapse_by <- vapply(specs, function(s) mean(m$collapse_frac[m$spec == s], na.rm = TRUE), numeric(1))
  anchor_role <- if (is.finite(collapse_by[["M1"]]) && collapse_by[["M1"]] >= collapse_by[["M0"]] - 0.05 && collapse_by[["M3"]] < 0.5 * collapse_by[["M0"]] && collapse_by[["M2"]] < collapse_by[["M1"]]) {
    paste0("Time-scale priors alone ($M_1$) did not prevent one-process solutions, whereas the specifications that include the spectral anchor ($M_2$, $M_3$) largely did; the lower ambiguity under $M_3$ thus comes largely from fixing the variance ratio of the two drives through the anchor, ",
           "and it should not be read as evidence that the apportionment between the drives is correct.")
  } else ""
  text <- anchor_role                                          # the per-specification values are in the table

  if (!is.null(sensitivity)) {
    s <- as.data.frame(sensitivity$summary)
    hyper <- s[s$kind == "hyperparameter", ]
    loo <- s[s$kind == "leave_one_out", ]
    worst <- hyper[which.max(pmax(hyper$med_abs_ln_tau_p, hyper$med_abs_ln_tau_s)), ]
    extra <- paste0("Across ", length(unique(sensitivity$sensitivity$SubjectID)), " recordings, varying one hyperparameter at a time changed the time constants by a median absolute log-ratio of at most ",
                   num(max(hyper$med_abs_ln_tau_p, hyper$med_abs_ln_tau_s)), " (largest for ", setting_label(worst$setting), ") and the lowest median correlation with the reference trajectories was ",
                   num(min(hyper$med_cor_x_p, hyper$med_cor_x_s, na.rm = TRUE)), ". Removing one prior component at a time, the largest change of the time constants occurred without ",
                   term_label(loo$setting[which.max(pmax(loo$med_abs_ln_tau_p, loo$med_abs_ln_tau_s))]), " (median absolute log-ratio ",
                   num(max(loo$med_abs_ln_tau_p, loo$med_abs_ln_tau_s)), ").")
    pp <- as.data.frame(sensitivity$prior_posterior)
    text <- paste0(text, " The $M_3$ estimates lay inside the 95% prior interval in ", pct(mean(pp$in_prior_tau_p)), " ($\\tau_p$) and ", pct(mean(pp$in_prior_tau_s)),
                   " ($\\tau_s$) of recordings, with a median posterior contraction of ", num(median(pp$contraction_tau_p, na.rm = TRUE)), " and ",
                   num(median(pp$contraction_tau_s, na.rm = TRUE)), " for $\\ln\\tau_p$ and $\\ln\\tau_s$, respectively.",
                   if (mean(pp$in_prior_tau_p) < 0.5) paste0(" For the vagal-like drive, most estimates therefore lay outside the 95% prior interval, so the data moved $\\tau_p$ away from the central region of its prior; this may reflect slower fluctuations being absorbed by the faster drive, and it limits the reading of $\\tau_p$ as a strictly vagal time constant.") else "")
    if ("Dataset" %in% names(pp) && uniqueN(pp$Dataset) > 1) {
      by_dataset <- as.data.table(pp)[, .(p = mean(in_prior_tau_p), s = mean(in_prior_tau_s)), by = Dataset]
      extra <- paste0(extra, " By dataset, $\\tau_p$ lay inside its 95% prior interval in ", join_and(sprintf("%s (%s)", pct(by_dataset$p), by_dataset$Dataset)), " of the recordings and $\\tau_s$ in ",
                      join_and(sprintf("%s (%s)", pct(by_dataset$s), by_dataset$Dataset)), ".")
    }
  }
  if (!is.null(sim_recovery)) {
    t <- as.data.frame(sim_recovery$table)
    err <- function(scenario, spec) median(abs(t$ln_err_tau_s[t$scenario == scenario & t$spec == spec & t$ok]), na.rm = TRUE)
    text <- paste0(text, " In simulations with known truth, when the true time constants lay outside the prior intervals, the median absolute log-error of $\\tau_s$ was ",
                   num(err("off_prior", "M0")), " for $M_0$ and ", num(err("off_prior", "M3")), " for $M_3$ (the $M_3$ error was ",
                   compare_phrase(err("off_prior", "M3"), err("off_prior", "M0")), " the $M_0$ error); when they lay inside, ",
                   num(err("in_prior", "M0")), " and ", num(err("in_prior", "M3")), ". These values delimit the claim: constraints that help when the truth is compatible with them can bias the estimates when it is not.")
  }
  if (part == "supplement") return(trimws(extra))
  text
})


# ---- Landscape sentence ------------------------------------------------------------------------------------------------------
#' "The mode of the occupied region moved by X baseline (or rest) SDs from <phase 1> to <phase 2>..."
landscape_sentence <- function(modes, reference = "baseline") guarded("landscape modes", {
  modes <- as.data.frame(modes)
  first <- modes[1, ]
  others <- modes[-1, ]
  moves <- vapply(seq_len(nrow(others)), function(i) sprintf("%s %s-SD units to the %s phase", num(others$shift_from_first[i], 1), reference, others$phase[i]), "")
  paste0("In the state plane, the strongest basin of the occupied region (the mode of the density of the filtered states) was displaced from its position in the ",
         first$phase, " phase by ", join_and(moves), " (Supplementary Figure S6); the number of basins detected was ",
         paste(sprintf("%d in the %s phase", modes$n_basins, modes$phase), collapse = ", "), ". The landscape is a descriptive map of where the states lie and is not a fitted potential of the model.")
})


# ---- Result 3: orthostatic perturbation -------------------------------------------------------------------------------------
#' part = "main": the orthostatic result; part = "supplement": arterial pressure, dose-response, comparison with the CWT moving spectrum and state-space landscape.
narrate_prcp <- function(prcp, part = c("main", "supplement")) guarded("results/prcp.rds", {
  part <- match.arg(part)
  s <- prcp$summary
  r <- prcp$response
  n <- s$n_detected
  both <- s$frac_both_expected
  pattern <- if (both >= NARRATIVE$agree_share) "consistent with" else if (both <= NARRATIVE$disagree_share) "not consistent with" else "only partly consistent with"
  lag <- s$lag_s_minus_p["med"]
  timing <- if (!is.finite(lag)) "" else if (abs(lag) < NARRATIVE$sync_s) "at about the same time as" else if (lag > 0) "after" else "before"
  tau_phrase <- compare_phrase(s$tau_s_resp["med"], s$tau_p_resp["med"])
  bound <- sprintf("%s of the $x_p$ and %s of the $x_s$ step-response fits reached a search bound", pct(s$frac_fit_at_bound_p), pct(s$frac_fit_at_bound_s))

  tilts <- if (is.null(s$n_events)) "" else paste0(" (", s$n_events, " tilts, each analyzed against its own baseline; a recording contributes the median over its tilts)")
  consistency <- if (is.null(s$frac_events_both_expected)) "" else paste0(" The same pattern was present in ", pct(s$frac_events_both_expected), " of the individual tilts.")
  apportion <- " Because heart rate rises during tilt, the net change of the two states is largely fixed by RR; what the spectral anchor and the priors determine is how that change is apportioned between the two drives."
  text <- paste0(
    "In ", n, " recordings with a defined perturbation", tilts, ", the vagal-like state decreased and the sympathetic-like state increased in ", count_of(round(both * n), n),
    ", which is ", pattern, " the expected pattern of vagal withdrawal with sympathetic engagement.", consistency, apportion, " ",
    "The sympathetic-like response began ",
    if (timing == "at about the same time as") paste0("at about the same time as the vagal-like response (median difference ", num(abs(lag), 1), " s; within ", NARRATIVE$sync_s, " s counted as simultaneous)")
    else paste0("a median of ", num(abs(lag), 1), " s ", timing, " the vagal-like response"), ", and its time constant was ", tau_phrase,
    " that of the vagal-like response (", num(s$tau_s_resp["med"], 1), " versus ", num(s$tau_p_resp["med"], 1), " s); ", bound, ", so these constants are approximate.")

  main <- text; text <- ""
  if (s$n_bp > 0) {
    xs <- s$xcorr[s$xcorr$state == "X_s", ]
    xp <- s$xcorr[s$xcorr$state == "X_p", ]
    text <- paste0(text, " In the ", s$n_bp, " recordings with arterial pressure, MAP changed by a median of ", num(s$d_map["med"], 1), " mmHg (pulse pressure ", num(s$d_pp["med"], 1),
                   " mmHg) with an onset ", num(s$lat_map["med"], 1), " s after the tilt onset; after removing slow trends, the peak lagged correlation of MAP with the states was ",
                   if (nrow(xp)) paste0(num(xp$peak_r, 2), " for $x_p$ (lag ", num(xp$lag_s, 0), " s; surrogate $p < 0.05$ in ", pct(xp$frac_p05), " of recordings)") else "not computed for $x_p$",
                   " and ", if (nrow(xs)) paste0(num(xs$peak_r, 2), " for $x_s$ (lag ", num(xs$lag_s, 0), " s; ", pct(xs$frac_p05), ")") else "not computed for $x_s$", ".",
                   if (all(abs(c(xp$peak_r, xs$peak_r)) < 0.3, na.rm = TRUE) && (nrow(xp) + nrow(xs)) > 0) " These correlations are weak in magnitude." else "",
                   { per <- prcp$xcorr
                     if (is.null(per) || !nrow(per)) "" else {
                       sign_text <- function(state, label) { r <- per$peak_r[per$state == state]; r <- r[is.finite(r)]
                         if (!length(r)) return(NULL); neg <- sum(r < 0)
                         if (neg == 0 || neg == length(r)) sprintf("%s was %s in all %d recordings", label, if (neg == 0) "positive" else "negative", length(r))
                         else sprintf("%s was negative in %d and positive in %d of %d recordings", label, neg, length(r) - neg, length(r)) }
                       signs <- c(sign_text("X_p", "the correlation with $x_p$"), sign_text("X_s", "that with $x_s$"))
                       mixed <- any(grepl("negative in [0-9]+ and positive", signs))
                       if (length(signs)) paste0(if (mixed) " The sign differed between recordings: " else " The sign was consistent within each state: ", paste(signs, collapse = "; "), ".") else "" } })
  }
  if (s$n_dose > 0) {
    rho <- s$dose_rho
    text <- paste0(text, " Over the tilt-angle plateaus (", s$n_dose, " recordings with at least four levels), the within-recording rank correlation with the angle was ",
                   paste(sprintf("%s for %s", num(rho$rho_median, 2), rho$variable), collapse = ", "), " (median).")
  }
  text <- paste0(text, " The sympathetic-like state left its baseline band ", num(s$lead_x_s["med"], 0), " s and the LF band power ", num(s$lead_LF["med"], 0),
                 " s after the tilt onset (median), and for the vagal-like state and the HF band power ", num(s$lead_x_p["med"], 0), " and ", num(s$lead_HF["med"], 0), " s, respectively.", cwt_timing_sentence(s))
  if (!is.null(prcp$landscape_modes)) text <- paste0(text, " ", landscape_sentence(prcp$landscape_modes))
  if (part == "supplement") return(trimws(text))
  paste0(main, " With ", n, " recordings and with events partly defined from the same RR series, this is a test of physiological coherence and not evidence that the components correspond to the activity of the respective autonomic branches.")
})


#' Comparison of the departure times of the states and of the moving band power (tilt analysis); wording follows the numbers.
cwt_timing_sentence <- function(s) {
  a <- s$lead_x_s["med"]; b <- s$lead_LF["med"]; c <- s$lead_x_p["med"]; d <- s$lead_HF["med"]
  if (!all(is.finite(c(a, b, c, d)))) return("")
  early_s <- a < b; early_p <- c < d
  paste0(" The sympathetic-like state therefore ", if (early_s) "left its baseline band earlier than" else "did not leave its baseline band earlier than", " the LF band power, and the vagal-like state ",
         if (early_p) "earlier than" else "not earlier than", " the HF band power",
         if (d < 0) " (a negative value for the band power means that it departed from its baseline before the tilt onset, which is consistent with the spreading of the moving window across the transition and not with an anticipatory response)" else "",
         "; the comparison with the moving spectrum is ", if (early_s && early_p) "in favor of the filter for both states." else "not uniformly in favor of the filter.")
}

#' Discussion version of the same comparison.
cwt_timing_discussion <- function(prcp) guarded("results/prcp.rds", {
  s <- prcp$summary; a <- s$lead_x_s["med"]; b <- s$lead_LF["med"]; c <- s$lead_x_p["med"]; d <- s$lead_HF["med"]
  if (!all(is.finite(c(a, b, c, d)))) return("the timing of the reconstructed responses relative to the CWT moving band power is described in the Results")
  paste0("the sympathetic-like state ", if (a < b) "left its baseline band earlier" else "did not leave its baseline band earlier", " than the LF band power, whereas for the vagal-like state the comparison with the HF band power was ",
         if (c < d) "also in favor of the filter" else "not in favor of the filter (see the Results)")
})

#' Clause on the innovation diagnostics: the wording follows the share of recordings whose KS distance is below its critical value.
diagnostics_clause <- function(frac_below, pacf_viol = NA, brief = FALSE, numbers = TRUE) {
  if (!is.finite(frac_below)) return("[PENDING]")
  share <- function(prefix, suffix) if (numbers) paste0(" (", prefix, pct(frac_below), suffix, ")") else ""
  if (brief) return(if (frac_below < 0.5) paste0("showed residual structure in most recordings", share("KS distance below its critical value in ", ""))
                    else paste0("were largely compatible with the assumed innovation distribution", share("KS distance below its critical value in ", "")))
  if (frac_below < 0.5) {
    paste0("indicated residual serial structure or a departure from the assumed innovation distribution in most recordings", share("the KS distance was below its critical value in ", " of them"))
  } else {
    paste0("were largely compatible with the assumed innovation distribution", share("the KS distance was below its critical value in ", " of recordings"))
  }
}

# ---- Result 3b: responses by type of maneuver (PRCP) -------------------------------------------------------------------
#' Kinds of maneuver in plural, for sentences ("slow tilts", "stand-ups").
kind_plural <- function(kind) ifelse(kind == "stand-up", "stand-ups", paste0(kind, "s"))

#' Responses by kind of maneuver (slow tilt, rapid tilt, stand-up): direction, size and paired within-recording contrasts.
#' Empty when there is nothing to stratify (no protocol annotations).
narrate_prcp_kinds <- function(prcp) guarded("results/prcp.rds", {
  k <- prcp$by_kind
  if (is.null(k) || nrow(k) < 2) return("")
  n_both <- round(k$frac_both_expected * k$n_events)
  both <- paste0(n_both, " of ", k$n_events, " ", kind_plural(k$kind), ifelse(k$n_events >= 10, paste0(" (", pct(k$frac_both_expected), ")"), ""))
  text <- paste0("The protocol annotations distinguish ", join_and(kind_plural(k$kind)), " (@tbl-prcp-kinds). Both states changed in the expected direction in ", join_and(both),
                 ". The median standardized change of the vagal-like state was ", join_and(num(k$z_d_x_p_med, 1)), " baseline SDs and that of the sympathetic-like state ",
                 join_and(num(k$z_d_x_s_med, 1)), ", in the same order.")
  contrasts <- prcp$kind_contrasts
  if (!is.null(contrasts) && nrow(contrasts) && any(is.finite(contrasts$p))) {
    variable_text <- c(z_d_x_p = "the change of $x_p$", z_d_x_s = "the change of $x_s$", lat_p = "the onset latency of $x_p$", lat_s = "the onset latency of $x_s$")
    n_raw <- sum(contrasts$p < 0.05, na.rm = TRUE); n_bh <- sum(contrasts$p_BH < 0.05, na.rm = TRUE)
    pairs <- range(contrasts$n_pairs[contrasts$n_pairs > 0])
    which <- contrasts[is.finite(p) & p < 0.05]
    text <- paste0(text, " Paired within-recording contrasts (rapid minus slow tilt, and stand-up minus tilts; ", nrow(contrasts), " exploratory contrasts on ", pairs[1],
                   if (pairs[2] != pairs[1]) paste0("-", pairs[2]), " pairs of recording medians) ",
                   if (n_raw == 0) "showed no difference at $p < 0.05$; with so few recordings this is not evidence that the responses are equivalent."
                   else paste0("had $p < 0.05$ in ", n_raw, " (", n_bh, " after Benjamini-Hochberg adjustment): ",
                               join_and(sprintf("%s for %s (median difference %s, $p$ = %s)", which$contrast, variable_text[which$variable], num(which$median_diff, 1), num(which$p, 3))), "."))
    if (any(grepl("^lat", which$variable))) text <- paste0(text, " Onset latencies are measured from the annotated time of the maneuver, which marks the start of the posture change and not a calibrated onset of the stimulus (a negative latency means that the response starts before the annotated time), so contrasts in latency between kinds are not interpreted physiologically.")
  }
  text
})

#' Limitation sentence on how many events the responses by maneuver rest on.
kinds_limit_sentence <- function(prcp) guarded("results/prcp.rds", {
  k <- prcp$by_kind
  if (is.null(k) || nrow(k) < 2) return("")
  per_recording <- if (is.null(prcp$response_events_all)) NA_integer_ else max(prcp$response_events_all[, .N, by = .(SubjectID, kind)]$N)
  paste0("The responses by type of maneuver rest on ", paste(k$n_events, kind_plural(k$kind), collapse = ", "),
         if (is.finite(per_recording)) paste0(" (at most ", per_recording, " events of a kind per recording)"), " and are descriptive.")
})


# ---- Result 4: Fantasia ---------------------------------------------------------------------------------------------------------------
#' part = "main": group contrasts and stability; part = "supplement": incremental information over HRV, respiration and baroreflex.
narrate_fantasia <- function(fantasia, hrv_results = NULL, part = c("main", "supplement")) guarded("results/fantasia.rds", {
  part <- match.arg(part)
  groups <- table(fantasia$table$Group)
  contrasts <- as.data.frame(fantasia$contrasts[is.finite(fantasia$contrasts$cliff_delta), ])
  describe <- function(family) {
    d <- contrasts[contrasts$family == family, ]
    big <- d[abs(d$cliff_delta) >= NARRATIVE$effect_size, ]
    if (!nrow(big)) return(paste0("none of the ", nrow(d), " ", family, " variables reached $|\\delta| \\geq ", NARRATIVE$effect_size, "$"))
    top <- big[which.max(abs(big$cliff_delta)), ]
    paste0(nrow(big), " of ", nrow(d), " ", family, " variables reached $|\\delta| \\geq ", NARRATIVE$effect_size, "$, the largest being ",
           sprintf("%s (%s in older adults, $\\delta$ = %s%s)", label_of(top$variable), ifelse(top$cliff_delta > 0, "higher", "lower"),
                   num(top$cliff_delta, 2), ifelse(top$p_BH < 0.05, "; adjusted $p < 0.05$", "")))
  }
  text <- paste0("Among ", groups[["young"]], " young and ", groups[["older"]], " older adults, ", describe("latent"), "; ", describe("HRV"), ". ",
                 "These latent summaries are not independent: the share of $x_p$ variance is a monotone function of the SD ratio of the two states, both depend strongly on the spectral anchor $\\lambda_P$, $\\ln(\\kappa_p/\\kappa_s)$ combines the two time constants, ",
                 "and $\\sigma^2_{sys}$ and the SDs of the states scale with the overall amplitude of variability, so that the count above overstates the number of distinct dimensions in which the groups differ. ")
  if (!is.null(fantasia$stability)) {
    st <- fantasia$stability
    text <- paste0(text, "Agreement between two non-overlapping windows of the same recording (ICC of the log values) had a median of ", num(median(st$icc_log[st$family == "latent"], na.rm = TRUE)),
                   " for the latent summaries and ", num(median(st$icc_log[st$family == "HRV"], na.rm = TRUE)), " for the conventional HRV metrics",
                   if (any(is.finite(st$icc_log))) paste0("; the lowest agreement was that of ", sym_tex(st$variable[which.min(st$icc_log)]), " (", num(min(st$icc_log, na.rm = TRUE)), "). ") else ". ")
  }
  main <- text                                           # the incremental information (below) is reported in the supplement
  inc <- fantasia$incremental
  text <- paste0("The leave-one-out AUC for separating the groups was ", num(inc$auc_hrv), " with HRV metrics, ", num(inc$auc_latent), " with latent summaries and ", num(inc$auc_both),
                 " with both, and HRV metrics predicted the latent summaries with a median leave-one-out $R^2$ of ", num(median(inc$r2_latent_from_hrv, na.rm = TRUE)), " (low values indicate variation that is not linearly predicted by these metrics, which may include estimation error as well as additional information). The AUC values are uncertain with these group sizes and were not compared formally. ")
  ph <- fantasia$physiology
  if (!is.null(ph) && "f_resp" %in% names(ph) && any(is.finite(ph$f_resp))) {
    both_ok <- is.finite(ph$resp_sd_x_p) & is.finite(ph$resp_sd_x_s)
    text <- paste0(text, "In the ", sum(is.finite(ph$f_resp)), " recordings with respiration, the respiratory frequency lay within the HF band in ", pct(mean(ph$resp_in_HF[is.finite(ph$f_resp)])),
                   ", and the SD of the vagal-like state within the respiratory band was larger than that of the sympathetic-like state in ", pct(mean(ph$resp_sd_x_p[both_ok] > ph$resp_sd_x_s[both_ok])),
                   " of recordings (median ratio ", num(median(ph$resp_sd_x_p[both_ok] / ph$resp_sd_x_s[both_ok]), 1), "), as expected from the faster time scale assigned to the vagal-like state. ")
  }
  if (!is.null(ph) && "alpha_LF" %in% names(ph) && any(is.finite(ph$alpha_LF))) {
    pc <- fantasia$physiology_contrasts
    row <- if (!is.null(pc)) pc[pc$variable == "alpha_LF", ] else NULL
    text <- paste0(text, "The baroreflex alpha index (LF) was ", num(median(ph$alpha_LF, na.rm = TRUE), 1), " ms/mmHg (median, ", sum(is.finite(ph$alpha_LF)), " recordings)",
                   if (!is.null(row) && nrow(row) && is.finite(row$cliff_delta)) paste0(" and ", ifelse(row$cliff_delta > 0, "higher", "lower"), " in older adults ($\\delta$ = ", num(row$cliff_delta, 2), ")") else "", ". ")
  }
  if (part == "supplement") return(trimws(text))
  paste0(main, "With these group sizes, all comparisons are exploratory, and no continuous age, sex or other covariates were available.")
})


# ---- Result 5: exercise ----------------------------------------------------------------------------------------------------------------
#' part = "main": responses and recovery times; part = "supplement": geometry of the state-plane path and its density landscape.
narrate_exercise <- function(exercise, hrv_results = NULL, part = c("main", "supplement")) guarded("results/exercise.rds", {
  part <- match.arg(part)
  s <- exercise$summary
  km <- exercise$km
  recovery_text <- if (!is.null(km) && nrow(km)) {
    join_and(sprintf("%s for the %s-like state (%s of recordings censored)",
                     ifelse(is.finite(km$median), paste0(num(km$median, 0), " s"), "not estimable (median not reached)"),
                     ifelse(km$state == "p", "vagal", "sympathetic"), pct(km$frac_censored)))
  } else "not estimable"
  recovery_note <- if (identical(as.character(s$recovery_definition), "frac90")) {
    " For a first-order relaxation, recovery to 90% of the excursion takes about $2.3\\,\\tau$, so these times are of the order implied by the estimated time constants, which are themselves shaped by their priors; they describe the fitted kinetics and are not independent measurements of recovery."
  } else " These recovery times are tied to the estimated time constants, which are shaped by their priors, and are not independent measurements of recovery."
  text <- paste0(
    "In ", s$n_usable, " of ", s$n, " exercise recordings the rest, exercise and recovery phases could be defined from heart rate. From rest to the end of exercise, the vagal-like state changed by ",
    median_range_text(exercise$rest_response$z_d_x_p, 1), " and the sympathetic-like state by ", median_range_text(exercise$rest_response$z_d_x_s, 1), " rest SDs, with a heart-rate increase of ",
    median_range_text(exercise$rest_response$d_hr, 0), " bpm (median [2.5th, 97.5th percentile] across recordings). Recovery times (", words_for(as.character(s$recovery_definition), RECOVERY_DEFINITION_TEXT), "; Kaplan-Meier medians) were ", recovery_text,
    ".", recovery_note)
  geometry <- paste0("The path in the state plane had a median length of ", median_range_text(exercise$geometry$path_length, 0), ", an enclosed area of ", median_range_text(exercise$geometry$loop_area, 0),
    " and a maximum energy $U$ of ", median_range_text(exercise$geometry$U_max, 0), " in rest-SD units.")
  if (!is.null(exercise$landscape_modes)) geometry <- paste0(geometry, " ", landscape_sentence(exercise$landscape_modes, "rest"))
  if (part == "supplement") return(geometry)
  if (!is.null(hrv_results) && !is.null(hrv_results$loo_r2)) {
    loo <- as.data.frame(hrv_results$loo_r2)
    r2 <- loo[loo$Dataset == "exercise", ]
    if (nrow(r2)) text <- paste0(text, " Conventional HRV metrics predicted the latent summaries of this cohort with a median leave-one-out $R^2$ of ", num(median(r2$loo_r2, na.rm = TRUE)), ".")
  }
  if (!is.null(hrv_results) && !is.null(hrv_results$response_association)) {
    association <- as.data.frame(hrv_results$response_association)
    ra <- association[association$Dataset == "exercise", ]
    if (nrow(ra)) {
      top <- ra[which.max(abs(ra$spearman)), ]
      text <- paste0(text, " The strongest rank correlation between a latent response and an HRV change was ", num(top$spearman, 2), " (", words_for(top$latent, RESPONSE_TEXT), " with the ", words_for(top$hrv, RESPONSE_TEXT), ", n = ", top$n, ").")
    }
  }
  paste0(text, " No subject characteristics or clinical and functional outcomes were available for this cohort, so the analysis is descriptive and no association with such variables is reported.")
})


# ---- Abstract ---------------------------------------------------------------------------------------------------------------------------
abstract_identifiability <- function(ident) guarded("results/identifiability.rds", {
  m <- as.data.frame(ident$metrics)
  a <- spec_median(m, "M0", "n_solutions"); b <- spec_median(m, "M3", "n_solutions")
  c0 <- mean(m$collapse_frac[m$spec == "M0"], na.rm = TRUE); c3 <- mean(m$collapse_frac[m$spec == "M3"], na.rm = TRUE)
  sprintf("multi-start fits reached a median of %s distinct solutions per recording without physiological information versus %s with it, and one process absorbed essentially all variability in %s versus %s of the fits",
          num(a, 1), num(b, 1), pct(c0), pct(c3))
})

abstract_prcp <- function(prcp) guarded("results/prcp.rds", {
  s <- prcp$summary
  sprintf("in %d tilt-table recordings, the vagal-like state decreased and the sympathetic-like state increased in %s", s$n_detected, count_of(round(s$frac_both_expected * s$n_detected), s$n_detected))
})

abstract_aging <- function(fantasia) guarded("results/fantasia.rds", {
  g <- table(fantasia$table$Group)
  latent <- fantasia$contrasts[is.finite(fantasia$contrasts$cliff_delta) & fantasia$contrasts$family == "latent", ]
  big <- latent[abs(latent$cliff_delta) >= NARRATIVE$effect_size, ]
  sprintf("between %d young and %d older adults, %d of %d partly redundant latent summaries reached $|\\delta| \\geq %s$", g[["young"]], g[["older"]], nrow(big), nrow(latent), NARRATIVE$effect_size)
})

abstract_exercise <- function(exercise) guarded("results/exercise.rds", {
  km <- exercise$km
  sprintf("in %d exercise recordings, median recovery times, which are closely tied to the fitted time constants, were %s",
          exercise$summary$n_usable,
          join_and(sprintf("%s (%s)", ifelse(is.finite(km$median), paste0(num(km$median, 0), " s"), "not estimable"),
                           ifelse(km$state == "p", "vagal-like", "sympathetic-like"))))
})


# ---- Discussion and conclusion -----------------------------------------------------------------------------------------------------------
sensitivity_sentence <- function(sensitivity) {
  if (is.null(sensitivity)) return("")
  guarded("results/sensitivity.rds", {
    s <- as.data.frame(sensitivity$summary)
    h <- s[s$kind == "hyperparameter", ]
    paste0("Varying the physiological hyperparameters one at a time changed the time constants by a median absolute log-ratio of at most ", num(max(h$med_abs_ln_tau_p, h$med_abs_ln_tau_s)),
           " and left the reconstructed trajectories correlated at least ", num(min(h$med_cor_x_p, h$med_cor_x_s, na.rm = TRUE)), " with the reference. ",
           "These correlations are insensitive to a change of scale and therefore do not address how a change is apportioned between the two drives, which was the most sensitive aspect (see the limitations).")
  })
}

#' What the three biological analyses showed, in words and without numbers (they are in the Abstract and the Results), for the Discussion and
#' the Conclusion. The wording follows the same data-dependent rules as the Results.
findings_summary <- function(prcp, fantasia, exercise, inline = FALSE) guarded("results", {
  both <- prcp$summary$frac_both_expected
  tilt <- if (both >= NARRATIVE$agree_share) "the vagal-like state decreased and the sympathetic-like state increased in most tilt-table recordings"
          else if (both <= NARRATIVE$disagree_share) "the vagal-like and sympathetic-like states did not change in the expected direction in most tilt-table recordings"
          else "the vagal-like and sympathetic-like states changed in the expected direction in only part of the tilt-table recordings"
  latent <- fantasia$contrasts[is.finite(fantasia$contrasts$cliff_delta) & fantasia$contrasts$family == "latent", ]
  share_big <- mean(abs(latent$cliff_delta) >= NARRATIVE$effect_size)
  aging <- if (share_big >= 0.5) "most of the partly redundant latent summaries differed between young and older adults"
           else if (share_big > 0) "some of the partly redundant latent summaries differed between young and older adults"
           else "none of the latent summaries differed between young and older adults beyond the effect-size threshold"
  km <- exercise$km
  recovery <- if (!is.null(km) && nrow(km) && any(is.finite(km$median))) "the exercise recordings showed a rest-exercise-recovery trajectory whose recovery times are closely tied to the fitted time constants"
              else "the recovery times of the exercise recordings could not be estimated"
  text <- paste0("In the tilt-table, aging and exercise analyses, ", tilt, "; ", aging, "; and ", recovery, ".")
  if (inline) paste0(tolower(substring(text, 1, 1)), substring(text, 2)) else text        # inline: continues a sentence
})

#' Discussion sentence on the orthostatic pattern; its wording follows the share of recordings with the expected direction.
prcp_discussion_sentence <- function(prcp) guarded("results/prcp.rds", {
  s <- prcp$summary
  both <- s$frac_both_expected
  pattern <- if (both >= NARRATIVE$agree_share) "consistent with" else if (both <= NARRATIVE$disagree_share) "not consistent with" else "only partly consistent with"
  paste0("During orthostatic tilt, the reconstructed vagal-like state decreased and the sympathetic-like state increased in ", count_of(round(both * s$n_detected), s$n_detected),
         " recordings (", pct(both), "), a pattern ", pattern, " the expected autonomic adaptation to orthostatic stress.")
})

prcp_pattern_clause <- function(prcp) guarded("results/prcp.rds", {
  both <- prcp$summary$frac_both_expected
  paste0("(", pct(both), " of the ", prcp$summary$n_detected, " recordings showed this pattern)")
})

hrv_information_sentence <- function(fantasia, hrv_results = NULL) guarded("results/fantasia.rds", {
  inc <- fantasia$incremental
  r2 <- median(inc$r2_latent_from_hrv, na.rm = TRUE)
  paste0("In the Fantasia recordings, HRV metrics predicted the latent summaries with a median leave-one-out $R^2$ of ", num(r2),
         if (is.finite(r2) && r2 < 0.5) ", so that a substantial part of their variation was not linearly recoverable from SDNN, RMSSD and mean heart rate" else
           ", so that most of their variation was linearly recoverable from SDNN, RMSSD and mean heart rate in this sample", ".")
})

epsilon_sentence <- function(checks) guarded("results/supp_checks.rds", {
  e <- checks$fitted_regime$eps
  paste0("in the fitted recordings, $\\varepsilon = s_{RR}/\\Delta t$ had a median of ", num(median(e), 4), " (maximum ", num(max(e), 4), "), corresponding to a bound of at most ",
         num(max(checks$fitted_regime$bound_cdf), 4), " on the CDF discrepancy")
})

prcp_n_text <- function(prcp) guarded("results/prcp.rds", sprintf("%d recordings", prcp$summary$n))

conclusion_identifiability <- function(ident) guarded("results/identifiability.rds", {
  m <- as.data.frame(ident$metrics)
  reduced <- spec_median(m, "M0", "n_solutions") > spec_median(m, "M3", "n_solutions") * (1 + NARRATIVE$similar_rel) ||
    mean(m$collapse_frac[m$spec == "M0"], na.rm = TRUE) > mean(m$collapse_frac[m$spec == "M3"], na.rm = TRUE) + NARRATIVE$similar_rel ||
    spec_median(m, "M0", "sd_ln_tau_s_equiv") > spec_median(m, "M3", "sd_ln_tau_s_equiv") * (1 + NARRATIVE$similar_rel)
  if (reduced) "In the recordings analyzed, adding the physiological information reduced the number of distinct latent solutions and the occurrence of one-process solutions relative to RR alone, at a cost in RR error."
  else "In the recordings analyzed, adding the physiological information did not materially change the ambiguity of the latent solutions relative to RR alone, which limits the claim that it is required."
})



# ---- Configuration used ------------------------------------------------------------------------------------------------------------------------
#' "600 s starting at 600 s" etc. from the configuration stored with the results.
window_text <- function(config_used) guarded("results/config_used.rds", {
  w <- config_used$windows$fantasia
  sprintf("%s s starting at %s s of the recording", w$duration_s, w$start_s)
})


# ---- Captions and configuration -----------------------------------------------------------------------------------------------------------
#' "recording 12734, largest standardized response": which recording a figure shows
#' (results/figure_info.rds, written by analysis/13_figures.R).
example_caption <- function(figure_info, figure = "fig7_prcp") {
  info <- figure_info[[figure]]
  if (is.null(info)) return("the showcased recording")
  id <- if (!is.null(info$prcp)) info$prcp else info$exercise
  rule <- if (!is.null(info$prcp)) info$prcp_rule else info$exercise_rule
  paste0("recording ", id, ", ", rule)
}

#' "600-1200 s and 1200-1800 s": the two windows of the stability analysis.
stability_window_text <- function(config_used) guarded("results/config_used.rds", {
  w <- config_used$windows$fantasia
  sprintf("%s\u2013%s s and %s\u2013%s s", w$stability_A[1], w$stability_A[2], w$stability_B[1], w$stability_B[2])
})
