# Tests of the text generators in R/report/manuscript_text.R: the wording must follow the outcome (both
# directions), degrade to "[PENDING" when results are missing, and never contain "NA" / "Inf".
#   Rscript tests/test_narrative.R
for (f in c("config.R")) source(f)
suppressPackageStartupMessages(library(data.table))
source("R/report/manuscript_text.R")
n_pass <- 0L
check <- function(name, cond) { if (!isTRUE(cond)) stop("FAILED: ", name, call. = FALSE); n_pass <<- n_pass + 1L; cat("ok   ", name, "\n") }
options(spark.narrative.debug = FALSE)

# Helpers ---------------------------------------------------------------------------------------------
check("compare_phrase: lower / similar / higher", identical(c(compare_phrase(1, 2), compare_phrase(1, 1.1), compare_phrase(3, 2)), c("lower than", "similar to", "higher than")))
check("join_and", join_and(c("a", "b", "c")) == "a, b and c" && join_and("a") == "a")
check("num: no negative zero, vectorized, non-finite becomes pending", num(-0.001, 1) == "0.0" && length(num(c(1, 2, NA))) == 3 && num(NA) == "[PENDING]")
check("txt catches Inf and errors", txt(Inf) == "[PENDING]" && txt(stop("x")) == "[PENDING]" && txt("ok") == "ok")
check("guarded: NULL input becomes pending", grepl("PENDING", narrate_prcp(NULL)) && grepl("PENDING", narrate_fantasia(NULL)) && grepl("PENDING", narrate_identifiability(NULL)))

# Identifiability: M0 more dispersed than M3 vs no difference ---------------------------------------------------
fake_metrics <- function(sd_m0, sd_m3, nsol_m0, nsol_m3, collapse_m0 = 0.3) {
  rows <- CJ(SubjectID = paste0("s", 1:6), spec = c("M0", "M1", "M2", "M3"))
  rows[, Dataset := rep(c("prcp", "fantasia", "exercise"), length.out = .N)]
  rows[, `:=`(n_solutions = ifelse(spec == "M0", nsol_m0, nsol_m3), n_equiv = 3, sd_ln_tau_s_equiv = ifelse(spec == "M0", sd_m0, sd_m3),
              sd_ln_tau_p_equiv = 0.1, share_p_range_equiv = 0.2, collapse_frac = ifelse(spec == "M0", collapse_m0, 0), best_rmse_ms = 1.0, frac_in_phys_range = 0.5)]
  list(metrics = rows)
}
ambiguous <- narrate_identifiability(fake_metrics(0.8, 0.1, 5, 1))
unambiguous <- narrate_identifiability(fake_metrics(0.1, 0.1, 1, 1, collapse_m0 = 0))
collapse_only <- narrate_identifiability(fake_metrics(0.1, 0.1, 1, 1, collapse_m0 = 0.4))
check("identifiability: ambiguity under M0 is reported as such", grepl("less well determined", ambiguous) && !grepl("did not differ materially", ambiguous))
check("identifiability: one-process solutions under M0 alone are reported as ambiguity", grepl("less well determined", collapse_only))
check("identifiability: no difference is reported as no difference", grepl("did not differ materially", unambiguous) && !grepl("less well determined", unambiguous))
check("conclusion follows the same rule", grepl("reduced the number of distinct latent solutions", conclusion_identifiability(fake_metrics(0.8, 0.1, 5, 1))) && grepl("did not materially change", conclusion_identifiability(fake_metrics(0.1, 0.1, 1, 1, collapse_m0 = 0))))

# PRCP: pattern shared by all recordings vs none ----------------------------------------------------------------------
fake_prcp <- function(frac_both) {
  mr <- function(m) c(med = m, lo = m - 1, hi = m + 1)
  list(summary = list(n = 8, n_detected = 8, frac_both_expected = frac_both, lag_s_minus_p = mr(10), tau_p_resp = mr(5), tau_s_resp = mr(30),
                      frac_fit_at_bound_p = 0, frac_fit_at_bound_s = 0, n_bp = 0, n_dose = 0, n_val = 0, lead_x_s = mr(5), lead_LF = mr(20), lead_x_p = mr(3), lead_HF = mr(10),
                      val_onset = mr(0), val_frac_within_15s = 1, xcorr = NULL, dose_rho = NULL), landscape_modes = NULL)
}
all_ok <- narrate_prcp(fake_prcp(1)); none_ok <- narrate_prcp(fake_prcp(0)); mixed <- narrate_prcp(fake_prcp(0.5))
check("PRCP: a shared pattern is called consistent", grepl("which is consistent with the expected", all_ok))
check("PRCP: an absent pattern is called not consistent", grepl("which is not consistent with the expected", none_ok))
check("PRCP: an intermediate share is called partly consistent", grepl("only partly consistent", mixed))
check("Discussion sentence on the tilt pattern follows the outcome", grepl("a pattern consistent with the expected", prcp_discussion_sentence(fake_prcp(1))) && grepl("a pattern not consistent with", prcp_discussion_sentence(fake_prcp(0))))
check("PRCP: a later sympathetic onset is described as after", grepl("a median of 10.0 s after the vagal-like", all_ok))

# Fantasia: strong group differences vs none -----------------------------------------------------------------------------------
fake_fantasia <- function(delta) {
  contrasts <- data.table(variable = c("Tau_P", "Tau_S", "SDNN", "RMSSD"), family = c("latent", "latent", "HRV", "HRV"), cliff_delta = delta, p_BH = ifelse(abs(delta) > 0.5, 0.01, 0.5))
  list(table = data.table(Group = factor(rep(c("young", "older"), each = 5), levels = c("young", "older"))), contrasts = contrasts, stability = NULL,
       incremental = list(auc_hrv = 0.7, auc_latent = 0.7, auc_both = 0.7, r2_latent_from_hrv = c(a = 0.2, b = 0.3)), physiology = NULL)
}
strong <- narrate_fantasia(fake_fantasia(c(-0.8, 0.6, -0.9, -0.7))); weak <- narrate_fantasia(fake_fantasia(c(0.05, 0.1, -0.1, 0.0)))
check("Fantasia: differences above the threshold are listed with direction", grepl("reached", strong) && grepl("lower in older adults", strong) && grepl("higher in older adults", strong))
check("Fantasia: no difference is reported as none", grepl("none of the 2 latent", weak) && !grepl("lower in older", weak))

# Missing optional inputs never break the text and no "NA" / "Inf" leaks ------------------------------------------------------------
all_text <- c(ambiguous, unambiguous, all_ok, none_ok, mixed, strong, weak)
# PRCP: responses by type of maneuver ----------------------------------------------------------------------------------
kinds_fake <- list(by_kind = data.table(kind = c("slow tilt", "rapid tilt", "stand-up"), n_events = c(20L, 20L, 18L), n_recordings = 10L, frac_both_expected = c(1, 0.95, 0.89),
                                        z_d_x_p_med = c(-1, -2, -3), z_d_x_s_med = c(5, 6, 7)),
                   kind_contrasts = data.table(contrast = c("rapid - slow tilt", "stand-up - tilts"), variable = "z_d_x_s", n_pairs = 10L, median_diff = c(1, 0.5), p = c(0.01, 0.4), p_BH = c(0.02, 0.4)),
                   response_events_all = data.table(SubjectID = c("a", "a", "b"), kind = c("slow tilt", "slow tilt", "stand-up")))
kinds_text <- narrate_prcp_kinds(kinds_fake)
check("PRCP by kind: the kinds, the shares in the expected direction and the contrasts below 0.05 are reported",
      grepl("slow tilts, rapid tilts and stand-ups", kinds_text) && grepl("19 of 20 rapid tilts \\(95%\\)", kinds_text) && grepl("had \\$p < 0.05\\$ in 1", kinds_text) && grepl("2 exploratory contrasts", kinds_text))
check("PRCP by kind: without protocol annotations (one kind or none) there is nothing to report", narrate_prcp_kinds(list(by_kind = NULL)) == "" && kinds_limit_sentence(list(by_kind = NULL)) == "")
check("PRCP by kind: the limitation sentence gives the events per kind and the maximum per recording", grepl("20 slow tilts, 20 rapid tilts, 18 stand-ups", kinds_limit_sentence(kinds_fake)) && grepl("at most 2 events", kinds_limit_sentence(kinds_fake)))
check("no 'NA', 'Inf' or 'NaN' tokens in the generated text", !any(grepl("(^|[^A-Za-z])(NA|Inf|NaN)([^A-Za-z]|$)", all_text)))
check("setting labels are readable", setting_label("tau_p_bounds=narrow") == "the $\\tau_p$ prior interval (narrow)" && term_label("M3 without anchor") == "the spectral anchor")
cat(sprintf("\n%d narrative checks passed.\n", n_pass))
