# ==============================================================================
# 06_fantasia.R -- Result 4: young vs older adults, stability, respiration and baroreflex.
#
# Reads   results/fit_fantasia.rds                    (from 01_fit_all.R) and the CSV files
#         (windows A and B are refitted for the test-retest analysis)
# Writes  results/fantasia.rds   table, contrasts, incremental, stability, physiology,
#                                physiology_contrasts, physiology_tests, physiology_assoc
# Method  Only group membership (young / older) is used as a subject characteristic: no
#         continuous age, sex or other covariates. Respiration and blood pressure exist only
#         for some subjects; their analyses use the subjects that have them.
#         The steps are in R/steps/fantasia.R.
# ==============================================================================
source("analysis/00_setup.R")

# ---- 1. Load -------------------------------------------------------------------------------
fits <- read_result("fit_fantasia", cfg)
dataset <- load_dataset("fantasia", cfg)

# ---- 2. Latent vs conventional HRV: contrasts and information ----------------------------------------
summary_table <- fantasia_summary_table(fits, dataset$rr, cfg)
contrasts <- fantasia_group_contrasts(summary_table, cfg)
incremental <- fantasia_incremental_info(summary_table)

# ---- 3. Test-retest stability across two windows ----------------------------------------------------------
stability <- fantasia_stability(fits$index$SubjectID, dataset, cfg)

# ---- 4. Respiration and blood-pressure context ----------------------------------------------------------------
physiology <- fantasia_physiology_table(fits, summary_table, cfg)
physiology_analysis <- fantasia_physiology_analysis(physiology, cfg)
log_msg("Physiology: respiration in %d, blood pressure in %d of %d recordings",
        sum(is.finite(physiology$f_resp)), sum(is.finite(physiology$sbp)), nrow(physiology))

# ---- 5. Save ----------------------------------------------------------------------------------------------------------
write_result(list(table = summary_table, contrasts = contrasts, incremental = incremental, stability = stability,
                  physiology = physiology, physiology_contrasts = physiology_analysis$contrasts,
                  physiology_tests = physiology_analysis$tests, physiology_assoc = physiology_analysis$associations),
             "fantasia", cfg)

log_msg("Fantasia: %d young, %d older", sum(summary_table$Group == "young"), sum(summary_table$Group == "older"))
print(contrasts[, .(variable, family, median1 = signif(median1, 3), median2 = signif(median2, 3),
                    cliff_delta = round(cliff_delta, 2), p = signif(p, 2), p_BH = signif(p_BH, 2))])
str(incremental)
shown <- intersect(c("SubjectID", "Group", "f_resp", "resp_sd_x_p", "resp_sd_x_s", "alpha_LF", "map"), names(physiology))
print(physiology[, ..shown])
str(physiology_analysis$tests)
if (!is.null(stability)) {
  print(stability[, .(variable, family, n, spearman = round(spearman, 2), icc_log = round(icc_log, 2))])
}
