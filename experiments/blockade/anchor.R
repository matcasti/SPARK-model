# ==============================================================================
# experiments/blockade/anchor.R -- what the empirical non-vagal share R(f) would do to the anchor lambda_P = V_p / V_s of the production recordings.
#   Rscript experiments/blockade/anchor.R         (needs _out/forms_lnrr_spontaneous.rds from forms.R and results/fit_*.rds of the pipeline)
# lambda = sum over the production range (0.04-0.40 Hz) of (1 - R(f)) P(f) / sum of R(f) P(f), P = band power of the recording's own heart-rate spectrum.
# Weights: production logistic (0.25 floor, crossover 0.10 Hz), the fitted free logistic and the Lorentzian-plus-floor of the blockade data (ln RR, spontaneous).
# ==============================================================================
source("experiments/blockade/lib.R")
fm <- readRDS("experiments/blockade/_out/forms_lnrr_spontaneous.rds"); p_lg <- fm$fits$logistic; p_lz <- fm$fits$lorentz_floor
slope <- -log(0.05 / 0.95) / (0.15 - 0.10)
R_prod <- function(f) 0.75 * (1 - 1 / (1 + exp(-slope * (f - 0.10))))
R_lg <- function(f) { lo <- plogis(p_lg[1]); hi <- plogis(p_lg[2]); hi + (lo - hi) / (1 + exp(exp(p_lg[4]) * (log(f) - p_lg[3]))) }
R_lz <- function(f) plogis(p_lz[1]) + plogis(p_lz[2]) / (1 + (f / exp(p_lz[3]))^2)
band <- BINS$f >= 0.04 & BINS$f <= 0.40
lam <- function(P, Rf) { R <- Rf(BINS$f[band]); sum((1 - R) * P[band]) / sum(R * P[band]) }
spectra <- function(tk) { sp <- segment_power(tk, 0, max(tk), "rate", trim = 0); if (is.null(sp)) NULL else sp$power }

# ---- 1. the Lancaster saline recordings (rest, spontaneous breathing, day 1) ----
seg <- readRDS("experiments/blockade/_out/segments.rds"); s1 <- seg[units == "rate" & cond == "saline" & day == 1 & breathing == "spontaneous"]
lan <- s1[, .(prod = lam(power, R_prod), logistic = lam(power, R_lg), lorentz = lam(power, R_lz)), by = subject]
cat("Lancaster saline spectra (heart rate, Hz), lambda_P under each weighting, per subject:\n"); print(lan[, lapply(.SD, function(x) round(x, 1))], row.names = FALSE)
cat(sprintf("median: production %.1f, free logistic %.0f, Lorentzian+floor %.0f\n", median(lan$prod), median(lan$logistic), median(lan$lorentz)))

# ---- 2. the 313 recordings of the pipeline, against the lambda_P stored in results/ ----
rows <- list()
for (ds in c("fantasia", "prcp", "exercise")) { f <- readRDS(sprintf("results/fit_%s.rds", ds)); sp <- as.data.table(f$spectral)
  for (id in unique(f$states$SubjectID)) { rr <- f$states[SubjectID == id, RR_Interval]; P <- spectra(cumsum(rr)); if (is.null(P)) next
    rows[[length(rows) + 1L]] <- data.table(dataset = ds, id = id, stored = sp[SubjectID == id, lambda_P], prod_bins = lam(P, R_prod), logistic = lam(P, R_lg), lorentz = lam(P, R_lz)) } }
r <- rbindlist(rows)
cat(sprintf("\n%d recordings. Bin approximation of the production partition against the stored lambda_P: Spearman %.3f, median ratio %.2f\n", nrow(r), cor(r$stored, r$prod_bins, method = "spearman"), median(r$prod_bins / r$stored)))
print(r[, .(n = .N, stored = round(median(stored), 2), prod_bins = round(median(prod_bins), 2), free_logistic = round(median(logistic)), lorentz_floor = round(median(lorentz)),
            logistic_p10 = round(quantile(logistic, .1)), logistic_p90 = round(quantile(logistic, .9))), by = dataset], row.names = FALSE)
cat(sprintf("Spearman between the production anchor and the blockade-derived anchor: free logistic %.2f, Lorentzian %.2f\n", cor(r$stored, r$logistic, method = "spearman"), cor(r$stored, r$lorentz, method = "spearman")))
saveRDS(list(lancaster = lan, recordings = r), "experiments/blockade/_out/anchor.rds")
