# ==============================================================================
# experiments/blockade/extras/profile_lambda.R -- does the data of the production recordings prefer the blockade-derived anchor? M3 refitted with V_p/V_s fixed at several values,
# on 60 recordings (20 per dataset). Profiled log-likelihood (LL_star) and tau_p relative to the production anchor.   Rscript experiments/blockade/extras/profile_lambda.R
# ==============================================================================
source("analysis/00_setup.R")
set.seed(7); grid <- c(3, 10, 30, 100, 300); rows <- list()
for (ds in c("fantasia", "prcp", "exercise")) { f <- readRDS(sprintf("results/fit_%s.rds", ds)); ids <- unique(f$states$SubjectID); ids <- sample(ids, min(20, length(ids)))
  res <- map_recordings(ids, function(id) { rr <- f$states[SubjectID == id, RR_Interval]; s <- prepare_series(rr); a0 <- s$priors$energy_ratio
    one <- function(lp, lab) { s2 <- s; s2$priors$energy_ratio <- lp; fit <- fit_variant(s2, model_spec("M3")); if (!isTRUE(fit$ok)) return(NULL)
      data.table(dataset = ds, id = id, lambda = lab, lp = lp, LL = fit$LL_star, tau_p = fit$tau_p, tau_s = fit$tau_s) }
    rbindlist(c(list(one(a0, "anchor")), lapply(grid, function(g) one(g, as.character(g))))) }, max(1L, cfg$cores - 2L))
  rows[[ds]] <- rbindlist(drop_failed(res, "recording")) }
r <- rbindlist(rows); r[, dLL := LL - LL[lambda == "anchor"], by = .(dataset, id)]
cat("\nChange of the profiled log-likelihood and tau_p relative to the production anchor (median over recordings; negative dLL = worse fit)\n")
print(r[lambda != "anchor", .(n = .N, dLL_median = round(median(dLL), 1), dLL_p10 = round(quantile(dLL, .1), 1), dLL_p90 = round(quantile(dLL, .9), 1),
                              share_better = round(mean(dLL > 2), 2), tau_p_median = round(median(tau_p), 1)), by = .(lambda = factor(lambda, c(as.character(grid))), dataset)][order(dataset, lambda)], row.names = FALSE)
cat("production anchor: tau_p median", round(median(r[lambda == "anchor", tau_p]), 1), "s\n")
saveRDS(r, "experiments/blockade/_out/profile_lambda.rds")
