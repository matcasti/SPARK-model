# ==============================================================================
# R/core/stats.R -- statistical helpers.
#   group_contrast()   older - young: medians, Hodges-Lehmann shift, Wilcoxon p, Cliff's delta (bootstrap CI)
#   cliffs_delta()     effect size with a percentile bootstrap interval
#   loo_r2(), loo_auc()   leave-one-out R^2 (OLS) and AUC (logistic) with small samples
#   km_median()        Kaplan-Meier median with right censoring
#   icc_consistency()  ICC(3,1);   paired_comparison()   paired Wilcoxon;   spearman_table()
# ==============================================================================

#' Cliff's delta with a percentile bootstrap CI (x vs y; positive when x > y).
cliffs_delta <- function(x, y, n_boot = 2000L, seed = 1L) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  delta_of <- function(a, b) mean(sign(outer(a, b, "-")))
  bootstrap <- with_seed(seed, replicate(n_boot, delta_of(sample(x, replace = TRUE), sample(y, replace = TRUE))))
  c(delta = delta_of(x, y), lo = unname(quantile(bootstrap, 0.025)), hi = unname(quantile(bootstrap, 0.975)))
}

#' Two-group contrast of one variable (second level minus first): medians, Hodges-Lehmann shift
#' with its CI, Wilcoxon p-value and Cliff's delta. Groups with fewer than 3 finite values return counts only.
group_contrast <- function(dt, var, group = "Group", lev = c("young", "older"), n_boot = 2000L) {
  first <- dt[[var]][dt[[group]] == lev[1]]
  second <- dt[[var]][dt[[group]] == lev[2]]
  first <- first[is.finite(first)]
  second <- second[is.finite(second)]
  if (length(first) < 3 || length(second) < 3) {
    return(data.table(variable = var, n1 = length(first), n2 = length(second)))
  }
  wilcoxon <- suppressWarnings(wilcox.test(second, first, conf.int = TRUE, exact = FALSE))
  delta <- cliffs_delta(second, first, n_boot)
  data.table(variable = var, n1 = length(first), n2 = length(second),
             median1 = median(first), median2 = median(second),
             hl_shift = unname(wilcoxon$estimate), hl_lo = wilcoxon$conf.int[1], hl_hi = wilcoxon$conf.int[2],
             p = wilcoxon$p.value, cliff_delta = delta[["delta"]], cliff_lo = delta[["lo"]], cliff_hi = delta[["hi"]])
}

#' Leave-one-out predictive R^2 of y on the predictors X (rows = recordings), by OLS.
loo_r2 <- function(y, X) {
  X <- as.data.frame(X)
  complete <- complete.cases(cbind(y, X))
  y <- y[complete]
  X <- X[complete, , drop = FALSE]
  n <- length(y)
  if (n < ncol(X) + 4) return(NA_real_)
  predicted <- vapply(seq_len(n), function(i) {
    model <- lm(y[-i] ~ ., data = X[-i, , drop = FALSE])
    as.numeric(predict(model, X[i, , drop = FALSE]))
  }, numeric(1))
  1 - sum((y - predicted)^2) / sum((y - mean(y))^2)
}

#' Leave-one-out AUC of a logistic model of a binary outcome on the (standardized) predictors X.
loo_auc <- function(y, X) {
  X <- as.data.frame(X)
  complete <- complete.cases(cbind(y, X))
  y <- y[complete]
  X <- X[complete, , drop = FALSE]
  n <- length(y)
  if (n < 8 || length(unique(y)) < 2) return(NA_real_)
  X <- as.data.frame(scale(X))
  predicted <- vapply(seq_len(n), function(i) {
    model <- suppressWarnings(glm(y[-i] ~ ., data = X[-i, , drop = FALSE], family = binomial(),
                                  control = glm.control(maxit = 50)))
    as.numeric(predict(model, X[i, , drop = FALSE], type = "response"))
  }, numeric(1))
  ranks <- rank(predicted)
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  (sum(ranks[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)       # Mann-Whitney form of the AUC
}

#' Kaplan-Meier median time with right censoring and its 95% CI (NA if the curve never reaches 0.5).
km_median <- function(time, event) {
  finite <- is.finite(time)
  if (sum(finite) < 3) return(c(median = NA_real_, lo = NA_real_, hi = NA_real_))
  fit <- survival::survfit(survival::Surv(time[finite], event[finite]) ~ 1)
  table <- suppressWarnings(summary(fit)$table)
  c(median = unname(table["median"]), lo = unname(table["0.95LCL"]), hi = unname(table["0.95UCL"]))
}


#' ICC(3,1): two-way mixed, consistency, single measurement (two windows as fixed raters).
icc_consistency <- function(a, b) {
  m <- cbind(a, b)
  n <- nrow(m); k <- 2; grand_mean <- mean(m)
  ms_rows <- k * sum((rowMeans(m) - grand_mean)^2) / (n - 1)
  residuals <- m - rowMeans(m) - matrix(colMeans(m), n, k, byrow = TRUE) + grand_mean
  ms_error <- sum(residuals^2) / ((n - 1) * (k - 1))
  (ms_rows - ms_error) / (ms_rows + (k - 1) * ms_error)
}

#' Paired comparison of two columns of `data` (rows with both finite): n, share where
#' `a > b`, median difference and the Wilcoxon signed-rank p-value; NULL with fewer than 4 pairs.
paired_comparison <- function(data, a, b) {
  both <- is.finite(data[[a]]) & is.finite(data[[b]])
  if (sum(both) < 4) return(NULL)
  list(
    n = sum(both), frac_a_greater = mean(data[[a]][both] > data[[b]][both]),
    median_diff = median(data[[a]][both] - data[[b]][both]),
    p = suppressWarnings(wilcox.test(data[[a]][both], data[[b]][both], paired = TRUE, exact = FALSE)$p.value)
  )
}

#' Spearman correlations of each `x_vars` column with each `y_vars` column (>= `min_n` complete pairs).
spearman_table <- function(data, x_vars, y_vars, x_label = "physio", y_label = "latent", min_n = 5) {
  rows <- lapply(x_vars, function(x) lapply(y_vars, function(y) {
    both <- is.finite(data[[x]]) & is.finite(data[[y]])
    if (sum(both) < min_n) return(NULL)
    out <- data.table(x, y, sum(both), suppressWarnings(cor(data[[x]][both], data[[y]][both], method = "spearman")))
    setnames(out, c(x_label, y_label, "n", "spearman"))
  }))
  rbindlist(unlist(rows, recursive = FALSE))
}
