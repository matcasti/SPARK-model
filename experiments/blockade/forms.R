# ==============================================================================
# experiments/blockade/forms.R -- which functional form describes the non-vagal share R(f) = P_atropine(f) / P_saline(f) (= 1 - vagal share) of the blockade data,
# against the production logistic; leave-one-subject-out (LOSO) error on log R. Run from the project root:
#   Rscript experiments/blockade/forms.R [lnrr|rate] [spontaneous|paced_0.25]
# R is the power ratio per subject and bin (high-dose atropine over saline of day 1); ln RR is scale free (atropine raises the heart rate from ~60 to ~105 bpm).
# ==============================================================================
source("experiments/blockade/lib.R")
args <- commandArgs(TRUE); u_sel <- if (length(args)) args[1] else "lnrr"; br_sel <- if (length(args) > 1) args[2] else "spontaneous"
seg <- readRDS("experiments/blockade/_out/segments.rds")[units == u_sel & breathing == br_sel]
sal <- seg[cond == "saline" & day == 1, .(subject, bin, f, sal = power)]; atr <- seg[cond == "atropine_high", .(subject, bin, atr = power)]
d <- merge(sal, atr, by = c("subject", "bin"))[, y := log(atr / sal)]            # per subject and bin: ln of the power ratio
subs <- sort(unique(d$subject)); in_range <- BINS$f >= 0.04 & BINS$f <= 0.40
cat(sprintf("units %s, %s breathing: %d subjects (%s), %d bins\n", u_sel, br_sel, length(subs), paste(subs, collapse = ","), nrow(BINS)))

slope <- -log(0.05 / 0.95) / (0.15 - 0.10); prodR <- function(f) pmax(0.75 * (1 - 1 / (1 + exp(-slope * (f - 0.10)))), 1e-3)   # production non-vagal weight (floored for the log)
forms <- list(
  prod     = list(npar = 0, par0 = NULL, fn = function(p, f) prodR(f)),
  constant = list(npar = 1, par0 = qlogis(0.05), fn = function(p, f) rep(plogis(p[1]), length(f))),
  logistic = list(npar = 4, par0 = c(qlogis(0.2), qlogis(0.01), log(0.05), 1.5), fn = function(p, f) { lo <- plogis(p[1]); hi <- plogis(p[2]); hi + (lo - hi) / (1 + exp(exp(p[4]) * (log(f) - p[3]))) }),
  lorentz_floor = list(npar = 3, par0 = c(qlogis(0.01), qlogis(0.2), log(0.03)), fn = function(p, f) plogis(p[1]) + plogis(p[2]) / (1 + (f / exp(p[3]))^2)),
  power_floor   = list(npar = 3, par0 = c(qlogis(0.01), qlogis(0.2), log(1.5)),  fn = function(p, f) plogis(p[1]) + plogis(p[2]) * (f / 0.05)^(-exp(p[3]))))
fit_form <- function(fm, dd) {                                                  # least squares on ln R over the bin means of the training subjects
  m <- dd[, .(y = mean(y)), by = .(bin, f)]
  if (fm$npar == 0) return(NULL)
  best <- NULL; for (jit in c(0, 0.5, -0.5)) { o <- tryCatch(optim(fm$par0 + jit * seq_along(fm$par0) %% 2, function(p) sum((log(fm$fn(p, m$f)) - m$y)^2), method = "BFGS", control = list(maxit = 2000)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  best$par
}
pred <- function(nm, p, f) { fm <- forms[[nm]]; if (fm$npar == 0) fm$fn(NULL, f) else fm$fn(p, f) }
rows <- list(); fits <- list()
for (nm in names(forms)) { p_all <- fit_form(forms[[nm]], d); fits[[nm]] <- p_all
  err <- sapply(subs, function(s) { p <- fit_form(forms[[nm]], d[subject != s]); h <- d[subject == s]; e <- (log(pred(nm, p, h$f)) - h$y)^2; c(all = mean(e), inrange = mean(e[h$f >= 0.04 & h$f <= 0.40])) })
  m <- d[, .(y = mean(y)), by = .(bin, f)]; r <- (log(pred(nm, p_all, m$f)) - m$y)^2
  rows[[nm]] <- data.table(form = nm, n_par = forms[[nm]]$npar, fit_rmse_all = sqrt(mean(r)), fit_rmse_0.04_0.40 = sqrt(mean(r[m$f >= 0.04 & m$f <= 0.40])),
                           loso_rmse_all = sqrt(mean(err["all", ])), loso_rmse_0.04_0.40 = sqrt(mean(err["inrange", ]))) }
res <- rbindlist(rows); cat("\nRMSE of ln R (a value of 1 is a factor 2.7 in the non-vagal share); LOSO = held-out subject\n"); print(res[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)], row.names = FALSE)
cat("\nFitted parameters (R = non-vagal share):\n")
for (nm in c("constant", "logistic", "lorentz_floor", "power_floor")) { p <- fits[[nm]]
  s <- switch(nm, constant = sprintf("R = %.3f", plogis(p[1])), logistic = sprintf("R_low = %.3f, R_high = %.4f, centre %.3f Hz, steepness %.2f per ln f", plogis(p[1]), plogis(p[2]), exp(p[3]), exp(p[4])),
              lorentz_floor = sprintf("floor %.4f + %.3f / (1 + (f / %.4f Hz)^2)   [corner = 1/(2 pi tau): tau = %.0f s]", plogis(p[1]), plogis(p[2]), exp(p[3]), 1 / (2 * pi * exp(p[3]))),
              power_floor = sprintf("floor %.4f + %.3f (f / 0.05)^-%.2f", plogis(p[1]), plogis(p[2]), exp(p[3])))
  cat(sprintf("  %-14s %s\n", nm, s)) }
m <- d[, .(R_emp = exp(mean(y)), lo = exp(quantile(y, .1)), hi = exp(quantile(y, .9))), by = .(f)][order(f)]
m[, `:=`(prod = prodR(f), best_logistic = forms$logistic$fn(fits$logistic, f), lorentz_floor = forms$lorentz_floor$fn(fits$lorentz_floor, f))]
cat("\nNon-vagal share R(f) by bin (geometric mean over subjects, P10-P90 across subjects) against the forms:\n"); print(m[, lapply(.SD, function(x) signif(x, 3))], row.names = FALSE)
saveRDS(list(res = res, fits = fits, m = m, d = d), sprintf("experiments/blockade/_out/forms_%s_%s.rds", u_sel, br_sel))
