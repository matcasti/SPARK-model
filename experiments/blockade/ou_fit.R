# ==============================================================================
# experiments/blockade/ou_fit.R -- the model's own spectrum fitted to the blockade spectra of each subject (Whittle likelihood on the periodogram of ln RR).
#   Rscript experiments/blockade/ou_fit.R [test]       "test" = recovery check on spectra drawn from the model
# Per subject and breathing condition, three conditions are fitted together:
#   saline    S = S_p(tau_p, V_p) + S_s(tau_s, V_s) + respiratory peak (A, f_r, w) + floor_sal
#   atropine  S = kappa S_s(tau_s, V_s) + floor_atr        (the vagal drive is blocked; kappa lets the slow component change by reflex readjustment)
#   double    S = floor_dbl                                (nothing left but the floor)
# S_i(f) = 4 V tau / (1 + (2 pi f tau)^2) is the spectrum of an OU drive with stationary variance V. The 5-min segments resolve f >= 0.011 Hz: only V_s / tau_s is
# identified, so tau_s is fixed on a grid and V_s (hence lambda_P = V_p / V_s) is reported for each value. Output: _out/ou_fit.rds
# ==============================================================================
source("experiments/blockade/lib.R")
args <- commandArgs(TRUE); test <- length(args) && args[1] == "test"
#' Negative joint Whittle log-likelihood. th (log scale): tau_p, V_p, V_s, A, w, floor_sal, kappa, floor_atr, floor_dbl; f_r is taken from `fr`. d = list(sal, atr, dbl) of pgram().
nll <- function(th, d, tau_s, fr, osc = TRUE, kfix = 1) {
  tau_p <- exp(th[1]); Vp <- exp(th[2]); Vs <- exp(th[3]); A <- if (osc) exp(th[4]) else 0; w <- exp(th[5]); c1 <- exp(th[6]); kap <- if (is.na(kfix)) exp(th[7]) else kfix; c2 <- exp(th[8]); c3 <- exp(th[9])
  S1 <- lor(d$sal$f, Vp, tau_p) + lor(d$sal$f, Vs, tau_s) + peak(d$sal$f, A, fr, w) + c1; S2 <- kap * lor(d$atr$f, Vs, tau_s) + c2; S3 <- rep(c3, length(d$dbl$f))
  sum(log(S1) + d$sal$I / S1) + sum(log(S2) + d$atr$I / S2) + (if (!is.null(d$dbl)) sum(log(S3) + d$dbl$I / S3) else 0)
}
LOW <- log(c(1, 1e-7, 1e-7, 1e-9, 0.004, 1e-7, 0.05, 1e-7, 1e-7)); UP <- log(c(40, 1e-1, 1e0, 1e-1, 0.12, 1e-2, 20, 1e-2, 1e-2))
fit_one <- function(d, tau_s, fr, osc = TRUE, kfix = 1) {
  best <- NULL; for (tp in c(3, 8, 20)) for (vp in c(1e-4, 1e-3)) { th0 <- log(c(tp, vp, vp / 3, vp / 3, 0.02, 1e-5, 1, 1e-5, 1e-5))
    o <- tryCatch(optim(th0, nll, d = d, tau_s = tau_s, fr = fr, osc = osc, kfix = kfix, method = "L-BFGS-B", lower = LOW, upper = UP, control = list(maxit = 500)), error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o }
  best
}

if (test) {                                                                       # recovery on spectra drawn from the model (exponential periodogram)
  set.seed(3); f <- seq(0.0107, 0.5, by = 1 / 280); tru <- c(tau_p = 5, V_p = 8e-4, V_s = 2e-4, A = 4e-4, w = 0.02, fl = 2e-6); tau_s <- 100; fr <- 0.25; res <- list()
  for (r in 1:20) { S1 <- lor(f, tru["V_p"], tru["tau_p"]) + lor(f, tru["V_s"], tau_s) + peak(f, tru["A"], fr, tru["w"]) + tru["fl"]; S2 <- 1.0 * lor(f, tru["V_s"], tau_s) + tru["fl"]
    d <- list(sal = list(f = f, I = S1 * rexp(length(f))), atr = list(f = f, I = S2 * rexp(length(f))), dbl = list(f = f, I = tru["fl"] * rexp(length(f))))
    o <- fit_one(d, tau_s, fr); res[[r]] <- exp(o$par[c(1, 2, 3, 4, 7)]); res[[r]][5] <- 1 }
  m <- do.call(rbind, res); cat("recovery from 20 simulated subjects (median [P10, P90] / truth)\n"); tr <- c(5, 8e-4, 2e-4, 4e-4, 1)
  for (i in 1:5) cat(sprintf("  %-6s %.2f [%.2f, %.2f]\n", c("tau_p", "V_p", "V_s", "A", "kappa")[i], median(m[, i]) / tr[i], quantile(m[, i], .1) / tr[i], quantile(m[, i], .9) / tr[i])); quit(save = "no") }

seg <- fread("experiments/blockade/_out/protocol.csv"); rows <- list(); TAU_S <- c(20, 50, 100, 200, 300)
for (br in c("spontaneous", "paced_0.25")) for (s in 1:7) {
  g <- function(cn, dy) { r <- seg[subject == s & cond == cn & day == dy & breathing == br]; if (!nrow(r)) return(NULL)
    tk <- clean_peaks(scan(sprintf("experiments/blockade/_out/rpeaks/rpeaks_Subject%d_day%d.csv", s, dy), quiet = TRUE)); pgram(tk, r$t0[1], r$t1[1]) }
  d <- list(sal = g("saline", 1), atr = g("atropine_high", 1), dbl = g("double", 1)); if (is.null(d$sal) || is.null(d$atr) || is.null(d$dbl)) next
  fr <- if (br == "paced_0.25") 0.25 else { f <- d$sal$f; m <- f >= 0.15 & f <= 0.40; sm <- as.numeric(stats::filter(d$sal$I[m], rep(1 / 5, 5), sides = 2)); sm[is.na(sm)] <- 0; f[m][which.max(sm)] }
  for (kf in c(1, NA)) for (ts_ in TAU_S) { o <- fit_one(d, ts_, fr, kfix = kf); o0 <- fit_one(d, ts_, fr, osc = FALSE, kfix = kf); p <- exp(o$par); if (!is.na(kf)) p[7] <- kf
    rows[[length(rows) + 1L]] <- data.table(breathing = br, subject = s, kappa_mode = if (is.na(kf)) "free" else "fixed 1", tau_s = ts_, f_r = fr, tau_p = p[1], V_p = p[2], V_s = p[3], A = p[4], w = p[5], kappa = p[7], floor_sal = p[6], floor_atr = p[8], floor_dbl = p[9],
      lambda_ou = p[2] / p[3], lambda_tot = (p[2] + p[4]) / p[3], nll = o$value, nll_no_osc = o0$value, tau_p_no_osc = exp(o0$par[1]), lambda_no_osc = exp(o0$par[2]) / exp(o0$par[3])) } }
r <- rbindlist(rows); saveRDS(r, "experiments/blockade/_out/ou_fit.rds")
cat("Fits by subject, breathing and fixed tau_s. Median over subjects (P10-P90):\n")
q <- function(x) sprintf("%.3g [%.3g, %.3g]", median(x), quantile(x, .1), quantile(x, .9))
print(r[, .(n = .N, tau_p = q(tau_p), lambda_ou = q(lambda_ou), lambda_total = q(lambda_tot), kappa = q(kappa), dNLL_osc = q(nll_no_osc - nll)), by = .(breathing, kappa_mode, tau_s)][order(breathing, kappa_mode, tau_s)], row.names = FALSE)
