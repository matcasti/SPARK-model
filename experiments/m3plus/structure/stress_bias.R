# EXPERIMENTAL. Do the variants inflate (or deflate) the sympathetic state during the stress? Simulated exercise recordings with RSA and known profiles: estimated vs true change
# of x_s and x_p from rest to the end of exercise (windows as in the pipeline, true onset and end), and the peak. Run from the project root:  Rscript experiments/m3plus/structure/stress_bias.R
source("analysis/00_setup.R"); source("experiments/m3plus/lib.R")
src <- readLines("experiments/m3plus/compare/compare_sim.R"); DUR <- 420; N <- 100L; n_cores <- max(1L, min(5L, cfg$cores - 2L))
eval(parse(text = src[grep("^osc_draw", src):grep("^one <- function", src)[1] - 1]))
SPEC <- list(M3 = list(), M3O = list(osc = TRUE), M3A = list(adapt = 25), M3OA = list(osc = TRUE, adapt = 25))
one <- function(job) {
  z <- build("exercise", job$seed); s <- z$s; series <- prepare_series(s$rr); dy <- series$dy; tt <- cumsum(dy); off <- s$t[1] - s$rr[1]; t_on <- z$true[1] - off; t_end <- z$true[2] - off
  rest <- tt >= cfg$exercise$rest_skip_s & tt < t_on - 5; endw <- tt >= t_end - cfg$exercise$end_window_s & tt <= t_end; during <- tt >= t_on & tt <= t_end
  tr_s <- s$x_s; tr_p <- s$x_p + s$z1; truth <- c(ds = mean(tr_s[endw]) - mean(tr_s[rest]), dp = mean(tr_p[endw]) - mean(tr_p[rest]), ps = max(tr_s[during]) - mean(tr_s[rest]))
  rbindlist(lapply(names(SPEC), function(m) { f <- fit_ext2(series, SPEC[[m]]); e <- state_estimates_ext2(f); xp <- e$xp; xs <- e$xs
    data.table(model = m, sim = job$sim, ds = mean(xs[endw]) - mean(xs[rest]), dp = mean(xp[endw]) - mean(xp[rest]), ps = max(xs[during]) - mean(xs[rest]), ds_true = truth[["ds"]], dp_true = truth[["dp"]], ps_true = truth[["ps"]],
               rest_sd_s = sd(xs[rest]), rest_sd_s_true = sd(tr_s[rest]), tau_s = f$tau_s, tau_p = f$tau_p) })) }
state_estimates_ext2 <- function(f) { r <- f$run(full = TRUE); list(xp = r$X_p + if (!is.null(f$p$osc)) r$X_o else 0, xs = r$X_s) }
jobs <- lapply(seq_len(N), function(i) list(sim = i, seed = id_seed(paste("compare exercise", i), cfg$seed)))
res <- rbindlist(drop_failed(map_recordings(jobs, one, n_cores), "simulation")); saveRDS(res, "experiments/m3plus/_out/stress_bias.rds")
res[, `:=`(r_s = ds / ds_true, r_p = dp / dp_true, r_peak = ps / ps_true, r_sd = rest_sd_s / rest_sd_s_true)]
cat(sprintf("\nSimulated exercise recordings: %d. Ratio estimated / true (median [IQR]); 1 = unbiased, > 1 inflated\n", uniqueN(res$sim)))
q <- function(x) sprintf("%.2f [%.2f, %.2f]", median(x), quantile(x, .25), quantile(x, .75))
print(res[, .(sympathetic_change = q(r_s), sympathetic_peak = q(r_peak), vagal_change = q(r_p), pct_symp_overestimated = round(100 * mean(r_s > 1.1), 1), pct_symp_underestimated = round(100 * mean(r_s < 0.9), 1), rest_sd_s_ratio = q(r_sd)), by = model][order(match(model, names(SPEC)))], row.names = FALSE)
cat("\nAbsolute error of the sympathetic change (Hz, median):\n"); print(res[, .(abs_err_s = round(median(abs(ds - ds_true)), 3), bias_s = round(median(ds - ds_true), 3), true_change_median = round(median(ds_true), 3)), by = model][order(match(model, names(SPEC)))], row.names = FALSE)
