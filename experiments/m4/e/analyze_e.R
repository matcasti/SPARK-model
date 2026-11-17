# ==============================================================================
# experiments/m4/e/analyze_e.R -- EXPERIMENTAL. Tuning (development) and the decision rule of README.md for variants assembled from the stored components.
#   Rscript experiments/m4/e/analyze_e.R dev [tag]    -> tunes q_osc, G, c_p, c_s on sim_dev<tag>, saves frozen<tag>.rds, evaluates on development data
#   Rscript experiments/m4/e/analyze_e.R conf [tag]   -> evaluates the frozen settings on sim_conf<tag>, real_conf<tag>, blockade<tag>
# A variant = list(q_osc, G, cal = list(filter = c(p, s), smooth = c(p, s))). M4c-w* = always oscillator, G = 9, no rescaling.
# ==============================================================================
args <- commandArgs(TRUE); phase <- args[1]; tag <- if (length(args) > 1) args[2] else ""
suppressMessages(source("analysis/00_setup.R"))
out <- "experiments/m4/e/_out"; f_in <- function(k) file.path(out, sprintf("%s_%s%s.rds", k, phase, tag))
COV_GRID <- seq(1, 2.5, by = 0.1); cc <- function(prefix, c) sprintf("%s%02d", prefix, round(10 * c))
pick <- function(lr_osc, gain, q_osc, G, have, rt = "R") { have <- strsplit(have, ",")[[1]]
  lin <- if (is.finite(lr_osc) && lr_osc >= q_osc && "O" %in% have) "O" else "B"
  if (is.finite(gain) && gain >= G) { r <- if (rt %in% c("RB", "QB")) rt else paste0(rt, lin); if (r %in% have) return(r) }; lin }
REGC <- c("RB", "RO", "QB", "QO", "JO", "QBo", "QJo"); rt_of <- function(v) if (is.null(v$rt)) "R" else v$rt
have_of <- function(rows, keys) rows[, .(have = paste(unique(comp), collapse = ",")), by = keys]
pick_lr <- function(lr_osc, gain, q_osc, G, have, cRB, cRO, cQB, cQO, qk = 13.28) { have <- strsplit(have, ",")[[1]]
  lin <- if (is.finite(lr_osc) && lr_osc >= q_osc && "O" %in% have) "O" else "B"
  if (!(is.finite(gain) && gain >= G) || !("RB" %in% have)) return(lin)
  ob <- if ("RO" %in% have && is.finite(cRO) && 2 * (cRB - cRO) >= q_osc) "O" else "B"; cr <- if (ob == "O") cRO else cRB; cq <- if (ob == "O") cQO else cQB
  if (paste0("Q", ob) %in% have && is.finite(cq) && 2 * (cr - cq) >= qk) paste0("Q", ob) else paste0("R", ob) }
#' E5: with the gate open, kinetics check K = 2 (cost_RB - cost_QB); K >= qk (mismatch) -> v$mm (RB or QB); otherwise JO if 2 (cost_RB - cost_JO) >= q_osc, else RB.
pick_tree <- function(lr_osc, gain, q_osc, G, have, cRB, cQB, cJO, qk, mm) { have <- strsplit(have, ",")[[1]]
  lin <- if (is.finite(lr_osc) && lr_osc >= q_osc && "O" %in% have) "O" else "B"
  if (!(is.finite(gain) && gain >= G) || !("RB" %in% have)) return(lin)
  if (is.finite(cQB) && 2 * (cRB - cQB) >= qk) return(if (mm %in% have) mm else "RB")
  if ("JO" %in% have && is.finite(cJO) && 2 * (cRB - cJO) >= q_osc) "JO" else "RB" }
#' Hybrid (M4e-h): states from pick_tree(mm = "RB"); parameters (tau, SDs) from QJo when the kinetics check signals a mismatch, else from the state component.
pick_par <- function(state, gain, G, have, cRB, cQB, qk) { have <- strsplit(have, ",")[[1]]
  if (is.finite(gain) && gain >= G && is.finite(cQB) && is.finite(cRB) && 2 * (cRB - cQB) >= qk && "QJo" %in% have) "QJo" else state }
#' st with a column `chosen` for variant v (rows: the component rows with their costs; keys: the record keys).
choose <- function(st, rows, keys, v) { st <- copy(st)
  if (!rt_of(v) %in% c("LR", "TREE", "HYB")) return(st[, chosen := mapply(pick, lr_osc, gain, v$q_osc, v$G, have, rt_of(v))])
  cw <- dcast(unique(rows[, c(keys, "comp", "cost"), with = FALSE]), as.formula(paste(paste(keys, collapse = " + "), "~ comp")), value.var = "cost")
  for (k in c("RB", "RO", "QB", "QO", "JO")) if (!k %in% names(cw)) cw[, (k) := NA_real_]
  setnames(cw, c("RB", "RO", "QB", "QO", "JO"), c("cRB", "cRO", "cQB", "cQO", "cJO")); z <- merge(st, cw[, c(keys, "cRB", "cRO", "cQB", "cQO", "cJO"), with = FALSE], by = keys, all.x = TRUE)
  qk <- if (is.null(v$qk)) 13.28 else v$qk
  if (rt_of(v) %in% c("TREE", "HYB")) { z[, chosen := mapply(pick_tree, lr_osc, gain, v$q_osc, v$G, have, cRB, cQB, cJO, qk, if (is.null(v$mm)) "RB" else v$mm)]
    if (rt_of(v) == "HYB") z[, chosen_par := mapply(pick_par, chosen, gain, v$G, have, cRB, cQB, qk)]; return(z[]) }
  z[, chosen := mapply(pick_lr, lr_osc, gain, v$q_osc, v$G, have, cRB, cRO, cQB, cQO, qk)][] }

# ---- tuning on development simulations -----------------------------------------------------------------------------------------------------
S <- readRDS(file.path(out, sprintf("sim_%s%s.rds", phase, tag)))
assemble_sim <- function(v) { st <- choose(S$stats, S$rows, c("scenario", "sim"), v)
  x <- merge(S$rows, st[, .(scenario, sim, chosen)], by = c("scenario", "sim"))[comp == chosen]
  if ("chosen_par" %in% names(st)) { pr <- merge(S$rows, st[, .(scenario, sim, chosen_par)], by = c("scenario", "sim"))[comp == chosen_par, .(scenario, sim, mode, tp2 = err_tau_p, ts2 = err_tau_s)]
    x <- merge(x, pr, by = c("scenario", "sim", "mode")); x[, `:=`(err_tau_p = tp2, err_tau_s = ts2)][, c("tp2", "ts2") := NULL] }
  for (md in c("filter", "smooth")) { cp <- if (is.null(v$cal)) 1 else v$cal[[md]][1]; cs <- if (is.null(v$cal)) 1 else v$cal[[md]][2]
    x[mode == md, `:=`(cover_p = ifelse(comp %in% REGC, get(cc("cp", cp)), cp10), cover_s = ifelse(comp %in% REGC, get(cc("cs", cs)), cs10))] }
  x }
if (phase == "dev") {
  st <- S$stats; q_osc <- max(11.3, quantile(st[scenario == "rest", lr_osc], 0.95, na.rm = TRUE))
  null <- st[scenario %in% c("rest", "rsa")]; null[, has_reg := grepl("RB", have)]
  cand <- sort(unique(c(9, null[is.finite(gain) & gain >= 9, gain + 1e-6]))); G <- cand[which(vapply(cand, function(g) mean(null$has_reg & null$gain >= g, na.rm = TRUE) <= 0.02, TRUE))[1]]
  tune_cal <- function(rt) tune_cal_v(list(q_osc = q_osc, G = G, cal = NULL, rt = rt))
  tune_cal_v <- function(v0) { v0$cal <- NULL; x <- assemble_sim(v0)[scenario %in% c("regime", "exercise", "two") & comp %in% REGC]
    cal <- lapply(c(filter = "filter", smooth = "smooth"), function(md) { y <- x[mode == md]
      c(p = COV_GRID[which(vapply(COV_GRID, function(c) mean(y[[cc("cp", c)]]) >= 0.95, TRUE))[1]], s = COV_GRID[which(vapply(COV_GRID, function(c) mean(y[[cc("cs", c)]]) >= 0.95, TRUE))[1]]) })
    lapply(cal, function(z) { z[is.na(z)] <- max(COV_GRID); z }) }
  frozen <- list(M4e = list(q_osc = q_osc, G = G, cal = tune_cal("R"), rt = "R"))
  if (any(S$rows$comp %in% c("QB", "QO"))) frozen$M4e2 <- list(q_osc = q_osc, G = G, cal = tune_cal("Q"), rt = "Q")
  if (any(S$rows$comp %in% c("QB", "QO")) && "cost" %in% names(S$rows)) frozen$M4e3 <- list(q_osc = q_osc, G = G, cal = tune_cal("LR"), rt = "LR")
  frozen$M4e4 <- list(q_osc = q_osc, G = G, cal = tune_cal("RB"), rt = "RB")
  if (any(S$rows$comp == "QJo")) frozen$M4eh <- { v <- list(q_osc = q_osc, G = G, rt = "HYB", mm = "RB"); v$cal <- tune_cal_v(v); v }
  if (any(S$rows$comp == "JO")) for (mm in c("RB", "QB")) frozen[[paste0("M4e5", if (mm == "RB") "a" else "b")]] <- { v <- list(q_osc = q_osc, G = G, rt = "TREE", mm = mm); v$cal <- tune_cal_v(v); v }
  saveRDS(frozen, file.path(out, sprintf("frozen%s.rds", tag)))
  for (nm in names(frozen)) { cal <- frozen[[nm]]$cal; cat(sprintf("\nTuned on development (%s): q_osc = %.2f, G = %.2f, cal filter = (%.1f, %.1f), smooth = (%.1f, %.1f)\n", nm, q_osc, G, cal$filter[1], cal$filter[2], cal$smooth[1], cal$smooth[2])) }
  cat("Oscillator kept (share):\n"); print(st[, .(osc = round(mean(lr_osc >= q_osc, na.rm = TRUE), 2), gate = round(mean(gain >= G, na.rm = TRUE), 2)), by = scenario])
} else frozen <- readRDS(file.path(out, sprintf("frozen%s.rds", tag)))
VARIANTS <- c(list(M4cw = list(q_osc = -Inf, G = 9, cal = NULL)), frozen)
if (exists("EXTRA_VARIANTS")) VARIANTS <- c(VARIANTS, EXTRA_VARIANTS)

# ---- simulation criteria (S1', S2, S3) --------------------------------------------------------------------------------------------------------
m3 <- S$rows[comp == "M3"][, `:=`(cover_p = cp10, cover_s = cs10)]
sim_summary <- function(x) x[, .(rmse_p = mean(rmse_p), rmse_s = mean(rmse_s), cover_p = mean(cover_p), cover_s = mean(cover_s), tau_p = median(err_tau_p), tau_s = median(err_tau_s)), by = .(scenario, mode)]
b <- sim_summary(m3); setnames(b, c("rmse_p", "rmse_s", "cover_p", "cover_s", "tau_p", "tau_s"), paste0("m3_", c("rmse_p", "rmse_s", "cover_p", "cover_s", "tau_p", "tau_s")))
SC <- c("rest", "rsa", "regime", "exercise", "two"); stationary <- c("rest", "rsa")
sim_eval <- function(nm, v) { s <- merge(sim_summary(assemble_sim(v)), b, by = c("scenario", "mode"))[, variant := nm]
  s[, lim := ifelse(scenario %in% stationary, 1.02, 0.90)]
  s[, `:=`(S1 = rmse_p <= lim * m3_rmse_p & rmse_s <= lim * m3_rmse_s, S2 = cover_p >= pmin(0.90, m3_cover_p - 0.05) & cover_s >= pmin(0.90, m3_cover_s - 0.05), S3 = tau_p <= m3_tau_p + 0.10 & tau_s <= m3_tau_s + 0.10)] }
SE <- rbindlist(lapply(names(VARIANTS), function(nm) sim_eval(nm, VARIANTS[[nm]])))
cat("\n== Simulation (filtered): RMSE ratio to M3, coverage, tau errors ==\n")
print(SE[mode == "filter", .(variant, scenario, rp = round(rmse_p / m3_rmse_p, 3), rs = round(rmse_s / m3_rmse_s, 3), cov_p = round(cover_p, 3), cov_s = round(cover_s, 3), m3cov_p = round(m3_cover_p, 3), m3cov_s = round(m3_cover_s, 3),
               tp = round(tau_p, 2), m3tp = round(m3_tau_p, 2), ts = round(tau_s, 2), m3ts = round(m3_tau_s, 2), S1, S2, S3)][order(variant, match(scenario, SC))], row.names = FALSE)
cat("\n== Simulation (smoothed): RMSE ratio, coverage ==\n")
print(SE[mode == "smooth", .(variant, scenario, rp = round(rmse_p / m3_rmse_p, 3), rs = round(rmse_s / m3_rmse_s, 3), cov_p = round(cover_p, 3), cov_s = round(cover_s, 3), S2)][order(variant, match(scenario, SC))], row.names = FALSE)

# ---- real-data criteria (R1-R4) ---------------------------------------------------------------------------------------------------------------
rf <- f_in("real")
if (file.exists(rf)) { R <- readRDS(rf); hv <- have_of(R$rows, c("ds", "id", "part")); stR <- merge(R$stats, hv, by = c("ds", "id", "part"))
  icc <- function(a, b) { x <- cbind(a, b); x <- x[stats::complete.cases(x), , drop = FALSE]; n <- nrow(x); gm <- mean(x); msr <- 2 * sum((rowMeans(x) - gm)^2) / (n - 1)
    mse <- sum((x - rowMeans(x) - matrix(colMeans(x) - gm, n, 2, byrow = TRUE))^2) / (n - 1); (msr - mse) / (msr + mse) }
  boot_lo <- function(x) { x <- x[is.finite(x)]; if (length(x) < 5) return(NA_real_); quantile(with_seed(7L, replicate(2000, mean(sample(x, replace = TRUE)))), 0.025) }
  wins <- if (phase == "dev") c("A", "B") else c("C", "D")
  real_eval <- function(nm, v) { s <- if (nm == "M3") copy(stR)[, chosen := "M3"] else choose(stR, R$rows, c("ds", "id", "part"), v)
    x <- merge(R$rows, s[, .(ds, id, part, chosen)], by = c("ds", "id", "part"))[comp == chosen]
    if ("chosen_par" %in% names(s)) { pr <- merge(R$rows, s[, .(ds, id, part, chosen_par)], by = c("ds", "id", "part"))[comp == chosen_par, .(ds, id, part, a1 = tau_p, a2 = tau_s, a3 = SD_p, a4 = SD_s)]
      x <- merge(x, pr, by = c("ds", "id", "part")); x[, `:=`(tau_p = a1, tau_s = a2, SD_p = a3, SD_s = a4)][, c("a1", "a2", "a3", "a4") := NULL] }
    f <- x[part == "full"]
    r1 <- f[, .(ks = sum(ks_ok, na.rm = TRUE), n = .N, pacf = mean(pacf, na.rm = TRUE))]
    h <- x[part == "head" & tail_bad < 0.1]; m3h <- R$rows[comp == "M3" & part == "head", .(ds, id, m3 = heldout)]; h <- merge(h, m3h, by = c("ds", "id")); d <- h$heldout - h$m3
    ex <- f[ds == "exercise" & is.finite(d_vagal)]; tw <- dcast(x[part %in% wins], id + ds ~ part, value.var = c("tau_p", "tau_s", "SD_p", "SD_s")); g <- function(k) tw[[paste0(k, "_", wins[1])]]; g2 <- function(k) tw[[paste0(k, "_", wins[2])]]
    data.table(variant = nm, ks = r1$ks, n = r1$n, pacf = round(r1$pacf, 2), held = if (nm == "M3") 0 else mean(d, na.rm = TRUE), held_lo = if (nm == "M3") NA_real_ else boot_lo(d),
               signs = round(100 * mean(ex$d_vagal < 0 & ex$d_symp > 0)), vfirst = round(100 * mean(ex$t50_vagal < ex$t50_symp, na.rm = TRUE)),
               icc_tp = icc(log(g("tau_p")), log(g2("tau_p"))), icc_ts = icc(log(g("tau_s")), log(g2("tau_s"))), icc_sdp = icc(log(g("SD_p")), log(g2("SD_p"))), icc_sds = icc(log(g("SD_s")), log(g2("SD_s"))),
               osc_share = mean(x$comp %in% c("O", "RO", "QO", "JO", "QJo")), gate_share = mean(x$comp %in% REGC),
               wsd_tp = sd(log(g("tau_p")) - log(g2("tau_p")), na.rm = TRUE) / sqrt(2), wsd_sdp = sd(log(g("SD_p")) - log(g2("SD_p")), na.rm = TRUE) / sqrt(2), wsd_sds = sd(log(g("SD_s")) - log(g2("SD_s")), na.rm = TRUE) / sqrt(2),
               age_r2 = { m <- (log(g("tau_p")) + log(g2("tau_p"))) / 2; grp <- grepl("y", tw$id); summary(lm(m ~ grp))$r.squared }) }
  RE <- rbindlist(lapply(c("M3", names(VARIANTS)), function(nm) real_eval(nm, VARIANTS[[nm]])))
  m3r <- RE[variant == "M3"]; RE[, `:=`(R1 = ks >= m3r$ks & pacf <= m3r$pacf, R2 = held > 0 & held_lo > 0, R3 = signs >= 95 & vfirst >= 95, R4 = icc_tp >= m3r$icc_tp - 0.05 & icc_sdp >= m3r$icc_sdp - 0.05 & icc_sds >= m3r$icc_sds - 0.05,
                R4p = wsd_tp <= 1.10 * m3r$wsd_tp & wsd_sdp <= 1.10 * m3r$wsd_sdp & wsd_sds <= 1.10 * m3r$wsd_sds & age_r2 >= m3r$age_r2 - 0.05)]
  cat("\n== Real recordings ==\n"); print(RE[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)], row.names = FALSE)
  cat("\nBy dataset (component shares, KS, held-out):\n")
  for (nm in names(VARIANTS)) { v <- VARIANTS[[nm]]; s <- choose(stR, R$rows, c("ds", "id", "part"), v); x <- merge(R$rows, s[, .(ds, id, part, chosen)], by = c("ds", "id", "part"))[comp == chosen]
    m3h <- R$rows[comp == "M3" & part == "head", .(ds, id, m3 = heldout, m3ks = NA)]; print(x[, .(variant = nm, osc = round(mean(comp %in% c("O", "RO", "QO", "JO", "QJo")), 2), gate = round(mean(comp %in% REGC), 2), ks_full = sum(ks_ok[part == "full"], na.rm = TRUE)), by = ds], row.names = FALSE) }
  print(R$rows[comp == "M3" & part == "full", .(variant = "M3", ks_full = sum(ks_ok, na.rm = TRUE)), by = ds], row.names = FALSE)
  saveRDS(list(sim = SE, real = RE), file.path(out, sprintf("eval_%s%s.rds", phase, tag)))
} else saveRDS(list(sim = SE), file.path(out, sprintf("eval_%s%s.rds", phase, tag)))

# ---- blockade (B'), confirmation only ---------------------------------------------------------------------------------------------------------
bf <- file.path(out, sprintf("blockade%s.rds", tag))
if (phase == "conf" && file.exists(bf)) { B <- readRDS(bf); hv <- B[, .(have = paste(unique(comp), collapse = ",")), by = seg]
  st_b <- merge(unique(B[, .(seg, lr_osc, gain)]), hv, by = "seg"); bst <- cfg$blockade
  bl <- function(nm, v) { s <- if (nm == "M3") copy(st_b)[, chosen := "M3"] else choose(st_b, B, "seg", v)
    if ("chosen_par" %in% names(s)) s[, chosen := chosen_par]; x <- merge(B, s[, .(seg, chosen)], by = "seg")[comp == chosen]; blockade_ratio_summary(blockade_ratios(x), bst)[, variant := nm] }
  BS <- rbindlist(lapply(c("M3", names(VARIANTS)), function(nm) bl(nm, VARIANTS[[nm]])))
  cat("\n== Blockade (ratio to saline) ==\n"); print(BS[, .(variant, breathing, cond, n, SD_p = sprintf("%.2f [%.2f, %.2f]", SD_p, SD_p_lo, SD_p_hi), SD_s = sprintf("%.2f [%.2f, %.2f]", SD_s, SD_s_lo, SD_s_hi))][order(cond, breathing)], row.names = FALSE)
  m3b <- BS[variant == "M3"]; far <- function(r) abs(log(r))
  BB <- BS[, { a <- .SD[cond == "atropine_high"]; p <- .SD[cond == "propranolol"]; ma <- m3b[cond == "atropine_high"]; mp <- m3b[cond == "propranolol"]
    .(Bprime = all(a$SD_p < 0.3) && all(far(p$SD_p[order(p$breathing)]) <= far(mp$SD_p[order(mp$breathing)]) + 1e-9) && all(far(a$SD_s[order(a$breathing)]) <= far(ma$SD_s[order(ma$breathing)]) + 1e-9)) }, by = variant]
  print(BB, row.names = FALSE); saveRDS(BS, file.path(out, sprintf("blockade_summary%s.rds", tag)))
}
