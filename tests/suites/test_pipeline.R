# Unit and smoke tests. Run from the project root:
#   Rscript tests/suites/test_pipeline.R            # unit checks (~1-2 min)
#   SPARK_E2E=1 Rscript tests/suites/test_pipeline.R # + end-to-end run of all scripts on synthetic data (~10 min)
Sys.setenv(SPARK_TEST_ROOT = Sys.getenv("SPARK_TEST_ROOT", "tests/_run"))
Sys.setenv(SPARK_CONFIG = normalizePath("tests/fixtures/config_test.R"))
source("analysis/00_setup.R"); source("pipeline.R")
n_pass <- 0L
check <- function(name, cond) { if (!isTRUE(cond)) stop("FAILED: ", name, call. = FALSE); n_pass <<- n_pass + 1L; log_msg("ok   %s", name) }

# Reproducible synthetic series with artifacts
set.seed(2); n <- 400; tm <- cumsum(rep(0.85, n))
dy <- 0.85 + 0.03 * sin(2 * pi * tm / 9) + 0.02 * sin(2 * pi * tm / 60) + rnorm(n, 0, 0.01); dy[c(100, 250)] <- c(1.6, 0.45)
ctl <- spark_control(); pri <- spark_priors()

# 1. compiled filter == engine filter
for (p in list(c(0.3, 0.006, 2.7, 0.1), c(0.5, 0.02, 1.2, 0.01), c(0.2, 1e-4, 5, 0.3))) {
  a <- run_gls_filter(dy, p[1], p[2], p[3], p[4], ctl); b <- filter_eval(dy, p[1], p[2], p[3], p[4], ctl, full = TRUE)
  check(sprintf("compiled filter equals engine filter (kp = %.1f, ks = %g)", p[1], p[2]),
        max(abs(a$LL_star - b$LL_star) / abs(a$LL_star), max(abs(a$states[, "X_p"] - b$X_p)), max(abs(a$states[, "X_s"] - b$X_s)), max(abs(a$implied_rr - b$implied_rr))) < 1e-9)
}
# 2. recorded reference values (taken before the M3 engine became the engine of record), the engine fit and the pipeline fits
gold <- readRDS("tests/fixtures/golden_m3.rds"); idx <- c(1, 50, 200, 400)
for (g in gold$filters) {
  a <- run_gls_filter(dy, g$par[1], g$par[2], g$par[3], g$par[4], ctl)
  check(sprintf("engine filter reproduces the recorded reference (kp = %.1f, ks = %g)", g$par[1], g$par[2]),
        abs(a$LL_star - g$LL) < 1e-8 * abs(g$LL) && abs(a$nu0 - g$nu0) < 1e-8 && abs(a$sig2 - g$sig2) < 1e-10 &&
        max(abs(a$states[idx, "X_p"] - g$xp), abs(a$states[idx, "X_s"] - g$xs), abs(a$implied_rr[idx] - g$rr)) < 1e-8)
}
ps <- prepare_series(dy, pri, ctl); sp <- model_spec("M3")
fe <- fit_model(dy, pri, ctl, verbose = FALSE); fv <- fit_variant(ps, sp); ff <- fit_model_fast(dy, pri, ctl)
check("engine fit_model reproduces the recorded M3 fit", max(abs(c(fe$params$ks, fe$params$kp, fe$params$nu0, fe$LL_star) / c(gold$fit$ks, gold$fit$kp, gold$fit$nu0, gold$fit$LL) - 1)) < 1e-4)
check("fit_variant(M3) reproduces the engine estimates", max(abs(c(fe$params$ks, fe$params$kp, fe$params$nu0) / c(fv$par$ks, fv$par$kp, fv$nu0) - 1)) < 1e-3)
check("fit_model_fast reproduces the engine estimates", max(abs(c(fe$params$ks, fe$params$kp, fe$params$nu0) / c(ff$params$ks, ff$params$kp, ff$params$nu0) - 1)) < 1e-3)
check("M3 optimum: the objective does not decrease in any coordinate direction", all(vapply(seq_along(fv$theta), function(i) { th <- fv$theta; th[i] <- th[i] - 0.05; th2 <- fv$theta; th2[i] <- th2[i] + 0.05
  variant_objective(th, ps$dy, sp, ps$priors, ps$control) >= fv$objective - 1e-6 && variant_objective(th2, ps$dy, sp, ps$priors, ps$control) >= fv$objective - 1e-6 }, logical(1))))
check("the engine default vagal floor fraction is 0.25 and the spectral anchor matches the record", spark_priors()$lower_lf_vagal_fraction == 0.25 && abs(ps$priors$energy_ratio - gold$energy_ratio) < 1e-8)
check("engine file is byte-identical to the recorded version", unname(tools::md5sum(cfg$paths$engine)) == trimws(readLines("tests/fixtures/engine.md5")[1]))
# 3. nested specifications behave as defined
fs <- fit_specs(ps, c("M0", "M1", "M2", "M3"), full = FALSE)
check("all specifications fit", all(vapply(fs, function(f) isTRUE(f$ok), TRUE)))
check("anchor specs fix lambda_P", all(abs(c(fs$M2$par$lp, fs$M3$par$lp) - ps$priors$energy_ratio) < 1e-12) && abs(fs$M0$par$lp - ps$priors$energy_ratio) > 1e-6)
check("kinetic ordering kp > ks in every spec", all(vapply(fs, function(f) f$par$kp > f$par$ks, TRUE)))
# 4. simulator: stationary moments and scale
x <- simulate_rr(900, kp = 1 / 4, ks = 1 / 150, lp = 2.5, sig2 = 9e-4, seed = 11)
check("simulator OU variances within 30% of target", abs(var(x$x_p_ou) / (2.5 * 9e-4) - 1) < 0.3 && abs(var(x$x_s_ou) / 9e-4 - 1) < 0.6)
check("simulated mean RR near 1/nu0", abs(mean(x$rr) * 1.2 - 1) < 0.08)
# 5. segmentation recovers known events; no event in rest
tl <- sim_tilt(5); e <- attr(tl, "events"); d <- detect_events(tl$rr, "tilt", cfg$segmentation)
check("tilt onset detected within 25 s", d$detected && abs(d$t_onset - e$t_onset) < 25)
ex <- sim_exercise(6); e <- attr(ex, "events"); d <- detect_events(ex$rr, "exercise", cfg$segmentation)
check("exercise profile onset/end within 15 s", abs(d$t_onset_profile - e$t_onset) < 15 && abs(d$t_end_profile - e$t_end) < 15)
# 5b. repeated tilts: angle, heart rate, matching, windows
rep_tilts <- sim_tilt_angle(5, 1500, "repeated"); truth_tilts <- attr(rep_tilts, "events")$tilts
angle_fun <- attr(rep_tilts, "angle"); grid_t <- seq(0, 1500, by = 1)
ang_ev <- tilt_events(grid_t + 10, angle_fun(grid_t), cfg$tilt)
check("tilt_events finds all 3 tilts, onset and return within 3 s", ang_ev$n_events == 3 && all(abs(ang_ev$events$t_onset - truth_tilts$t_onset) < 3) && all(abs(ang_ev$events$t_return - truth_tilts$t_return) < 3))
check("tilt_events keeps the first tilt in t_onset / t_return", ang_ev$t_onset == ang_ev$events$t_onset[1] && ang_ev$t_return == ang_ev$events$t_return[1])
hr_ev <- detect_events(rep_tilts$rr, "tilt", cfg$segmentation)
check("detect_events finds the 3 repeated tilts from heart rate alone (onsets within 30 s)", hr_ev$n_events == 3 && all(abs(hr_ev$events$t_onset - truth_tilts$t_onset) < 30))
check("exercise mode still returns one event", detect_events(sim_exercise(6)$rr, "exercise", cfg$segmentation)$n_events == 1)
dip <- pair_perturbations(c(100, 150, 200, 300), c(1.0, 1.4, 1.3, 1.5, 1.0))
check("pair_perturbations: a dip inside a perturbation does not end it", nrow(dip) == 1 && dip$t_onset == 100 && dip$t_recovery == 300)
two <- pair_perturbations(c(100, 200, 400, 500), c(1.0, 1.4, 1.0, 1.4, 1.0))
unfinished <- pair_perturbations(c(100, 200, 400), c(1.0, 1.4, 1.0, 1.4))
check("pair_perturbations: two separate perturbations and an unfinished last one", nrow(two) == 2 && nrow(unfinished) == 2 && is.na(unfinished$t_recovery[2]))
mt <- match_events(c(100, 300, 500), c(110, 480, 310, 900), 60)
check("match_events: one-to-one, nearest first, none beyond the tolerance", identical(mt$detected, c(1L, 3L, 2L)) && is.na(match_events(c(100), c(300), 60)$detected) && sum(!is.na(match_events(c(100, 105), 110, 60)$detected)) == 1)
win <- tilt_windows(c(300, 700, 880), c(500, 800, 1000), 0, 1200, cfg$tilt)
check("tilt_windows: baseline starts a washout after the previous return; short baselines are not usable", win$t_base_from[2] == 620 && win$baseline_s[2] == 80 && !win$usable[3] && all(win$usable[1:2]))
ph <- assign_tilt_phases(seq(0, 1200, by = 5), win)
check("assign_tilt_phases: no beat belongs to two events and tilt beats are inside the tilts", all(ph[phase == "tilt", event] %in% 1:2) && sum(!is.na(ph$phase)) > 0 && !any(ph[event == 3, !is.na(phase)]))
dep <- first_departure(1:200, c(rep(0, 100), rep(5, 100)) + rep(c(-.1, .1), 100), c(10, 90), t_max = 150)
check("first_departure respects t_max", dep >= 100 && dep <= 101 && is.na(first_departure(1:200, c(rep(0, 100), rep(5, 100)) + rep(c(-.1, .1), 100), c(10, 90), t_max = 50)))
per_event <- data.table(SubjectID = c("a", "a", "a", "b"), event = c(1, 2, 3, 1), d_x_p = c(-1, -3, 2, -2), d_x_s = c(1, 2, 3, -1), tau_p_resp = c(2, 4, 6, 8))
by_recording <- prcp_by_recording(per_event)
check("prcp_by_recording: median over a recording's tilts, expected pattern recomputed", by_recording[SubjectID == "a", d_x_p] == -1 && by_recording[SubjectID == "a", n_events] == 3 && by_recording[SubjectID == "a", expected_p] && !by_recording[SubjectID == "b", expected_s])
check("no event in a long stationary recording", !detect_events(sim_rest(3, "young", 600)$rr, "tilt", cfg$segmentation)$detected)

# 5c. protocol annotations of PRCP (.anI): parsing, protected beats, kinds, events with stand-ups, responses by kind
prep <- new.env(parent = globalenv()); assign("PREP_AUTORUN", FALSE, prep); sys.source("tools/prep-data.R", prep)      # parsing functions only
labels_real <- c("Start slow tilt up", "End slow tilt up", "Start slow tilt down", "End slow tilt down", "Start fast tilt up", "End fast tilt up", "Start fast tilt down ",
                 "End fast tilt down", "Start slow tilt up", "End slow tilt up", "Start slow tilt down ", "End slow tilt back down", "Stand-up from supine", "Back to supine",
                 "Stand-up from supine", "Back to supine", "Start rapid tilt up", "End rapid tilt up", "Start rapid tilt down", "End rapid tilt down")
times_real <- c(300.7, 352.1, 545.5, 596.3, 904.1, 906.4, 1093.5, 1095.3, 1362.8, 1414.4, 1600.5, 1650.4, 2077.6, 2237.1, 2529.9, 2696.9, 2980.2, 2981.7, 3167.4, 3169.2)
mv <- prep$parse_maneuvers("12734", times_real, labels_real)
check("parse_maneuvers: the 20 annotations of a real PRCP record give 2 slow tilts, 2 rapid tilts and 2 stand-ups with their times, none unrecognized",
      nrow(mv) == 6 && identical(mv$kind, c("slow tilt", "rapid tilt", "slow tilt", "stand-up", "stand-up", "rapid tilt")) && length(attr(mv, "unparsed")) == 0 &&
        all(abs(mv$onset_s[c(1, 4, 6)] - c(300.7, 2077.6, 2980.2)) < 1e-9) && abs(mv$return_s[5] - 2696.9) < 1e-9 && abs(mv$return_end_s[1] - 596.3) < 1e-9 && all(is.na(mv$onset_end_s[4:5])))
mv2 <- prep$parse_maneuvers("x", c(10, 20, 30, 40, 50, 60, 70), c("START FAST TILT UP", "end fast tilt up", " Start  fast tilt down", "End fast tilt down", "Stand up from supine", "Back to supine", "Sit down"))
check("parse_maneuvers: case, spacing and fast/rapid wording are tolerated, and an unknown label is reported",
      identical(mv2$kind, c("rapid tilt", "stand-up")) && identical(attr(mv2, "unparsed"), "Sit down"))
check("protected_beats: beats within 5 s before to 30 s after a maneuver boundary",
      identical(prep$protected_beats(c(94, 96, 104, 130, 131, 170, 196, 199), data.table(onset_s = 100, onset_end_s = NA_real_, return_s = 200, return_end_s = NA_real_), 5, 30),
                c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)))
check("label_tilt_kinds: angle-derived onsets take the kind of the annotated tilt that contains them; a stand-up is never matched",
      identical(label_tilt_kinds(c(305, 905, 2200, 5000), mv, 5), c("slow tilt", "rapid tilt", "tilt", "tilt")) && identical(label_tilt_kinds(305, NULL), "tilt"))
two_subjects <- data.table(id = c("a", "a", "b"), kind = c("slow tilt", "stand-up", "slow tilt"), onset_s = c(1, 2, 3), onset_end_s = NA_real_, return_s = c(4, 5, 6), return_end_s = NA_real_, label = "")
check("annotated_maneuvers: only the rows of the requested recording (a column is called id)", nrow(annotated_maneuvers(two_subjects, "a")) == 2 && nrow(annotated_maneuvers(two_subjects, "b")) == 1 && is.null(annotated_maneuvers(two_subjects, "z")))
sim_mixed <- sim_tilt_angle(11, 1500, "mixed")
ann_mixed <- cbind(id = "s1", attr(sim_mixed, "events")$maneuvers, label = "")
states_mixed <- data.table(Time = sim_mixed$t, RR_Interval = sim_mixed$rr, angle = attr(sim_mixed, "angle")(sim_mixed$t))
ev_mixed <- prcp_define_events(states_mixed, "s1", cfg, annotated_maneuvers(ann_mixed, "s1"))
labels_12815 <- c("Initiation of slow tilt up", "Conclusion of slow tilt up", "Initiation of slow tilt down", "Conclusion of slow tilt down", "Blood pressure signal fades", "Switched off Finapres",
                   "Switched Finapres on", "Blood pressure artifacts", "Switched off Finapres", "Switched Finapres on", "Initiation of rapid tilt up", "Conclusion of rapid tilt up",
                   "Initiation of rapid tilt down", "Initiation of rapid tilt down", "Blood pressure signal fades", "Switched off Finapres", "Switched Finapres on", "Initiation of slow tilt up",
                   "Conclusion of slow tilt up", "Initiation of slow tilt down", "Conclusion of slow tilt down", "Initiation of rapid tilt up", "Conclusion of rapid tilt up",
                   "Initiation of rapid tilt down", "Conclusion of rapid tilt down", "Stand up", "Transition back to supine", "Stand up", "Transition back to supin")
mv3 <- prep$parse_maneuvers("12815", seq(100, by = 60, length.out = length(labels_12815)), labels_12815)
check("parse_maneuvers: the wording of the other PRCP records (Initiation of / Conclusion of, Stand up, Transition back to supine; one 'Conclusion' missing) gives 2 slow tilts, 2 rapid tilts and 2 stand-ups",
      identical(mv3$kind, c("slow tilt", "rapid tilt", "slow tilt", "rapid tilt", "stand-up", "stand-up")) && length(attr(mv3, "unparsed")) == 0 && !all(is.na(mv3$return_s[5:6])))
check("parse_maneuvers: Finapres losses become intervals without arterial pressure (fade or switch-off until switch-on), and an artifact label is only counted",
      nrow(attr(mv3, "bp_gaps")) == 3 && attr(mv3, "n_bp_flags") == 1L && all(is.finite(attr(mv3, "bp_gaps")$to_s)) && attr(mv3, "bp_gaps")$from_s[1] == 100 + 4 * 60)
labels_12821 <- c("## time resolution: 250", "Initiate slow tilt up", "Conclude slow tilt up", "Initiate slow tilt down", "Conclude slow tilt down", "Initiate rapid tilt up", "Conclude rapid tilt up",
                  "Initiate rapid tilt down", "Conclude rapid tilt down", "Stand up", "Transition back to supine", "Initiate rapid tilt up", "Conclude rapid tilt up", "Initiate rapid tilt down",
                  "Conclude rapid tilt down", "Stand up", "Transition back to supine", "Initiate slow tilt up", "Initiate slow tilt up", "Initiate slow tilt down", "Conclude slow tilt down")
mv4 <- prep$parse_maneuvers("12821", seq(10, by = 60, length.out = length(labels_12821)), labels_12821)
check("parse_maneuvers: 'Initiate / Conclude', a comment line and a duplicated start (the first one stays without a return) are handled",
      identical(mv4$kind, c("slow tilt", "rapid tilt", "stand-up", "rapid tilt", "stand-up", "slow tilt", "slow tilt")) && is.na(mv4$return_s[6]) && !is.na(mv4$return_s[7]) && length(attr(mv4, "unparsed")) == 0)
ann_extra <- rbind(ann_mixed, data.table(id = "s1", kind = "slow tilt", onset_s = 50, onset_end_s = NA_real_, return_s = NA_real_, return_end_s = NA_real_, label = ""))
ev_noangle <- prcp_define_events(states_mixed[, .(Time, RR_Interval)], "s1", cfg, annotated_maneuvers(ann_extra, "s1"))
check("prcp_define_events: without an angle the annotated tilts are used, and a tilt with no annotated return is dropped", ev_noangle$recording$source == "annotation" && nrow(ev_noangle$events) == 3 && all(is.finite(ev_noangle$events$t_return)))
check("id_seed: distinct seeds for 5000 consecutive indices and for typical recording identifiers, stable, and dependent on the base seed",
      length(unique(id_seed(1:5000, 1))) == 5000 && length(unique(id_seed(c(sprintf("%03d", 1:266), sprintf("f1o%02d", 1:10), sprintf("f1y%02d", 1:10), "12726", "12734", "13960"), 1))) == 266 + 23 &&
        identical(id_seed("a", 1), id_seed("a", 1)) && id_seed("a", 1) != id_seed("a", 2))
duplicated_benchmark <- data.table(seed = c(1, 1), real = c(1, 1), state = "p", end_method = "peak", definition = "frac90")
check("benchmark: repeated (seed, realization) pairs are refused instead of producing a meaningless test-retest correlation",
      inherits(try(assert_unique_realizations(duplicated_benchmark), silent = TRUE), "try-error") && isTRUE(assert_unique_realizations(duplicated_benchmark[1])))
check("fmt_fixed: negative numbers use a true minus sign and keep three significant digits", fmt_fixed(-0.0793) == "\u22120.0793" && fmt_fixed(165.2) == "165" && fmt_fixed(0.000370) == "0.000370")
check("prcp_define_events: slow tilt and rapid tilt from the angle, the stand-up (angle unchanged) from the annotations, in time order",
      identical(ev_mixed$events$kind, c("slow tilt", "stand-up", "rapid tilt")) && identical(ev_mixed$events$source, c("angle", "annotation", "angle")) && ev_mixed$recording$n_stand == 1L && !is.unsorted(ev_mixed$events$t_onset))
check("prcp_define_events: the baseline of the rapid tilt starts a washout after the end of the stand-up",
      ev_mixed$windows$t_base_from[3] >= ev_mixed$events$t_return[2] + cfg$tilt$washout_s - 1e-9)
cfg_off <- cfg; cfg_off$tilt$use_annotations <- FALSE
ev_off <- prcp_define_events(states_mixed, "s1", cfg_off, annotated_maneuvers(ann_mixed, "s1"))
check("prcp_define_events: with use_annotations = FALSE only the two angle-derived tilts remain, unlabelled", nrow(ev_off$events) == 2 && all(ev_off$events$kind == "tilt") && all(ev_off$events$source == "angle"))
fake <- data.table(SubjectID = rep(letters[1:5], 2), kind = rep(c("slow tilt", "rapid tilt"), each = 5), expected_p = TRUE, expected_s = TRUE, z_d_x_p = -2, z_d_x_s = rep(c(8, 9), each = 5) + (1:5) / 10,
                   lat_p = 0, lat_s = rep(c(5, 3), each = 5), tau_p_resp = 3, tau_s_resp = 40)
bk <- prcp_kind_summary(fake); kc <- prcp_kind_contrasts(fake)
check("prcp_kind_summary / prcp_kind_contrasts: counts per kind, and the paired difference rapid - slow of the recording medians",
      identical(bk$kind, c("slow tilt", "rapid tilt")) && all(bk$n_events == 5) && all(bk$n_recordings == 5) && abs(bk$z_d_x_s_med[2] - 9.3) < 1e-9 &&
        kc[contrast == "rapid - slow tilt" & variable == "z_d_x_s", n_pairs] == 5 && abs(kc[contrast == "rapid - slow tilt" & variable == "z_d_x_s", median_diff] - 1) < 1e-9)
# 5c. notation of figures and tables (R/report/symbols.R)
check("tex(): italic variable, upright subscript, italic Greek", tex("$x_p$") == "<i>x</i><sub>p</sub>" && tex("$\\tau_s$") == "<i>\u03c4</i><sub>s</sub>")
check("tex(): index subscripts stay italic, labels upright, hat inside the italic run", tex("$t_{k-1}$") == "<i>t</i><sub><i>k</i>\u22121</sub>" && grepl("<i>\u03c3\u0302</i>", tex("$\\hat\\sigma^2_{sys}$"), fixed = TRUE))
check("tex(): plain text is escaped (no stray markdown), newline and **bold** work", tex("a_b *c* <d>") == "a&#95;b &#42;c&#42; &lt;d&gt;" && tex("a\nb") == "a<br>b" && tex("**a**") == "<b>a</b>")
check("tex(): an unsupported command is an error, not a silent wrong label", inherits(try(tex("$\\notacommand$"), silent = TRUE), "try-error"))
check("sym(): every variable of the group contrasts and the HRV comparators has a label", all(c(FANTASIA_LATENT_VARIABLES, FANTASIA_HRV_VARIABLES) %in% names(SYMBOL_TEX)))
check("sym_tex(): unknown names stay readable; spec and setting labels are LaTeX", sym_tex("some_new_variable") == "some new variable" && spec_tex("M3") == "$M_3$" && grepl("\\tau_p", setting_axis_label("tau_p_bounds=narrow"), fixed = TRUE))
check("table labels contain no pipe (it would split a Quarto pipe table)", !any(grepl("|", c(unname(SYMBOL_TEX), unname(SCENARIO_LABELS), setting_axis_label(c("crossover_hz=0.08", "M3 without tau"))), fixed = TRUE)))
check("figure_font() returns an installed family or the default", { f <- figure_font(cfg); identical(f, "") || f %in% systemfonts::system_fonts()$family })
# 6. recovery definitions on a plateau followed by an exponential decay with known tau
set.seed(4); t0 <- seq(0.5, 600, by = 0.5)
x0 <- ifelse(t0 < 100, 0, ifelse(t0 < 200, 0.5, 0.5 * exp(-(t0 - 200) / 40))) + rnorm(length(t0), 0, 0.01)
check("recovery_metrics returns NULL without a rest window", is.null(recovery_metrics(t0, x0, c(1000, 1100), 200, cfg$exercise)))
rm <- recovery_metrics(t0, x0, c(10, 90), 200, cfg$exercise)
check("exponential-constant recovery recovers tau = 40 s", abs(rm[definition == "expfit", time] - 40) < 6)
check("90% recovery time near tau*ln(10) = 92 s", abs(rm[definition == "frac90", time] - 92) < 10)
check("a state without excursion is flagged", recovery_metrics(t0, rnorm(length(t0), 0, 0.01), c(10, 90), 200, cfg$exercise)[1, status] == "no_excursion")
# 7. HRV
sine <- 0.8 + 0.020 * sin(2 * pi * 0.1 * cumsum(rep(0.8, 750)))
check("LF band power of a 20-ms sine is ~200 ms^2", abs(band_power(sine, c(0.04, 0.15)) / 200 - 1) < 0.05)
# 8. IO: CSV loader, unit guard, join artifacts, covariates
tmp <- tempfile(); dir.create(tmp)
f <- file.path(tmp, "ms.csv"); fwrite(data.table(time = cumsum(rep(0.8, 100)), rri = 1000 * rep(0.8, 100)), f)
check("RR in ms is rejected", inherits(try(read_recording(f, cfg), silent = TRUE), "try-error"))
check("placeholder paths are rejected", inherits(try(check_data_path("data/<X>", "prcp", cfg), silent = TRUE), "try-error"))
check("Fantasia groups from names", identical(as.character(fantasia_group(c("f1y03", "f2o10"), cfg)), c("young", "older")))
sp <- c(rep(0.8, 60), 3.5, rep(0.8, 60)); y <- clean_sample(sp)
check("clean_sample replaces an isolated outlier and keeps the rest", abs(y[61] - 0.8) < 0.05 && max(abs(y[-61] - sp[-61])) < 1e-9)
set.seed(8); sim <- sim_rest(21, "young", 300, m_p = rsa_drive(0.25, 0.04))
dtmp <- file.path(tmp, "d"); dir.create(dtmp)
write_fantasia_csv(sim, "fx", dtmp, has_bp = TRUE, has_resp = FALSE, seed = 5, legacy = TRUE); write_tmst_csv(sim, "tx", dtmp)
write_fantasia_csv(sim, "fn", dtmp, has_bp = TRUE, has_resp = FALSE, seed = 5)
fx <- read_recording(file.path(dtmp, "fx.csv"), cfg, "auto"); tx <- read_recording(file.path(dtmp, "tx.csv"), cfg, "auto")
check("legacy rolling-join file: auto picks dtime (rri not aligned with time)", fx$qc$rr_source == "dtime" && fx$qc$time_rri_consistency < 0.5)
check("TMST file: auto keeps rri (diff(time) == rri)", tx$qc$rr_source == "rri" && tx$qc$time_rri_consistency > 0.99)
check("dtime reconstruction equals the true RR except where clean_sample() overwrites beats (250 Hz quantization, up to ~15% of beats)", { k <- seq_len(min(length(fx$rr), length(sim$rr)) - 5); median(abs(fx$rr[k] - sim$rr[k])) < 0.003 && mean(abs(fx$rr[k] - sim$rr[k]) < 0.010) > 0.85 && cor(fx$rr[k], sim$rr[k]) > 0.95 })
fn <- read_recording(file.path(dtmp, "fn.csv"), cfg, "auto")
check("current-format file: auto keeps rri (aligned with time) and rri equals the true RR (median error < 3 ms)", fn$qc$rr_source == "rri" && fn$qc$time_rri_consistency > 0.95 && median(abs(fn$rr[seq_len(300)] - sim$rr[seq_len(300)])) < 0.003)
check("covariates: BP kept, absent respiration dropped", all(c("systolic", "diastolic", "map") %in% names(fx$cov)) && !"resp" %in% names(fx$cov) && fx$qc$has_bp && !fx$qc$has_resp)
check("TMST file has no covariates", ncol(tx$cov) == 0)
w <- window_rec(fx, start_s = 60, duration_s = 100)
check("window_rec cuts rr, time and covariates consistently", length(w$rr) == length(w$time) && nrow(w$cov) == length(w$rr) && all(w$time > 60 & w$time <= 160))
# 9. covariate analyses against known truth
set.seed(3); tt <- seq(0, 540, by = 0.85); ang <- ifelse(tt < 120, 0, ifelse(tt < 180, 30, ifelse(tt < 240, 60, ifelse(tt < 300, 90, 0)))) + rnorm(length(tt), 0, 0.8)
ev <- tilt_events(tt + 500, ang, cfg$tilt)
check("tilt_events recovers onset (120 s) and return (300 s)", ev$detected && abs(ev$t_onset - 120) < 3 && abs(ev$t_return - 300) < 3)
lv <- angle_levels(tt + 500, ang, cfg$tilt)
check("angle_levels finds the 4 plateau levels", all(c(0, 30, 60, 90) %in% lv$level))
nn <- length(tt); fr <- 0.27; resp <- sin(2 * pi * fr * tt) + rnorm(nn, 0, 0.2)
check("respiratory frequency recovered (0.27 Hz)", abs(resp_frequency(tt, resp)[["f_resp"]] - fr) < 0.01)
xcp <- 0.5 * sin(2 * pi * fr * tt + 0.6) + rnorm(nn, 0, 0.1); xin <- rnorm(nn, 0, 0.1) + 0.5 * sin(2 * pi * 0.05 * tt)
cp <- resp_coupling(tt, xcp, resp, fr); ci <- resp_coupling(tt, xin, resp, fr)
check("respiratory-band SD: coupled ~0.35, uncoupled much smaller", abs(cp[["band_sd"]] - 0.35) < 0.08 && ci[["band_sd"]] < 0.5 * cp[["band_sd"]])
check("respiratory coherence higher for the coupled series", cp[["coh"]] > 0.8 && ci[["coh"]] < 0.5)
t2 <- seq(0, 900, by = 0.85); s0 <- function(u) 120 + 4*sin(2*pi*0.10*u) + 3*sin(2*pi*0.07*u + 1) + 3*sin(2*pi*0.12*u + 2) + 2*sin(2*pi*0.05*u + 3)
set.seed(9); sbp <- s0(t2) + rnorm(length(t2), 0, 0.3); rr2 <- 0.85 + 8e-3 * (s0(t2 - 2) - 120) + rnorm(length(t2), 0, 2e-3)
b <- brs_alpha(t2, sbp, rr2, cfg$physio$brs_bands, cfg$physio$brs_min_coh)
check("baroreflex alpha (LF) recovers the true gain of 8 ms/mmHg within 15%", abs(b[["alpha_LF"]] / 8 - 1) < 0.15)
xx <- as.numeric(stats::filter(rnorm(length(t2)), 0.9, "recursive")); yy <- c(rep(0, 6), head(xx, -6)) + rnorm(length(t2), 0, 0.3)
xc <- xcorr_surrogate(t2, xx, yy, 30, n_surr = 50)
check("lagged correlation recovers a 5-s lag with a small surrogate p-value", abs(xc$lag_s - 5) <= 1 && xc$p < 0.05)
# 10. synthetic data writers produce loadable files in all three datasets (isolated temp folder)
tr <- write_synthetic_data(file.path(tempfile(), "data"), n_prcp = 4, n_ex = 2, n_fant = 3, seed = 3, prcp_duration = 400, ex_duration = 500, fant_duration = 300)
cfgt <- cfg; for (ds in cfg$datasets) cfgt$paths[[ds]] <- tr$dirs[[ds]]
dl <- lapply(cfg$datasets, load_dataset, cfg = cfgt); names(dl) <- cfg$datasets
check("synthetic CSV files load through load_dataset()", all(vapply(dl, function(d) nrow(d$index) > 0, TRUE)))
check("missing covariates are detected per recording (PRCP: one without BP, one without angle)", sum(!dl$prcp$index$win_systolic) == 1 && sum(!dl$prcp$index$win_angle) == 1)
check("Fantasia: one recording without BP and one without respiration", sum(!dl$fantasia$index$win_systolic) == 1 && sum(!dl$fantasia$index$win_resp) == 1)
check("rr_source: files in the current format are read as rri in all three datasets", all(unlist(lapply(dl, function(d) d$index$rr_source == "rri"))))

# 11. State-space landscape, flow field and example selection
set.seed(1)
cluster_points <- data.table(u = c(rnorm(600, -2, 0.6), rnorm(400, 2, 0.6)), v = c(rnorm(600, 1, 0.6), rnorm(400, -1, 0.6)), weight = 1)
limits <- landscape_limits(cluster_points$u, cluster_points$v)
two_basins <- density_landscape(cluster_points, limits, n_grid = 80)
check("landscape finds both basins, at (-2, 1) and (2, -1)", nrow(two_basins$basins) >= 2 &&
        any(abs(two_basins$basins$u + 2) < 0.4 & abs(two_basins$basins$v - 1) < 0.4) && any(abs(two_basins$basins$u - 2) < 0.4 & abs(two_basins$basins$v + 1) < 0.4))
cluster_points[, weight := c(rep(0.2 / 600, 600), rep(1 / 400, 400))]
check("point weights decide which basin is strongest", {
  b <- density_landscape(cluster_points, limits, n_grid = 80)$basins[1]; abs(b$u - 2) < 0.4 && abs(b$v + 1) < 0.4 })
check("the potential is 0 at the strongest basin and capped at potential_max", min(two_basins$grid$potential) == 0 && max(two_basins$grid$potential) <= 6)
set.seed(2)
ou_paths <- rbindlist(lapply(1:6, function(i) { x <- matrix(NA, 400, 2); x[1, ] <- 3 * runif(2, -1, 1)
  for (k in 2:400) x[k, ] <- x[k - 1, ] + 0.1 * (c(1, -1) - x[k - 1, ]) + rnorm(2, 0, 0.15); data.table(id = i, t = 1:400, u = x[, 1], v = x[, 2]) }))
flow <- flow_field(ou_paths, list(u = c(-3, 3), v = c(-3, 3)), n_bins = 10, min_n = 15)
check("flow arrows point towards the attractor at (1, -1)", mean((sign(flow$du) == sign(1 - flow$u) | abs(1 - flow$u) < 0.4) & (sign(flow$dv) == sign(-1 - flow$v) | abs(-1 - flow$v) < 0.4)) > 0.9)
scores <- c(a = 3, b = 9, c = 5, d = 7)
check("choose_example: largest, median and a given id", choose_example(scores, "largest")$id == "b" && choose_example(scores, "median")$id == "c" && choose_example(scores, "d")$id == "d")
check("choose_example: an unknown id warns and falls back to the largest", {
  r <- withCallingHandlers(choose_example(scores, "zzz"), warning = function(w) invokeRestart("muffleWarning")); r$id == "b" })
check("figure settings are in config.R", all(c("format", "dpi", "base_size", "example", "palette", "landscape") %in% names(cfg$figures)))

if (nzchar(Sys.getenv("SPARK_E2E"))) {
  log_msg("End-to-end run of all scripts on synthetic data...")
  system2("Rscript", "tests/fixtures/make_synthetic_data.R")
  check("run_all.R completes", system2("Rscript", "tools/run_all.R") == 0)
  need <- unlist(lapply(PIPELINE, `[[`, "writes"))                       # pipeline.R
  check("all result files exist", all(file.exists(file.path(cfg$paths$results, paste0(need, ".rds")))))
  check("all figures exist", all(paste0(names(FIGURE_REGISTRY), ".png") %in% list.files(cfg$paths$figures)))
}
log_msg("%d checks passed.", n_pass)
