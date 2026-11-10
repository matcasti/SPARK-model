# Candidate: two OU drives in ln(RR), exact per-beat filter, blockade layer

Experimental, 2026-11-09. Spec: *Hoja de implementación OU bivariado en log-RR con capa de bloqueo* (Matías Castillo-Aguilar). Nothing here touches the production engine, the pipeline, the manuscript or the app (`v1.0.0` is unchanged).

Done: roadmap steps 1 (core, unit tests) and 2 (recovery on synthetic data). Not done: step 0 and steps 3-8 (see *Pending*).

## Files

| File | Content |
|---|---|
| `lib.R` | `ou_transition`, `obs_noise`, `kf_run` (double scalar Kalman filter, level `c` profiled), `profile_level`, `rts_smooth`, `ou_priors`, `neg_log_post`, `fit_candidate`, `innovations`. `h = (1 - b_p, -(1 - b_s))` is already an argument, so the blockade layer needs no new filter |
| `test_core.R` | unit tests of step 1: `Rscript experiments/ou_logrr/test_core.R` |
| `sim.R` | generator G1, scenarios E1/E2, recovery: `Rscript experiments/ou_logrr/sim.R E2 60 variant=wide_oracle` |
| `_out/` | `recovery_*.rds` and `summary_*.txt` (git-ignored) |

Pure R: a fit takes 0.4 s (10 min of RR) to 2 s (exercise, 940 beats), so no C++ is needed.

## The spec against this repository's engine (rows marked "verificar")

The sheet's "motor actual" is not `engine/model_engine.R`. Differences found:

- **Prior on tau_s:** the engine default is log-normal 100-300 s (`spark_priors`, not overridden in `config.R`); the sheet says 60-100 s. The reference value for record 133 (tau_S near 68.6 s) cannot come from the repo defaults (not checked by running).
- **Noise model:** the sheet's M1 (scale depending on X_p, gamma, 15-s low-pass, Lomax prior) is not in this repository. Rate-dependent noise was tried in `experiments/m3plus` (gamma negative on real data, not supported).
- **Optimizer:** the engine runs BFGS from three starting points (since 2026-11-04); the sheet says one.
- **Matches:** phase-domain terminal-threshold observation, lambda_P = V_p/V_s fixed from the Morlet CWT, sigma^2 and nu0 profiled by GLS, MAP with BFGS.
- `01_model_engine_v2.R` does not exist in this project nor in `2026-08 RRi-IC` or `2026-09 STATE-theory`. Records 133 (1192 beats) and 008 (1881 beats) exist in `data/tmst`.

## Step 1: core (14 checks, all pass)

Transition of d1 + d2 against two steps (1e-19); filter likelihood against the dense multivariate normal, N = 200 (1e-12); `c = A/B` against numerical optimisation (2e-9); covariance positive definite at every beat; smoother variance not above the filtered one and last point identical; skipping a beat's update against merging its interval (1e-17); smoother mean and variance against the dense Gaussian posterior (1e-15, an extra check not in the spec).

## Step 2: recovery with known truth (G1, sigma_t = 1 ms, tau_p 4.5 s, tau_s 75 s, V_p = V_s = 0.00125, c = ln 0.85)

Median [P10, P90] of estimate / truth. The thresholds of the spec are proposals ("recalibrate after the first run").

**E1, stationary.**

| | 10 min (100 records) | 30 min (60 records) | Proposed band |
|---|---|---|---|
| tau_p | 0.96 [0.84, 1.10] | 0.97 [0.83, 1.09] | 0.67-1.5 |
| tau_s | 1.03 [1.01, 1.04] | 1.02 [0.99, 1.06] | 0.67-1.5 |
| V_p | 0.99 [0.84, 1.17] | 0.97 [0.87, 1.13] | 0.75-1.33 |
| V_s | **0.73** [0.19, 1.62] | 0.86 [0.48, 1.29] | 0.75-1.33 |
| 95% coverage x_p / x_s | 0.92 / **0.83** | 0.93 / 0.90 | 0.90-0.98 |

Innovations pass (mean 0.01, variance 1.00, Ljung-Box p > 0.05 in 96%). V_s is poorly determined at 10 min because tau_s = 75 s gives about four slow cycles: an information limit that shrinks with record length, not a bug (the core is exact).

**E2, exercise (300 s rest, 120 s exercise with a fast vagal step to -0.25 and a linear sympathetic ramp to +0.30, 300 s recovery), 60 records per variant.** Change from rest to the end of exercise, estimated / true (true: x_p -0.25, x_s +0.29):

| Variant | x_p change | x_s change | V_p | V_s | cover x_p / x_s |
|---|---|---|---|---|---|
| **Spec as written** (tau priors strict, V free) | **0.03** | **1.84** | 0.36 | 11.7 | 0.57 / 0.10 |
| tau_p prior wide (3-60 s), V free | unstable | unstable | 13 | 0.015 | 0.05 / 0.09 |
| strict tau, V_p/V_s fixed at truth | 0.74 | 1.22 | 2.4 | 2.4 | 0.69 / 0.53 |
| **wide tau_p, V_p/V_s fixed at truth** | 1.09 | 0.92 | 6.8 | 6.8 | 0.94 / 1.00 |
| wide tau_p, ratio x 3 too high | 1.61 | 0.45 | 10.5 | 3.5 | 0.97 / 0.86 |
| wide tau_p, ratio / 3 too low | 0.53 | 1.41 | 3.2 | 9.7 | 0.77 / 0.87 |

## What step 2 says

1. **The spec as written fails E2 (the criterion to advance).** The vagal plateau is lost (3% of the true change) and the sympathetic state absorbs the whole net change (1.84x); V_s is inflated 12x and V_p drops to 0.36x, i.e. the collapse of P that the sheet lists as risk 1, here in simulation with truth. Mechanism: a zero-mean OU with tau_p of 3-6 s returns to zero within seconds and cannot hold a plateau; the sum of the two changes is fixed by ln RR, the split is not. The "correct shape" of step 3 (P falls and stays on a plateau) is out of reach under that prior.
2. **Only a wide tau_p together with an accurate ratio recovers the split** (last-but-two row). Widening tau_p alone swaps the branches (non-identifiability, as in the manuscript); fixing the ratio alone leaves the plateau biased.
3. **The split then follows the ratio:** an error of 3x in V_p/V_s moves the vagal change to 1.6x or 0.5x of the truth. The production anchor misses lambda_P by a median factor of about 3 in simulation (`project-anchor-experiments`), so the same dependence on the anchor is expected on real data.
4. The "efector" variant (tau_p 0.5-2 s) was not run on E2; it should make the plateau worse, not better.

## Pending / decisions

- **Which engine is the baseline of step 0?** The sheet's "motor simple actual" is not in this repository. Options: run `engine/model_engine.R` (M3) on 133 and 008 as the baseline, or point to the file the sheet means.
- Step 3 (real records 133 and 008, `tau_prior = "spark"`, `anchor_ratio = TRUE`) needs the CWT anchor from `extract_spectral_priors`; given 1-3 above it would likely show the same plateau problem, so the spec may need a wider tau_p prior first. Not run.
- Not implemented: `diagnose()` (three criteria), `fit_blockade()` and the nested blockade models (steps 6 and 8; the Lancaster "Blocker data" would feed step 8), generator G2 and scenarios E3/E4 (steps 6-7).
- Caveats: one generator (G1) and one parameter set; 60-100 records per cell; the E2 offsets (time constants 5, 20, 50 s) are my choice of the spec's "fast step and slow ramp".

---

# Comparison against production M3 (2026-11-09, decision: does it become the default?)

Scripts in `compare/`: `lib.R` (variants, units, scoring), `sim_compare.R` (known truth), `real_compare.R` (313 real recordings). Outputs in `_out/compare_*.rds`.

## Variants compared

| | Candidate | Priors on tau | V_p / V_s |
|---|---|---|---|
| **M3** | production (`fit_variant(model_spec("M3"))` / `fit_model_fast`) | tau_p 3-6, tau_s 100-300 s | production CWT anchor |
| F1 | spec + anchor (roadmap step 3) | 3-6 / 60-100 s | anchor |
| F2 | spec, V free (step 4) | 3-6 / 60-100 s | free |
| F3 | "effector" (step 5) | 0.5-2 / 10-20 s | anchor |
| F4 | wide tau_p (my addition, from the exercise recovery) | 3-60 / 60-100 s | anchor |
| F5 | M3 priors with the candidate's observation model (isolates the observation layer) | 3-6 / 100-300 s | anchor |

Units: production states are in Hz, the candidate's in ln RR; the candidate is converted with x_Hz = exp(-c) x (first-order), and M3 to ln RR with 1/nu0 for the candidate's own generator. Scores are unit-free (RMSE / SD of truth, correlation, 95% coverage). Variances of the candidate include the uncertainty of the profiled level (as M3 does for nu0); a first run without it under-covered, which was an omission of mine, not of the spec.

## Decision rule (written before the scenarios other than `stationary` and the real data were looked at)

The burden of proof is on the candidate: M3 is the default unless a variant is non-inferior to M3 on everything below and clearly better on something that matters.

Simulation, filtered states (what production reports), every scenario of both generators:
- S1 point accuracy: mean standardised RMSE not above M3's by more than 10% (relative) and correlation not lower by more than 0.03, in both states.
- S2 calibration: 95% coverage within 0.90-0.98 or not lower than M3's by more than 0.03.
- S3 exercise and tilt: change of each state from rest to the end of the event, estimated / true, within 0.67-1.5; vagal reaches half of its change before the sympathetic one as often as in the truth.
- S4 time constants: median |ln error| of tau_p, tau_s not above M3's by more than 0.1.

Real data (313 recordings):
- R1 admissibility: convergence 100%; share of tau_p and tau_s inside the production prior ranges not lower than M3's by more than 10 points.
- R2 diagnostics: innovation KS pass (p > 0.05) and PACF violations (lags 1-20) not worse than M3's.
- R3 structure: exercise, vagal faster than sympathetic in at least as many recordings as M3 minus 5 points; PRCP expected signs (x_p down, x_s up) not lower than M3's minus 5 points.
- R4 reliability and known groups (Fantasia): test-retest ICC of SD_p, SD_s, tau_p, not lower than M3's by more than 0.05; young vs older AUC of SD_p not lower by more than 0.03.
- R5 agreement: correlation of the filtered states with M3's (descriptive: a variant that disagrees with M3 must be better on S1-S4, not only different).

Added value that could justify a variant: calibrated absolute units for blockade studies (not available in M3). That is a capability, assessed separately: a variant can be the module for blockade data without being the default.

## Results

Simulation (6 scenario groups; 30-120 series each; filtered states; `verdict.R` applies the rule). Cells that fail, by criterion (a group is one of: stationary x4, rsa, exercise, tilt, G1_E1, G1_E2; S3 applies to exercise, tilt, G1_E2):

| | S1 accuracy | S2 calibration | S3 exercise/tilt | S4 time constants | Total |
|---|---|---|---|---|---|
| F1 spec + anchor | 0 | 5 | 0 | 4 | 9 |
| F2 spec, V free | 3 | 6 | 3 | 4 | 16 |
| F3 effector | 5 | 6 | 3 | 4 | 18 |
| F4 wide tau_p | 1 | 1 | 1 | 5 | 8 |
| **F5 M3 priors, candidate observation** | **0** | 5 | **0** | 1 | **6** |

Real data (313 recordings; Fantasia reliability on 40). `ok` = rule met:

| | F1 | F2 | F3 | F4 | F5 | M3 |
|---|---|---|---|---|---|---|
| R1 converged | 100% | 100% | 100% | 100% | 100% | 100% |
| R1 tau_s inside the production prior | 17% (own prior 60-100 s) | 19% | 1% | 7% | 93% | 99% |
| R2 innovation KS pass / PACF violations | 21% / 9.7 ok | 28% / 7.5 ok | 23% / 9.6 ok | 26% / 8.2 ok | 20% / 9.6 ok | 22% / 9.3 |
| R3 exercise: vagal reaches half first | 100% ok | 99% ok | **47%** | **69%** | 98% ok | 99% |
| R3 PRCP: both expected signs | 98% ok | 98% ok | 100% ok | 98% ok | 98% ok | 98% |
| R4 ICC SD_p / SD_s (windows A, B) | .95 / .92 ok | **.44 / .01** | .95 / .92 ok | .94 / .89 ok | .95 / .91 ok | .95 / .90 |
| R4 ICC tau_p | **.75** | **.81** | **.73** | **.74** | **.75** | .88 |
| R4 AUC young vs older, SD_p | .97 ok | .97 ok | .97 ok | .94 ok | .96 ok | .93 |
| R5 correlation with M3 (x_p / x_s) | .98 / .99 | .63 / .94 | .91 / .94 | .93 / .96 | .99 / 1.00 | — |

Exercise cohort, change from rest to the end of exercise in rest SDs (median; x_p / x_s) and time to half (s): M3 -8.4 / 37.7, 13 / 57; F1 -6.2 / 28.4, 9 / 47; F5 -5.8 / 31.7, 9 / 51; F2 -1.1 / 16.2 (vagal response lost); F3 -11.1 / 15.2, vagal 19 s after 26 s (order reversed in half); F4 -10.1 / 18.8, 17 / 33.

## Verdict: not the default; M3 stays

No variant meets the rule. The ones that reproduce M3's point estimates (F1, F5; F4 in simulation) are not better than M3 on anything measured and are worse on calibration and on the reliability of tau_p:

- **The spec as written does not pass.** F1 (anchored) matches M3's states (correlation 0.98-0.99, S1 and S3 pass) but under-covers (S2 fails in 5 of 6 groups: 95% bands cover 0.83-0.89 against M3's 0.87-0.97), misses the time constants (tau_s prior 60-100 s against the 120-250 s of the production generator) and loses reliability of tau_p (ICC 0.75 against 0.88). F2 (free V, the spec's step 4) collapses the vagal response (-1.1 rest SDs against -8.4) and is unreliable (ICC of SD_s 0.01).
- **Variations on the spec do not repair it.** F3 (effector) reverses the order of the responses in half of the recordings; F4 (wide tau_p) cuts the vagal-first share to 69% and fails accuracy in exercise.
- **F5 isolates the observation layer:** with M3's priors and anchor the candidate's ln RR observation gives the same states as M3 (correlation 0.99-1.00, accuracy and exercise/tilt criteria pass) with no gain, a lower coverage and lower tau_p reliability (0.75 against 0.88). The log-RR observation is therefore not what improves M3; the remaining differences come from fixing the noise (sigma_t = 1 ms instead of M3's estimated jitter) and from the priors.
- **What the candidate adds is a capability, not accuracy:** a level `c` in ln RR units and the `h = (1 - b_p, -(1 - b_s))` hook for blockade data (steps 6-8, not done). F5 is the natural starting point for that module (it does not change what production reports); it would need steps 6-8 on the Lancaster "Blocker data" (or equivalent) before any claim.

Caveats: the thresholds are mine (S3 amended after seeing that M3 itself misses the band in tilt); two simulation generators and one parameter draw each; Fantasia reliability is 40 recordings; 313 real recordings come from three datasets with the production windows. **Not tested:** a variant with the timing noise sigma_t estimated under M3's jitter prior (the sheet fixes it on purpose, but it is the most likely cause of the differences in calibration, KS on PRCP and tau_p reliability), an oscillator, blockade steps 6-8, and the sheet's generator G2 / scenarios E3-E4.
