# M4e: M4c-w with three changes, development and confirmation

Experimental, started 2026-11-17. It follows `../README.md`: no prototype passed, and the closest was the post hoc M4c-w. The goal is a candidate that can replace M3 under stationary and non-stationary conditions, for both the parameters and the latent states. Nothing here touches the engine, the pipeline, the manuscript or the app.

## Pre-registration (written before any M4e run)

**Candidate M4e** = M4c-w with three changes. M4c-w is the multi-start base, the production wavelet anchor and the respiratory oscillator, plus the two-regime HMM with reciprocity prior behind a change-point gate.

- **C1, gated oscillator.** The oscillator is kept only if LR = 2 (cost without - cost with oscillator) >= q_osc. q_osc is the 95th percentile of LR in the development simulations without an oscillator (rest), with a floor of 11.3, the 0.99 quantile of chi-squared with 3 df (frequency search plus two parameters).
- **C2, regime gate tuned on rest.** The gate threshold G is the smallest value >= 9 that opens in at most 2% of the development rest and RSA simulations.
- **C3, recalibrated intervals.** When the gate is open, the posterior SDs of the two states are multiplied by factors c_p and c_s. Each factor is the smallest value on the grid 1.0-2.5 that brings the pooled 95% coverage of the development regime, exercise and two-regime simulations to >= 0.95, estimated separately for filtered and smoothed states. A closed gate is not rescaled.

q_osc, G, c_p and c_s are estimated **only** on the development data and then frozen.

**Data split.**

| | Development (tuning and sub-experiments) | Confirmation (run once per frozen candidate) |
|---|---|---|
| Simulation | seeds `m4e-dev`, 40 per scenario | seeds `m4e-conf`, 40 per scenario, never used for tuning |
| Exercise | the 60 recordings of `../real.R` | a disjoint random sample of 60 of the remaining 203 |
| Fantasia | full, head, windows A (600-1200 s) and B (1200-1800 s) | full, head, windows C (3600-4200 s) and D (4200-4800 s) |
| PRCP | the 10 recordings, for development only (too few for confirmation) | |
| Blockade | none | the 81 segments |

**Decision rule.** Defined in `../README.md` (S1-S3, R1-R4, B), with two amendments fixed now, before any M4e result, because the original versions could not be passed by M3 itself or were decided by noise:

- **S1'.** Mean RMSE of both states <= 1.02 x M3 in the stationary scenarios (rest, RSA), and <= 0.90 x M3 in the non-stationary scenarios (regime, exercise, two). The 2% is a non-inferiority margin; the original rule failed M4c-w on a difference of 0.4%.
- **B'.** Atropine SD_p ratio < 0.3 in both breathing modes, and the propranolol SD_p ratio and the atropine SD_s ratio no farther from 1 than M3's in either mode. M3 itself fails the original propranolol condition (0.54 with paced breathing).

S2, S3, R1 (pooled), R2, R3 and R4 are unchanged. A candidate is **strong** if it passes every criterion on the confirmation data.

**Sub-experiments.** If the frozen M4e fails on development, the failing criterion is diagnosed on development data only. A change is then specified, registered in the log below before it is run, and tested on development. The confirmation data are used only for candidates frozen in the log.

## Log

(entries are added in order; each says what was changed, why, and whether it was decided before seeing the result)

**E1 (2026-11-17 00:47), M4e as pre-registered, development simulations.**

- Tuned values:
  - q_osc = 11.3, the chi-squared floor. The oscillator is kept in 0% of rest and 100% of RSA simulations.
  - G = 9. The gate opens in 0% of rest/RSA, 98-100% of regime/exercise and 80% of two-regime simulations.
  - c_p = c_s = 2.5, the top of the grid. Even at 2.5 the pooled 95% coverage is short of 0.95.
- Outcome: passes S1'-S3 in rest, RSA, regime and two-regime. Fails exercise: vagal RMSE 0.947 x M3, against the 0.90 required, and tau_s error 0.48 against M3 + 0.10 = 0.44. Exercise coverage after rescaling is 0.84 / 0.83.
- Diagnosis: in the exercise generator (`draw_exercise_pars`) the set-point responses have their own kinetics: vagal on 6-12 s and off 12-60 s, sympathetic on 15-30 s and off 25-100 s. The regime HMM forces the response to a regime change to follow the OU rates (tau_s about 170 s), so the sympathetic response is far too slow. This misfit is also why the intervals cannot be calibrated by a constant.

**E2 (registered before implementation), regime response kinetics.**

- The state of each drive is u_b + m_b. u_b is the OU deviation, with the same rates and noise as before. m_b is a deterministic set-point response that relaxes towards mu_b(s) at a rate that depends on the regime being entered: r_on,b towards the active regime and r_off,b back to the reference. With r = kappa_b this is exactly the current regime model.
- Priors: ln r ~ N(ln kappa_b, 1.5^2), with bounds 1/600 to 2 s^-1.
- Starts: r = kappa, and physiological values (1/8, 1/30, 1/20, 1/60 s^-1).
- Filter: GPB1 with regime-specific transition matrices (`reg4_filter.cpp`). The smoother is RTS on the collapsed sequence, with the prior-weighted transition matrix (an approximation).
- Gates: q_osc and G unchanged. c_p and c_s are re-tuned on development by the same rule.
- Filter check: with rate = kappa it reproduces the old filter exactly (same log-likelihood and filtered states). In one exercise simulation it recovered the sympathetic on/off constants (16.6 / 26.7 s against the true 19.3 / 27.9 s).

**E2 result (development simulations).** Not better on average. In exercise the vagal RMSE was 0.972 x M3 and the sympathetic 0.790, against 0.947 and 0.668 for E1; tau_s improved to 0.17. Component by component:

- In exercise the regime models are better in the median (vagal RO 0.84 x M3), but a few catastrophic fits (up to 4.7 x M3) dominate the mean.
- RB, the regimes without the oscillator, is better than RO in exercise: 0.84 / 0.61 against 0.96 / 0.67.

The problem is therefore selection and robustness between components, not structure.

**E1 result (development real recordings).**

- M4e passes R1, R2 and R3, and fails R4: Fantasia ICC of tau_p 0.66 against 0.88.
- Diagnosis on windows A/B:
  - With the oscillator, the within-subject SD of ln tau_p is almost unchanged (0.16 -> 0.18). The between-subject SD falls (0.45 -> 0.29).
  - The age-group separation is kept: R^2 of ln tau_p on young/old 0.26 with the oscillator, 0.25 for M3.
  - Within age groups, the extra between-subject spread of M3's tau_p correlates with the strength of the respiratory oscillation (r = -0.39 with the oscillator LR). It does not correlate with the respiratory frequency (r = 0.03).
  - So the ICC of M3's tau_p rewards the stability of the respiratory amplitude that tau_p absorbs. It does not measure the precision of tau_p.

**Amendment R4' (post hoc, from development data only; confirmation windows C/D are unseen).** Two conditions, both reported alongside the original R4:

- (i) the within-subject SD of ln tau_p, ln SD_p and ln SD_s between the two windows is <= 1.10 x M3's;
- (ii) the R^2 of ln tau_p on the Fantasia age group is >= M3's - 0.05.

This amendment is post hoc and is labelled as such wherever it is used.

**E3 (registered before looking at the component costs), likelihood-ratio selection between components.** Same components as E2. Selection rules:

- Linear stage: oscillator if LR_osc >= 11.3 (unchanged).
- Regime gate: G = 9 (unchanged).
- With the gate open, the regime model is chosen within {RB, RO, QB, QO} by likelihood ratios of the regime fits themselves:
  - the oscillator if 2 (cost_RB - cost_RO) >= 11.3;
  - the response kinetics (Q) if 2 (cost_R - cost_Q) >= 13.28 (chi-squared, 4 df, 0.99).
- c_p and c_s are re-tuned on development by the same rule.

**E3 result (development).**

- Likelihood-ratio selection does not fix exercise: vagal 1.015 x M3, sympathetic 0.769.
- Component analysis on development:
  - An oracle choosing the best component per simulation would reach a vagal RMSE of 0.62 x M3 in exercise, so the information is in the components.
  - RB alone (regimes without the oscillator) reaches 0.84 / 0.61 in exercise, against 0.96 / 0.67 for RO.
- Diagnosis on the six worst exercise simulations:
  - The oscillator estimated by the stationary linear fit is strongly biased under exercise: lo = 0.02-0.09 against the true 0.3. Held fixed in the regime model, it degrades the states, while the regime likelihood still prefers it.
  - The vagal set-point shift is over-assigned (about 1.7 x the true one), because the sympathetic response is forced to the slow OU rate.
  - One simulation has a sign error: vagal +0.42 with the sympathetic shift at 0, which the reciprocity prior does not penalize.

**E4 (post hoc, suggested by the E3 component analysis on development; registered before its evaluation), regimes without the oscillator.** Linear stage as M4e (gated oscillator). With the gate open, the regime model is RB (no oscillator). c_p and c_s are re-tuned by the same rule.

**E4 result (development).**

- Simulation: exercise passes S1' (0.838 / 0.608). S3 fails: tau_p in regime and two-regime (0.39 / 0.38 against M3 + 0.10), and tau_s in exercise (0.55).
- Real data: vagal-first falls to 93% (R3), and the tau_p ICC to 0.49 (R4), because windows mix models with and without the oscillator.
- Trading one criterion for another by choosing components is not the way forward.

**E5 (registered before its run; the kinetics statistic was inspected on the E3 development costs), a model check on the response kinetics, plus a jointly estimated oscillator.**

- New component JO: the regime model with response rates = kappa and the oscillator estimated jointly (alpha and lo free, omega from the linear fit). Checked on eight development cases:
  - regime scenario: tau_p closer to the truth than RB, states equal or better;
  - exercise: still worse than RB.
- Kinetics check K = 2 (cost_RB - cost_QB). It is a likelihood-ratio test of "set-point responses follow the OU rates", with threshold 13.28 (chi-squared, 4 df, 0.99).
  - On E3 development costs it is >= 13.28 in 88% of exercise simulations, 0% of regime, 3% of two-regime and 0% of rest/RSA.
  - On development real data it is >= 13.28 in 71% of exercise, 72% of PRCP and 9% of Fantasia records with the gate open.
- Rule: linear stage as M4e. With the gate open:
  - K >= 13.28 (misspecified kinetics) -> RB (variant M4e5a) or QB (M4e5b);
  - otherwise JO if 2 (cost_RB - cost_JO) >= 11.3, else RB.
- c_p and c_s are re-tuned by the same rule.

**E5 result (development).**

- M4e5a:
  - passes S1' in every scenario (exercise 0.846 / 0.610) and S2;
  - fails S3 on tau_s in exercise (0.55);
  - on real data it passes R1 and R2, fails R3 (vagal-first 93%; RB in real exercise), and fails R4' marginally (within-subject SD of tau_p 0.179 against the limit of 0.176).
- M4e5b (QB under mismatch) is worse for the states.

**E6 (registered before its run), ordered response kinetics.** Q models with vagal on-rate >= sympathetic on-rate (the analogue of tau_p < tau_s), without (QBo) and with (QJo) the joint oscillator. Result on development simulations:

- The ordering does not fix the states: exercise vagal 0.98 (QBo) and 1.08 (QJo) x M3.
- QJo has the most accurate parameters in exercise: tau_p 0.12, tau_s 0.16.
- In exercise the regime models therefore trade states against parameters. With free kinetics, the split of the RR step between the two drives is weakly identified from the RR alone; RB's fixed kinetics regularize it.
- On real exercise data, R3 is met by JO (95 / 97) and RO (97 / 95), not by RB (95 / 93; one recording below the threshold).

**Freeze (2026-11-17 01:50; before any confirmation result was read).** To stop iterating on the development data, two candidates are frozen, with settings in `_out/frozen_o.rds`: q_osc = 11.3, G = 9, qk = 13.28, c_p = c_s = 2.5 (filtered and smoothed, regime components only).

- **M4e**, the pre-registered design: gated oscillator; regimes RO/RB behind the gate; rescaled intervals.
- **M4e-h**, the result of E1-E6, with states and parameters from different components:
  - Linear stage: oscillator if LR >= 11.3.
  - Gate open, kinetics check K = 2 (cost_RB - cost_QB) >= 13.28: states from RB; parameters (tau_p, tau_s, SD_p, SD_s) from QJo.
  - Gate open, no mismatch: JO if its oscillator LR >= 11.3, else RB, for both states and parameters.
  - Rationale: under a detected kinetic mismatch, the constrained model regularizes the weakly identified split of the state, while the model whose kinetics pass the check estimates the parameters.
  - Development results: passes S1', S2 and S3 in every simulated scenario. Real development data with QJo were not run before the freeze.

Confirmation: seeds `m4e-conf`; 60 new exercise recordings; Fantasia C/D; the 81 blockade segments. Criteria S1', S2, S3, R1, R2, R3, R4 (and the post hoc R4'), B'. Each candidate is reported as is.

## Confirmation (2026-11-17 02:00; frozen settings, unseen data)

**Simulation** (seeds `m4e-conf`, 40 per scenario, filtered RMSE ratio to M3):

| Scenario | M4e vagal / symp | M4e-h vagal / symp | M4e-h tau_p, tau_s error (M3) |
|---|---|---|---|
| rest | 1.000 / 1.000 | 1.000 / 1.000 | 0.08, 0.19 (0.08, 0.19) |
| RSA | 0.981 / 0.987 | 0.981 / 0.987 | 0.11, 0.20 (0.35, 0.19) |
| regime | 0.385 / 0.524 | 0.380 / 0.522 | 0.14, 0.10 (0.21, 0.45) |
| exercise | **1.112** / 0.692 | **0.976** / 0.620 | 0.11, 0.28 (1.14, 0.45) |
| two | 0.663 / 0.826 | 0.663 / 0.838 | 0.16, 0.20 (0.13, 0.38) |

- Coverage: S2 is met by both, with filtered coverage 0.79-1.00 in the non-stationary scenarios against M3's 0.56-0.77.
- S3 is met by M4e-h. M4e fails tau_s in exercise (0.58).
- Both fail S1' in exercise, on the vagal state. For M4e-h's state component (RB), the per-simulation ratio has a median of 0.77 and is better than M3 in 72% of the simulations, but 5% of the fits are catastrophic (> 1.5 x M3) and drive the mean to 0.97. In development the same rule gave 0.85, so the exercise vagal result is also unstable between seed sets.

**Real recordings** (60 new exercise recordings, Fantasia full / head / C / D; 100 per part):

| | M3 | M4e | M4e-h |
|---|---|---|---|
| KS pass (full) | 32 | 45 | 46 |
| PACF violations | 11.07 | 7.43 | 7.52 |
| held-out gain [lower 95%] | — | 0.031 [0.004] | 0.026 [-0.003] |
| exercise signs / vagal-first | 98 / 98 | 98 / 98 | 98 / 97 |
| ICC tau_p / SD_p / SD_s | 0.86 / 0.95 / 0.90 | 0.65 / 0.92 / 0.88 | 0.61 / 0.92 / 0.86 |
| within-subject SD tau_p / SD_p / SD_s | 0.162 / 0.144 / 0.198 | 0.151 / 0.166 / 0.217 | 0.155 / 0.175 / 0.250 |
| R^2 of tau_p on age group | 0.19 | 0.23 | 0.35 |

**Blockade** (atropine SD_p ratio; atropine SD_s; propranolol SD_p, spontaneous / paced):

| | atropine SD_p | atropine SD_s | propranolol SD_p |
|---|---|---|---|
| M3 | 0.16 / 0.19 | 0.26 / 0.38 | 1.42 / 0.54 |
| M4e | 0.13 / 0.12 | 0.30 / 0.27 | 1.07 / 0.46 |
| M4e-h | 0.12 / 0.11 | 0.23 / 0.22 | 1.07 / 0.46 |

**Decision.**

| Criterion | M4e | M4e-h |
|---|---|---|
| S1' | fails (exercise vagal 1.11) | fails (exercise vagal 0.98) |
| S2 | passes | passes |
| S3 | fails (exercise tau_s) | passes |
| R1 | passes | passes |
| R2 | passes | fails (lower bound -0.003) |
| R3 | passes | passes |
| R4 | fails | fails |
| R4' (post hoc) | fails (within-subject SD of SD_p 0.166 > 0.158) | fails (SD_p, SD_s) |
| B' | fails (paced propranolol 0.46 vs 0.54; paced atropine SD_s farther from 1) | fails (same, plus atropine SD_s) |

**Neither frozen candidate is a strong replacement for M3.** The confirmation data (simulation seeds, exercise sample, Fantasia C/D, blockade) have now been used; a new candidate would need fresh confirmation data.

## Reading

1. **What replicates on unseen data** (both candidates):
   - stationary records equal to M3, and the RSA bias of tau_p removed (0.35 -> 0.11);
   - set-point regimes with responses at the drives' own rates (regime, two): states 34-62% better, with calibrated intervals;
   - the tau_p bias in exercise removed (1.14 -> 0.11-0.22);
   - better residual diagnostics (KS 45-46 vs 32 of 100, fewer PACF violations);
   - a forward-prediction gain of about +0.03 nats per beat (significant for M4e);
   - the spontaneous-breathing propranolol ratio closer to 1 (1.07 vs 1.42).
2. **What does not replicate: the vagal state in exercise.** When the set-point responses have their own kinetics, splitting an RR step into a vagal and a sympathetic part is weakly identified from the RR alone:
   - fixed kinetics (RB) regularize the split but bias tau_s;
   - free kinetics (Q) estimate the parameters well but split the step badly;
   - 5% of the regime fits fail catastrophically, including sign errors that the reciprocity prior does not exclude.
   No selection rule fixed this robustly in six development rounds (E1-E6).
3. **The tau_p reliability trade-off is real.** With the oscillator, the within-subject precision of tau_p is equal or better (0.151-0.155 vs 0.162) and the age signal is stronger (R^2 0.23-0.35 vs 0.19). The ICC falls (0.61-0.65 vs 0.86) because M3's between-subject spread partly reflects the respiratory amplitude. The stationary amplitude SD_p is somewhat less stable with the M4 models (within-subject SD 0.166-0.175 vs 0.144).
4. **What would be needed.** Information that is not in the RR step itself:
   - a vagal modulation of the respiratory oscillation amplitude (RSA falls with vagal withdrawal). This is a real physiological coupling, but the simulator has no such coupling, so a fair test needs data with measured respiration or blockade;
   - or measured inputs (workload, tilt angle), as M4d-I showed.
   Either is a new model and needs new confirmation data. Iterating further on the present data would be fitting the evaluation.
