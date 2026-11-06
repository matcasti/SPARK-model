# Latent regime (Markov) states on top of M3: prototype and comparison

Experimental, overnight 2026-11-05/06. Nothing here touches production M3 (engine, pipeline, manuscript and app are unchanged).

## Question

Can latent discrete regimes (a hidden Markov chain over beats) explain autonomic transitions better than the step functions of the events model, and do it at least as well as M3 in state reconstruction, informativeness, physiological realism and residual structure?

## Model (`reg_filter.cpp`, `lib.R`)

M3 with regime-dependent set-points: each drive relaxes toward `mu_b(s_k)`, where the regime `s_k` is constant within a beat and follows a Markov chain. Regime 0 is the reference (`mu = 0`, i.e. M3 itself). Covariances do not depend on the regime, so regimes only shift the predicted mean; the Kalman filter is run per regime and collapsed after each beat (GPB1). `nu0` and `sigma²` are explicit parameters; priors on tau_p, tau_s, separation, jitter as in M3, plus `mu ~ N(0, 0.8²)` and a persistence prior on the stay probabilities.

Variants compared (all fitted by BFGS in whitened coordinates, initialised from a change-point fit):

| Variant | Description |
|---|---|
| M3 | production model |
| M3E_det | M3 plus step events estimated from the same RR (the earlier comparison) |
| CP | one rectangular regime `[t0, t1]` found by grid + Nelder–Mead, amplitudes by GLS |
| HMM2 | two regimes (reference + one), persistence prior 7 (sd 0.7) |
| HMM3 | three regimes (overfits and chatters; dropped) |
| **HMM2g** | HMM2 behind a gate: the regime model is used only if the change-point search gains at least 9 log-likelihood units; otherwise M3 |

Tuning (`run_tune2.R`, `run_sim2.R`): gate 9, persistence prior 7 with sd 0.7 was the best compromise across simulations and PRCP. Stricter priors (9, 11) lost the real transitions; looser ones chatter.

## Results

Simulations (`_out/sim.rds`, `_out/tune2_sim.rds`; 40-60 per scenario, truth known). RMSE of the filtered states relative to M3 (negative is better), gated HMM2:

| Scenario | Vagal | Sympathetic | Onsets found (true) |
|---|---|---|---|
| One regime | -45% | -24% | 1.02 (1) |
| Two regimes | -28% | -20% | 1.81 (2) |
| Exercise-like step | +4% | -29% | 1.03 (1) |
| Rest (no regime) | 0% (gate closed) | 0% | 0 |

Without the gate, rest worsens (vagal +55%, sympathetic +54%) because the model finds regimes in noise; the gate removes it (0% gate opening in simulated rest). In the exercise-like step, M3E_det (which has the true step structure) is better than HMM2 in vagal RMSE (-23% vs +4%).

Real data (150 recordings; `_out/real.rds`, `analyze_real.R`), HMM2g unless stated:

| Dataset | Gate opens | KS pass % (M3 -> HMM2g) | PACF violations | tau_p median (M3 -> HMM2g) | dLL vs M3 (median) |
|---|---|---|---|---|---|
| Exercise (100) | 100% | 13 -> 25 | 14.1 -> 9.7 | 16.5 -> 5.3 s | +71.5 |
| PRCP (10) | 90% | 20 -> 0 | 14.1 -> 11.3 | 7.1 -> 3.4 s | +75 |
| Fantasia (40) | 20% | 60 -> 60 | 6.8 -> 6.7 | 3.9 -> 3.9 s | 0 |

- KS in PRCP gets worse (20% to 0% of 10 recordings) even though dLL and PACF improve; with 10 recordings this is weak evidence either way.
- Realism (exercise, rest-to-end change): vagal faster than sympathetic in 93% of recordings (M3 98%, M3E_det 94%), expected signs in 94% (M3 98%, M3E_det 75%). The M3 split is preserved better than by M3E_det or CP.
- Held-out forward score (fit first 70%, score last 30%), HMM2 minus M3 in nats per beat: exercise -0.106 (better in 49%), PRCP +0.019 (better in 90%), Fantasia ~0.
- PRCP annotated tilt changes: HMM2 finds 21% of onsets and 19% of all annotated changes within 30 s, with 71 detections of which 48 have no annotation nearby; M3E_det finds 50% / 28% with 50 detections. CP finds 5%.

Figure: `_out/figure.png` (exercise recording and PRCP 12734: heart rate, regime posterior, vagal and sympathetic states under M3 and HMM2). In PRCP the posterior follows the tilt blocks closely; in exercise it follows the ramp and adds some flicker at rest.

## Verdict

**HMM2g is the winner of the regime variants**, and it meets "at least equal to M3" on the criteria that can be checked:

- State reconstruction: -20% to -45% RMSE in regime scenarios, equal to M3 at rest (by construction of the gate).
- Informativeness and residual structure: KS and PACF clearly better in exercise, tau_p moves toward its prior (5 s instead of 16 s), dLL +71 with AIC better in 99% of exercise recordings.
- Realism: vagal-fast, sympathetic-slow ordering and expected signs preserved.

It is **not** a candidate for production yet:

1. It does not beat the simpler M3E_det where the transition is a single step (exercise simulation, PRCP timing).
2. Held-out forward prediction in exercise is worse than M3 (-0.106 nats/beat), so part of the in-sample gain is fit, not prediction.
3. The gate opens in 20% of Fantasia recordings, none in simulated rest; the HMM detects many changes that the PRCP annotations do not list (they may be real, but they cannot be validated).
4. GPB1 is an approximation; the regime posterior chatters at rest in exercise.
5. Another layer (gate, persistence prior, K) means more to defend to reviewers.

## Round 2

See `v2/README.md`: reciprocity prior (no regime may push both drives the same way), smoothed regime probabilities and an oscillator variant. The recommended model is now the reciprocity-prior HMM with gate (better than M3 or equal in nearly every criterion, with the exceptions listed there).

## Next steps if pursued

- Replace GPB1 with a proper smoother over regimes (or an exact two-regime forward-backward on a discretised drive) and check chatter at rest.
- Held-out evaluation of the regime model with parameters fitted on the whole of a different recording (population `mu`), which is what the hierarchical idea needs.
- Compare against ground-truth transitions in a dataset with known autonomic manoeuvres (more than the 10 PRCP recordings).

## Files

`reg_filter.cpp`, `lib.R` (engine, fitting, change-point search, simulation); `run_sim.R`, `run_sim2.R`, `run_tune2.R` (simulations and tuning); `run_real.R <stay_prior> <gate>` and `analyze_real.R` (real data); `figure.R`; outputs in `_out/`.
