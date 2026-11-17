# M4 prototypes: a unified evaluation against M3

Experimental, 2026-11-16. Nothing here touches the engine, the pipeline, the manuscript or the app. Builds on the extended filter of `experiments/m3plus` (oscillator, RTS smoother, inputs), the regime model of `experiments/regimes/v2` (HMMr: two regimes, reciprocity prior, change-point gate) and the feature anchor of `experiments/features`.

## Prototypes (a nested ladder)

| Model | Definition |
|---|---|
| M3 | production specification (extended filter, three production starts) |
| M4a | M3 + multi-start (production starts + 12 dispersed starts, tau_p 2-20 s x tau_s 60-600 s) + damped respiratory oscillator in the vagal pathway (frequency on a 0.15-0.40 Hz grid, two damping starts); states reported filtered and RTS-smoothed |
| M4b | M4a + probabilistic anchor: ln lambda_P ~ N(ln m, 0.5^2), m = scale-free feature anchor (entropy, SD1/SD2, DFA, spectral ratios; elastic net trained on the 1000 simulations of `experiments/features`, target = the OU part lambda_P); MAP over a 5-point grid m x e^{-1, -0.5, 0, 0.5, 1} |
| M4c | M4b + latent regimes of the set-points: change-point gate (gain >= 9 log-likelihood units), then the 2-regime HMM with reciprocity prior and the M4b oscillator; closed gate = M4b |
| M4c-w (post hoc) | M4a + the same regimes as M4c, with the production wavelet anchor and the M4a oscillator (closed gate = M4a). Added after the first simulation run showed that the feature anchor of M4b fails on non-stationary records; it separates the effect of the regimes from that of the anchor. Not part of the pre-specified ladder |
| M4d-H | M4c with empirical-Bayes time-scale priors: mean and between-recording SD (floor 0.15) of ln tau_p, ln tau_s from weak-prior fits (95% intervals 0.5-60 s, 10-2000 s) of the OTHER recordings of the same dataset (or simulation scenario) |
| M4d-I | M4b + known inputs (rectangular events with GLS amplitudes) where they are independent of the RR: true events in simulation, annotated tilts in PRCP; not applicable to exercise (events would come from the same RR) or Fantasia |

Not included: the adaptive noise scale (a two-stage rescaling, not a likelihood; the regime filter does not take it) and the rejected ideas of `experiments/m3plus` (rate-dependent noise, Student-t, random-walk set-points).

## Evaluation (`sim.R`, `real.R`, `blockade.R`)

- **Simulation with truth** (40 per scenario; generator `sim_regime` of `experiments/regimes` with the M3 kinetics, 420 s): rest (model, no oscillator), RSA (rest + oscillator), regime (one rectangular regime + oscillator), exercise (asymmetric on/off profiles + oscillator), two (two tilt-like regimes, 640 s, + oscillator). Normalized RMSE and 95% coverage of the vagal-like (x_p + oscillator) and sympathetic-like states, filtered and smoothed; |ln error| of tau_p, tau_s.
- **Real recordings**: 60 exercise (same sample as `regimes/v2`), 40 Fantasia, 10 PRCP. In-sample: KS pass, PACF violations, tau_p (median, share inside 3-6 s); held-out forward score (fit on the first 70% of beats, filter the whole record with the head's parameters, mean log predictive density of the last 30%; regimes inferred online; degenerate tails excluded as in `regimes/v2`); exercise realism (rest -> end of exercise: vagal-like falls and sympathetic-like rises, vagal-like reaches half its change first); Fantasia test-retest ICC(3,1) of ln tau_p, ln SD_p, ln SD_s between windows A and B.
- **Blockade** (Lancaster, 81 segments): ratios to saline of SD_p (vagal-like total: OU + oscillator) and SD_s, as in the paper.

**Decision rule (fixed before the run).** A prototype improves on M3 if it passes all of:
S1 mean RMSE of both states <= M3 in every simulated scenario; S2 95% coverage >= min(0.90, M3 - 0.05); S3 median |ln error| of tau_p and tau_s <= M3 + 0.10;
R1 KS pass and PACF violations no worse than M3; R2 mean held-out gain > 0 with the bootstrap 95% lower bound > 0 (all recordings);
R3 exercise realism: expected signs and vagal-first in >= 95% of recordings; R4 Fantasia ICC of tau_p, SD_p, SD_s >= M3 - 0.05;
B blockade: atropine SD_p ratio < 0.3, atropine SD_s ratio closer to 1 than M3's, propranolol SD_p ratio >= 0.8.

## Files

`lib.R` (the ladder), `sim.R`, `real.R`, `blockade.R`; outputs in `_out/` (git-ignored). Follow-up: `e/` (M4e, development and confirmation); `f/` (M4f design, RSA amplitude coupled to the vagal level; stage-0 gate not met).

## Results

Run of 2026-11-16/17: 200 simulations, 110 real recordings (300 jobs: full, head, Fantasia A/B), 81 blockade segments. Logs in `_out/*_log.txt`.

### Simulation (filtered states; mean normalized RMSE vagal / sympathetic)

| Model | rest | RSA | regime | exercise | two |
|---|---|---|---|---|---|
| M3 | 0.500 / 1.62 | 0.522 / 1.29 | 0.680 / 1.32 | 0.784 / 0.81 | 0.958 / 1.65 |
| M4a | 0.499 / 1.62 | 0.514 / 1.27 | 0.656 / 1.28 | 0.781 / 0.83 | 0.929 / 1.61 |
| M4b | 0.484 / 1.54 | 0.508 / 1.24 | 0.951 / 1.83 | 1.095 / 0.72 | 1.075 / 1.88 |
| M4c | 0.487 / 1.54 | 0.508 / 1.24 | 0.381 / 0.84 | 0.506 / 0.47 | 1.000 / 1.76 |
| M4c-w | 0.502 / 1.61 | 0.514 / 1.27 | **0.245 / 0.61** | 0.680 / 0.56 | **0.619 / 1.16** |
| M4d-H | 0.482 / 1.52 | 0.510 / 1.24 | 0.490 / 1.07 | 1.296 / 0.89 | 1.051 / 1.94 |
| M4d-I | — | — | 0.280 / 0.71 | 0.559 / 0.59 | 0.412 / 1.07 |

Median |ln error| of tau_p (M3 -> best): exercise 1.09 -> 0.09 (M4c), 0.15 (M4c-w); RSA 0.40 -> 0.11 (every oscillator model); M4a alone worsens it in exercise (1.21) and regime (0.44); M4d-H 1.88 in exercise. tau_s: M4c and M4c-w 0.18-0.23 in the non-stationary scenarios (M3 0.36-0.46); M4d-H 0.68-0.76 everywhere.
Coverage of the 95% intervals: M4b, M4c and M4d-H cover x_s in only 0.80 of the beats at rest (M3 0.95), because the feature anchor moves lambda_P; M4c-w covers x_p 0.57 and x_s 0.46 in exercise (M3 0.66 and 0.60): better point estimates, too-narrow intervals.
The feature anchor (M4b) fails out of its training distribution: it was trained on stationary simulations, and in regime/exercise/two it worsens the vagal RMSE by 12-43% and is worse than M3 in 85-100% of the simulations.

### Real recordings

| Model | KS pass (ex / Fa / PRCP) | PACF viol. (ex / Fa / PRCP) | held-out gain vs M3, nats/beat [boot. lower 95%] | exercise signs / vagal-first, % | Fantasia ICC tau_p / tau_s / SD_p / SD_s |
|---|---|---|---|---|---|
| M3 | 18 / 60 / 20 | 14.2 / 6.8 / 14.1 | — | 97 / 97 | 0.88 / 0.31 / 0.95 / 0.90 |
| M4a | 20 / 65 / 0 | 13.1 / 5.2 / 13.7 | 0.034 [0.022] | 97 / 97 | 0.62 / 0.39 / 0.95 / 0.88 |
| M4b | 22 / 70 / 0 | 10.6 / 5.2 / 10.9 | 0.035 [0.020] | 95 / 93 | 0.65 / 0.27 / 0.95 / 0.82 |
| M4c | 18 / 72 / 0 | 8.2 / 5.2 / 10.4 | 0.040 [0.020] | 95 / 88 | 0.68 / 0.24 / 0.93 / 0.81 |
| M4c-w | 18 / 68 / 0 | 8.4 / 5.0 / 10.3 | 0.034 [0.014] | 97 / 95 | 0.67 / 0.15 / 0.94 / 0.87 |
| M4d-H | 18 / 70 / 0 | 8.2 / 4.8 / 9.2 | 0.038 [0.015] | 98 / 82 | 0.68 / 0.26 / 0.93 / 0.78 |
| M4d-I (PRCP only) | — / — / 0 | — / — / 13.6 | not scoreable | — | — |

- The held-out gain is almost the same for M4a to M4d-H (0.034-0.040), so it comes mainly from the oscillator: Fantasia +0.06 (90% of the recordings better), PRCP +0.06-0.08, exercise +0.01.
- The regime gate opens in 98-100% of the exercise recordings, 18-20% of Fantasia and 90% of PRCP, and in 14% of the blockade segments.
- With the feature anchor, the median lambda_P in exercise falls from 1.52 to 0.18, which is the same failure as in simulation.
- The oscillator frees tau_p from the respiratory band: the median in Fantasia rises from 3.9 to 4.8 s and the share inside 3-6 s from 55% to 70-80%, but the test-retest ICC of tau_p falls from 0.88 to 0.62-0.68. The stability of the M3 tau_p is therefore partly the stability of the respiratory peak that it absorbs.
- M4d-I is not scoreable on the held-out tail. Tilts that start after the head have no amplitude estimated from the head, so their extrapolation is degenerate in 4 of the 10 records. It also fails the KS check in PRCP.

### Blockade (ratio to saline; geometric mean, n = 6-7 volunteers)

| Model | atropine SD_p (sp / paced) | atropine SD_s (sp / paced) | propranolol SD_p (sp / paced) |
|---|---|---|---|
| M3 | 0.16 / 0.19 | 0.26 / 0.38 | 1.42 / 0.54 |
| M4a | 0.16 / 0.17 | 0.38 / 0.38 | 1.08 / 0.45 |
| M4b | 0.13 / 0.11 | 0.42 / 0.51 | 1.09 / 0.38 |
| M4c | 0.12 / 0.09 | 0.41 / 0.39 | 1.09 / 0.38 |
| M4c-w | 0.13 / 0.12 | 0.30 / 0.27 | 1.08 / 0.46 |
| M4d-H | 0.15 / 0.12 | 0.51 / 0.54 | 1.11 / 0.39 |

Every model keeps atropine SD_p < 0.3 and gives a spontaneous-breathing propranolol ratio closer to 1 than M3 (1.08-1.11 vs 1.42). The atropine SD_s ratio is closer to 1 than M3's with the feature anchor (M4b, M4c, M4d-H: 0.39-0.54); M4c-w is mixed (0.30 vs 0.26 spontaneous, 0.27 vs 0.38 paced). The paced propranolol ratio is below 0.8 for every model, M3 included, with intervals of about 0.1-1.1.

### Decision rule

| Criterion | M4a | M4b | M4c | M4c-w (post hoc) | M4d-H | M4d-I |
|---|---|---|---|---|---|---|
| S1 RMSE <= M3, all scenarios | no (exercise x_s +2%) | no | no (two) | no, by 0.4% (rest x_p 0.502 vs 0.500) | no | yes (3 applicable scenarios) |
| S2 coverage | no, by 0.002 (exercise x_s) | no | no (rest x_s 0.80) | no (exercise 0.57 / 0.46) | no | yes |
| S3 tau errors <= M3 + 0.10 | no (tau_p exercise, regime) | no (tau_s exercise) | yes | yes | no | no (tau_s exercise 0.58) |
| R1 KS / PACF (pooled) | yes | yes | yes | yes | yes | no (PRCP KS) |
| R2 held-out gain, lower bound > 0 | yes | yes | yes | yes | yes | not scoreable |
| R3 exercise >= 95% | yes | no (93) | no (88) | yes (97 / 95) | no (82) | n/a |
| R4 Fantasia ICC | no (tau_p 0.62) | no | no | no (tau_p 0.67) | no | n/a |
| B blockade | no (paced propranolol; M3 also fails) | no | no | no (paced propranolol; M3 also fails) | no (paced propranolol) | n/a |

On R1: pooled, KS passes in 38-41 of 110 recordings with M4 and 37 with M3, and PACF violations are fewer in every dataset. In PRCP, however, KS falls from 2 of 10 to 0 of 10 with every M4.

**No prototype passes the pre-specified rule; nothing is adopted.**

## Reading

1. **Latent regimes are the real gain for non-stationary records.** M4c and M4c-w reduce the state RMSE by 30-60% in the regime, exercise and two-regime scenarios. They also correct the largest M3 bias, tau_p in exercise (|ln error| 1.09 -> 0.09-0.15). M4c-w passes S3, R1, R2 and R3. It fails S1 by 0.4% at rest, S2 in exercise (intervals too narrow), and R4.
2. **The feature anchor is harmful outside its training distribution.** It was trained on stationary simulations. On exercise and regime records it pulls lambda_P down by an order of magnitude, and it drags M4c down with it. Any learned anchor would need non-stationary training data, or a gate that applies it only to stationary records.
3. **The oscillator has two sides.** It improves the forward prediction (Fantasia, PRCP) and removes the tau_p bias under RSA. However, it barely changes the states (0-4%) and costs test-retest stability of tau_p (ICC 0.88 -> 0.62). For the paper, the reliability of tau_p reported for M3 rests in part on the respiratory peak.
4. **Empirical-Bayes priors on tau (M4d-H) hurt.** Pooling time constants across heterogeneous recordings biases tau_s (|ln error| about 0.7) and degrades exercise. Known inputs (M4d-I) give the best states where true events exist, but they are not general and cannot be scored forward.
5. **Next candidate, if pursued.** M4c-w with three changes, to be fixed in advance and tested on new seeds:
   - intervals recalibrated in the regime filter;
   - the oscillator only when its likelihood gain is significant, to keep the tau_p stability;
   - the regime gate tuned on rest, which opens in 2% of rest and 20% of Fantasia.

## Fixes during the run

- `real.R`: invalid data.table syntax in the PRCP subset (base subsetting used).
- `blockade.R`: in stage 2, `protocol$subject[j]` inside `W[...]` resolved `j` to the column `j` of W, so the other-subject set was empty and all 81 M4d-H refits failed. The local variable was renamed and the blockade script was re-run.
