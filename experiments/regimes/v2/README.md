# Regime model, round 2: toward "at least M3 on everything"

Experimental, 2026-11-06. Production M3 is unchanged. Starts from the gated 2-regime HMM of round 1 (`../README.md`) and tries to remove its weaknesses. Figure: `_out/figure.png`.

## What was wrong in round 1, and what it turned out to be

- **Held-out exercise mean was negative (-0.106 nats/beat).** Not the model: three recordings dominate it, and one of them (104) has 560 s of constant 15 bpm at the end (sensor dropout). With degenerate tails excluded the median is ~0 and the mean is -0.025.
- **Vagal state worse than M3 in the exercise-like simulation (+7%, 30% of simulations worse, worst 7.9x).** The change-point initialisation sometimes lands on unphysiological targets (vagal and sympathetic amplitudes of the same sign, cancelling each other in the rate: e.g. mu_p = -1.5, mu_s = -1.5). Nothing in the model excluded it.
- **Regime posterior flickers at rest and onsets are late** (filtered probabilities only see the past).

## What was added

1. **Reciprocity prior** (`recip = 20` in `fit_regime`): hinge penalty `20 * max(mu_p * mu_s, 0)`, i.e. a regime may not push both drives in the same direction (it allows either reciprocal pattern: withdrawal + activation, or the opposite as in sleep). The change-point initialisation is sign-corrected.
2. **Smoothed regime probabilities** (Kim's backward recursion on the filtered probabilities, `regime_smooth`), used for onset timing and plots; states can use GPB1 + RTS (`reg2_filter`, `smooth = TRUE`).
3. **Oscillator variant** (HMMO): the M3O damped oscillator (frequency, damping and amplitude from the M3O fit) inside the regime filter.

## Variants that did not help (kept in `rounds/`)

| Idea | Result |
|---|---|
| Tight prior on the targets (sd 0.4, bound 1.5) | Exercise vagal +45% vs M3 (the true amplitudes are larger) |
| Anchor the time scales to the M3 fit | Exercise vagal +10 to +18%: M3's tau_p is inflated by the step |
| Graded chain 0-1-2 (staircase for ramps) | Worse everywhere (exercise vagal +62%) |
| Recovery regime with its own slower drive rates (`reg3_filter.cpp`) | Log-likelihood gain 0 to 2 nats, no RMSE gain |
| Three starts for the HMM fit | No change (fits already reach the same optimum, 0.1 s each) |
| Kim regime smoother feeding a soft-input state smoother | No better than GPB1 + RTS for the states (but the smoothed regime probability itself is much better for timing) |
| Second-stage gate on the fitted log-likelihood gain | No benefit: the change-point gate already never opens in simulated rest; it only costs power |

## Results

Simulations (40 per scenario, truth known, respiratory oscillator in the simulator; change-point gate 9). RMSE relative to M3 (%; negative is better), mean over simulations, filtered estimates:

| Scenario | Vagal | Sympathetic | 95% coverage vagal / sympathetic (M3) | Simulations worse than M3 (vagal) |
|---|---|---|---|---|
| One regime | -62 | -49 | 0.93 / 0.78 (0.67 / 0.80) | 0% |
| Two regimes | -38 | -21 | 0.76 / 0.70 (0.49 / 0.78) | 2% |
| Exercise-like (slow, asymmetric responses) | -16 | -40 | 0.59 / 0.49 (0.70 / 0.62) | 25% |
| Rest (no regime) | 0 | 0 | gate closed in 100% | 0% |

Smoothed (RTS): -71/-56 (one regime), -42/-22 (two), -17/-31 (exercise). Round 1 (agnostic prior) in the same exercise scenario: +7 vagal, -23 symptomatic. The oscillator does not improve state recovery (HMMO exercise vagal -1%; M3O vs M3 -3% at best). Gate opening in simulated rest: 0%.

Real recordings (60 exercise, 20 Fantasia, 10 PRCP). HMM0 = round 1, HMMr = reciprocity prior, HMMOr = reciprocity + oscillator; regime models gated (closed gate = M3, or M3O for HMMOr):

| | M3 | M3O | HMM0 | **HMMr** | HMMOr |
|---|---|---|---|---|---|
| Exercise: KS pass % | 18 | 20 | 20 | **23** | 18 |
| Exercise: PACF violations | 14.2 | 13.2 | 9.3 | **8.8** | 8.4 |
| Exercise: log-likelihood vs M3 (median) | 0 | +34 | +66 | **+72** | +92 |
| Exercise: tau_p median (s) | 17.9 | 19.4 | 5.3 | **5.0** | 6.2 |
| Exercise: vagal faster than sympathetic % | 97 | 97 | 90 | **93** | 95 |
| Exercise: expected signs % | 97 | 97 | 88 | **95** | 97 |
| Exercise held-out vs M3 (nats/beat, mean) | 0 | +0.010 | -0.025 | **+0.003** | +0.007 |
| PRCP held-out vs M3 (mean) | 0 | +0.056 | +0.019 | **+0.021** | +0.069 |
| PRCP PACF violations | 14.1 | 13.7 | 11.3 | **11.0** | 10.2 |
| Fantasia: gate opens in-sample | | | 25% | **25%** | 25% |

Annotated PRCP onsets found within 30 s (62 in total): M3 0, round 1 HMM0 13, HMMr filtered 8, **HMMr smoothed 18 (24 within 60 s)**; M3E_det (round 1) 31. Smoothed regime probabilities remove the rest flicker and line up with the annotated tilts (figure).

## Verdict

**Recommended regime model: HMMr (2 regimes, reciprocity prior, change-point gate), with Kim-smoothed regime probabilities and GPB1 + RTS states for retrospective use.** Against M3 it is better or equal on state reconstruction (-16% to -62% vagal, -21% to -49% sympathetic in every scenario with a regime; identical at rest), informativeness (regime posterior; 18/62 annotated tilts timed within 30 s where M3 has no notion of one), residual structure (PACF 14.2 to 8.8, KS 18% to 23%), physiological realism (tau_p 5.0 s instead of 17.9 s, in line with the events model; ordering and signs 93-95% vs 97%) and held-out prediction (parity to slightly better; PRCP +0.021).

It is not better in everything:

1. Exercise-like recoveries slower than the drive time scale: 25% of the simulated recordings have a worse vagal state than M3 (worst 1.7x) and 95% coverage drops to 0.59/0.49 (M3 0.70/0.62).
2. Vagal-faster and sign checks are 93-95% against 97% for M3.
3. Against M3O (not M3) the held-out gain in exercise is ~0 (-0.007); the oscillator gives the predictive gain, the regimes the structure. HMMOr (both) is the best predictor (PRCP +0.069 over M3) but its states are not better than HMMr's.
4. PRCP KS pass is 0% (M3 20%) while PACF and likelihood improve; the PRCP detector still finds fewer onsets (18/62) than M3E_det (31/62).
5. The gate opens in 25% of Fantasia recordings (no harm measured, but not zero).

If this were to be pursued for production, the model to carry is HMMr; the open problem is the exercise recovery shape (state-dependent time scales did not help, so it likely needs a deterministic recovery profile as in the events model).

## Round 2b: hybrid (HYB) = regime model dates the transitions, events model estimates the states

HYB: fit the gated reciprocity-prior HMM, take the episodes where the Kim-smoothed probability of the non-reference regime exceeds 0.5 for at least 30 beats (at most 6-8 longest) as events, and fit M3E (rectangular inputs with GLS amplitudes, as in the pipeline's events model). Scripts: `hybrid_sim.R`, `hybrid_real.R`; code in `lib2.R` (`regime_episodes`).

Simulations (40 per scenario, filtered estimates; RMSE vs M3, % negative better; M3E_true uses the true events and is the oracle):

| Scenario | M3E_det | M3E_true (oracle) | HMMr | **HYB** |
|---|---|---|---|---|
| One regime, vagal / symp | -30 / -14 | -63 / -50 | -62 / -49 | **-62 / -50** |
| Two regimes | -38 / -11 | -63 / -32 | -38 / -21 | -36 / -17 |
| Exercise-like | -30 / -42 | -31 / -43 | -16 / -40 | -25 / -41 |
| Rest | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| Simulations worse than M3 (vagal), exercise | 5% | 2% | 25% | **10%** |
| 95% coverage vagal / symp, exercise | 0.67 / 0.63 | 0.68 / 0.64 | 0.59 / 0.49 | **0.64 / 0.63** |
| 95% coverage vagal / symp, one regime | 0.78 / 0.80 | 0.95 / 0.97 | 0.93 / 0.78 | **0.95 / 0.97** |

Smoothed estimates give the same ordering (HYB: -72/-59 one regime, -29/-35 exercise, -40/-18 two).

Real recordings (in-sample, smoothed states; HYB only where the gate opens: 57/60 exercise, 8/10 PRCP):

| | M3 | M3E_det | HMMr | **HYB** |
|---|---|---|---|---|
| Exercise: KS pass % | 18 | 23 | 23 | **26** |
| Exercise: PACF violations | 14.2 | 9.4 | 8.8 | 8.9 |
| Exercise: log-likelihood gain over M3 / minus parameter count | 0 | 55 / 53 | **72 / 68** | 63 / 60 |
| Exercise: tau_p (s) | 17.9 | 6.6 | 5.0 | 6.6 |
| Exercise: expected signs % | 97 (filtered) | 70 | 93 | 91 |
| PRCP: KS pass % | 20 | 11 | 0 | **12** |
| PRCP: log-likelihood gain / minus parameters | 0 | 100 / 88 | 83 / 79 | **121 / 109** |
| PRCP: onsets within 30 s (of 62) / detections without annotation | 0 / 0 | 31 / 17 | 18 / 18 | 16 / **12** |

The smoothed M3 states break the "vagal faster than sympathetic" check (17%); the regime and events models keep it at 91-93%. The M3 figure in the table (97%) uses filtered states.

**HYB is the best all-rounder so far.** It reaches the oracle in the one-regime scenario with calibrated intervals, closes most of the exercise gap of HMMr (worse-than-M3 vagal 25% to 10%, coverage 0.59/0.49 to 0.64/0.63), and has the best KS in exercise and the best likelihood in PRCP. It gives up some likelihood to HMMr in exercise (63 vs 72) and two onsets in PRCP, and it does not match M3E_det in exercise-like recordings (-25 vs -30 vagal), where the detector is built for that profile. It uses the whole record to date the events, so no forward (held-out) score has been computed yet; that, and the calibration of the gate, are the next runs.

## Round 2c: forward score of the hybrid and calibration of the gate

**Forward (held-out) score** (`hybrid_held.R`; everything detected and fitted on the first 70%, last 30% scored with the head's parameters; an episode still active at the head's end is assumed to persist; recordings with a degenerate tail excluded). Nats per beat minus M3:

| | M3E_det | HMMr (regime inferred online) | HYB |
|---|---|---|---|
| Exercise (60): mean / median / trimmed mean | -0.022 / +0.005 / -0.004 | **+0.003 / +0.010 / +0.001** | -0.073 / +0.005 / -0.011 |
| Exercise: better than M3 / worst | 55% / -0.53 | 57% / -0.36 | 52% / -1.87 |
| PRCP (10): mean | -0.016 (10% better) | **+0.021 (90% better)** | -0.016 (0% better) |

Both events-based models lose to M3 out of sample (their per-recording scores correlate at 0.92) and only the regime model with online inference wins. Fixed events cannot anticipate the changes that happen in the tail (PRCP alternates every few minutes), and the parameters fitted with events (shorter tau_p, smaller sigma^2) are overconfident on unmodelled changes. So HYB is a retrospective estimator, not a predictor.

**Gate calibration** (`gate_calib.R`; 20 Fantasia recordings, 6 parametric bootstraps each from their M3 fit and from their M3O fit): the gain of the change-point search under the stationary null has a 95th percentile of 8.5 (M3) and 8.9 (M3O); the current gate (9) therefore has a false-open rate of about 3-5% under the null. The observed gain in real rest is larger (median 5.8 vs 4.8; 90th percentile 10.7 vs 7.7; 30% of these 20 and 20% of 40 exceed 9), so the Fantasia openings are not a calibration artefact of the gate or of the respiratory oscillator: real rest is less stationary than the model. The regimes found there are small (round 1: log-likelihood gain median 3-4 nats, no change in the metrics when gated). Decision: keep 9 (about 5% under the null); a stricter gate of 13 would cut Fantasia openings to about 5% at a cost in power (PRCP median gain 15).

## Round 2d: repeated manoeuvres as a test-retest (PRCP)

The exercise cohort has one record per person and no session labels, so there is no between-session retest. PRCP has 6-7 annotated tilts or stand-ups per recording (60 usable blocks in 10 people), which gives a within-recording retest (`prcp_reliability.R`). Index per block: mean state in [onset + 30 s, return] minus mean state in the supine window [onset - 90 s, onset - 10 s]. SNR = mean/SD of the block effect (sign flipped for the vagal state, so positive is the expected direction); ICC = one-way, between recordings.

| Model (states smoothed) | Vagal: SNR / expected sign / ICC | Sympathetic: SNR / expected sign / ICC |
|---|---|---|
| Heart rate (reference) | 1.88 / - / 0.83 | |
| **M3** | 1.31 / 93% / **0.86** | **1.78 / 100% / 0.81** |
| M3E_det | 1.73 / 96% / 0.74 | 0.77 / 78% / 0.41 |
| HMMr | 1.32 / 92% / **0.88** | 1.39 / 93% / 0.55 |
| HYB | 1.05 / 93% / 0.71 | 1.28 / 92% / 0.51 |

(Filtered states give the same picture.) **The regime models do not make the block-level tilt indices more reliable than M3's.** M3's sympathetic index is already close to the heart-rate change itself (SNR 1.78 vs 1.88, ICC 0.81) and has the expected sign in every block; the regime models take part of that change into the regime target and leave a smaller, noisier state deviation (sympathetic ICC 0.55). The vagal index is equally reliable in M3 and HMMr (0.86 vs 0.88). So the benefit of the regime models is in structure (timing, calibration, simulation accuracy, residuals), not in the reliability of "state change during a manoeuvre" as an index. A reliability test of an index built from the regime parameters themselves (target amplitudes, regime duration) is the fair next step for them. Caveats: 10 people, three kinds of manoeuvre mixed.

## Round 3: frozen validation on untouched exercise recordings

`validation/` holds a protocol written before the run (`protocol.md`, with md5 hashes of the frozen code), the run on the 109 exercise recordings no earlier analysis used, and the report (`report.md`). Pre-specified endpoints, paired bootstrap, Holm adjustment.

**Result: superiority in residual structure is established** (PACF violations -5.1 per recording [-5.8, -4.5]; KS pass +13.8 pp [+4.6, +22.9]; HMMg), and HMMg is non-inferior in physiological ordering and signs. **The forward-score endpoint fails** (mean -0.036 nats/beat [-0.094, 0.002] against a margin of -0.01; median 0, better in 50%, driven by two recordings: -2.53 and -0.91), so the protocol's verdict is "not established" for HMMg; HYB additionally fails the expected-signs margin (-4.6 pp [-9.2, -0.9]). Details and the exploratory (post hoc) sensitivity are in `validation/report.md`.

## Where this leaves the candidates

| Use | Model |
|---|---|
| Online or forward prediction, event timing | HMMr (reciprocity prior, gate, online regime inference) |
| Retrospective states with calibrated intervals | HYB (smoothed regime episodes -> M3E); best KS and calibration, worse forward score |
| Exercise with the profile detector already in the pipeline | M3E_det remains the strongest (-30 vagal) in exercise-like recordings |
| Block-level amplitude of a manoeuvre as a person index | M3 (as reliable or more than any regime model, round 2d) |

## Files

`lib2.R` (loads `../lib.R`; regime smoothing; recovery-regime and chain fits), `reg2_filter.cpp` (oscillator, RTS, soft-input smoother), `reg3_filter.cpp` (regime-specific rates), `sim.R` (simulations), `real.R` (real recordings), `summarize.R [sim|real]`, `figure.R`, `rounds/` (ablation scripts for the table above), `_out/` (results, figure). Fit options added to `../lib.R`: `mu_sd`, `mu_max`, `anchor_sd`, `osc`, `recip`, `multi`.
