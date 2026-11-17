# M4f — RSA amplitude coupled to the vagal level (design)

Status: design and stage-0 check only (2026-11-17). Nothing is implemented. Nothing is adopted.

## Why

M4e confirmed one structural limit: the RR series alone does not identify how an RR step splits between the vagal and the sympathetic drive (exercise vagal state; about 5% catastrophic fits).

A second observable could separate the two drives. Vagal withdrawal should shrink the respiratory oscillation of the heart rate, and a sympathetic change should not. If that holds, an RR step with a falling RSA amplitude is vagal, and one with an unchanged amplitude is sympathetic.

## Model (M4f)

M4f starts from the M4e RB component: the two-regime HMM (reg4 filter) with the respiratory oscillator. The rate is ν0 − x_p + x_s − z1, so a lower vagal level μ_p means a faster rate.

- **Coupling (C1, primary).** The oscillator's stationary variance depends on the regime through the vagal shift:
  lo(s) = lo0 · exp(2γ · (μ_p(s) − μ_p(0))), with γ ≥ 0.
  The RSA amplitude ratio between regimes is exp(γ Δμ_p). There is no coupling to μ_s.
- **Frequency.** Each regime gets its own oscillator frequency on the O grid (respiration changes with exercise).
- **Calibration of γ.** γ comes from the Lancaster blockade segments, which were already used in the M4e confirmation and so are calibration data here. The prior is γ ~ N(γ̂, 1) truncated at 0. γ cannot be estimated freely for each recording, because it is confounded with the split it is meant to resolve.
- **Continuous variant (C2, secondary).** The oscillator noise is scaled by exp(2γ x̂_p) using the predicted vagal state (a plug-in, approximate filter). It is kept only if C1 passes.

## Simulator and data

**Simulator.** sim_regime gets a coupling option: lo(t) = lo0 · exp(2γ_true · m_p(t)). The truth includes γ_true ∈ {0, 1.5, 3}.

- γ_true = 0 is the current simulator. Here the coupled model is misspecified, and it must not lose to M3 or M4e.
- An M4f gain at γ_true > 0 is expected by construction. On its own it is not evidence.

**Fresh data.** The M4e confirmation data are used up.

| | Development | Confirmation |
|---|---|---|
| Simulation seeds | "m4f-dev" | "m4f-conf" |
| Exercise | the M4e dev sample (seed 21) | 60 of the 143 never-used recordings (seed 23) |
| Fantasia windows | A, B | E (4800–5400 s), F (5400–6000 s); all recordings ≥ 6233 s |
| Blockade | calibration only (no test) | — |

**Rules.** S1', S2, S3, R1–R4 and R4' (now pre-registered), on γ_true = 0 and γ_true = 3 separately. B' becomes descriptive only. Real exercise has no ground truth, so the inferred vagal share of the heart-rate rise is reported, not scored.

## Stage 0: does the identifying assumption hold in the data? (post hoc, development data only)

This check was run before any implementation. It used the respiratory-band envelope of the instantaneous heart rate in Hz (FFT band-pass + Hilbert transform, 20 s detrend).

The gate is post hoc: it was defined after the blockade and first exercise numbers below had been seen.

**Blockade.** Amplitude ratio and rate change against saline on the same day, with n per breathing mode:

| Condition | Amplitude ratio (spont / paced) | Δrate (Hz) | n (spont / paced) |
|---|---|---|---|
| Atropine high | 0.10 / 0.09 | +0.77 | 6 / 7 |
| Double | 0.10 / 0.13 | +0.47 to +0.52 | 6 / 7 |
| Propranolol | 0.97 / 1.02 | −0.06 | 7 / 7 |
| Atropine low | 0.70 / 0.74 | about 0 | 5 / 5 |

- Atropine high and double (vagal removed): the amplitude collapses.
- Propranolol (sympathetic removed): the amplitude is unchanged.
- Atropine low: the amplitude drops by about 30% with no rate change. A level-only coupling cannot produce this.
- Implied γ̂ ≈ −ln(0.095) / 0.77 ≈ **3.1 Hz⁻¹**. This assumes the whole atropine shift is vagal and that the relation is log-linear up to full blockade.

**Step-test exercise (M4e dev sample, 60 recordings).** Resting heart rate is 67 bpm.

| Window | ΔHR | Amplitude ratio, median [IQR] (0.15–0.5 Hz) | Amplitude ratio, median (0.15–0.8 Hz) | Share < 0.5 |
|---|---|---|---|---|
| Early (10–40 s) | +16 bpm | 0.90 [0.69–1.36] | 0.96 | 0.12–0.15 |
| End (last 30 s) | +37 bpm | 0.79 [0.55–1.25] | 0.89 | 0.15–0.20 |

- The amplitude does not track the heart-rate change (Spearman ρ between −0.07 and 0.09).
- The respiratory peak moves 0.22 → 0.17 Hz early in exercise and returns to 0.23–0.25 Hz at the end.

**Reading.**

- The assumption holds under pharmacology but fails in the step test.
- With γ̂ ≈ 3, the end-of-exercise ratios imply a vagal shift of 0.04–0.08 Hz. That is 6–13% of the heart-rate rise, against about 40% in the simulator's exercise profile and the usual physiological reading of a rise from 67 to 104 bpm.
- So M4f would move real exercise steps to the sympathetic drive. The reason would be the assumption, not the data. Two candidate explanations cannot be told apart without a ground truth:
  - larger tidal volume or mechanical (non-vagal) RSA during stepping;
  - a genuinely small vagal share.
- In simulation, M4f wins only when the truth is built with the coupling. That result would be circular.

**Gate G0 (post hoc): not met.** G0 required a median end-of-exercise ratio ≤ 0.6, about what γ̂ = 3 predicts for the simulator's typical vagal share (exp(−3 × 0.24) ≈ 0.49). The observed values are 0.79–0.89.

## Options

1. **Stop the M4 line here (recommended).** Keep M3. M4/M4e/M4f document why: M3 is not beaten on parameters and states at the same time, and the vagal/sympathetic split needs information the RR does not carry.
2. **Build M4f as a drug/posture model, not an M3 replacement.** Scope it to conditions where the coupling is supported, such as blockade or tilt. No fresh blockade data are public (7 Lancaster subjects, all used).
3. **Build M4f anyway and report its exercise split as a finding.** This would be exploratory and could not be confirmed.

Script: `stage0.R` (`Rscript experiments/m4/f/stage0.R`, about 1 min, no outputs written).
