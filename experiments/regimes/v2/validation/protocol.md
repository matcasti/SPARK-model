# Frozen validation protocol: regime models against M3 on untouched exercise recordings

Written before the validation run (2026-11-06). Code and hyperparameters are frozen (md5 hashes in `hashes.txt`, recorded before the first run). Any change after the run is a deviation and is listed at the end of the report.

## Data

The 109 exercise recordings (`data/tmst`) that no earlier analysis of the regime models used: not in the round-1 sample (seed 5, 100 recordings), the round-2 real-data sample (seed 21, 60), nor the screening sample (seed 11 plus six named cases, 40). Ids in `ids.txt`. Exercise only: Fantasia and PRCP were used for tuning or for the gate calibration and have no clean subset.

## Models (frozen)

| Name | Definition |
|---|---|
| M3 | production model (`fit_ext` without extensions) |
| HMMg | 2-regime HMM, reciprocity prior (`recip = 20`), persistence prior 7 (sd 0.7), mu ~ N(0, 0.8^2); used only if the change-point gain is at least 9, otherwise M3. States and innovations of the forward filter |
| HYB | gated; events = episodes of the Kim-smoothed non-reference regime probability above 0.5 for at least 30 beats (at most 8), then M3E (GLS amplitudes) |
| M3E_det | M3E with events from the pipeline's detector |

Fit failures and closed gates count as M3 for that recording (intention to treat).

## Endpoints (paired over recordings; HMMg and HYB each against M3)

| | Endpoint | Rule |
|---|---|---|
| E1 | Forward score (fit first 70%, last 30% scored; nats per beat), mean difference; recordings whose tail has more than 10% of beats of 2 s or more excluded; HMMg only | non-inferiority margin -0.01; superiority if the lower bound is above 0 |
| E2 | PACF violations of the standardised innovations (mean difference) | superiority if the upper bound is below 0 |
| E3 | KS pass (difference in proportions) | superiority if the lower bound is above 0 |
| E4 | Expected signs of the rest-to-end-of-exercise change (vagal down, sympathetic up; difference in proportions) | non-inferiority margin -0.05 |
| E5 | Vagal state reaches half of its change before the sympathetic state (difference in proportions) | non-inferiority margin -0.05 |

Inference: percentile bootstrap over recordings (10,000 resamples), 95% two-sided intervals; one-sided bootstrap p-values (one per endpoint) adjusted by Holm across the endpoints of each model (five for HMMg, four for HYB).

Verdict per model: **superior** if E2 or E3 superiority holds (Holm p < 0.05) and every non-inferiority endpoint holds; **non-inferior** if only the non-inferiority endpoints hold; **not established** otherwise.

## Secondary (descriptive)

Log-likelihood gain over M3 (explicit), tau_p, gate opening, fit failures, run time, distribution of the forward difference (median, trimmed mean, worst), exercise-onset realism (t50 vagal and sympathetic).

## Not claimed

No claim outside exercise recordings, none about index reliability for manoeuvres (round 2d: M3 is as reliable or better there), and none about HYB's forward prediction (events use the whole record).
