# Frozen validation: result (109 untouched exercise recordings, 2026-11-06)

Run with the code recorded in `hashes.txt` (the md5 of `run.R` was checked again after the run); analysis exactly as in `protocol.md` (`analyze.R`, output `../_out/validation_endpoints.csv`). No recording failed; the gate opened on the full record in 100% of recordings and on the head in 100%.

## Primary endpoints (difference to M3, 95% bootstrap interval, Holm-adjusted one-sided p within model)

| Model | Endpoint | Rule | Estimate [95% CI] | Margin | Holm p | Met |
|---|---|---|---|---|---|---|
| HMMg | E1 forward score, nats/beat | non-inferiority | -0.036 [-0.094, 0.002] | -0.01 | 0.85 | **no** |
| HMMg | E2 PACF violations | superiority | **-5.13** [-5.75, -4.51] | 0 | <0.0001 | yes |
| HMMg | E3 KS pass | superiority | **+13.8 pp** [+4.6, +22.9] | 0 | 0.003 | yes |
| HMMg | E4 expected signs | non-inferiority | -0.9 pp [-2.8, 0.0] | -5 pp | 0.002 | yes |
| HMMg | E5 vagal before sympathetic | non-inferiority | -0.9 pp [-2.8, 0.0] | -5 pp | 0.0004 | yes |
| HYB | E2 PACF violations | superiority | **-4.95** [-5.63, -4.29] | 0 | <0.0001 | yes |
| HYB | E3 KS pass | superiority | **+12.8 pp** [+3.7, +22.0] | 0 | 0.014 | yes |
| HYB | E4 expected signs | non-inferiority | -4.6 pp [-9.2, -0.9] | -5 pp | 0.38 | **no** |
| HYB | E5 vagal before sympathetic | non-inferiority | -1.8 pp [-4.6, 0.0] | -5 pp | 0.026 | yes |

**Verdict by the protocol: "not established" for both models.** HMMg fails only E1 (forward score); HYB fails only E4 (expected signs).

## What the result says

- **Residual structure: superiority is established, out of sample of the model's development.** In 109 recordings never used before, both regime models reduce PACF violations by about 5 per recording (13.7 to 8.6) and raise the KS pass rate from 16% to 29%. This is the same magnitude as in the development samples.
- **Physiological realism: non-inferior for HMMg** (signs and ordering within 1 pp of M3). HYB loses 4.6 pp in expected signs (interval reaches -9.2 pp), which fails the 5 pp margin.
- **Forward prediction: not established for HMMg.** Median difference to M3 is 0 nats per beat, the 10% trimmed mean is -0.003 and it is better in 50% of recordings, but the mean is -0.036 because of a heavy left tail: 163 (-2.53; 7% of its tail beats are 2 s or longer, below the exclusion threshold) and 164 (-0.91), then -0.33 and smaller. The pre-specified mean-based endpoint therefore fails. Exploratory and post hoc (not part of the verdict): without the two worst recordings the mean is -0.013 [-0.036, 0.005]; without the three worst, -0.002 [-0.013, 0.011].
- **Secondary:** log-likelihood gain over M3 is +68 (HMMg), +66 (HYB), +58 (M3E_det) per recording; tau_p 4.8, 5.9, 6.2 against 16.9 s for M3; run time about 8 s per recording for all models.

## Consequences

1. The claim that is supported by pre-specified, untouched data: **the regime models explain the exercise recordings' innovation structure better than M3 (PACF, KS), with HMMg keeping the physiological ordering and signs.**
2. The claim of at-least-equal forward prediction is **not** supported as specified: the average loses to a few recordings in which the regime model is overconfident out of sample. The next change to test is a protection against that tail (for example inflating the variance in the forward pass by the regime uncertainty, or shrinking sigma^2 towards M3's), evaluated on a fresh sample because these 109 are now used.
3. HYB should be reported as a retrospective estimator; its sign check is below M3 by a margin that cannot be excluded as meaningful.
