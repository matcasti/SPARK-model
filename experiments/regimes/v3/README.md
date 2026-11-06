# Round 4: protecting the forward tail and an index from the regime parameters

Experimental, 2026-11-06. Goal: remove the two obstacles to making the gated regime model (HMMg) the default instead of M3: the forward-score tail (frozen validation, `../v2/validation/report.md`) and the lower reliability of its person-level indices (`../v2/README.md`, round 2d). Code changes after the validation freeze are additive (`eps`, `kappa` in `reg2_filter`/`reg_run`; defaults reproduce the frozen results, and the K = 1 regime filter equals M3's to 1e-13).

## 1. Forward tail

**Diagnosis.** The losing recordings (163, 164 in validation; 104, 099 in development) lose almost everything in a few beats: in 163, 93% of the tail loss is in the worst 10% of beats (minimum -127 nats in one beat), around stretches of 20 bpm or abrupt jumps that the artifact gate does not flag. The regime model's sigma^2 is 3 to 12 times smaller than M3's (the regimes absorb the non-stationarity), so an outlier beat costs it that many times more.

**Candidates, designed only on development data** (154 exercise recordings used before the validation, plus 10 PRCP; `tail_dev.R`). Forward score minus M3, mean / median / worst:

| Candidate | Exercise (153) | PRCP (10) |
|---|---|---|
| none (as validated) | -0.016 / +0.007 / -2.35 | +0.021 / +0.015 / 0.00 |
| sigma^2 inflated x2, x4, x8, x16 | +0.023 / -0.069; -0.12 / -0.28; -0.37; -0.66 | worse (-0.06 to -0.78) |
| sigma^2 floor at 0.25, 0.5, 1 times M3's | -0.08, -0.10, -0.26 (worst -5.2 to -5.9) | +0.020, -0.004, -0.100 |
| **contamination density (eps 0.01, kappa 25 or 100), same for both models** | **+0.005 / +0.008 / -0.35** | **+0.021 / +0.017 / 0.00** |

Inflating or flooring the variance is wrong (the model is not overconfident on average, only on outlier beats). A contamination density bounds the cost of an outlier beat (`logdens_rob` in `reg2_filter`) and brings the worst loss from -2.35 to -0.35. Both models are scored with the same density, so M3 gets the same protection.

**On the 109 validation recordings** (`tail_validate.R`; not independent, they had been seen under the plain density, and their worst tails motivated this check): plain density -0.036 [-0.094, 0.002] (as in the report); robust density **-0.002 [-0.016, 0.010]**, median +0.002, trimmed mean +0.001, better than M3 in 52%, worst -0.43 (was -2.53). The point estimate is at parity, but the lower bound (-0.016) is below the -0.01 margin (one-sided p = 0.12), so **non-inferiority is not established** under the robust density either; the data cannot exclude a loss of 0.01 to 0.016 nats per beat. The robust-density endpoint was declared after seeing the plain result, so even a pass would have been exploratory.

## 2. Index from the regime parameters (split-half reliability, PRCP)

Each recording is cut between the third and fourth manoeuvre; each half is fitted independently (`prcp_halves.R`). Correlation between halves across people (n = 9-10; wide intervals):

| Index | Pearson | Spearman |
|---|---|---|
| M3 block effect, vagal | 0.93 | 0.88 |
| M3 block effect, sympathetic | 0.88 | 0.89 |
| Regime model block effect, vagal | 0.84 | 0.72 |
| Regime model block effect, sympathetic | 0.34 | 0.18 |
| Regime target mu_p (vagal) | 0.08 | 0.38 |
| Regime target mu_s (sympathetic) | 0.74 | 0.00 |
| Rate change (-mu_p + mu_s) | 0.31 | 0.30 |
| Fraction of time in the regime | 0.73 | 0.55 |

**None of the regime-based indices is as reliable as M3's block effects.** mu_s has a high Pearson correlation but a null rank correlation (driven by one or two people), mu_p and the rate change are unreliable, and the fraction of time in the regime mostly reflects how many manoeuvres a half contains. The gate opened in both halves for 9 of 10 people.

## Conclusion

Neither obstacle is removed. The tail loss is mostly outlier beats and is neutralised by a robust density (parity, but non-inferiority not formally shown); the person-level indices of the regime model are less reliable than M3's, with the regime parameters being worse than its states. Given the stated priorities (informative, robust latent states, fit not worse than M3), **the regime model does not displace M3 as the default for person-level latent states**; it remains the better choice when the question is fit quality, residual structure or detecting and dating transitions.

## Files

`tail_dev.R`, `tail_validate.R`, `prcp_halves.R`, `_out/` (results).
