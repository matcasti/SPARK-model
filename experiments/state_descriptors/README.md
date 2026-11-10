# STATE descriptors computed on the M3 states

Experimental, 2026-11-10. Nothing here touches the pipeline, the manuscript or the app (v1.0.0 stays frozen). Question: the STATE theory (Castillo-Aguilar, Mabe-Castro, Núñez-Espinosa, bioRxiv 10.64898/2026.11.02.763675, section 3.3) defines perturbation-response descriptors on the latent state; can they be computed from SPARK's reconstruction, with its uncertainty, and are they trustworthy?

## Descriptors (`lib.R`; specification fixed before looking at the data)

State A_t = (x_p, x_s); rest region R_rest = disc around the rest mean in coordinates whitened by the stationary SDs of the two drives (V_p = lambda_P sigma^2, V_s = sigma^2); distance d_t = Euclidean norm, smoothed 5 s.

- **Dmax**: largest d inside the perturbation window.
- **Trec(delta, tau)**: first delay after the end of the perturbation at which d stays <= delta for tau = 15 s (STATE: earliest delay with a dwell time).
- **Preturn(H, tau)**: Trec <= H (H = 300 s exercise, 100 s tilt; NA if the follow-up is shorter than H and there was no return).
- **d_H**: distance at the horizon.
- delta = sqrt(-2 ln 0.001) = 3.717, the exact 0.999 quantile of the stationary distance of two independent OU drives (the model-based version of STATE's rule "calibrate delta from the resting regime alone", Appendix B.1).
- **Ptransition(H)**: the persistent regions other than rest are the regimes of the second-stage HMM of `experiments/hmm_states/` (aware emission with the rest-calibrated variance, number of regimes K = 1..3 by BIC). Post-perturbation state = the regime occupying the last 30 s before t_end + H; transition = it differs from the resting regime (the majority regime of the rest window); also the mean posterior probability of not being in the resting regime. Simulations add a third scenario, **shift** (the displacement of both drives persists after the perturbation: a true regime shift, true Preturn = 0).

Two estimators. *Plug-in*: the filtered mean path, scale and delta from its own rest window (the literal rule). *Draws*: 100 paths from the filter's marginal posterior N(x_t, v_t), AR(1) errors with the OU correlation of each drive (approximation: the joint filtering law is not stored), summarised by median and 90% interval, Preturn as the share of draws that returned.

## Finding 1: the scale must come from the resting regime, not from the full fit (`sim.R`, 60 simulated recordings per scenario, truth known)

With the stationary variances of the full-recording M3 fit, the estimated variance is **5.9x (V_p) and 7.1x (V_s) the true one in exercise-like recordings** (1.7x and 2.1x in tilt-like ones): sigma^2 absorbs the perturbation. Dmax then comes out 0.39 of the true value and Trec 28 s too short. A scale taken from the variance of the resting window alone is worse (the slow drive barely moves in 100 s). **M3 refitted on the resting beats alone** (the STATE rule) brings the variance ratio to 1.4-1.8.

| simulations, rest-refit scale | exercise | tilt |
|---|---|---|
| Dmax, true / estimated (median) | 17.5 / 14.0 | 5.0 / 4.4 |
| Dmax, Spearman with truth | 0.48 | 0.74 |
| 90% interval covers the true Dmax | 27% | 77% |
| Trec, true / estimated (median, s) | 69 / 96 | 3 / 5 |
| Trec, Spearman with truth | 0.68 | 0.84 |
| Trec, bias (median, s) | +30 | 0 |
| 90% interval covers the true Trec | 68% | 92% |
| share of rest samples above delta (nominal 0.1%) | 6.7% | 0.09% |

Reading: the descriptors are recoverable only in rank and order of magnitude; in strong perturbations (exercise) the recovery time is overestimated by about 30 s because the resting scale is underestimated (rest above delta 6.7% of the time) and the amplitude is attenuated; the intervals under-cover. The plug-in estimator is worse (Trec 174 s against 69 s true, 10% of recordings never return).

## Finding 2: Ptransition in simulations (60 recordings per scenario, truth known; H = 300 s exercise and shift, 120 s tilt)

| scenario | true state at the horizon | regimes found (K >= 2) | Ptransition (share of recordings) | Preturn from draws (rest scale) | Preturn plug-in |
|---|---|---|---|---|---|
| exercise (recovers) | returned in 100% | 100% | **0.00** | 0.99 | 0.85 |
| tilt (recovers) | returned in 100% | 93% | **0.12** | 0.99 | 0.68 |
| shift (persistent displacement) | not returned in 100% | 100% | **1.00** | 0.00 | 0.00 |

The regime-based Ptransition separates a persistent shift from a recovery without false alarms after exercise (12% after tilt, the HMM's occasional extra regime); the draws-based Preturn agrees; the plug-in Preturn is biased downward in recovering recordings (0.85 and 0.68 against 1.00).

## Finding 3: real recordings (`real.R`; scale and time constants from M3 refitted on the resting beats)

**Exercise, 263 recordings** (rest 342 s, follow-up after exercise 407 s median; 255 with follow-up >= 300 s):
- Dmax 18 [10, 34] stationary SDs (median [10th, 90th percentile]); Trec 173 s [80, 311]; share of draws returned within 300 s, median 0.93 [0.17, 1.00], so recordings differ a lot; plug-in returned 56% of 260. Rest exceedance of delta 13% (the simulations say this inflates Trec).
- Sensitivity: Trec (plug-in) hardly moves with delta quantile (0.99/0.999) or dwell (5-30 s): 199-215 s. Preturn depends almost entirely on H (0.10 at 120 s, 0.56 at 300 s, 0.89 at 600 s): **H must be reported with the descriptor**, as STATE says.
- Agreement with the pipeline's recovery times (frac90): Spearman 0.72 (x_p), 0.63 (x_s). Consistent, not the same quantity.
- Redundancy with static baseline HRV (H1/H2-like, descriptive; no outcome exists): Trec is nearly independent of rest HR, SDNN, RMSSD (R^2 of the ranks 0.07 on the draws, 0.03 plug-in), Dmax is not (R^2 0.35 / 0.28, mostly rest SDs: larger baseline variability means a smaller distance in whitened units).

**PRCP, 60 tilts and stand-ups in 10 people** (H = 100 s): Dmax about 7.3 for all three kinds; Trec 3 / 24 / 30 s (slow tilt / rapid tilt / stand-up); draws say all return, the plug-in 70%. Reproducibility (H3-like, ICC(1) across the two repeats of each kind, 10 people, so imprecise): Dmax draws 0.00 / 0.81 / 0.04 and plug-in 0.38 / 0.76 / 0.63; Trec draws 0.43 / 0.55 / -0.45; with the kind removed, Dmax 0.13 (draws) and 0.55 (plug-in), Trec 0.42.

**Ptransition, real** (HMM regimes, `real.R`): exercise, BIC chooses K = 3 in 261 of 263 recordings (rest, exercise, recovery-like regimes; K = 2 in 2); Ptransition at H = 300 s is **0.31** (of 255 with follow-up; mean posterior probability 0.31). Agreement with the distance-based return (draws: more than half of the draws returned): 195 of 255 recordings agree (76%, Cohen kappa about 0.41): of the 80 flagged as transitions, 40 had returned by the distance criterion (the slow drive can sit in a recovery-like regime while the distance is already within delta, or the reverse). PRCP: transition 15% (slow tilt), 10% (rapid tilt), 15% (stand-up), the same as the false-alarm rate in the tilt simulations (12%): no evidence of a persistent shift within 100 s of the end of a tilt.

## What it supports, what it does not

- Supports: the descriptors are computable with uncertainty on the same states; Trec carries something static baseline HRV does not; the dependence on H, delta and tau is large and has to be stated (the STATE message).
- Does not support: Dmax in absolute terms (attenuated), a precise Trec (interval under-covers; +30 s bias in strong perturbations), reproducibility (0.0-0.8 depending on kind and descriptor; n = 10 people).
- Ptransition: separates a persistent shift (100%) from a recovery (0% after exercise, 12% after tilt) in simulations; in real exercise 31% are in a non-rest regime at 300 s, which agrees only moderately (kappa about 0.41) with the distance criterion; the HMM's K = 3 almost always (the recovery phase gets its own regime), so a transition is "another regime than rest", not necessarily a pathological one.
- Open: calibrated HMM gate per recording length; posterior draws from the smoother with the joint covariance; a protocol with repeated exercise tests; outcome data for H4.

## Reproduce

```bash
Rscript experiments/state_descriptors/sim.R test   # unit checks of the descriptor code
Rscript experiments/state_descriptors/sim.R        # ~5 min
Rscript experiments/state_descriptors/real.R       # ~3 min; needs results/exercise.rds, prcp.rds and the fit_*.rds
```
