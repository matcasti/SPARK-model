# Second-stage hidden Markov model on the estimated states (uncertainty-aware regime detection)

Experimental, 2026-11-10. Nothing here touches the pipeline, the manuscript or the app. Question: if production M3 is left as it is and a hidden Markov chain is fitted to the estimated states afterwards, respecting their uncertainty, can it detect regime changes (rest -> tilt, rest -> exercise -> recovery) in a model-driven way, without false regimes at rest? Different from `experiments/regimes/` (there the regimes live inside the filter and M3 is refitted).

## Model (`lib.R`)

Stage 1: production M3 states x_p, x_s (filtered) and their estimation variances v_t (`X_p_var`, `X_s_var`), averaged on a 5-s grid. Stage 2: K-regime HMM, EM with a sticky Dirichlet prior on transitions, forward-backward, Viterbi. Two emissions:

- **aware (model-driven):** y_t | regime k ~ N(mu_k, V + v_t), per drive. V is the stationary variance of the OU drive from the M3 fit (V_p = lambda_P sigma^2, V_s = sigma^2), so only the regime means are free, and a state the filter knows poorly weighs little.
- **naive:** free diagonal variance per regime, uncertainty ignored (a standard Gaussian HMM).

Detection statistic: LR = 2 (logL_K=2 - logL_K=1); the gate is the 95th percentile under the null (stationary M3 simulations, no regimes: 100 recordings, `sim.R`): **aware 72, naive 199** (null medians 8.5 and 119). Number of regimes by BIC (K = 1, 2, 3). Baseline regime = the one occupying the first 15 bins (75 s). Approximations: bins are treated as independent given the regime (the slow drive is strongly autocorrelated, hence the high gate), diagonal uncertainty only (no covariance between the drives, no parameter uncertainty), and M3 shrinks the states towards zero (regime shifts are underestimated).

## Simulations with known events (`sim.R`; 60 recordings each, the generator of the M3O comparison)

| Scenario | Mode | detected | onset within 30 s | onset error (median) | end error | shift x_p / x_s, estimated / true | BIC picks K = 2 |
|---|---|---|---|---|---|---|---|
| exercise | **aware** | 100% | 100% | +12.8 s | +35 s | 1.07 / 0.60 | 0% (K = 3: rest, exercise, recovery) |
| exercise | naive | 100% | 100% | +5.7 s | +124 s | 0.59 / 0.54 | 0% |
| tilt | **aware** | 67% | 50% | +20 s | +57 s | 0.39 / 1.08 | **100%** |
| tilt | naive | 62% | 58% | +14 s | +75 s | 0.34 / 1.00 | 7% |

The aware model is slower to call the onset (the states respond with their own time constants) but returns to baseline more accurately and chooses the number of regimes correctly; the naive one chatters (K = 3 in 84% of the null recordings, 7% right in tilt).

## Real data (313 recordings; `real.R`, example `_out/hmm_examples.png`)

| | aware | naive |
|---|---|---|
| **Fantasia, rest (40)**: recordings above the null gate | **0%** | 10% |
| BIC chooses 2 or more regimes | 12% | 100% |
| transitions per hour | 6 | 30 |
| **PRCP (10)**: above the gate | 100% | 100% |
| AUC of "left the baseline regime" against tilt angle > 20 deg (K = 2) | **0.88** | 0.82 |
| tilt onsets with a transition within 30 s | 48% | 72% |
| **Exercise (263)**: above the gate | 100% | 100% |
| BIC chooses K = 3 | 95% | 100% |
| onset within 30 s of the heart-rate detector | 62% (median +22 s) | 68% (+6.5 s) |
| second transition within 60 s of its end of exercise | **57%** | 34% |

## Reading

1. **At rest, the uncertainty-aware second stage sees no regimes** (0 of 40 above the gate, 12% with BIC >= 2), whereas the standard Gaussian HMM finds regimes in every resting recording: respecting the uncertainty and the model's own variance is what makes the detection usable.
2. **It follows the perturbations**: tilt periods are separated with AUC 0.88, the exercise cohort is split into rest / exercise / recovery in 95% of the recordings, and the recovery returns to the rest regime (see the example). The price is a lag of 13-22 s in the onset and moderate hit rates (48% of tilt onsets, 62% of exercise onsets within 30 s), about double those of the in-filter regimes of `experiments/regimes/` (HMM2g: 21% of PRCP onsets) at a fraction of the cost.
3. The method is retrospective (forward-backward over the whole record, a baseline regime at the start) and descriptive: it dates changes in the reconstructed states, it does not test them against any independent event.

Not done: a calibrated gate per recording length, group-specific baselines, covariance between the drives, the smoothed states of M3+ instead of the filtered ones, a regime-level summary (duration, shift in rest SDs) as a latent index, comparison with a change-point detector on the RR itself.
