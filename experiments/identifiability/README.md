# Where can the vagal / sympathetic split be recovered? (identifiability boundary by simulation)

Experimental, 2026-11-10. Nothing here touches the pipeline, the manuscript or the app. Question: for which time-scale ratios, variance ratios and recording lengths does each specification of the paper (M0 structural constraints only, M1 + time-scale priors, M2 + spectral anchor, M3 = both) recover the two latent drives? Links the paper's identifiability result with the theorem of the STATE preprint (isotropic dynamics are not separable; distinct time scales or coupling restore observability).

`sim.R` (990 simulated recordings: 33 cells x 30, the production generator, tau_p = 4 s, sigma^2 and lambda_R as in the pipeline, 300 / 600 / 1800 s; time-scale axis tau_s/tau_p = 1.5 ... 100 at lambda_P = 2; variance axis lambda_P = 0.2 ... 20 at tau_s/tau_p = 30) -> `_out/sim.rds`; `summary.R` -> tables and `_out/identifiability.png`.

## Results (median over 30 simulations per cell; correlation of the filtered state with the truth)

M3, time-scale axis (lambda_P = 2):

| tau_s / tau_p | 5 min | 10 min | 30 min | one-process fits, M0 |
|---|---|---|---|---|
| 1.5-2 | 0.38-0.45 | 0.40-0.43 | 0.40-0.43 | 67-83% |
| 4 | 0.48 | 0.48 | 0.55 | 40-57% |
| 10 | 0.68 | 0.65 | 0.66 | 3-37% |
| 30 | 0.62 | 0.71 | 0.79 | 0-57% |
| 100 | 0.55 | 0.54 | 0.80 | 7-53% |

- Below a ratio of 10, the unconstrained model collapses to one drive in 66% of fits (89% with time-scale priors only, 46% with the anchor only); with ratio >= 10 and >= 10 min, 17% (M0), 23% (M1), 8% (M2). **M3 never collapses (0%)**: the priors and the anchor guard against collapse, but below a ratio of 10 they return a structured answer that is poorly correlated with the truth (0.44 against 0.71).
- Variance axis (tau_s/tau_p = 30, M3): the weaker drive is lost. lambda_P = 20 (vagal dominant): correlation of the sympathetic-like state 0.25 (5 min) - 0.38 (30 min); lambda_P = 0.2: vagal-like 0.57-0.63, sympathetic-like 0.81-0.88. The dominant drive is recovered at 0.85-0.90.
- The production setting (tau_p about 4 s, tau_s 100-300 s, 10 min) lies where the correlation is about 0.7; recording length matters mainly at large ratios (0.62 -> 0.79 from 5 to 30 min at ratio 30).

## Reading

The split is identifiable to a moderate degree when the time scales are separated by an order of magnitude or more and the recording is at least 10 min; the constraints do not create identifiability but prevent the collapse. Ready-to-add text: `docs/REVIEWER_FAQ.md`, K3. Caveats: one generator (the model's own), tau_p fixed at 4 s, 30 simulations per cell, noise level of the pipeline, no oscillator in the truth.
