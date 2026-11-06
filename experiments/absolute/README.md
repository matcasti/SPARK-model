# Absolute autonomic levels from the fitted states and parameters

Experimental, 2026-11-06. Question: can the fitted M3 states and parameters give an absolute (comparable across people and times) sympathetic and parasympathetic level? Follows the 2026-11-05 analysis (`experiments/softanchor`): nu0, mu_p and mu_s enter the likelihood only through the combination nu0 - mu_p + mu_s, so levels need external information. Data: Fantasia (40 resting recordings, M3 fitted on the whole record and on four 30-minute windows) and 100 exercise recordings (M3E fits). Scripts: `run.R`, `analyze.R`; results in `_out/`.

## What is stable and comparable without extra assumptions

| Index (physical units) | ICC across four windows | Young vs older (whole record) |
|---|---|---|
| Mean heart rate | 0.99 | 61.8 vs 57.9 bpm, d = -0.52 |
| Vagal-like variability SD_p = sqrt(lambda_P sigma^2) | 0.86 | 0.077 vs 0.026 Hz, **d = -2.06, AUC 0.93** |
| Sympathetic-like variability SD_s = sqrt(sigma^2) | 0.83 | 0.044 vs 0.019 Hz, d = -1.41, AUC 0.85 |
| tau_p | 0.75 | 3.2 vs 5.4 s, d = +1.53 |
| tau_s | 0.21 | 165 vs 167 s, d = -0.10 |

The variability amplitudes and tau_p are stable within a person across half-hour windows and separate young from older adults clearly (with the expected direction: older adults have less variability and a slower vagal process). tau_s is not stable at this window length. These are measures of autonomic variability, not of tonic level; and the ratio SD_p / SD_s is the spectral anchor (sqrt of lambda_P), so their split inherits the anchor's accuracy.

## Net tone (sympathetic minus vagal), with an external anchor

Net tone = mean rate - intrinsic rate, with the intrinsic rate from the age formula (118.1 - 0.57 age, bpm; group midpoints 27.5 and 76.5 years because Fantasia records only the group). Result: -40.7 bpm in young adults, -16.6 in older adults (d = 2.49), the expected direction. But this separation is by construction (the formula differs by 26 bpm between the groups while the groups' heart rates are similar), and within a group the index is just the heart rate shifted by a constant (ICC 0.99): **the states and parameters add nothing to the net tone**, which needs the age and carries the formula's error (about 5-10 bpm per person).

## Separating vagal from sympathetic level: not supported

A closure is needed (tone proportional to variability). Tests on Fantasia:

- **Two coefficients** (mu_p = beta_p SD_p, mu_s = beta_s SD_s) are not identified: beta_s = -0.8 (SE 4.6), beta_p = 6.9 (SE 3.0), SD_p and SD_s correlate 0.91 (the anchor ties them).
- **Common coefficient of variation:** c = 5.7% (young 5.8%, older 4.2%), the same 4-6% as on 2026-11-05. It reproduces the group contrast in net tone (uncentred R^2 0.76) but not person-to-person differences: within each group the correlation between predicted and observed net tone is -0.20 (older) and -0.14 (young). It is a group-level fit with one parameter on two groups.
- **Implied levels are implausibly large:** young adults median vagal tone 81 bpm and sympathetic tone 46 bpm (older 27 and 20); implied intrinsic rate 93 and 62 bpm against 102 and 74 from the formula; in 25% of young adults the implied vagal tone exceeds the intrinsic rate (negative rate). Autonomic blockade studies give resting vagal tone of the order of 30-50 bpm and sympathetic tone of the order of 10-20 bpm in young supine adults (general knowledge, not from these data), so the closure overstates both, mostly the sympathetic one.
- **Exercise cohort (coherence bound):** implied intrinsic rate median 79 bpm (10th-90th percentile 64-109; 75% within 70-125 bpm); vagal withdrawal over implied resting vagal tone has median 0.08 and no recording exceeds 1. The check passes, but it is weak because the implied tone is large.

## Conclusion

The fitted states and parameters give **comparable absolute variability amplitudes (in Hz or bpm), time scales and heart rate**, stable across time and discriminating between age groups; they do **not** give absolute vagal and sympathetic tone levels. An age-based intrinsic rate turns the heart rate into a net-tone index (assumption: the formula, which contributes the between-group difference by itself), and the split into two levels rests on a closure that fails the person-level test and produces implausible magnitudes. A defensible split needs data with pharmacological blockade (atropine, propranolol) to calibrate it.

Suggested reporting if pursued: SD_p, SD_s (bpm-equivalent), tau_p and heart rate as the comparable quantities; net tone only as an explicitly age-conditional index; no absolute tone claim for either branch.
