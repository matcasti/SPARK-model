# Soft anchor and anchor alternatives

Experimental. Nothing here touches the pipeline, the manuscript or the app. The earlier results (2026-11-05: soft anchors of width 0.15-1.2 around the own anchor, no anchor, oracle anchor; `run.R`, `sim.R`) are in the project memory and in `docs/REVIEWER_FAQ.md` B7. This file records the blockade-informed variant of 2026-11-10.

## A soft anchor centred on the blockade value (`blockade_prior.R`)

lambda_P free with a log-normal prior, centre and width in ln: **const_s1.4** = ln 3.5 and 1.4 (the structural estimate of `experiments/blockade/ou_fit.R` for tau_s 100-300 s and the between-subject spread; the recording's own anchor ignored); **const_s0.7** = ln 3.5, 0.7; **anchor2x_s0.6** = twice the own anchor, 0.6; **mix_s1.0** = geometric mean of the own anchor and 3.5, 1.0; **hard** = production; **free** = no anchor.

Simulations (400 recordings, the pipeline's four scenarios, truth known), median:

| Variant | abs ln error of lambda_P | cor x_p / x_s | RMSE x_s | coverage x_p / x_s |
|---|---|---|---|---|
| hard (production) | 1.12 | 0.862 / 0.585 | 1.275 | 0.93 / 0.90 |
| const_s1.4 | **0.64** | 0.864 / **0.601** | **1.164** | 0.91 / 0.85 |
| const_s0.7 | 0.61 | 0.872 / 0.577 | 1.156 | 0.92 / 0.84 |
| anchor2x_s0.6 | 0.72 | 0.870 / 0.578 | 1.168 | 0.92 / 0.87 |
| mix_s1.0 | 0.63 | 0.863 / 0.593 | 1.197 | 0.93 / 0.89 |
| free | 2.40 | 0.853 / 0.582 | 1.215 | 0.88 / 0.59 |

In simulation, any centre near 3 reduces the error of lambda_P from a factor 3 to about 1.9 and the RMSE of x_s by about 8%, at the cost of 3-6 points of coverage of x_s. This repeats the earlier finding that the own anchor is no better than a constant (here the constant is in the range of the simulated scenarios: by scenario, in_prior 0.72 / 1.20 for const_s1.4 against 0.61 / 1.37 for hard; vagal-dominant 0.30 / 1.31 against 1.90 / 1.74, symp-dominant 0.68 against 1.47, but off_prior worse, 1.30 against 0.76).

**Real data (313 recordings) does not follow the prior.** Median lambda_P ends at 0.04-0.11 (production anchor 1.49) whatever the centre, tau_p at 4.4-5.2 s (61-76% inside 3-6 s against 11%), KS pass 29-31% (22%), and the profiled log-likelihood is better by +27 to +35 nats (87-89% of recordings better by more than 2): the likelihood wants a sympathetic-dominated split (the slow drive taking the variance), which no prior of width 0.7-1.4 around 3.5 holds, and an improvement of fit of this kind is the one the earlier experiments identified as not more information. Test-retest of the reliability (Fantasia, windows A and B, ICC of the log values): amplitudes SD_p 0.95 (hard) against 0.92-0.95; **SD_s 0.90 against 0.54 for const_s1.4 and 0.72 for mix_s1.0** (0.87 for const_s0.7 and anchor2x_s0.6); tau_p 0.88 against 0.85-0.88; tau_s 0.31 against 0.13-0.39. Young-vs-older AUC of SD_p 0.93-0.94 in all.

**Reading.** A constant or blockade-informed centre improves the recovery of the split in simulation but the production recordings pull lambda_P to values 30-100 times lower than the centre, the reliability of SD_s falls with the wide prior, and the fit gain is not a gain in information. Not for production. The structural estimate would be best used as the centre of a narrow prior (0.6) or, better, replaced by an oscillator-based model (see `experiments/blockade/README.md`, step 2) before being used as an anchor.
