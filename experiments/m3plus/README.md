# M3 extensions (oscillator, smoother, inputs, ...): code, and the go / no-go for the respiratory oscillator

Experimental. Nothing here touches the pipeline, the manuscript or the app. The earlier results (2026-11-05, in the project memory and in `docs/REVIEWER_FAQ.md` B8) are not repeated; this file records the decision exercise of 2026-11-10.

## Does M3 with the respiratory oscillator (M3O) become the default? (`compare/go_nogo.R`)

The decision rule used for the log-RR candidate (written in the script header), applied to M3O against production M3, on filtered states (what production reports), with the stored simulations (80 per scenario) and the stored real fits (313), plus new simulations for exercise and tilt (60 each) and the Fantasia test-retest (windows A and B).

| Criterion | Result |
|---|---|
| S1 accuracy (RMSE, correlation) | **passes** in all four scenario groups (stationary, rsa, exercise, tilt); e.g. rsa vagal RMSE 0.53 -> 0.50 |
| S2 calibration (95% coverage) | **fails** in exercise (vagal / sympathetic 0.70 / 0.62 -> 0.62 / 0.56) and tilt (0.83 / 0.91 -> 0.82 / 0.87); passes in stationary and rsa |
| S3 exercise / tilt, change rest -> end, estimated / true | passes: exercise vagal 1.12 -> 1.21, sympathetic 0.90 -> 0.84; tilt vagal 0.47 -> 0.60, sympathetic 1.70 -> 1.58 (closer to 1); vagal drive reaches half of its change first in 100% (truth 92-93%) |
| S4 time constants | fails in exercise (error of tau not below M3's + 0.1); passes elsewhere |
| R1 real diagnostics | passes: KS pass 21.7% -> 26.2%, PACF violations 13.0 -> 11.9; tau_p inside 3-6 s 11% -> 10% (median tau_p 15.1 -> 16.7 s) |
| R2 forward score (held-out nats per beat) | **M3O better**: +0.0205 [0.0006, 0.0347], better in 65% of 313 recordings (Fantasia +0.064, PRCP +0.056, exercise +0.013) |
| R3 structure | passes: PRCP tilts with both expected signs 100% in both; exercise responses unchanged (z_p -6.3 -> -6.5, z_s 19.8 -> 19.9) |
| R4 Fantasia test-retest and groups | amplitudes pass (ICC SD_p 0.95 -> 0.93, SD_s 0.90 -> 0.88; AUC young-vs-older SD_p 0.93 -> 0.94); **tau_p ICC falls 0.88 -> 0.62** (tau_s 0.31 -> 0.39) |

**Reading.** M3O is better where the model is stationary and where the criterion is predictive (diagnostics, held-out score, rsa), equal on the physiological summaries, and worse on two things that matter for the paper: the coverage of its bands in exercise and tilt, and the test-retest reliability of tau_p. By the rule it is not a default today; it is a well-supported optional extension (the held-out score is the strongest evidence, a lower bound of +0.0006 for the mean difference). Not tested here: its effect on the exercise recovery times and on the summaries of the manuscript (the oscillator moves about 10% of the net response from sympathetic to vagal, see the project memory), a calibrated variant, a respiratory signal as input (PRCP and Fantasia have one for some recordings).
