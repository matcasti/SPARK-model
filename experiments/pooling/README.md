# Empirical-Bayes partial pooling of the time constants

Experimental, 2026-11-10. Nothing here touches the pipeline, the manuscript or the app. Question: if the prior of tau_p and tau_s is estimated from the population instead of fixed (3-6 s, 100-300 s), do the latent summaries become more reliable? `run.R`: Fantasia, 40 subjects, windows A and B (600-1200 s and 1200-1800 s); M3 refitted under (i) production priors, (ii) weak priors (tau_p 0.5-60 s, tau_s 10-2000 s), (iii) empirical-Bayes priors N(mu, sd_between) for ln tau_p and ln tau_s estimated on the other half of the subjects (two folds; sd_between from the covariance of the two windows of the weak fits). Output `_out/fits.rds`, `_out/run.log`.

| ICC(3,1) of the log values between windows | production | weak | empirical Bayes |
|---|---|---|---|
| SD_p | 0.95 | 0.91 | 0.91 |
| SD_s | 0.90 | 0.88 | 0.88 |
| tau_p | 0.88 | 0.54 | 0.70 |
| tau_s | 0.31 | 0.48 | 0.54 |

Young-vs-older AUC (young larger = above 0.5): tau_p 0.14 / 0.23 / 0.14, tau_s 0.48 / 0.27 / 0.17; share of tau_p inside 3-6 s 55% / 7% / 15%, tau_s inside 100-300 s 100% / 48% / 25%. Population values from the weak fits: tau_p 3.2-4.0 s, tau_s 72-81 s.

**Reading.** The amplitudes keep their reliability whatever the prior. The production tau_p reliability (0.88) is partly the prior compressing the values (SD of ln tau_p 0.46 against 1.04 with weak priors); the production tau_s (ICC 0.31) is the prior, and with a population prior it becomes both more reliable (0.54) and informative about group (AUC 0.17: older adults have longer tau_s). The prior would have to be group-specific, and it is estimated from the same population, so this is a direction (hierarchical model), not a drop-in change. Ready-to-add text: `docs/REVIEWER_FAQ.md`, K4. Caveats: 40 subjects, one dataset, two folds, between-subject variance from a two-window covariance.
