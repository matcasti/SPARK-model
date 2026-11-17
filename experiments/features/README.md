# Nonlinear and geometric HRV features as information for the anchor and the latent decomposition

Experimental, 2026-11-16. Nothing here touches the engine, the pipeline, the manuscript or the app.

**Questions.** (1) Do information domains other than the spectrum (entropy and its derivatives; Poincaré geometry and asymmetry; detrended fluctuation analysis) help to identify the spectral anchor lambda_P = V_p / V_s? (2) Do they help to identify the other parameters (tau_p, tau_s) and the latent states?

## Design (fixed before the first run)

**Theory.** Under the SPARK generative model the RR series is, to first order, a linear Gaussian process: its distribution is fixed by its second-order structure (spectrum / autocovariance), and the Kalman likelihood already uses all of it. For such a process sample entropy, permutation entropy, multiscale entropy, SD1/SD2 and DFA exponents are functionals of the spectrum, and the Poincaré asymmetry indices (Guzik, Porta) are 50% in expectation (time reversibility). By the data-processing inequality no feature can carry more information about the parameters than the likelihood. So in model-generated data the features can only help as *estimators* (a learned map from summaries to lambda_P), and any gain must be compared with a learned map from spectral summaries alone; information beyond the spectrum can only come from real data that depart from the model, and whether it bears on the vagal/sympathetic split can only be judged against an external manipulation (blockade).

**Feature sets** (`lib.R`), all from the RR series alone:
- S, spectral and time domain: wavelet anchor, its LF and HF energies, periodogram LF, HF, VLF (0.01-0.04 Hz), LF/HF, total power, mean RR, SDNN, RMSSD (log scale).
- E, entropy: sample entropy (m = 2, r = 0.2 SD), permutation entropy (orders 3 and 4, delays 1 and 2), multiscale sample entropy (scales 1-8: values, area, slope).
- G, geometric: SD1, SD2, SD1/SD2, CSI, CVI, Guzik and Porta asymmetry indices.
- F, fractal: DFA alpha1 (4-16 beats) and alpha2 (16-64 beats).

**Simulations** (`sim.R`): 1000 recordings of 420 s, half from the model (`simulate_rr`) and half with an unmodelled respiratory oscillator in the vagal channel (`simulate_rr_osc`). Parameters drawn independently: tau_p U(2.5, 10) s, tau_s logU(60, 400) s, lambda_P logU(0.25, 12), lambda_R logU(0.005, 0.08), **sigma^2 logU(1e-4, 1e-2) and nu0 U(0.8, 1.8) Hz** (the earlier anchor experiment used a fixed sigma^2, which lets total power reveal lambda_P trivially). Oscillator: f U(0.15, 0.40) Hz, damping 0.18 s^-1, variance logU(0.1, 1.5) x sigma^2.
- Q1: elastic net (glmnet, alpha 0.5, features + squares, 5-fold CV stratified by generator) predicting ln lambda_P from each feature set: const, wavelet anchor as is, S, E, G, F, E+G+F, S+E, S+G, S+F, S+E+G+F. Score: out-of-fold median |ln error| and R^2, by generator.
- Q2a: M3 refitted with lambda_P fixed at the out-of-fold prediction of S and of S+E+G+F, against the production anchor, a constant and the oracle (true lambda_P). Score: |ln error| of tau_p, tau_s; correlation, normalized RMSE and 95% coverage of x_p, x_s.
- Q2b: direct prediction of ln tau_p, ln tau_s and the true variance share from S, from S+E+G+F and from the M3 estimates plus each set (stacking): do the features add information about the time constants beyond M3?

**Real data** (`real.R`):
- R1, is there information beyond the spectrum? Features E, G, F of the 313 production analysis windows against 19 IAAFT surrogates each (same amplitude distribution and periodogram, randomized phases). Share of recordings with the real value outside the surrogate range (two-sided rank test, alpha 0.10 for 19 surrogates). Caveat: IAAFT assumes stationarity, so non-stationary recordings (exercise) can reject for that reason alone.
- R2, external test with the Lancaster blockade segments (7 volunteers): the sim-trained S and S+E+G+F anchors applied to every segment, M3 refitted, and the ratios to saline of SD_p and SD_s computed as in the paper. Criterion fixed in advance: a better anchor keeps SD_p/saline under atropine below 0.3 **and** moves SD_s/saline under atropine towards 1 (the production anchor gives 0.26; atropine does not block sympathetic outflow), while not raising SD_p under propranolol alone.

**Predictions.** P1: in simulation S+E+G+F improves on S by less than 10% (relative median |ln error|). P2: learned anchors beat the production anchor, but mostly through the training prior (compare with the constant). P3: gains in tau and states stay below the oracle's. P4: E/G/F differ from surrogates in a sizable share of real recordings. P5: no prediction on the blockade criterion.

## Files

| Step | File |
|---|---|
| Features (R + `features.cpp`: sample entropy, permutation entropy) | `lib.R` |
| Simulations, Q1 and Q2 | `sim.R` |
| Surrogates (R1) and blockade (R2) | `real.R` |
| Full second-order set and M3 likelihood profile (added after the first run) | `second_order.R` |

Outputs in `_out/` (git-ignored). Run order: `sim.R`, `real.R`, `second_order.R` (~15 min in total on 6 cores).

## Deviations from the design (found after the first run)

1. **Amplitude leak.** Even with sigma^2 randomized, amplitude features (total power, SDNN, SD1, SD2, band powers) correlate with ln lambda_eff (Spearman 0.43-0.61) because total variance = (1 + lambda) sigma^2 and the two are drawn independently: high power "reveals" high lambda through the simulation prior, not through the recording. Scale-free sets were added (`Sf`: wavelet anchor, LF/HF, relative band powers, RMSSD/SDNN; `Gf`: SD1/SD2, CSI, Guzik, Porta; E and F are scale-free already) and are the fair comparison. The pre-specified S and S+E+G+F are reported but not interpreted.
2. **IAAFT miscalibration for SD1/SD2-type features.** In model-generated (linear) series the surrogate test rejected SD1, SD2, SD1/SD2, CSI and CVI in 84% (model) and 60% (oscillator) of recordings: IAAFT alters the lag-1 structure enough for these. They are excluded from R1. The other features reject in 4-20% of simulated series (nominal 10%).
3. **Decisive checks added** (`second_order.R`): a full scale-free second-order set (A: autocorrelations at lags 1-20 + relative log power in 12 bins 0.01-0.5 Hz) and the profiled M3 log-likelihood over a 9-point grid of lambda_P 0.25-12 (L).

## Results (2026-11-16)

### Q1. Anchor: out-of-fold median |ln error| of ln lambda_eff (R^2), 1000 simulations

| Set | all | model | oscillator |
|---|---|---|---|
| production wavelet anchor | 1.087 (-0.47) | 1.032 | 1.170 |
| constant (training median) | 0.843 (0.00) | 1.017 | 0.734 |
| Sf, scale-free spectral | 0.746 (0.09) | 0.922 | 0.639 |
| A, full second order | 0.611 (0.34) | 0.693 | 0.496 |
| L, M3 likelihood profile | 0.658 (0.28) | 0.800 | 0.552 |
| E + Gf + F (entropy, geometry, DFA) | 0.550 (0.44) | 0.625 | 0.475 |
| Sf + E + Gf + F | 0.530 (0.45) | 0.616 | 0.477 |
| Sf + A + E + Gf + F | 0.512 (0.50) | 0.552 | 0.488 |
| L + E + Gf + F | 0.541 (0.40) | 0.626 | 0.481 |
| L + Sf + A + E + Gf + F | 0.535 (0.49) | 0.575 | 0.493 |
| (leaky) S / S + E + G + F | 0.475 / 0.459 | | |

- The non-spectral domains **do** carry information about lambda_P that the band ratios, a 32-number second-order description and the M3 likelihood profile do not extract with this learner: E+Gf+F alone beats all of them, and adding E+Gf+F to Sf cuts the error by 29%, to A by 12%, to L by 18%.
- The signal is in short-term correlation measures: the largest standardized coefficients are DFA alpha1 (squared), SD1/SD2, RMSSD/SDNN, sample entropy and permutation entropy (order 4). For a Gaussian process these are functionals of the autocorrelation, so the gain is best read as a *better summary* of the beat-scale correlation structure (and of the mild non-Gaussianity of integrate-and-fire intervals), not as proof of information outside second order.
- The M3 likelihood is nearly flat in lambda_P (median range 3.8 log-likelihood units over a 48-fold range; its argmax misses by e^0.88), which is why a learned summary can do better: the learner borrows the joint prior of the simulation (amortized inference). The gain is therefore conditional on the simulated parameter ranges.
- The production anchor is worse than a constant (as in `experiments/anchors`).

### Q2a. M3 refitted with each anchor (1000 simulations; medians, coverage = mean)

| Anchor | abs ln err lambda | abs ln err tau_p | abs ln err tau_s | cor x_s | RMSE x_s | coverage x_s | abs share error |
|---|---|---|---|---|---|---|---|
| wavelet (production) | 1.087 | 0.354 | 0.446 | 0.615 | 1.246 | 0.853 | 0.099 |
| constant | 0.843 | 0.344 | 0.455 | 0.627 | 1.216 | 0.865 | 0.081 |
| Sf | 0.746 | 0.348 | 0.452 | 0.629 | 1.203 | 0.864 | 0.074 |
| Sf + E + Gf + F | 0.530 | 0.355 | 0.455 | 0.629 | 1.176 | 0.881 | 0.067 |
| oracle (true lambda) | 0 | 0.357 | 0.451 | 0.637 | 1.143 | 0.904 | 0.055 |

x_p is insensitive to the anchor (cor 0.847-0.853, coverage 0.92-0.93). The best feature anchor recovers about two thirds of the oracle's gain in x_s RMSE, coverage and variance share. **The time constants do not move with the anchor at all**, not even with the true lambda.

### Q2b. Direct prediction of the time constants and the share (out-of-fold; in brackets: the 223 simulations whose truth lies inside the production priors)

- tau_p: M3 0.354 (0.190), Sf+E+Gf+F 0.203 (0.228), M3 + Sf+E+Gf+F 0.196 (0.202). The overall gain of the learners is prior mismatch (truth drawn over 2.5-10 s, M3 prior 3-6 s); inside the prior the features add nothing to M3.
- tau_s: nothing predicts it (R^2 <= 0.04 for every set; M3 0.446, best 0.443). In 420-s recordings tau_s is not identified by any domain.
- variance share: M3 0.099 (0.093), M3 + Sf+E+Gf+F 0.070 (0.060), the same gain as through a better anchor.

### R1. Information beyond the spectrum in real recordings (share outside 19 IAAFT surrogates; nominal 0.10)

| Feature | Fantasia (40) | PRCP (10) | exercise (263) | simulated, model / oscillator |
|---|---|---|---|---|
| sample entropy | 0.58 | 0.6 | 0.62 | 0.18 / 0.12 |
| permutation entropy, order 4 | 0.60 | 0.8 | 0.54 | 0.14 / 0.14 |
| MSE area | 0.69 | 0.8 | 0.39 | 0.12 / 0.06 |
| Guzik asymmetry | 0.62 | 0.9 | 0.65 | 0.04 / 0.18 |
| Porta asymmetry | 0.52 | 0.9 | 0.55 | 0.04 / 0.10 |
| DFA alpha1 / alpha2 | 0.22 / 0.05 | 0.1 / 0.1 | 0.17 / 0.21 | 0.12 / 0.04-0.12 |

Real RR series carry nonlinear, time-irreversible structure that a linear Gaussian process (and SPARK) does not represent, including the resting Fantasia recordings (where non-stationarity is least likely to drive the rejection). DFA, a second-order measure, is close to nominal.

### R2. Blockade (Lancaster, sim-trained anchors; ratio to saline, geometric mean [95% CI], spontaneous breathing)

| Anchor | atropine SD_p | atropine SD_s | double SD_s | propranolol SD_p | anchor saline -> atropine |
|---|---|---|---|---|---|
| wavelet (production) | 0.16 [0.11, 0.21] | 0.26 [0.16, 0.40] | 0.13 | 1.42 [0.96, 2.06] | 4.57 -> 1.02 |
| Sf | 0.15 [0.11, 0.20] | 0.16 [0.12, 0.20] | 0.10 | 1.45 | 1.64 -> 1.31 |
| Sf + E + Gf + F | 0.13 [0.09, 0.19] | 0.38 [0.22, 0.76] | 0.11 | 1.32 [0.90, 1.96] | 3.64 -> 0.79 |

Paced breathing, atropine SD_s: 0.36 (wavelet), 0.21 (Sf), 0.55 [0.31, 1.02] (Sf+E+Gf+F). The pre-set criterion is met only partly: with the feature anchor SD_p stays below 0.3 and propranolol does not raise it, and SD_s under atropine moves towards 1 (0.26 -> 0.38 spontaneous, 0.36 -> 0.55 paced) but stays well below it and the intervals overlap with the production anchor's. The scale-free spectral anchor alone is worse than the production anchor on this criterion.

### Predictions

- P1 (< 10% gain over the spectrum): **not supported** for the fair, scale-free comparison (29% over Sf, 12% over A); supported only for the leaky S (3%).
- P2 (learned anchors beat the production anchor, mostly through the prior): beat it, and beat the constant too, so not only the prior; part of the gain is the simulated joint prior (amortization).
- P3 (gains below the oracle's): supported; tau_p and tau_s do not change even with the oracle.
- P4 (real recordings differ from surrogates): supported for entropy and asymmetry (0.39-0.9 vs 0.04-0.20 in simulations).

## Reading

1. **Q1, anchor: yes, modestly and conditionally.** Entropy, Poincaré-ratio and DFA summaries estimate lambda_P better than the production wavelet anchor (median error factor e^0.53 = 1.7 vs e^1.09 = 3.0) and better than spectral ratios or the likelihood profile. The information is about the beat-scale correlation structure; a learned anchor depends on a simulation prior that has to be defended, and on real data it only partly passes the blockade test (7 subjects).
2. **Q2, other parameters and states: little.** A better anchor improves the sympathetic-like state and the variance share (about two thirds of the oracle's gain: RMSE x_s -6%, coverage 0.85 -> 0.88, share error -32%), but the time constants are not identified better by any domain: tau_p gains are prior mismatch and tau_s is not recoverable in 420 s from anything.
3. **The non-spectral structure of real RR is outside the model.** Entropy and asymmetry differ from linear surrogates in most real recordings; that structure (nonlinear, time-irreversible) cannot inform parameters of a linear Gaussian latent model except through a learned map, and points instead to model extensions (asymmetric or state-dependent dynamics) if it is to be used.
4. **Nothing to adopt for the submitted paper.** A learned anchor would need a defended training prior, an external validation larger than 7 subjects, and longer records to say anything about tau_s. Candidate for future work: an amortized (simulation-based) anchor from scale-free short-term features, validated on blockade data with more subjects.

## Caveats

420-s simulations only (longer records may change tau_s and the anchor results); elastic net with squares as the only learner (a more flexible learner could change the ranking of A vs E/G/F); lambda_R in real data (1e-4 to 2.3) is far wider than simulated (0.005-0.08); IAAFT assumes stationarity; blockade n = 6-7 subjects.
