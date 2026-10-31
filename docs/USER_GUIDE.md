# SPARK user guide: comparing people, groups and conditions

For researchers and clinical researchers who use SPARK (this repository or the web app) and want to compare recordings. Research use only: there are no
reference values and no clinical thresholds. The numbers below come from the manuscript, the supplementary material and `results/`.

## 1. Read this first

SPARK returns many quantities. They are **not equally trustworthy for comparisons**:

- some describe the recording and are stable (amplitude and baseline rate);
- some are strongly shaped by the priors (the slow time constant);
- the split of variability between the two drives depends on the spectral anchor.

The two latent drives are *vagal-like* and *sympathetic-like* reconstructions, conditional on the model's assumptions. They are not measurements of nerve activity.
A comparison is only as good as the quantity it uses, so choose the quantity first (section 2).

## 2. Which quantities to compare

Reliability = test–retest ICC of the log values between two non-overlapping 10-minute windows of the same resting recording (Fantasia, 20 young + 20 older
adults, same session; the pooled sample is heterogeneous, which raises ICC). "Window-to-window change" is the median relative change between those windows:
a within-person difference smaller than this is not interpretable.

| Quantity (app label) | Reliability (ICC) | Window-to-window change | How it is determined | Use |
|---|---|---|---|---|
| Baseline rate ν₀ (Overview card "Baseline rate ν₀", = ν₀ × 60; not the mean heart rate) | 0.99 | 1.5% | data | Comparable |
| Overall amplitude σ²_sys (parameter report) | 0.90 | 44% | data | Comparable (SDNN and RMSSD describe the same thing) |
| SD of x_p and of x_s (per-beat CSV) | 0.91 / 0.82 | 17% / 34% | data (amplitude) | Comparable |
| V_P / V_S (λ_P) | 0.81 | 23% | **is the spectral anchor** (from the recording's wavelet spectrum) | Comparable as a spectral summary |
| τ Parasympathetic (τ_p) | 0.88 | 15% | partly data; 81% of estimates fall outside the prior | With caution: a "fast time constant" that can absorb baroreflex-range oscillations; not a vagal constant |
| Split between drives (share of x_p variance, SD ratio, κ_P/κ_S, σ_P/σ_S) | 0.8–0.9 | variable | **conditional on the anchor and priors** (κ_P > κ_S is enforced) | Within one design only, after the sensitivity check (section 4) |
| **τ Sympathetic (τ_S)** | **0.31** | 4.7% | **97% inside the prior; mostly the prior** | **Do not compare between people** |
| Areal velocity, mean excursion index, balance σ, total drive, log balance index (Overview and Autonomic tabs) | not evaluated | not evaluated | descriptive displays | **Do not use for comparisons**: the paper does not support interpreting energy, entropy production, dissipation or allostatic load |
| Fit RMSE, gated beats, log-likelihood | – | – | quality of the fit | Quality control, not an outcome (log-likelihood only between fits of the same series) |

## 3. Three designs, three rules

**A. Different people or groups** (e.g. young vs older, patients vs controls)
- Same condition, posture and recording length for everyone (at least 10 minutes for the spectral anchor to be meaningful).
- Same settings for everyone: the defaults, with the priors fixed. Never tune priors or the anchor per person or per group.
- Compare the "Comparable" quantities of section 2. Report effect sizes with intervals, and mark the analysis as exploratory if n is small (in the paper, 20 + 20 gives only exploratory results).
- Compare SDNN and RMSSD alongside: for overall variability they carry the same group information (δ about −0.7 to −0.8 in the aging analysis). SPARK adds information when the question concerns the dynamics.
- Time constants must be read against their priors: use the *Prior and posterior summary* card. A posterior spread close to the prior spread means the prior is doing the work.

**B. The same person over time** (e.g. before and after an intervention)
- Use the same condition and the same time of day. Record two windows at baseline if possible: their difference is the person's own variability.
- A change is interpretable only if it clearly exceeds the window-to-window change in section 2 (about 15–20% for τ_p and SD of x_p, 40% for σ²_sys).
- What is supported: stability inside one session. **Stability across days or weeks has not been measured**, so do not claim longitudinal monitoring.

**C. The same person in different conditions (the strongest design)**
- This is how the paper uses the model (tilt, exercise): each event against its own baseline.
- Define rest, stressor and recovery windows (automatic detection or draw them on the tachogram; the *Phase summary* card lists the means per phase).
- Standardized change = (mean of the state in the stressor window − mean in the baseline window) / SD of the state in the baseline window; compute it from the per-beat CSV.
- Compare the order and delay of the two drives, and the response size; treat fitted step-response constants as approximate (70% of x_p and 40% of x_s fits reach a search bound in the PRCP data).
- Check that the pattern is not a drift artifact: split the baseline into two halves and apply the same contrast (in PRCP, the expected sign pattern appeared in 35% of baseline halves against 98% of tilts).

## 4. Before you trust a comparison

1. **Fit quality.** Gated beats below 5% (green chip). Cleaning on or off (off by default) must be the same for all recordings being compared.
2. **Residual diagnostics will often "fail".** In the paper, the KS test passed in only 22% of recordings and the innovation PACF shows structure (respiratory and baroreflex components are not modelled). Record the result, do not use it to exclude recordings, and read the states as conditional on this.
3. **Where the prior dominates.** Look at *Prior vs posterior* and the τ_P–τ_S regions; if τ sits at the prior, do not compare it.
4. **Sensitivity to the anchor.** Repeat the key comparison with the spectral crossover at 0.08 and 0.12 Hz and the vagal floor at 0.1 and 0.5 (Settings → Model). Direction and time course were robust in the paper (trajectory correlation ≥ 0.99), but the vagal variance share moved by up to 0.23. If your conclusion depends on the split between the drives and does not survive this, report it as conditional.
5. **Same software state.** Export the session JSON or the R script that reproduces the fit, and report the app version and settings.

## 5. What to report (minimum set)

- Recording length, condition and posture; cleaning on/off; settings (default or which ones changed).
- The comparable quantities used (section 2), with their intervals.
- For designs B and C: the baseline window used and the within-person variability.
- Fit quality: gated beats and the diagnostic chips.
- The sensitivity check, if the conclusion involves the split between the drives.
- The wording "vagal-like / sympathetic-like reconstructed states", not "vagal / sympathetic activity".

## 6. Things SPARK does not give you

- Reference values or thresholds; a clinical diagnosis; validation against nerve activity or pharmacological blockade.
- Independent estimates of the slow (sympathetic) time constant.
- Smoothed estimates: the states are causal (filtered), so they lag sudden changes. In the Live section, states are provisional until about 5 minutes of recording.
- Results from populations with overt cardiac autonomic disease (not tested).

## 7. Where to find the details

- Manuscript: *Identification structure* (what is fixed by assumption), *Identifiability* results (what RR alone cannot determine), *Limitations*.
- Supplement: *Spectral anchoring*, *Identifiability and sensitivity analyses: details*, *Additional analyses of the Fantasia recordings* (stability).
- Reproducing a fit: `README.md` of this repository, and the exported R script of the app.
