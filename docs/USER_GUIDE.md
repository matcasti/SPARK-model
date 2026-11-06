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
a within-person difference smaller than this is not interpretable. How to word each quantity physiologically is in section 6.

| Quantity (app label) | Reliability (ICC) | Window-to-window change | How it is determined | Use |
|---|---|---|---|---|
| Baseline rate ν₀ (Overview card "Baseline rate ν₀", = ν₀ × 60; not the mean heart rate) | 0.99 | 1.5% | data | Comparable |
| Overall amplitude σ²_sys (parameter report) | 0.90 | 44% | data | Comparable (SDNN and RMSSD describe the same thing) |
| SD of x_p and of x_s (per-beat CSV) | 0.91 / 0.82 | 17% / 34% | data (amplitude) | Comparable |
| Stationary SD of the drives, SD_S = √σ²_sys and SD_P = √(λ_P σ²_sys) (parameter table, Hz; × 60 = bpm equivalent) | 0.90 / 0.95 | 22% / 15% | data (amplitude); their ratio is the anchor | Comparable; the variability amplitudes that separate young from older adults most clearly (Cliff's δ −0.71 and −0.86); not tonic levels |
| Net tone (sympathetic minus vagal; needs the age in Settings → Data → Person) = mean heart rate − intrinsic rate of the age formula | as the heart rate (0.99) | as the heart rate | **assumption**: population formula, individual error about ±8 bpm | Only within an age-conditional design; adds nothing to the heart rate within one age and is not a measure of either branch |
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
- The wording "vagal-like / sympathetic-like reconstructed states", not "vagal / sympathetic activity" (section 6 gives the wording quantity by quantity).

## 6. How to word each quantity physiologically

Section 2 says *which* quantities to compare; this section says *what each one means* and *how to write it up*, so that a number is not given a physiological meaning that the model cannot support. The rule behind every row:

> Describe the **reconstructed state or the fitted parameter**, never the **nerve or the branch**. Words such as *tone*, *activity*, *outflow*, *reactivity* and *balance* belong to measured physiology. SPARK returns model-conditional reconstructions from RR intervals alone, with no validation against nerve activity or pharmacological blockade.

### 6.1 Quantity by quantity

| Quantity | What it means physiologically (in the model) | Do report as | Do not report as |
|---|---|---|---|
| Baseline rate ν₀ (× 60 = bpm equivalent) | The pacing rate the model assigns to the sinoatrial node when both drives sit at their baseline: the overall rate level of the recording, net of both branches. It is not the observed mean heart rate. | "Baseline pacing rate (model estimate) of [value] bpm-equivalent", next to the observed mean heart rate | "Intrinsic heart rate"; "resting heart rate" if the person was not at rest; "vagal tone" or "sympathetic tone"; "equal to the mean heart rate" |
| Mean heart rate (60 / mean RR) | The observed average rate in the window, the result of both branches, the intrinsic rate and the condition. | "Mean heart rate" | A measure of either branch |
| SD_P and SD_S (Hz; × 60 = bpm equivalent), SD of x_p and of x_s | How large the fluctuations of each reconstructed drive are around its own baseline: how much heart-rate modulation the model attributes to the fast (vagal-like) and to the slow (sympathetic-like) component. | "Variability amplitude of the vagal-like (fast) drive, SD_P = [value] Hz ([value] bpm equivalent)"; "lower in group A than in group B (exploratory; Cliff's δ [value], 95% CI [..])" | "Vagal tone" or "vagal activity"; "sympathetic outflow"; "[value] bpm of vagal tone"; reading SD_P as a nerve firing rate or as a level |
| Overall amplitude σ²_sys | The overall variability of RR that the model has to explain (global scale of both drives). Carries the same information as SDNN and RMSSD. | "Overall modelled variability (σ²_sys), in line with SDNN and RMSSD" | "Total autonomic activity", "autonomic function" or "autonomic reserve" |
| V_P / V_S (λ_P) | How the model splits the variability between the fast and the slow component. It is **set by the recording's wavelet spectrum** (the spectral anchor), not estimated independently. | "Spectral partition ratio (fast-to-slow variance), derived from the wavelet spectrum of the recording" | "Sympathovagal balance", "autonomic balance", "an LF/HF-type marker of sympathetic dominance", an independent finding about the person |
| τ Parasympathetic (τ_p) | The time scale of the fast drive: how quickly it relaxes back to baseline. Partly prior-informed (81% of estimates fall outside the prior interval) and able to absorb baroreflex-range oscillations. | "Fast time constant (τ_p), read against its prior"; "higher in group A (δ [value]); [x]% of estimates outside the prior interval" | "Vagal time constant", "parasympathetic reactivation speed", "vagal recovery rate" |
| τ Sympathetic (τ_S) | The time scale of the slow drive. 97% of estimates fall inside the prior interval, so it mostly returns the prior. | "Not interpreted; prior-dominated (reported for completeness)" | Any between-person or between-group comparison; "sympathetic time constant" or "sympathetic recovery time" |
| Split between drives (share of x_p variance, SD ratio, κ_P/κ_S, σ_P/σ_S) | How variability is allocated between the fast and the slow component. Depends on the spectral anchor and the priors, and κ_P > κ_S is enforced. | "Allocation of variability between the fast and the slow component, conditional on the spectral anchor ([sensitivity check result])", within one design | "Vagal dominance" or "sympathetic dominance"; "autonomic balance"; "the vagal share of autonomic activity" |
| Net tone (needs age) | The mean heart rate minus the intrinsic rate of a population formula for that age (118.1 − 0.57 × age; individual error about ±8 bpm). | "Heart rate relative to the age-predicted intrinsic rate (population formula, an assumption)", only in age-conditional designs | "Net sympathetic tone"; "sympathovagal balance"; "this person's measured intrinsic rate"; any individual-level reading |
| Standardized change of a state in a perturbation (design C) | Direction and relative size of the change of a reconstructed state from its own baseline, in SDs of that baseline. | "The vagal-like state decreased and the sympathetic-like state increased (standardized change against the own baseline), the expected pattern for vagal withdrawal with sympathetic engagement" | "Vagal withdrawal of [x]%" or "sympathetic activation of [y]%"; comparing the size of the two states (x_s can change by tens of baseline SDs because its baseline SD is small); "proves sympathetic activation" |
| Order and delay of responses (onset latency, relative onset lag, step-response constant) | Relative timing of the fast and slow reconstructed components around an event. Approximate: many step-response fits reach a search bound, and events detected from the same RR series are not independent evidence. | "The sympathetic-like response began [value] s after the vagal-like response (median; [x]% of fits at a search bound)" | "Neural latency" or "conduction delay"; "sympathetic latency of [value] s" as a reaction time; absolute latencies when the event was detected from the RR series itself |
| Recovery time (90% of the excursion; Kaplan–Meier median) | How long the fitted state takes to return to 90% of its excursion. About 2.3 × τ, so it reflects the fitted kinetics (shaped by the priors); censored if not reached in the window. | "Model-based recovery time of the vagal-like / sympathetic-like state (Kaplan–Meier median; [n] censored)" | "Vagal reactivation time", "autonomic recovery time", "parasympathetic recovery"; equivalence to heart-rate-recovery indices; ignoring censoring |
| Path geometry in the (x_p, x_s) plane (length, enclosed area, maximum quadratic energy E) | A descriptive shape of the response in rest-SD units. | "Descriptive geometry of the response (path length, enclosed area, maximum distance from rest)" | "Energy expenditure"; "dissipation"; "entropy production"; "allostatic load"; "adaptive capacity" |
| Areal velocity, mean excursion index, balance σ, total drive, log balance index (app displays) | Descriptive displays with no validated physiological meaning. | Not reported | Any physiological quantity |
| 95% bands of the states | Uncertainty of the filtered state given the fitted parameters; it leaves out the uncertainty of the parameters themselves. | "Filtered states with 95% intervals, conditional on the maximum a posteriori parameters" | "Confidence interval of the person's autonomic activity" |
| Fit RMSE, gated beats, KS and PACF diagnostics, log-likelihood | Quality of the in-sample fit and of the noise assumptions. | "In-sample filtered fit (RMSE [value] ms); gated beats [x]%; residual diagnostics [result]" | "Prediction accuracy"; "validation of the states"; "the model is right because the fit is good"; a reason to exclude a recording |

### 6.2 Vocabulary

| Say | Instead of |
|---|---|
| vagal-like / sympathetic-like state, drive or component | vagal / sympathetic activity, tone, outflow |
| reconstructed, model-based, model-conditional | measured, observed, recorded |
| consistent with the expected pattern | demonstrates, confirms, proves |
| conditional on the spectral anchor and the priors | independent of the settings |
| lower / higher in group A (exploratory, effect size with interval) | impaired, preserved or healthy autonomic function |
| baseline pacing rate | intrinsic heart rate (unless the age formula is stated as an assumption) |
| fast and slow time constant | vagal and sympathetic time constant |

### 6.3 Example sentences to adapt

Replace the bracketed values with your own results and delete any clause you cannot support.

- **Methods.** "Latent autonomic states were reconstructed from RR intervals with SPARK ([version], default settings: spectral crossover 0.10 Hz, vagal floor 0.25; cleaning [on/off]). The two states are vagal-like and sympathetic-like reconstructions, defined by their kinetics and spectral allocation; they are not measurements of neural activity. Gated beats were [x]% of the recording."
- **Different groups (design A).** "Compared with [group B], [group A] had a lower variability amplitude of the vagal-like drive (SD_P, median [value] versus [value] Hz; Cliff's δ = [value], 95% CI [..]; exploratory, n = [..] per group), with a similar baseline pacing rate. SDNN and RMSSD showed the same direction."
- **Same person over time (design B).** "After [intervention], the fast time constant changed by [x]%, against a window-to-window change of [y]% between two baseline windows of the same person; the change [exceeds / does not exceed] this within-person variability."
- **Perturbation (design C).** "During [tilt], the vagal-like state decreased ([value] baseline SD) and the sympathetic-like state increased ([value] baseline SD) in [n] of [N] recordings, the expected pattern of vagal withdrawal with sympathetic engagement; the same contrast between two halves of the baseline gave it in [x]%. The sympathetic-like response began [value] s after the vagal-like response ([x]% of fits at a search bound). How the change is split between the two drives is conditional on the spectral anchor ([sensitivity result])."
- **Limitation.** "The states are model-based reconstructions conditional on the stated assumptions. They have not been validated against microneurography or pharmacological blockade, and RR intervals alone do not identify the tonic level of either branch."

### 6.4 Check your text before submitting

1. Does any sentence say *tone*, *activity*, *outflow* or *balance* of a branch for a SPARK quantity? Rewrite it as a reconstructed state or a parameter.
2. Does any comparison use τ_S, or the split between drives without the sensitivity check? Remove it or mark it as conditional.
3. Does any sentence give the size of a change to one branch, or compare the sizes of the two states? Report the direction and the timing instead.
4. Is every time constant read against its prior, and is every recovery time reported with its censoring?
5. Is the heart rate (or SDNN and RMSSD) reported next to any SPARK quantity that it largely reproduces?

## 7. Things SPARK does not give you

- Reference values or thresholds; a clinical diagnosis; validation against nerve activity or pharmacological blockade.
- Independent estimates of the slow (sympathetic) time constant.
- Absolute vagal or sympathetic tone: ν₀, μ_P and μ_S enter the likelihood only through ν₀ − μ_P + μ_S, so RR intervals alone do not identify the tonic level of either branch (see `experiments/absolute/`: a tone proportional to the variability reproduces a group contrast but not individual differences and gives implausible magnitudes). The comparable quantities are the variability amplitudes, the time constants and the heart rate; a split into two levels would need pharmacological blockade.
- Smoothed estimates: the states are filtered, not smoothed (the parameters and the spectral anchor use the whole recording, so the procedure is retrospective rather than strictly causal), and they lag sudden changes. In the Live section, states are provisional until about 5 minutes of recording.
- Results from populations with overt cardiac autonomic disease (not tested).

## 8. Where to find the details

- Manuscript: *Model assumptions and identification structure* (what is fixed by assumption), *Identifiability* results (what RR alone cannot determine), *Limitations*.
- Supplement: *Spectral anchoring*, *Identifiability and sensitivity analyses: details*, *Additional analyses of the Fantasia recordings* (stability).
- Reproducing a fit: `README.md` of this repository, and the exported R script of the app.
