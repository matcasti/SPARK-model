# Empirical vagal / sympathetic shares from pharmacological blockade (Lancaster "Blocker data")

Experimental, 2026-11-10. Nothing here touches the engine, the pipeline, the manuscript or the app. Question: how do the shares of RR-interval power per frequency bin removed by atropine and by propranolol compare with the production logistic partition (the anchor lambda_P), and does another functional form describe them better?

## Data and preparation

`BlockerData.zip` (doi 10.17635/lancaster/researchdata/280, CC BY-NC, 514.8 MB, downloaded by hand because the server blocks scripted downloads; **not in git**: `_data/` is ignored). 7 healthy volunteers (6 men, 29 +- 2 y), 14 sessions (VCU 1998, Hoag, Cooke, Eckberg): day 1 saline, low-dose atropine, high-dose atropine (0.04 mg/kg), propranolol on top (0.2 mg/kg; double blockade); day 2 saline, then propranolol alone. ECG at 500 Hz, BP, CO2, nerve signal; no respiration channel. Each condition has 5 min of spontaneous breathing and 5 min at 0.25 Hz.

| Step | File |
|---|---|
| Segment table (81 segments of 300 s) read from the free-text logs, checked by hand | `prep/protocol.py` -> `_out/protocol.csv` |
| R peaks from the ECG, read from the zip (the other signals are never extracted) | `prep/ecg_to_rr.py` -> `_out/rpeaks/` |
| Spectra per segment and bin (smoothed periodogram of heart rate or ln RR, 12 bins 0.01-0.40 Hz), double detections removed | `lib.R` |
| Pooled power ratios with bootstrap over subjects | `analysis.R` |
| Functional forms with leave-one-subject-out error | `forms.R` |
| Anchor lambda_P implied on the 310 production recordings | `anchor.R` |
| The model's own spectrum fitted to the blockade spectra (Whittle) | `ou_fit.R` |
| M3 refitted with lambda_P fixed (likelihood profile), figures | `extras/` |

All 81 segments are usable. Heart rate (bpm, spontaneous): saline 47-87, high-dose atropine 93-117, double 81-99, propranolol alone 44-71. Subject 4 has no spontaneous segment under atropine or double blockade; subject 2 has no low-dose atropine. Ratios use the saline of the same day (day 1 for atropine and double; day 2 for propranolol alone); saline day 2 over day 1 ranges 0.65-1.45 per bin: that is the day-to-day noise floor.

## What the data show (ln RR, spontaneous breathing, 6 subjects; pooled ratio of power, 95% bootstrap CI)

Non-vagal share R(f) = P(atropine) / P(saline). Production assigns the sympathetic branch the weight 0.75 (1 - L(f)).

| f (Hz) | R(f), data | double blockade / saline | sympathetic residual (atropine - double) / saline, heart rate | production sympathetic weight |
|---|---|---|---|---|
| 0.014 | 0.076 [0.046, 0.117] | 0.023 | 0.19 | 0.745 |
| 0.035 | 0.076 [0.023, 0.138] | 0.037 | 0.17 | 0.734 |
| 0.049 | 0.047 [0.017, 0.093] | 0.026 | 0.10 | 0.715 |
| 0.069 | 0.011 [0.004, 0.023] | 0.010 | 0.017 | 0.644 |
| 0.089 | 0.002 [0.001, 0.005] | 0.003 | 0.002 | 0.488 |
| 0.112 | 0.001 [0.000, 0.002] | 0.001 | 0.000 | 0.250 |
| 0.224 | 0.002 [0.001, 0.004] | 0.004 | -0.007 | 0.001 |
| 0.274 | 0.008 [0.004, 0.017] | 0.007 | 0.006 | 0.000 |

- The vagal share is 0.92-0.95 at 0.014-0.035 Hz, about 0.95 at 0.05 Hz, 0.99 from 0.07 Hz up. The sympathetic share left under vagal blockade is 0.10-0.19 below 0.05 Hz and zero above 0.09 Hz. The production logistic puts the transition at 0.10 Hz and gives the sympathetic branch 0.74-0.64 of the power between 0.01 and 0.07 Hz: 10-60 times what the blockade leaves.
- Atropine alone leaves 0.1-8% of the power (by bin); with propranolol on top 0.1-4%. Propranolol alone does not remove power: it raises it (ratio to saline 1.8-2.6 at 0.09-0.17 Hz, ln RR): the effects are not additive.
- Low-dose atropine is not a partial version of the high dose (subject 1: heart rate falls, power at 0.11-0.14 Hz rises 24-30 times); it is excluded from the shares and not analysed further.

## Functional forms (log of the power ratio; RMSE in ln units, 1 = factor 2.7)

| Form | parameters | leave-one-subject-out RMSE, four settings (ln RR / heart rate, spontaneous / paced 0.25 Hz) |
|---|---|---|
| Production logistic (fixed) | 0 | 3.4-3.9 / 2.7-3.2 |
| Constant | 1 | 2.1-2.3 |
| **Free logistic in ln f** | 4 | **1.3-1.8** |
| Lorentzian + floor (R = i + a / (1 + (f/fc)^2)) | 3 | 1.4-1.9 |
| Power law + floor | 3 | 1.4-1.9 |

Free logistic (ln RR, spontaneous): R falls from 0.043 (low f) to 0.0015 (high f), centre 0.057 Hz, steepness 10 per ln f; across the four settings the centre is 0.052-0.057 Hz, the steepness 8-10, R_high 0.0015-0.005, R_low 0.04-0.26 (depends on units: in heart rate the atropine-induced tachycardia inflates the ratio, ln RR is scale free). Lorentzian + floor: corner 0.007-0.011 Hz (tau about 15-22 s), floor 0.001-0.004. With 6-7 subjects the three free forms are not separable (their differences are within the between-subject spread); the production logistic and the constant are clearly worse. None reproduces the small residual peak at the respiratory frequency (0.22-0.35 Hz, R about 0.005-0.02).

## What it does to the anchor (`anchor.R`; figure `_out/blockade_shares.png`)

The bin version of the production partition reproduces the stored lambda_P (Spearman 0.94, median ratio 0.99 on 310 recordings). With the blockade-derived weights lambda = sum (1 - R) P / sum R P:

| | production | free logistic | Lorentzian + floor |
|---|---|---|---|
| Lancaster saline spectra (6 subjects), median | 3.9 | 289 | 337 |
| Fantasia (40), median | 2.1 | 129 | 247 |
| PRCP (10) | 1.7 | 113 | 229 |
| Exercise (260) | 1.4 | 84 | 203 |

The rank order of the recordings is kept (Spearman 0.85-0.92 with the production anchor), the level is about 60-140 times higher. Refitting M3 with lambda_P fixed at 30-300 (60 recordings, `extras/profile_lambda.R`) lowers the profiled log-likelihood in PRCP (median -24 to -26) and in exercise (-1.3 to -1.6, P10 -20), is neutral in Fantasia (+0.7), and moves tau_p up (exercise 17 s at lambda 3 -> 24 s at 300; PRCP 8 -> 12 s). The likelihood of the production recordings does not support a blockade-derived lambda_P of the order of 100.

## Reading

1. **The production logistic is not supported by blockade at rest**: its crossover (0.10 Hz) and its 0.25 vagal floor are at odds with a vagal share of 0.92-0.99 over the whole 0.01-0.4 Hz range. A form with a lower transition (about 0.05 Hz) and a saturating vagal share describes the data far better; with this many subjects, 3-4 parameters are the most that is justified, and a Lorentzian + floor has the advantage of a reading in terms of the model (a slow component plus a floor).
2. **The share of band power is not V_p / V_s.** The anchor equates the ratio of the in-band energies with the ratio of the stationary variances of two OU processes; the slow process puts most of its variance below the band, so a very high in-band vagal share does not mean a variance ratio of the same size, and the likelihood says so (PRCP). Using the shares directly as a default is not supported.
3. A defensible use would be the model's own spectrum: fit the two-Lorentzian OU spectrum (plus a respiratory oscillator) to the blockade spectra, with the atropine condition constraining the slow component, and take V_p / V_s and the time scales from there. The 5-min segments resolve frequencies above about 0.01 Hz, not a tau_s of 100-300 s, so the slow component cannot be pinned down with this dataset.

## Caveats

6-7 young, mostly male subjects, supine at rest, a single dose, no respiration signal; atropine raises the heart rate from about 60 to 105 bpm (hence ln RR); blockade is not additive and incomplete blockade is possible; one simple detector for the R peaks (subject 5 day 2 had 73 RR below 0.3 s, cleaned as double detections); spectra from 300 s segments (0.0033 Hz resolution, lowest bin 0.01-0.02 Hz); propranolol-alone and low-dose atropine only described; the nerve signal (MSNA) was not used. Not tested: the direct fit of the OU spectrum, a default derived from more subjects (HMS-MIT-FFMS and the other sets need author requests). The data are CC BY-NC: a derived default shipped in the app would need that checked.


---

# Step 2: the model's own spectrum fitted to the blockade spectra (`ou_fit.R`)

The band share above is not the ratio of variances of two OU drives, so the OU spectrum itself was fitted. Per subject, three conditions together (Whittle likelihood on the periodogram of ln RR, 0.011-0.5 Hz): saline = vagal OU (tau_p, V_p) + slow OU (tau_s, V_s) + respiratory peak + floor; atropine = the slow OU + floor; double blockade = floor only. S(f) = 4 V tau / (1 + (2 pi f tau)^2). The slow component is assumed unchanged by atropine (kappa = 1); with kappa free the fit runs into its bounds (not identified), so it is fixed. Recovery on spectra drawn from the model (20 simulated subjects, median [P10, P90] / truth): tau_p 0.88 [0.64, 1.25], V_p 0.92 [0.76, 1.09], V_s 1.00 [0.90, 1.25], respiratory variance 0.95 [0.82, 1.21].

**tau_s is not identified.** The 5-min segments only reach 0.011 Hz, where a slow OU is a tail proportional to V_s / tau_s; the likelihood is flat in tau_s (median change <= 0.5 log-likelihood units from 20 to 300 s). So V_s grows with tau_s and lambda_P = V_p / V_s is reported as a function of it.

| tau_s fixed (s) | lambda_P vagal OU only, paced / spontaneous | lambda_P counting the respiratory peak as vagal, paced / spontaneous |
|---|---|---|
| 20 | 20 / 44 | 28 / 68 |
| 50 | 8.3 / 18 | 12 / 29 |
| 100 | 4.2 / 8.9 | 5.8 / 14 |
| 200 | 2.1 / 4.4 | 2.9 / 7.3 |
| 300 | 1.4 / 3.0 | 1.9 / 4.8 |

Medians over 7 (paced) and 6 (spontaneous) subjects; the between-subject range is wide (spontaneous, tau_s = 100: 2.0 to 99 across the six subjects; P10-P90 of lambda 2-61).

- **tau_p is supported by blockade.** 3.2 s [1.8, 5.2] (paced) and 4.3 s [2.3, 23] (spontaneous; one subject at the bound), inside the production prior (3-6 s). Without the respiratory peak the fit puts tau_p at 1.7 s [1.2, 2.1] and the log-likelihood is worse by 37-77 units (3 more parameters): the peak is needed, as in the M3 oscillator experiments. The M3 estimates of tau_p on the production recordings (median about 14 s in the manuscript) are not what the blockade spectra show.
- **The level of the production anchor is plausible; 100 is not.** For the production prior of tau_s (100-300 s) the structural lambda_P is 1.4-9 (vagal OU only 1.4-8.9; with the respiratory power 1.9-14): the production median (1.4-2.1) sits at the low end, which corresponds to tau_s of 200-300 s for the OU alone. The band-share value (about 100) is an artefact of comparing in-band energies, as the likelihood profile already indicated.
- **The fit is not perfect** (`_out/ou_fit_spectra.png`): the observed saline spectrum has a trough at 0.03-0.08 Hz that the sum of one vagal and one slow OU overpredicts (by about 2-5 times in the pooled plot, part of it the 0.56 bias of the mean of log periodogram values). Power below 0.05 Hz is largely removed by atropine (R = 0.05-0.08), so at rest the slow part of the spectrum belongs to the vagal pathway; the slow OU of the model is "sympathetic-like" only by name and time scale, not by what the blockade removes.

## What this suggests

1. Keep the order of magnitude of the production anchor; the structural estimate (about 2-9 for tau_s 100-300 s, wide between people) argues, if anything, for a modest increase and for reporting the anchor with an uncertainty, not for the band-share value.
2. Adding a respiratory oscillator to M3 is supported by the blockade spectra (tau_p 3-4 s with it, 1.7 s without it, 14 s in the production fits).
3. The label "sympathetic-like" for the slow drive is conditional on rest: what the blockade removes at slow frequencies is vagal. Not a change for the submitted manuscript (it states that the latent states are model-conditional reconstructions), but a point to keep in mind for the discussion.
4. A default from these data needs more subjects (HMS-MIT-FFMS and the other sets by request) and records long enough to see the 100-300 s scale.
