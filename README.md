# SPARK analysis

Reproduces every result, figure and supplementary check of the manuscript from RR-interval recordings.
The estimation engine `engine/model_engine.R` is used **unmodified** (byte-identical to the original; checked by the tests).

## Quick start

```bash
Rscript tools/checks/check_env.R          # packages, C++ compiler, Python wfdb, Quarto: what is installed and what is missing
Rscript tools/download-data.R      # once: Fantasia and PRCP from PhysioNet into data-raw/ (checksum-verified; `make download`)
Rscript tools/prep-data.R          # once: raw WFDB / text files in data-raw/ -> data/<dataset>/<id>.csv (+ prep_qc.csv, prcp_events.csv)
Rscript tools/run_all.R            # the whole analysis: results/*.rds and manuscript/figures/
Rscript tests/run_tests.R          # checks (--quick: a few minutes; --e2e: also the full pipeline on synthetic data)
quarto render manuscript           # manuscript.qmd and supplementary.qmd read results/ and manuscript/figures/
```

`make help` lists the same commands as targets. Settings (data folders, windows, workload, figures) are in `config.R`; the steps, and what each reads and writes, in `pipeline.R`.

## Workflow

```
data-raw/  --tools/prep-data.R-->  data/*.csv  --analysis/01 ... 13-->  results/*.rds + manuscript/figures/  --quarto-->  manuscript
```

| Command | What it does |
| :-- | :-- |
| `Rscript tools/run_all.R --list` | the 13 steps with what each reads and writes (no run) |
| `Rscript tools/run_all.R --status` | which steps are `fresh`, `stale` (an input is newer than an output) or `missing` |
| `Rscript tools/run_all.R --stale` | runs only the steps that are stale or missing |
| `Rscript tools/run_all.R --only=05,13` | only those steps (also `--from=05`, `--to=09`, `--skip=02`, `--dry-run`) |
| `make submission` | renders the LaTeX route (profile `latex`, keeping the `.tex`) and gathers `submission/`: `manuscript.tex` and `supplementary.tex` (pdfLaTeX), `Figure_N.pdf` / `FigureS_N.pdf`, `highlights.docx`, `manuscript.pdf` and `supplementary.pdf` (both from Typst) |

Each step is a script of `analysis/` that can also be run on its own once the results it reads exist (`Rscript analysis/datasets/05_prcp.R`). After each step `run_all.R` checks
that the outputs declared in `pipeline.R` exist and were written by that run. A run leaves in `results/` its own record:

- `run_log.csv`: status, seconds and time of the last run of every step;
- `provenance.txt`: git commit and whether code changed since it, R and package versions, seed, md5 of `config.R`, the engine and the data files;
- `MANIFEST.csv`: size and md5 of every file of `results/` and `manuscript/figures/`;
- `sessionInfo.txt`.

To compare two runs (for example before and after a change of the code): `Rscript tools/checks/compare_results.R <results_a> <results_b>`.

## Layout

```
config.R      every setting, in numbered sections (what is analysed and how)
pipeline.R    the 13 steps as data: script, reads, writes (what runs, in which order, from what to what)
engine/       model_engine.R (the model, unmodified), fast_filter.cpp (verified C++ port of its filter)
analysis/     one short script per result, numbered in run order (00_setup.R is sourced by all); each reads, calls R/steps, saves
  model/        01-04: fit, recovery under known truth, identifiability, sensitivity
  datasets/     05, 06, 09: PRCP, Fantasia, exercise cohort
  validation/   07, 08, 10-12: benchmarks, event detection, HRV comparators, checks
R/            functions only, in layers (load.R loads them in this order):
  core/         logging, results store, CSV loader, statistics
  model/        model specifications M0-M3, fits, identifiability, simulator, synthetic recordings
  signals/      events, physiology, HRV, state summaries, landscapes
  steps/        the per-recording steps of each analysis script
  report/       notation, plots, figures, manuscript text
tools/        prep-data.R, download-data.R, run_all.R, make_docs.R, make_submission.R
  checks/       check_env.R, figure checks, compare_results.R, benchmark_cleaning_gate.R
tests/        run_tests.R
  suites/       the test suites
  fixtures/     test config, recorded reference values, synthetic-data generators
docs/         DATAFLOW.md: scripts, reads/writes and a function index (generated: Rscript tools/make_docs.R)
manuscript/   manuscript.qmd, supplementary.qmd, model_code.R, figures/
  templates/    Typst and LaTeX templates
  support/      citations.bib, highlights.md
data-raw/  data/  results/     inputs and outputs
```

A layer may use the layers before it, not the ones after it. `analysis/<group>/NN_*.R` scripts have the same four parts, **Load**, **Analyse every recording**, **Combine and summarise**,
**Save**, and call named functions of `R/steps/`, so a script reads as a recipe. Start with `docs/DATAFLOW.md`.

## Data preparation

The repository holds code and the one non-public input, the exercise cohort (`data-raw/tmst/`, RR series without identifiers). The public datasets and every output are not versioned:
`make download` (`Rscript tools/download-data.R`) fetches Fantasia (doi:10.13026/C2RG61) and PRCP (doi:10.13026/C2PC7R) from PhysioNet into `data-raw/fantasia/` and `data-raw/prcp/`, verifying
each file against PhysioNet's SHA-256 list (re-run it to resume; `--check` only verifies), and `make data` does that and the preparation. `data/`, `results/*.rds`, the figures and
`manuscript/_site/` are rebuilt by `make data`, `make run` and `make render`; `results/MANIFEST.csv`, `provenance.txt` and `sessionInfo.txt` are the versioned record of a run.
To prepare by hand, put the raw WFDB / text files in `data-raw/` (`fantasia/`, `prcp/`, `tmst/`) and run `Rscript tools/prep-data.R`. It needs Python's `wfdb` through `reticulate`
(`python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`; reticulate uses `.venv` automatically) and writes `data/<dataset>/<id>.csv` plus `data/prep_qc.csv`.
CSV columns: `time` (s, end of the interval), `rri` (s), and when available `systolic`, `diastolic`, `map`, `resp`, `angle`.
Check `config.R` afterwards (data folders, Fantasia group regexes `y<number>` / `o<number>`, analysis windows, workload).

## Requirements

`Rscript tools/checks/check_env.R` lists them. Required: `data.table`, `ggplot2`, `patchwork`, `ggtext`, `biwavelet`, `numDeriv`, `MASS`, `Rcpp`, `survival`, `Matrix`, `scales`, `viridisLite`;
recommended: `ragg`, `systemfonts`, `svglite` (anti-aliased PNG with full Unicode: Greek letters, sub- and superscripts); Python `wfdb` and Quarto only for the preparation and the manuscript. Use a UTF-8 locale.
A C++ compiler is optional: without it the multi-start analyses use the engine's R filter (same results, slower). On macOS the compiler must be Apple's clang: `.Renviron` puts `/usr/bin`
first in `PATH`, because a Homebrew LLVM ahead of it mixes its libc++ headers with the SDK and the build fails.

## Figures and manuscript text

- **Figures.** `Rscript analysis/13_figures.R` regenerates all figures from `results/` in seconds. Time series that share an axis are facets of one plot and plots are combined with patchwork, so panels stay aligned (`R/report/plots.R`). Their look is set in `config.R`, section 6
  (`cfg$figures`): file format and dpi, font and size, every palette, and the density-landscape settings (colour map, smoothing, flow arrows, basins).
- **Which recording to showcase.** `cfg$figures$example$prcp` and `$exercise` take `"largest"` (largest standardized response, the default),
  `"median"` or a recording id such as `"prcp02"`. The choice is recorded in `results/figure_info.rds` and quoted in the figure captions.
- **State-space landscapes.** `R/signals/landscape.R` turns the filtered states into a density, a quasi-potential U = -ln(p / p_max), basins and a flow field
  (supplementary figure S6). It is a descriptive map of where the states lie, not a fitted potential.
- **Vector figures (journal submission).** Every figure is also written as PDF (`cfg$figures$vector_format`, Cairo, fonts embedded) next to its PNG.
  The figures are checked at the journal width (`journal_width_in`, 190 mm): `13_figures.R` reports those whose text would be smaller than `min_text_pt` there.
- **Figure 1 (conceptual figure).** It is not generated by the pipeline: design it separately and save it as `manuscript/figures/fig0_schematic.png`
  (the manuscript includes that file; keep the PDF or vector original for submission).
- **Notation.** Symbols in figures are written LaTeX-style and rendered by ggtext: `tex("$x_p$ and $\\tau_s$ (s)")` gives italic variables, upright subscripts and Greek letters
  (`R/report/symbols.R`; unsupported commands raise an error). Variable names of the result tables map to labels through `SYMBOL_TEX`, `sym()` (figures) and
  `sym_tex()` (LaTeX for `knitr::kable` headers).
- **Repeated tilts.** A PRCP recording can hold several tilts. `tilt_events()` (angle) and `detect_events()` (heart rate) return every one; each is analysed against its own baseline (`tilt_windows()`, settings in `config.R`, `cfg$tilt`) and the between-recording statistics use the median over a recording's tilts.
- **Manuscript text.** Descriptive sentences of the Results, Discussion, abstract and conclusion are generated from the results by
  `R/report/manuscript_text.R` (rules and thresholds in `NARRATIVE`); the wording follows the outcome. Only the literature-citation placeholders remain manual.

## Model specifications (`R/model/model_variants.R`)

M0 structural constraints only, lambda_P free; M1 + time-scale and separation priors, lambda_P free; M2 + spectral anchor only;
M3 = time-scale + separation + anchor + jitter prior: the reported model (every analysis, table and figure uses M3; M0-M2 are the nested specifications of the identifiability analysis).
The engine (`engine/model_engine.R`) is the M3 model: its defaults are the pipeline's, notably `lower_lf_vagal_fraction = 0.25`; the manuscript copies it to
`manuscript/model_code.R` (the deposited code). `tests/suites/test_pipeline.R` checks its byte content (`tests/fixtures/engine.md5`), recorded reference values
(`tests/fixtures/golden_m3.rds`) and the compiled filter (`engine/fast_filter.cpp`) against the R filter.

## Tests

`Rscript tests/run_tests.R` runs the suites below (`--quick`: the first three only; `--e2e`: also the full pipeline on synthetic data).

| Suite | What it checks |
| :-- | :-- |
| `tests/suites/test_workflow.R` | structure and bookkeeping: R/ is laid out in the layers of `R/load.R` and no name is defined twice; `pipeline.R` lists every script and its reads / writes match the `read_result()` / `write_result()` calls in the code; `docs/DATAFLOW.md` is what the code generates; files named in the docs exist; the engine is byte-identical to `tests/fixtures/engine.md5` and to `manuscript/model_code.R` |
| `tests/suites/test_narrative.R` | the text generators: wording follows the outcome in both directions, degrades to [PENDING], no NA/Inf leaks |
| `tests/suites/test_pipeline.R` | unit checks (C++ filter vs engine, recorded reference values `golden_m3.rds`, objective equality, simulator, events, recovery definitions, HRV, covariate estimators, loader, landscapes, example selection); with `SPARK_E2E=1` also every script and figure on synthetic data (~10 min) |
| `tests/suites/test_prep.R` | runs the real `tools/prep-data.R` on synthetic WFDB records with known truth (needs Python `wfdb`) |
| `tests/suites/test_edge_cases.R` | scripts on datasets with absent covariates (no blood pressure, single-step tilt, no angle, no events) and the R chunks of the manuscript |
| `tools/checks/compare_results.R <dir_a> <dir_b>` | not a suite: numerical comparison (tolerance 1e-8) of two `results/` folders, used to verify refactors |

`tools/checks/benchmark_cleaning_gate.R` is the benchmark behind the cleaning-gate setting of `tools/prep-data.R`.
Synthetic data only exercise the code; they are not results.

## Limits to keep in mind

- Tilt angle exists only for PRCP subjects that have it; elsewhere (and for the exercise cohort) events are detected from heart rate, so
  timing relative to them is not independent of the RR series. The RR-based detection is validated against the angle where both exist.
- Blood pressure and respiration exist only for some subjects; the model does not use them (they are external comparators), and
  every analysis reports how many recordings entered.
- No subject covariates except young / older in Fantasia, and no clinical outcomes: no association analyses are implemented.
- The recovery-definition benchmark and the event-detection validation rely on simulated profiles.

## Figure fonts and symbols
Symbols in figure labels (Greek letters, true minus, operators) are drawn from a font that must contain them. `figure_font()` (R/report/symbols.R) takes the first installed family (`cfg$figures$font`, then Arial, Helvetica, Liberation Sans, DejaVu Sans, ...) that has every symbol of `figure_glyphs()`, skips families that lack some and says so; `Rscript tools/checks/check_figure_fonts.R` lists what each candidate lacks. PNG figures are drawn with ragg and vector PDFs with Cairo, where a missing glyph shows as an empty box, so a complete font matters most for the PDFs.

`Rscript tools/checks/check_figure_rendering.R` checks, in the R installation at hand, that symbols in axis titles, tick labels (discrete and continuous), facet strips, legends and titles are drawn as formatted text (also inside a patchwork) and not as literal markup such as `<i>x</i><sub>p</sub>`; it prints the package versions and, for any failure, which workarounds work. `theme_spark()` sets the markdown text element for every position of the axis text and the strips.

## PRCP protocol annotations
The PRCP records contain, besides the beat annotations, protocol annotations (`.anI`: start and end of every slow and rapid tilt, stand-up from supine, back to supine). `tools/prep-data.R` parses them into `data/prcp_events.csv` (one row per maneuver: `id, kind, onset_s, onset_end_s, return_s, return_end_s, label`; kinds `slow tilt`, `rapid tilt`, `stand-up`) and reports labels it does not recognize. `analysis/datasets/05_prcp.R` then
- keeps the tilts from the tilt angle (labelled slow or rapid from the annotations), uses the annotated times for recordings without a usable angle, and adds the stand-ups, which the table angle cannot show, as events;
- analyses every maneuver against its own baseline and lets every maneuver, not only tilts, space the next baseline (`tilt_windows()`);
- pools the tilts as before, and reports the responses by kind (`by_kind`, `kind_contrasts`; Table `tbl-prcp-kinds`, supplementary figure S7).

`cfg$tilt$use_annotations = FALSE` (or no `prcp_events.csv`) restores the earlier behaviour: events from the tilt angle and heart rate only. In the preparation, RR intervals between 0.3 and 2 s within 5 s before to 30 s after a maneuver boundary are not replaced by the cleaning gate (`PREP_SETTINGS$events`; set `extension = NULL` to disable the annotations altogether). After changing any of these, rerun `tools/prep-data.R`, `analysis/model/01_fit_all.R` and `analysis/datasets/05_prcp.R` and `13_figures.R`.

## License and citation

The code is released under the MIT license (`LICENSE`); `CITATION.cff` gives the citation metadata. The PRCP and Fantasia recordings are not redistributed here: they are downloaded from PhysioNet under their own license (Open Data Commons Attribution 1.0). The exercise-cohort RR series in `data-raw/tmst` come from the open-access study cited in the manuscript. `manuscript/support/elsevier-vancouver.csl` is the Elsevier numeric style (parent of the style of *Chaos, Solitons & Fractals*; CC BY-SA 3.0, Citation Style Language project).
