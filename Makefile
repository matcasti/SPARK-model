# Shortcuts for the commands of the README. Run from the project root: `make help`.
# UTF-8 for R and Quarto whatever the shell locale (see .Renviron)
export LANG := en_US.UTF-8
export LC_ALL := en_US.UTF-8

.PHONY: help env download data prep run status stale figures test test-quick test-e2e docs render submission

help:             ## list the targets
	@grep -E '^[a-z0-9-]+:.*##' $(MAKEFILE_LIST) | sed -E 's/:[^#]*## /\t/'
env:              ## check packages, compiler, Python wfdb, Quarto
	Rscript tools/check_env.R
download:         ## public datasets (Fantasia, PRCP) from PhysioNet into data-raw/, checksum-verified
	Rscript tools/download-data.R
data:             ## download + prep: everything the analysis reads
	$(MAKE) download
	$(MAKE) prep
prep:             ## raw data (data-raw/) -> CSV recordings (data/)
	Rscript tools/prep-data.R
run:              ## the whole analysis pipeline
	Rscript tools/run_all.R
status:           ## which steps are fresh, stale or missing
	Rscript tools/run_all.R --status
stale:            ## run only the steps that are stale or missing
	Rscript tools/run_all.R --stale
figures:          ## regenerate the figures from results/
	Rscript analysis/13_figures.R
test-quick:       ## structure, text generators and unit checks
	Rscript tests/run_tests.R --quick
test:             ## every suite except the end-to-end run
	Rscript tests/run_tests.R
test-e2e:         ## all suites + the full pipeline on synthetic data (~10 min)
	Rscript tests/run_tests.R --e2e
docs:             ## regenerate docs/DATAFLOW.md
	Rscript tools/make_docs.R
render:           ## manuscript and supplement
	quarto render manuscript
submission:       ## files for the journal in submission/ (tex, vector figures, highlights)
	quarto render manuscript/manuscript.qmd --profile latex --to pdf -M keep-tex:true
	Rscript tools/make_submission.R
