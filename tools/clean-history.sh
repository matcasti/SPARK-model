#!/usr/bin/env bash
# One-off: writes a copy of this repository's history without generated outputs and public raw data (not pushed anywhere).
#   1. commit your work here first (the copy is made from the committed history)
#   2. bash tools/clean-history.sh ../SPARK-clean
#   3. inspect the copy; if it is right, force-push from it (the command is printed at the end)
# Needs git-filter-repo (python3 -m pip install git-filter-repo). Delete this script once it has been used.
set -euo pipefail
out="${1:?usage: bash tools/clean-history.sh <new-empty-directory>}"
[ -e "$out" ] && { echo "$out already exists"; exit 1; }
remote="$(git remote get-url origin)"
git clone -q --no-local . "$out"
cd "$out"
du -sh .git | sed 's/^/before: /'
git filter-repo --force --invert-paths \
  --path manuscript/_site --path manuscript/figures --path manuscript/df_population_summary.rds --path data \
  --path-glob 'data-raw/prcp/*.dat' --path-glob 'data-raw/prcp/*.hea' --path-glob 'data-raw/prcp/*.anI' --path-glob 'data-raw/prcp/*.wabp' --path-glob 'data-raw/prcp/*.wqrs' \
  --path-glob 'data-raw/fantasia/*.dat' --path-glob 'data-raw/fantasia/*.hea' --path-glob 'data-raw/fantasia/*.ecg' \
  --path-glob '*.DS_Store' --path-glob 'results/*.rds' --path results/run_log.csv --path results/population_summary.csv
git reflog expire --expire=now --all
git gc --prune=now --aggressive -q
du -sh .git | sed 's/^/after:  /'
git remote add origin "$remote"
echo "Check the copy ($(pwd)): git log --oneline | head; git ls-files | wc -l"
echo "Then publish it (rewrites the remote history; collaborators must re-clone):  git push --force origin main"
