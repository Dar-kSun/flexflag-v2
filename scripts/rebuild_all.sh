#!/usr/bin/env bash
# Rebuild every dataset and result from the download cache, in order.
#   bash scripts/rebuild_all.sh        (from the repo root; ~1 h with a warm cache)
# Stops at the first failing step. Log: data/cache/rebuild.log
set -euo pipefail
export PYTHONUTF8=1
PY=.venv/Scripts/python.exe; [ -x "$PY" ] || PY=.venv/bin/python
LOG=data/cache/rebuild.log
: > "$LOG"
run() { echo "=== $* ($(date +%T))" | tee -a "$LOG"; "$PY" -u "$@" >> "$LOG" 2>&1; }

# Caches that depend on the labels and must be recomputed.
rm -f results/external/pocket.csv results/robustness/per_protein.csv results/robustness/pockets.json

run scripts/01_build_dataset.py --workers 24
run scripts/02_train_eval.py
run scripts/03_pocket.py
run scripts/04_pocket_model.py
run scripts/05_external_dataset.py --workers 24
run scripts/06_external_eval.py
run scripts/07_export_rule.py
run scripts/08_robustness.py
run scripts/09_retest.py --workers 24
run scripts/10_p2rank.py --threads 12
run scripts/11_esm_embed.py
run scripts/12_esm_eval.py
run scripts/13_multipair.py --workers 24
run -m pytest -q
echo "=== all done ($(date +%T))" | tee -a "$LOG"
