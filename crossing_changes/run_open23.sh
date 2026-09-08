#!/usr/bin/env bash
# Second half of run_all.sh: rebuild results/open23 (candidate crossings of the open alternating [2,3] knots) from
# the current reference table.  About 15 minutes on 8 cores.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python3}"; N="${N:-8}"
O=../results/open23/work; rm -rf "$O"; mkdir -p "$O"
for i in $(seq 0 $((N-1))); do "$PY" apply833.py "$i/$N" & done; wait
"$PY" rank833.py > ../results/open23/rank833.txt
for i in $(seq 0 $((N-1))); do "$PY" flype727.py "$i/$N" & done; wait
cat "$O"/flype727_*.jsonl > ../results/open23/flype_orbit_check.jsonl
"$PY" merge_tiers.py > ../results/open23/merge_tiers.txt
"$PY" resolve_unres.py > ../results/open23/resolve_unres.txt
cat "$O"/rows833_*.jsonl | gzip > ../results/open23/rows833.jsonl.gz
echo "run_open23: done"
