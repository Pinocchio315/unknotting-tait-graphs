#!/usr/bin/env bash
# Rebuild everything in results/crossing_changes and results/open23 from the reference table
# (results/u_table_2026-09-07.json, written by consolidate_results.py).  About 40 minutes on 8 cores.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PY:-python3}"; N="${N:-8}"
W=../results/crossing_changes/work; O=../results/open23/work
rm -rf "$W" "$O"; mkdir -p "$W" "$O"
# 1. the dataset of single crossing changes in the KnotInfo diagrams of the alternating knots with known u <= 4
for i in $(seq 0 $((N-1))); do "$PY" crossing_dataset.py "$W/rows_$i.jsonl" "$i/$N" & done; wait
for i in $(seq 0 $((N-1))); do "$PY" augment.py "$i/$N" & done; wait          # signatures, effective resistances
"$PY" identify_rest.py > "$W/identify_rest.log"                                # composites and remaining identifications
"$PY" relabel.py > ../results/crossing_changes/relabel.txt                     # -> dataset_v2.json, bj_diagram_check.json
"$PY" analyze_v2.py > ../results/crossing_changes/analyze_v2.txt               # the signature criterion and the predictor
gzip -f ../results/crossing_changes/dataset_v2.json
# 2. the open alternating [2,3] knots: candidate crossings, child knots, flype-orbit validation
for i in $(seq 0 $((N-1))); do "$PY" apply833.py "$i/$N" & done; wait
"$PY" rank833.py > ../results/open23/rank833.txt
for i in $(seq 0 $((N-1))); do "$PY" flype727.py "$i/$N" & done; wait
cat "$O"/flype727_*.jsonl > ../results/open23/flype_orbit_check.jsonl
"$PY" merge_tiers.py > ../results/open23/merge_tiers.txt
"$PY" resolve_unres.py > ../results/open23/resolve_unres.txt
gzip -f -c "$O"/rows833_*.jsonl > ../results/open23/rows833.jsonl.gz 2>/dev/null || cat "$O"/rows833_*.jsonl | gzip > ../results/open23/rows833.jsonl.gz
echo "run_all: done"
