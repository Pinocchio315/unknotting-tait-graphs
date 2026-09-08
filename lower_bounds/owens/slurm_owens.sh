#!/bin/bash
#SBATCH --job-name=owens-u2
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/owens_%A_%a.out
# Owens u=2 obstruction sweep (|sigma| = 4 alternating [2,3] knots) as a Slurm array.
#   mkdir -p logs && PROJ=$(pwd) sbatch --export=ALL,PROJ=$(pwd) --array=0-7 slurm_owens.sh 13a
# arg1 = knot-name prefix (13a / 12a / empty for all alternating [2,3] sigma4 rows)
set -euo pipefail
PREFIX=${1:-13a}
PY=${PY:-$HOME/unknot-venv/bin/python}
N=${SLURM_ARRAY_TASK_COUNT:-1}; I=${SLURM_ARRAY_TASK_ID:-0}; W=${SLURM_CPUS_PER_TASK:-4}
ROOT="${PROJ:-${SLURM_SUBMIT_DIR:-$PWD}}"
[ -f "$ROOT/run_sweep.py" ] || { echo "ERROR: owens folder not found (PROJ='${PROJ:-}')" >&2; exit 2; }
cd "$ROOT"; mkdir -p logs
exec "$PY" run_sweep.py --out "results_${PREFIX:-all}_s4_shard${I}.jsonl" --prefix "$PREFIX" --shard "$I/$N" --workers "$W"
