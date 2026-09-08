#!/usr/bin/env bash
#SBATCH --job-name=greene-sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --output=greene_%A_%a.out
#SBATCH --array=0-15
#
# Correction-term obstructions for knots whose double branched cover is an L-space (paper, Section 3).
# Sixteen array tasks each take one sixteenth of the list and write their own JSONL file, so the tasks
# never contend for a file and any task can be resubmitted on its own.  A finished knot is skipped when
# a task is rerun, so a time limit costs only the unfinished knots.
#
#   sbatch slurm_greene.sh                                  # the 2,003 open targets
#   sbatch --export=ALL,SET=controls slurm_greene.sh        # the 1,506 controls; every one must not obstruct
#   sbatch --export=ALL,TEST=rank4,ARRAY_SIZE=4 --array=0-3 slurm_greene.sh    # the slow rank-four knots
#   sbatch --export=ALL,TPK=1800 slurm_greene.sh            # abandon a knot after 30 minutes
#   bash slurm_greene.sh --dry-run                          # print the resolved command only
#
# Knobs (all optional): PY python interpreter with numpy and sympy; SET targets|controls;
# TEST u1|rank2|rank3|rank4; TPK seconds per knot; MAXDET skip larger determinants; WORKERS concurrent
# knots inside one task (leave at 1 with --cpus-per-task=1); OUTDIR result directory; ARRAY_SIZE must
# match the --array range; PKG path to this directory.
set -euo pipefail
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"

# Under Slurm the batch script runs from a spool copy, so locate the package explicitly: this script's
# own directory, the submission directory (or a lower_bounds/greene below it), or the submitted path.
if [ -z "${PKG:-}" ]; then
    if [ -f "$SRC_DIR/greene_sweep.py" ]; then
        PKG="$SRC_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/greene_sweep.py" ]; then
        PKG="$SLURM_SUBMIT_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/lower_bounds/greene/greene_sweep.py" ]; then
        PKG="$SLURM_SUBMIT_DIR/lower_bounds/greene"
    elif [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        SUBMITTED="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^[[:space:]]*Command=//p' | head -n 1)"
        SUBMITTED="${SUBMITTED%% *}"
        if [ -n "$SUBMITTED" ] && [ -f "$SUBMITTED" ]; then
            PKG="$(cd "$(dirname "$SUBMITTED")" && pwd)"
        fi
    fi
fi
if [ -z "${PKG:-}" ] || [ ! -f "$PKG/greene_sweep.py" ]; then
    echo "cannot locate greene_sweep.py: submit from this directory or set PKG=/path/to/lower_bounds/greene" >&2
    exit 1
fi

SET="${SET:-targets}"
TARGETS="${TARGETS:-$PKG/data/greene_targets.json}"
ARRAY_SIZE="${ARRAY_SIZE:-16}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
WORKERS="${WORKERS:-1}"
TPK="${TPK:-10800}"
OUTDIR="${OUTDIR:-$PKG/runs/$SET}"
mkdir -p "$OUTDIR"

CMD=("$PY" "$PKG/greene_sweep.py" --targets "$TARGETS" --shard "$TASK/$ARRAY_SIZE"
     --workers "$WORKERS" --time-per-knot "$TPK"
     --out "$OUTDIR/results_${SET}_${TASK}_${ARRAY_SIZE}.jsonl")
[ "$SET" = controls ] && CMD+=(--controls)
[ -n "${TEST:-}" ] && CMD+=(--test "$TEST")
[ -n "${MAXDET:-}" ] && CMD+=(--max-det "$MAXDET")

if [ "${1:-}" = "--dry-run" ]; then
    printf '%q ' "${CMD[@]}"; echo
    exit 0
fi

echo "package $PKG, set $SET, shard $TASK/$ARRAY_SIZE, ${WORKERS} worker(s), ${TPK}s per knot"
cd "$PKG"
exec env PYTHONUNBUFFERED=1 "${CMD[@]}"
