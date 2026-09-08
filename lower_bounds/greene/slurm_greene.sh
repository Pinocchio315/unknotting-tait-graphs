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
# One submission runs everything: each task first sweeps its share of the controls, and only if none of
# them is obstructed does it go on to the targets.  Resubmitting the same command resumes, because a knot
# already recorded in the output file is skipped.
#
#   sbatch slurm_greene.sh                                  # controls, then targets
#   sbatch --export=ALL,SET=targets slurm_greene.sh         # targets only
#   sbatch --export=ALL,SET=controls slurm_greene.sh        # controls only
#   sbatch --export=ALL,SKIP=rank4 slurm_greene.sh          # leave the slow rank-four knots out, then
#   sbatch --export=ALL,TEST=rank4,ARRAY_SIZE=4 --array=0-3 slurm_greene.sh    # run them on their own
#   sbatch --export=ALL,TPK=1800 slurm_greene.sh            # abandon a knot after 30 minutes
#   bash slurm_greene.sh --dry-run                          # print the resolved commands only
#
# Knobs (all optional): PY python interpreter with numpy and sympy; SET all|targets|controls;
# TEST one test kind; SKIP comma-separated test kinds to leave out; TPK seconds per knot; MAXDET skip
# larger determinants; WORKERS concurrent knots inside one task (leave at 1 with --cpus-per-task=1);
# OUTDIR result directory; ARRAY_SIZE must match the --array range; PKG path to this directory.
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

if [ ! -f "${TARGETS:-$PKG/data/greene_targets.json}" ]; then
    echo "missing input file: ${TARGETS:-$PKG/data/greene_targets.json}" >&2
    echo "copy data/greene_targets.json next to this script (it is part of the package built by" >&2
    echo "make_server_package.py), or point TARGETS=/path/to/greene_targets.json at it" >&2
    exit 1
fi

SET="${SET:-all}"
TARGETS="${TARGETS:-$PKG/data/greene_targets.json}"
ARRAY_SIZE="${ARRAY_SIZE:-16}"
TASK="${SLURM_ARRAY_TASK_ID:-0}"
WORKERS="${WORKERS:-1}"
TPK="${TPK:-10800}"
OUTDIR="${OUTDIR:-$PKG/runs}"
case "$SET" in
    all)      SETS="controls targets" ;;
    targets)  SETS="targets" ;;
    controls) SETS="controls" ;;
    *) echo "SET must be all, targets or controls" >&2; exit 1 ;;
esac

build_command() {
    local which="$1"
    CMD=("$PY" "$PKG/greene_sweep.py" --targets "$TARGETS" --shard "$TASK/$ARRAY_SIZE"
         --workers "$WORKERS" --time-per-knot "$TPK"
         --out "$OUTDIR/$which/results_${which}_${TASK}_${ARRAY_SIZE}.jsonl")
    # Plain [ ... ] && ... would make this function return 1 whenever the test is false, which under
    # set -e aborts the job, so the optional flags are added with if blocks.
    if [ "$which" = controls ]; then CMD+=(--controls); fi
    if [ -n "${TEST:-}" ];    then CMD+=(--test "$TEST"); fi
    if [ -n "${SKIP:-}" ];    then CMD+=(--skip-test "$SKIP"); fi
    if [ -n "${MAXDET:-}" ];  then CMD+=(--max-det "$MAXDET"); fi
}

if [ "${1:-}" = "--dry-run" ]; then
    for which in $SETS; do build_command "$which"; printf '%q ' "${CMD[@]}"; echo; done
    exit 0
fi

cd "$PKG"
for which in $SETS; do
    mkdir -p "$OUTDIR/$which"
    build_command "$which"
    echo "package $PKG, set $which, shard $TASK/$ARRAY_SIZE, ${WORKERS} worker(s), ${TPK}s per knot"
    # The controls run first; a nonzero exit means an obstructed control, and set -e stops the task
    # before any target is computed with an implementation that is known to be wrong.
    env PYTHONUNBUFFERED=1 "${CMD[@]}"
done
