#!/usr/bin/env bash
#SBATCH --job-name=unknot34-x3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=48:00:00
#SBATCH --output=unknot34_x3_%j.out

# Exhaustive +3-crossing campaign for the 224 knots whose current range is [3,4].
#
# Slurm:  sbatch /path/to/server_code/slurm_enum.sh
# Direct: bash   /path/to/server_code/slurm_enum.sh
# Check:  bash   /path/to/server_code/slurm_enum.sh --dry-run
#
# The defaults remove the old caps on pass routes, determinant hits and states.
# A knot can still be incomplete when its fair-share time expires; that fact is
# recorded as stop_reason="knot_time" and is never reported as exhaustion.
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"

if [ -z "${PKG:-}" ]; then
    if [ -f "$SRC_DIR/enum_search.py" ]; then
        PKG="$SRC_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/enum_search.py" ]; then
        PKG="$SLURM_SUBMIT_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/server_code/enum_search.py" ]; then
        PKG="$SLURM_SUBMIT_DIR/server_code"
    elif [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        SUBMITTED_SCRIPT="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
            | sed -n 's/^[[:space:]]*Command=//p' | head -n 1)"
        SUBMITTED_SCRIPT="${SUBMITTED_SCRIPT%% *}"
        if [ -n "$SUBMITTED_SCRIPT" ] && [ -f "$SUBMITTED_SCRIPT" ]; then
            PKG="$(cd "$(dirname "$SUBMITTED_SCRIPT")" && pwd)"
        fi
    fi
fi

if [ -z "${PKG:-}" ] || [ ! -f "$PKG/enum_search.py" ]; then
    echo "cannot locate server_code/enum_search.py" >&2
    echo "submit the original script by absolute path, or set PKG=/path/to/server_code" >&2
    exit 1
fi
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "Python interpreter not found: $PY" >&2
    exit 1
fi

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
OUTDIR="${OUTDIR:-$SUBMIT_DIR/enum34_x3_runs}"
mkdir -p "$OUTDIR"

# Keep a buffer inside the scheduler's 48-hour allocation so active workers can
# finish writing JSON certificates and partial-state summaries.
WALL_SECONDS="${WALL_SECONDS:-172800}"
RESERVE_SECONDS="${RESERVE_SECONDS:-900}"
if [ "$WALL_SECONDS" -le "$RESERVE_SECONDS" ]; then
    echo "WALL_SECONDS must be greater than RESERVE_SECONDS" >&2
    exit 1
fi
RUN_SECONDS=$((WALL_SECONDS - RESERVE_SECONDS))

if [ -n "${WORKERS:-}" ]; then
    NWORKERS="$WORKERS"
elif [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then
    NWORKERS="$SLURM_CPUS_PER_TASK"
else
    NWORKERS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
    DIRECT_WORKER_CAP="${DIRECT_WORKER_CAP:-16}"
    if [ "$NWORKERS" -gt "$DIRECT_WORKER_CAP" ]; then
        NWORKERS="$DIRECT_WORKER_CAP"
    fi
fi
if [ "$NWORKERS" -lt 1 ]; then
    echo "WORKERS must be at least 1" >&2
    exit 1
fi

# Array submission remains supported, although it is not needed: each array
# member receives a disjoint subset and then schedules that subset dynamically.
TASK="${SLURM_ARRAY_TASK_ID:-0}"
NTASKS="${SLURM_ARRAY_TASK_COUNT:-1}"
SHARD="$TASK/$NTASKS"
if [ "$TASK" -ge 224 ]; then
    echo "array index $TASK has no target (there are 224 knots)"
    exit 0
fi
TARGETS_THIS_JOB=$(((223 - TASK) / NTASKS + 1))
if [ "$NWORKERS" -gt "$TARGETS_THIS_JOB" ]; then
    NWORKERS="$TARGETS_THIS_JOB"
fi

# Fair-share cap: 90% of the usable CPU time divided by the maximum number of
# knots assigned to one worker.  Fast knots return their unused time naturally.
KNOTS_PER_WORKER=$(((TARGETS_THIS_JOB + NWORKERS - 1) / NWORKERS))
TPK="${TPK:-$((RUN_SECONDS * 9 / 10 / KNOTS_PER_WORKER))}"

MAXC="${MAXC:-16}"
EXTRA="${EXTRA:-3}"
ROUTES="${ROUTES:-1,2}"
PASS_EXTRA="${PASS_EXTRA:-3}"
MAX_ROUTES="${MAX_ROUTES:-0}"
HITS_PER_STATE="${HITS_PER_STATE:-0}"
MAX_STATES="${MAX_STATES:-0}"
RESULT_TAG="${RESULT_TAG:-34_x3_${TASK}_${NTASKS}}"

echo "[3,4] campaign: targets=$TARGETS_THIS_JOB workers=$NWORKERS shard=$SHARD"
echo "budget=+${EXTRA} (absolute cap $MAXC), wall=${RUN_SECONDS}s, time/knot=${TPK}s"
echo "routes=$ROUTES pass_routes=$MAX_ROUTES determinant_hits=$HITS_PER_STATE states=$MAX_STATES"
echo "output=$OUTDIR"

CMD=("$PY" "$PKG/enum_search.py" \
    --shard "$SHARD" \
    --k-class 3 \
    --max-crossings "$MAXC" \
    --extra-crossings "$EXTRA" \
    --time-per-knot "$TPK" \
    --wall-time "$RUN_SECONDS" \
    --workers "$NWORKERS" \
    --routes "$ROUTES" \
    --pass-extra "$PASS_EXTRA" \
    --max-routes "$MAX_ROUTES" \
    --hits-per-state "$HITS_PER_STATE" \
    --max-states "$MAX_STATES" \
    --result-tag "$RESULT_TAG" \
    --out-dir "$OUTDIR")

if [ "${1:-}" = "--dry-run" ]; then
    printf 'command:'
    printf ' %q' "${CMD[@]}"
    printf '\n'
    exit 0
fi

exec env PYTHONUNBUFFERED=1 "${CMD[@]}"
