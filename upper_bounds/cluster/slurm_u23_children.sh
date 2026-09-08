#!/usr/bin/env bash
#SBATCH --job-name=u23-children
#SBATCH --array=0-111
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=u23_children_%A_%a.out

# u = 1 witness search for the child knots of the open alternating [2,3] knots
# (data/targets_u23_children.json, 112 non-alternating knots with u in [1,2] or [1,3], k = 1).
# A witness (one crossing change to the unknot in some diagram of the closure) proves u(child) = 1
# and hence u = 2 for every parent listed in the target entry.  One knot per array task
# (--shard i/112 with i = SLURM_ARRAY_TASK_ID); the search is single-threaded.
#
# Knobs (environment variables):
#   EXTRA      crossing budget relative to the knot (default 4)
#   MAXC       absolute cap (default 17)
#   R2         1 adds Reidemeister II additions to the move set (much larger closure)
#   PASS_EXTRA pass moves may add this many crossings (default 3)
#   TPK        seconds per knot, 0 = unlimited (default 0; the Slurm wall time is the cap)
#   OUTDIR     result directory (default single_runs/u23_children_x<EXTRA>[_r2])
#
#   sbatch slurm_u23_children.sh                                   # +4 crossings, pass/R3/flype moves
#   sbatch --export=ALL,EXTRA=3,R2=1 slurm_u23_children.sh          # +3 crossings with R2 additions
#   sbatch --array=0-19 slurm_u23_children.sh                       # only the 20 highest-priority children
#   bash slurm_u23_children.sh --dry-run                            # print the command
set -euo pipefail
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PY:-python3}"
# Under Slurm the batch script runs from a spool copy, so locate the package like slurm_enum.sh does:
# the script's own directory, the submission directory (or its server_code/ subdirectory), or the
# directory of the submitted script; PKG=/path/to/package overrides everything.
if [ -z "${PKG:-}" ]; then
    if [ -f "$SRC_DIR/enum_search.py" ]; then
        PKG="$SRC_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/enum_search.py" ]; then
        PKG="$SLURM_SUBMIT_DIR"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/server_code/enum_search.py" ]; then
        PKG="$SLURM_SUBMIT_DIR/server_code"
    elif [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
        SUBMITTED_SCRIPT="$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null | sed -n 's/^[[:space:]]*Command=//p' | head -n 1)"
        SUBMITTED_SCRIPT="${SUBMITTED_SCRIPT%% *}"
        if [ -n "$SUBMITTED_SCRIPT" ] && [ -f "$SUBMITTED_SCRIPT" ]; then
            PKG="$(cd "$(dirname "$SUBMITTED_SCRIPT")" && pwd)"
        fi
    fi
fi
if [ -z "${PKG:-}" ] || [ ! -f "$PKG/enum_search.py" ]; then
    echo "cannot locate enum_search.py: submit from the package directory or set PKG=/path/to/package" >&2
    exit 1
fi
TARGETS="${TARGETS:-$PKG/data/targets_u23_children.json}"
KNOWN="${KNOWN:-$PKG/data/known_u_codes_noDKT.json}"
EXTRA="${EXTRA:-4}"
MAXC="${MAXC:-17}"
PASS_EXTRA="${PASS_EXTRA:-3}"
TPK="${TPK:-0}"
NT="$("$PY" -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" "$TARGETS")"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
TAG="${TAG:-u23_children_x${EXTRA}${R2:+_r2}}"
OUTDIR="${OUTDIR:-${SLURM_SUBMIT_DIR:-$PWD}/single_runs/$TAG}"
mkdir -p "$OUTDIR"
CMD=("$PY" "$PKG/enum_search.py" --targets "$TARGETS" --known "$KNOWN" --k-class 1 --shard "$IDX/$NT"
     --max-crossings "$MAXC" --extra-crossings "$EXTRA" --time-per-knot "$TPK"
     --routes 1 --pass-extra "$PASS_EXTRA" --max-routes 0 --hits-per-state 0 --max-states 0
     --result-tag "${TAG}_${IDX}" --out-dir "$OUTDIR")
[ -n "${R2:-}" ] && CMD+=(--r2)
if [ "${1:-}" = "--dry-run" ]; then printf 'command:'; printf ' %q' "${CMD[@]}"; printf '\n'; exit 0; fi
echo "child search: task $IDX of $NT targets=$TARGETS budget=+$EXTRA (cap $MAXC) r2=${R2:-0} out=$OUTDIR"
exec env PYTHONUNBUFFERED=1 "${CMD[@]}"
