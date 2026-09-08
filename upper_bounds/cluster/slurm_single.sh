#!/usr/bin/env bash
#SBATCH --job-name=unknot-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=unknot_single_%j.out

# One knot, one configuration, one core: the level-by-level search of enum_search.py is
# single-threaded, so several configurations (budgets, --r2, route sets) are submitted as
# separate jobs.  Knobs (environment variables):
#   TARGETS   target list, default data/targets_11a_14.json
#   KNOWN     partner list, default data/known_u_codes_noDKT.json (KnotInfo values that come
#             from the Jones-only DKT pipeline removed; data/dkt_104.json lists them)
#   KCLASS    2 (u <= 2 witnesses) or 3 (u <= 3 witnesses)
#   EXTRA     crossing budget relative to the knot (default 4)
#   MAXC      absolute cap (default 20)
#   ROUTES    1,2 (k = 2) or 1,2,3 (k = 3)
#   R2        1 adds Reidemeister II additions to the move set (much larger closure)
#   PASS_EXTRA pass moves may add this many crossings (default 3)
#   TPK       seconds per knot, 0 = unlimited (default 0; the Slurm wall time is the cap)
#   OUTDIR    result directory
#
#   sbatch --export=ALL,TARGETS=data/targets_11a_14.json,KCLASS=2,EXTRA=4 slurm_single.sh
#   sbatch --export=ALL,TARGETS=data/targets_11a_14.json,KCLASS=2,EXTRA=3,R2=1 slurm_single.sh
#   sbatch --export=ALL,TARGETS=data/controls.json,NAMES=13a_647,KCLASS=3,EXTRA=3,ROUTES=1,2,3 slurm_single.sh
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
TARGETS="${TARGETS:-$PKG/data/targets_11a_14.json}"
KNOWN="${KNOWN:-$PKG/data/known_u_codes_noDKT.json}"
KCLASS="${KCLASS:-2}"
EXTRA="${EXTRA:-4}"
MAXC="${MAXC:-20}"
ROUTES="${ROUTES:-1,2}"
PASS_EXTRA="${PASS_EXTRA:-3}"
TPK="${TPK:-0}"
TAG="${TAG:-single_${KCLASS}_x${EXTRA}${R2:+_r2}}"
OUTDIR="${OUTDIR:-${SLURM_SUBMIT_DIR:-$PWD}/single_runs/$TAG}"
mkdir -p "$OUTDIR"
CMD=("$PY" "$PKG/enum_search.py" --targets "$TARGETS" --known "$KNOWN" --k-class "$KCLASS"
     --max-crossings "$MAXC" --extra-crossings "$EXTRA" --time-per-knot "$TPK"
     --routes "$ROUTES" --pass-extra "$PASS_EXTRA" --max-routes 0 --hits-per-state 0 --max-states 0
     --result-tag "$TAG" --out-dir "$OUTDIR")
[ -n "${NAMES:-}" ] && CMD+=(--names $NAMES)
[ -n "${R2:-}" ] && CMD+=(--r2)
if [ "${1:-}" = "--dry-run" ]; then printf 'command:'; printf ' %q' "${CMD[@]}"; printf '\n'; exit 0; fi
echo "single-knot run: targets=$TARGETS known=$KNOWN k=$KCLASS budget=+$EXTRA (cap $MAXC) routes=$ROUTES r2=${R2:-0} out=$OUTDIR"
exec env PYTHONUNBUFFERED=1 "${CMD[@]}"
