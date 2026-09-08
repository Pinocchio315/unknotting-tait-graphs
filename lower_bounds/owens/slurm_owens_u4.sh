#!/usr/bin/env bash
#SBATCH --job-name=owens-u4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=owens_u4_%A_%a.out
#SBATCH --array=0-114
# Rank-4 Owens-type obstruction (u = 4 excluded for alternating knots with |sigma| = 8), one knot per
# array task: 24 open [4,x] targets followed by the 91 controls with known u = 4 (every control must PASS).
#   sbatch slurm_owens_u4.sh                      # needs numpy + database_knotinfo in $PY
# Results: one JSON line per task in owens_u4_<jobid>_<task>.out
set -euo pipefail
PY="${PY:-python3}"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAME=$("$PY" - "$SLURM_ARRAY_TASK_ID" <<'PYEOF'
import sys, re, database_knotinfo as dk
rows=[r for r in dk.link_list() if str(r.get('crossing_number','')).strip().isdigit()]
def rng(u):
    m=re.findall(r'\d+',str(u)); return (int(m[0]),int(m[-1])) if m else None
alt=lambda r: str(r.get('alternating','')).upper().startswith('Y')
s8=[r for r in sorted(rows,key=lambda r:int(r['determinant'])) if alt(r) and abs(int(r['signature']))==8]
targets=[r['name'] for r in s8 if rng(r['unknotting_number'])[0]==4 and rng(r['unknotting_number'])[1]>4]
controls=[r['name'] for r in s8 if rng(r['unknotting_number'])==(4,4)]
names=targets+controls
print(names[int(sys.argv[1])])
PYEOF
)
cd "$DIR"
exec env PYTHONUNBUFFERED=1 "$PY" owens_u4.py "$NAME"
