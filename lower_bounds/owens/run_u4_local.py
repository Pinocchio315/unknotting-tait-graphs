#!/usr/bin/env python
"""Run the rank-4 test (owens_u4.py) locally for every alternating knot with |sigma| = 8 that is not yet in
the partial result file: the open [4, x] targets first, then the u = 4 controls (which must all PASS).
    python run_u4_local.py --workers 6 --out results_u4_local_2026-09-07.jsonl
    python run_u4_local.py --workers 3 --targets-only        # the open [4, x] knots only, no controls
"""
import argparse, json, os, re, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=6); ap.add_argument('--out', default=os.path.join(HERE, 'results_u4_local_2026-09-07.jsonl'))
ap.add_argument('--skip', default=os.path.join(HERE, 'results_u4_partial_2026-09-07.json')); ap.add_argument('--targets-only', action='store_true'); ap.add_argument('--only', default='', help='comma-separated knot names: run only these'); args = ap.parse_args()
import database_knotinfo as dk
rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
def rng(u):
    m = re.findall(r'\d+', str(u)); return (int(m[0]), int(m[-1])) if m else None
alt = lambda r: str(r.get('alternating', '')).upper().startswith('Y')
s8 = [r for r in sorted(rows, key=lambda r: int(r['determinant'])) if alt(r) and abs(int(r['signature'])) == 8]
targets = [r['name'] for r in s8 if rng(r['unknotting_number'])[0] == 4 and rng(r['unknotting_number'])[1] > 4]
controls = [r['name'] for r in s8 if rng(r['unknotting_number']) == (4, 4)]
done = set()
if os.path.exists(args.skip):
    d = json.load(open(args.skip)); done |= set(d['targets']) | set(d['controls'])
if os.path.exists(args.out):
    for l in open(args.out):
        try: done.add(json.loads(l)['name'])
        except Exception: pass
todo = [n for n in (targets if args.targets_only else targets + controls) if n not in done]
if args.only: todo = [n for n in todo if n in set(args.only.split(','))]
print(f'{len(targets)} targets, {len(controls)} controls, {len(done)} done, {len(todo)} to run', flush=True)
def run(nm):
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, 'owens_u4.py'), nm], capture_output=True, text=True)
    line = [l for l in p.stdout.splitlines() if l.startswith('{')]
    rec = json.loads(line[-1]) if line else {'name': nm, 'verdict': 'ERROR', 'stderr': p.stderr[-500:]}
    rec['role'] = 'target' if nm in targets else 'control'
    with open(args.out, 'a') as h: h.write(json.dumps(rec) + '\n')
    print(f"{nm} {rec.get('verdict')} {round(time.time() - t0)}s", flush=True)
with ThreadPoolExecutor(args.workers) as ex: list(ex.map(run, todo))
print('done', flush=True)
