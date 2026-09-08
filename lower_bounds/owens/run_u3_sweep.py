"""Sweep the rank-3 Owens-type obstruction over alternating knots with |sigma| = 6.
    python run_u3_sweep.py validate out.jsonl   # knots with known u = 3: every verdict must be PASS
    python run_u3_sweep.py apply out.jsonl      # gapped rows with lower bound 3: OBSTRUCTED => u >= 4
"""
import sys, json, time, ast, re
import database_knotinfo as dk
from owens_u3 import obstruct_u3
mode, out = sys.argv[1], sys.argv[2]
shard = sys.argv[3] if len(sys.argv) > 3 else '0/1'
si, sn = (int(x) for x in shard.split('/'))
rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
def rng(u):
    m = re.findall(r'\d+', str(u)); return (int(m[0]), int(m[-1]))
sel = []
for r in rows:
    if not str(r.get('alternating', '')).upper().startswith('Y') or abs(int(r['signature'])) != 6:
        continue
    lo, hi = rng(r['unknotting_number'])
    if mode == 'validate' and lo == hi == 3: sel.append(r)
    if mode == 'apply' and lo == 3 and hi > 3: sel.append(r)
import glob, os
done = set()
for f in glob.glob(os.path.join(os.path.dirname(out) or '.', mode + '_u3*.jsonl')):
    for line in open(f):
        try: done.add(json.loads(line)['name'])
        except ValueError: pass
sel = sel[si::sn]
print(f'{mode}: {len(sel)} knots, {len(done)} done', flush=True)
with open(out, 'a') as h:
    for r in sel:
        if r['name'] in done: continue
        t0 = time.time()
        try:
            res = obstruct_u3(ast.literal_eval(r['pd_notation']), int(r['signature']), int(r['determinant']))
        except Exception as e:
            res = {'verdict': 'ERROR', 'reason': repr(e)[:200]}
        rec = {'name': r['name'], 'u': str(r['unknotting_number']).strip(), 'sigma': int(r['signature']), 'det': int(r['determinant']), **res, 'seconds': round(time.time() - t0, 1)}
        h.write(json.dumps(rec) + '\n'); h.flush()
        print(f"{rec['name']} u={rec['u']} det={rec['det']} -> {rec['verdict']} ({rec.get('candidates','-')} forms, {rec['seconds']}s)", flush=True)
