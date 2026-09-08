from paths import *
from candidates import candidate_status
import os
import json, glob, collections, ast
from fractions import Fraction
import numpy as np, database_knotinfo as dk
from sklearn.ensemble import HistGradientBoostingClassifier
OUT = OPEN23
S = os.path.join(OPEN23, 'work')
FEATS = ['sign', 'sig_rel', 'dsig', 'sig_drop', 'R_gt_half', 'Rf', 'ratio', 'deg_min', 'deg_max', 'face_min', 'face_max', 'twist_par', 'twist_ser', 'twist_max_diagram', 'reduced_crossings', 'crossings', 'det', 'det_changed', 'u', 'sigma', 'n_black', 'n_white']
def derive(r):
    r['Rb'] = Fraction(r['R_black']); r['Rf'] = float(r['Rb']); r['dsig'] = abs(r['sigma_c']) - abs(r['sigma'])
    r['sig_drop'] = int(r['dsig'] == -2); r['sig_rel'] = r['sign'] * (1 if r['sigma'] > 0 else -1 if r['sigma'] < 0 else 0)
    r['R_gt_half'] = int(r['Rb'] > Fraction(1, 2)); r['ratio'] = r['det_changed'] / r['det']
import gzip
_p = os.path.join(CC, 'dataset_v2.json')
train = [r for r in (json.load(open(_p)) if os.path.exists(_p) else json.load(gzip.open(_p + '.gz', 'rt'))) if r['label'] is not None]
for r in train: derive(r)
X = np.array([[float(r[f]) for f in FEATS] for r in train]); y = np.array([r['label'] == 'good' for r in train])
model = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05).fit(X, y)
T = [json.loads(l) for f in sorted(glob.glob(S + '/rows833_*.jsonl')) for l in open(f)]
for r in T: derive(r)
XT = np.array([[float(r[f]) for f in FEATS] for r in T]); p = model.predict_proba(XT)[:, 1]
for r, s in zip(T, p): r['score'] = float(s)
# Under u(K)=2 AND the Bernhard–Jablan property, a minimal-diagram neighbor
# must have u=1. Machine-learning scores only order the surviving candidates.
stat = collections.Counter(); byk = collections.defaultdict(list); witnesses = []
for r in T:
    ur = r['u_range']
    r['status'] = candidate_status(r['result'], r['result_how'], ur, r['sigma_c'])
    if r['status'] == 'WITNESS': witnesses.append(r)
    stat[r['status']] += 1; byk[r['knot']].append(r)
print('crossing status:', dict(stat)); print('witnesses (would give u=2):', [(w['knot'], w['result']) for w in witnesses])
# per-knot summary
knots = []
for nm, rs in byk.items():
    cands = sorted([r for r in rs if r['status'].startswith('candidate')], key=lambda r: -r['score'])
    knots.append({'knot': nm, 'P_u3': rs[0]['P_u3'], 'crossings': len(rs), 'sigma_dropping': sum(r['sig_drop'] for r in rs), 'sigma': rs[0]['sigma'], 'sigma_ok': sum(abs(r['sigma_c']) <= 2 for r in rs), 'n_candidates': len(cands),
                  'candidates': [{'crossing': r['crossing'], 'result': r['result'], 'u_range': r['u_range'], 'status': r['status'], 'score': round(r['score'], 3), 'reduced_crossings': r['reduced_crossings'], 'det_changed': r['det_changed']} for r in cands],
                  'excluded_results': sorted({r['result'] for r in rs if r['status'] == 'excluded:u>=2'}, key=str)})
knots.sort(key=lambda k: (-max([c['score'] for c in k['candidates']], default=-1)))
for k in knots: k['n_heuristic_exclusions'] = sum(1 for r in byk[k['knot']] if r['status'] == 'excluded:u>=2(jones)')
zero = [k for k in knots if k['n_candidates'] == 0]
zero_rig = [k for k in zero if k['n_heuristic_exclusions'] == 0]
print(f'\nknots with NO candidate crossing in the KnotInfo diagram: {len(zero)} (rigorous: {len(zero_rig)}; the rest rely on Jones-identified composites)')
print('  ', ', '.join(f"{k['knot']}(P3={(k['P_u3'] or 0.0):.2f})" for k in sorted(zero, key=lambda k: -(k['P_u3'] or 0.0))[:60]))
dist = collections.Counter(min(k['n_candidates'], 6) for k in knots); print('candidate-count distribution (capped at 6):', dict(sorted(dist.items())))
# aggregate by child knot K_c
child = collections.defaultdict(lambda: {'parents': [], 'scores': []})
for r in T:
    if r['status'].startswith('candidate') and r['result']:
        child[r['result']]['parents'].append(r['knot']); child[r['result']]['scores'].append(r['score']); child[r['result']]['u_range'] = r['u_range']; child[r['result']]['how'] = r['result_how']
ch = [{'child': c, 'u_range': d.get('u_range'), 'n_parents': len(set(d['parents'])), 'parents': sorted(set(d['parents'])), 'max_score': round(max(d['scores']), 3), 'mean_score': round(float(np.mean(d['scores'])), 3), 'how': d.get('how')} for c, d in child.items()]
ch.sort(key=lambda d: (-d['n_parents'], -d['max_score']))
print(f'\ndistinct candidate child knots K_c: {len(ch)}; with >=2 parents: {sum(1 for d in ch if d["n_parents"] >= 2)}')
print('top 25 children by (number of parent knots, max score):')
for d in ch[:25]: print(f"   {d['child']:12s} range={d['u_range']} parents={d['n_parents']:2d} max={d['max_score']:.3f} mean={d['mean_score']:.3f}  {', '.join(d['parents'][:6])}")
print('\ntop 25 (knot, crossing) by score:')
top = sorted([r for r in T if r['status'].startswith('candidate')], key=lambda r: -r['score'])
for r in top[:25]: print(f"   {r['knot']:9s} c{r['crossing']:2d} -> {str(r['result']):12s} {r['u_range']} score={r['score']:.3f} P3={(r['P_u3'] or 0.0):.2f} red={r['reduced_crossings']}")
# server target file for the children (u = 1 witness search), priority order
ki = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
tg = []
for d in ch:
    nm = d['child']
    if nm not in ki: continue
    r = ki[nm]
    tg.append({'name': nm, 'range': d['u_range'], 'k': 1, 'crossings': int(r['crossing_number']), 'det': int(r['determinant']), 'alternating': r['alternating'], 'n_parents': d['n_parents'], 'parents': d['parents'], 'max_score': d['max_score'], 'pd': ast.literal_eval(r['pd_notation'])})
json.dump(tg, open(os.path.join(ROOT, 'upper_bounds', 'cluster', 'data', 'targets_u23_children.json'), 'w'))
json.dump({'knots': knots, 'children': ch, 'status_counts': dict(stat), 'witnesses': [(w['knot'], w['result']) for w in witnesses]}, open(OUT + '/priority_u23.json', 'w'), indent=1)
with open(OUT + '/priority_u23.md', 'w') as f:
    f.write('# Search priority for the 833 Owens-PASS alternating [2,3] knots (hypothesis u = 2)\n\n')
    f.write(f'Crossing status counts: {dict(stat)}\n\nKnots whose KnotInfo diagram has no candidate crossing ({len(zero)}): ' + ', '.join(k['knot'] for k in sorted(zero, key=lambda k: -(k['P_u3'] or 0.0))) + '\n\n')
    f.write('## Child knots K_c to attack first (proving u(K_c) = 1 gives u = 2 for every parent)\n\n| K_c | u-range | parents | max score | parent knots |\n|---|---|---|---|---|\n')
    for d in ch[:60]: f.write(f"| {d['child']} | {d['u_range']} | {d['n_parents']} | {d['max_score']} | {', '.join(d['parents'][:8])}{' …' if len(d['parents']) > 8 else ''} |\n")
    f.write('\n## Per-knot top candidate\n\n| knot | P(u=3) | candidates | best crossing | K_c | u-range | score |\n|---|---|---|---|---|---|---|\n')
    for k in knots:
        c = k['candidates'][0] if k['candidates'] else None
        f.write(f"| {k['knot']} | {(k['P_u3'] or 0.0):.2f} | {k['n_candidates']}/{k['sigma_dropping']}/{k['crossings']} | {c['crossing'] if c else '-'} | {c['result'] if c else '-'} | {c['u_range'] if c else '-'} | {c['score'] if c else '-'} |\n")
print('\nwritten:', OUT + '/priority_u23.{json,md}', 'and server_code/data/targets_u23_children.json with', len(tg), 'targets')
