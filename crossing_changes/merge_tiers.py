from paths import *
import os
import json, glob, collections, ast, statistics
import database_knotinfo as dk
OUT = OPEN23
pri = json.load(open(OUT + '/priority_u23.json')); knots = {k['knot']: k for k in pri['knots']}
fl = {}
for l in open(OUT + '/flype_orbit_check.jsonl'): r = json.loads(l); fl[r['knot']] = r
withc = [k for k in pri['knots'] if k['n_candidates'] > 0]
zero = [k for k in pri['knots'] if k['n_candidates'] == 0]
zero_rig = [k for k in zero if k['n_heuristic_exclusions'] == 0]
# validation: flype orbit of the zero-candidate knots
real = {nm: [c for c in fl[nm]['candidates'] if c['status'] == 'candidate'] for nm in fl}
art = {nm: [c for c in fl[nm]['candidates'] if c['status'] != 'candidate'] for nm in fl}
print('zero-candidate knots (KnotInfo diagram):', len(zero), 'rigorous', len(zero_rig))
print('flype-orbit validation on them: knots with a genuine candidate in another minimal diagram:', sum(1 for nm in real if real[nm]), '; with only unidentified/ambiguous results in other diagrams:', sum(1 for nm in fl if not real[nm] and art[nm]))
for nm in [nm for nm in real if real[nm]][:5]: print('   genuine:', nm, real[nm][:2])
print('mean number of minimal diagrams:', round(statistics.mean(r['diagrams'] for r in fl.values()), 2))
sig4 = [k for k in pri['knots'] if abs(k['sigma']) == 4]; sig2 = [k for k in pri['knots'] if abs(k['sigma']) == 2]
print('|sigma|=4:', len(sig4), 'of which zero-candidate', sum(1 for k in sig4 if k['n_candidates'] == 0), '; |sigma|=2:', len(sig2), 'zero-candidate', sum(1 for k in sig2 if k['n_candidates'] == 0))
p3 = lambda L: round(statistics.mean((x['P_u3'] or 0.0) for x in L), 3) if L else None
print('mean P(u=3) of the earlier ML ranking: with candidates', p3(withc), 'zero-candidate', p3(zero))
child = collections.defaultdict(lambda: {'parents': set(), 'max_score': 0.0, 'u_range': None})
for k in withc:
    for c in k['candidates']:
        if c['result'] and c['status'] == 'candidate':
            d = child[c['result']]; d['parents'].add(k['knot']); d['max_score'] = max(d['max_score'], c['score']); d['u_range'] = c['u_range']
ki = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
ch = sorted([{'child': c, 'u_range': d['u_range'], 'n_parents': len(d['parents']), 'parents': sorted(d['parents']), 'max_score': round(d['max_score'], 3)} for c, d in child.items()], key=lambda d: (-d['n_parents'], -d['max_score']))
covered = {p for d in ch for p in d['parents']}
print('distinct child knots:', len(ch), '; parents covered:', len(covered), '; knots with only unresolved candidates:', len([k for k in withc if k['knot'] not in covered]))
for d in ch[:12]: print(f"   {d['child']:9s} {d['u_range']} parents={d['n_parents']:2d} maxscore={d['max_score']:.3f}  {', '.join(d['parents'][:7])}")
unres = [(k['knot'], c) for k in withc for c in k['candidates'] if c['status'] != 'candidate']
tg = []
for d in ch:
    if d['child'] in ki:
        r = ki[d['child']]; tg.append({'name': d['child'], 'range': d['u_range'], 'k': 1, 'crossings': int(r['crossing_number']), 'det': int(r['determinant']), 'alternating': r['alternating'], 'n_parents': d['n_parents'], 'parents': d['parents'], 'max_score': d['max_score'], 'pd': ast.literal_eval(r['pd_notation'])})
json.dump(tg, open(os.path.join(ROOT, 'upper_bounds', 'cluster', 'data', 'targets_u23_children.json'), 'w'))
fin = {'criterion': 'u(K)=2 together with the Bernhard-Jablan property needs a crossing c of a minimal diagram with u(K_c)=1; crossing changes commute with flypes for prime alternating knots, so one minimal diagram suffices; necessary: |sigma(K_c)| <= 2; certified exclusions: known u(K_c) >= 2, verified nontrivial composites (Scharlemann); Jones-only identifications remain heuristic',
       'knots_with_candidates': sorted(k['knot'] for k in withc), 'zero_candidate_knots': sorted([(k['knot'], round((k['P_u3'] or 0.0), 3), k['sigma'], k['n_heuristic_exclusions']) for k in zero], key=lambda x: -x[1]),
       'zero_candidate_rigorous_count': len(zero_rig), 'flype_validation': {'checked': len(fl), 'genuine_candidates_elsewhere': sum(1 for nm in real if real[nm]), 'artifacts_only': sum(1 for nm in fl if not real[nm] and art[nm])},
       'children': ch, 'unresolved_candidates': unres, 'knots_with_only_unresolved_candidates': sorted(k['knot'] for k in withc if k['knot'] not in covered)}
json.dump(fin, open(OUT + '/priority_u23_final.json', 'w'), indent=1)
with open(OUT + '/priority_u23.md', 'w') as f:
    f.write('# Search priority for the 833 open alternating [2,3] knots (hypothesis u = 2)\n\n')
    f.write('If u(K) = 2, the Bernhard–Jablan conjecture gives a crossing c of a minimal diagram with u(K_c) = 1. Crossing changes commute with flypes, so the multiset of change results is the same for every minimal diagram and the KnotInfo diagram suffices (validated: no zero-candidate knot gains a genuine candidate in its flype orbit). A crossing is a *candidate* if |sigma(K_c)| <= 2 (necessary for u(K_c) = 1) and the known u-range of K_c contains 1. Composites are excluded by Scharlemann. Identification: canonical code, SnapPy isometry against all KnotInfo knots with the same det and Jones polynomial, block decomposition; Jones-only composite identification only for non-hyperbolic results (flagged heuristic).\n\n')
    f.write(f'- {len(withc)} knots have candidate crossings; their u = 2 hinges on {len(ch)} child knots with u in [1,2] (or [1,3])\n- {len(zero)} knots have NO candidate crossing ({len(zero_rig)} rigorous, the rest use Jones-identified composites): u = 2 would contradict the Bernhard–Jablan conjecture\n- |sigma| = 4: {len(sig4)} knots ({sum(1 for k in sig4 if k["n_candidates"] == 0)} without candidates); |sigma| = 2: {len(sig2)} knots ({sum(1 for k in sig2 if k["n_candidates"] == 0)} without candidates)\n\n')
    f.write(f'Mean P(u=3) of the earlier ML ranking: with candidates {p3(withc)}, without {p3(zero)}\n\n')
    f.write('## Child knots to attack (u(K_c) = 1 proves u = 2 for all parents; u(K_c) >= 2 removes the candidates)\n\n| K_c | u-range | parents | max score | parent knots |\n|---|---|---|---|---|\n')
    for d in ch: f.write(f"| {d['child']} | {d['u_range']} | {d['n_parents']} | {d['max_score']} | {', '.join(d['parents'][:10])}{' …' if len(d['parents']) > 10 else ''} |\n")
    f.write('\n## Knots without any candidate crossing, by P(u=3)\n\n')
    f.write(', '.join(f"{nm} ({p:.2f}{'' if h == 0 else ', jones'})" for nm, p, s, h in fin['zero_candidate_knots']) + '\n')
    f.write('\n## Unresolved candidates\n\n' + '\n'.join(f"- {nm}: {c}" for nm, c in unres) + '\n')
print('targets written:', len(tg))
