#!/usr/bin/env python
"""For the signature-sharp alternating knots with u = 2 and their sigma-dropping crossings: which known
obstruction explains u(K_c) = 2 when the crossing does not lower u?  Reads the data set and the electrical
Lickorish verdicts; writes results/crossing_changes/residual_breakdown_2026-09-07.json (counts used in the paper).
    python residual_breakdown.py
"""
import json, gzip, os, re, collections
import database_knotinfo as dk
from paths import *
rows = json.load(gzip.open(os.path.join(CC, 'dataset_v2.json.gz'), 'rt'))
L = {(r['knot'], r['crossing']): r for r in json.load(open(os.path.join(CC, 'lickorish_electrical_2026-09-07.json')))}
ki = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
cc = json.load(open(os.path.join(RESULTS, 'cyclic_cover', 'cyclic_cover_bound_2026-09-07.json')))
mont = json.load(open(os.path.join(RESULTS, 'montesinos', 'montesinos_u1_2026-09-07.json')))
mont_obstructed = {nm for part in ('validate', 'apply') for nm, v in mont[part].items() if v['verdict'] == 'OBSTRUCTED'}
def num(s):
    m = re.findall(r'-?\d+', str(s)); return int(m[0]) if m else None
def explain(nm):
    r = ki.get(nm)
    if r is None: return 'composite' if '#' in str(nm) else 'not tabulated'
    if r['alternating'].upper().startswith('Y'): return 'alternating (McCoy)'
    if (num(r['unknotting_number_algebraic']) or 0) >= 2 or (num(r['nakanishi_index']) or 0) >= 2 or cc.get(nm, {}).get('bound', 0) >= 2: return 'Alexander module / cyclic covers'
    if (num(r['smooth_four_genus']) or 0) >= 2 or abs(num(r['ozsvath_szabo_tau_invariant']) or 0) >= 2: return 'four-genus / tau'
    if nm in mont_obstructed: return 'Montesinos correction terms'
    return 'other (correction terms of Owens-Strle)'
for r in rows:
    r['sharp'] = 2 * r['u'] == abs(r['sigma']); r['dsig'] = abs(r['sigma_c']) - abs(r['sigma'])
out = {}
u1 = [r for r in rows if r.get('u_result') == 1]
out['rows'] = len(rows); out['u1_rows'] = len(u1); out['u1_rows_lickorish_pass'] = sum(1 for r in u1 if L[(r['knot'], r['crossing'])]['lickorish'] == 'pass')
out['cyclic_rows'] = sum(1 for l in L.values() if l.get('cyclic')); out['x_generates_rows'] = sum(1 for l in L.values() if l.get('x_generates'))
S2 = [r for r in rows if r['sharp'] and r['u'] == 2 and r['dsig'] == -2 and r['label'] in ('good', 'neutral')]
out['sharp2_rows'] = len(S2); out['sharp2_good'] = sum(1 for r in S2 if r['label'] == 'good')
neutral = [r for r in S2 if r['label'] == 'neutral']
out['sharp2_neutral'] = len(neutral)
out['sharp2_neutral_lickorish_fail'] = sum(1 for r in neutral if L[(r['knot'], r['crossing'])]['lickorish'] != 'pass')
rest = [r for r in neutral if L[(r['knot'], r['crossing'])]['lickorish'] == 'pass']
out['sharp2_neutral_lickorish_pass'] = len(rest); out['sharp2_neutral_pass_children'] = len({r['result'] for r in rest})
bre = collections.Counter(explain(r['result']) for r in rest); out['breakdown_rows'] = dict(bre)
out['breakdown_children'] = dict(collections.Counter(explain(nm) for nm in {r['result'] for r in rest}))
out['other_children'] = sorted({r['result'] for r in rest if explain(r['result']).startswith('other')})
json.dump(out, open(os.path.join(CC, 'residual_breakdown_2026-09-07.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
