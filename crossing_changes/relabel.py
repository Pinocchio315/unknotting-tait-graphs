from paths import *
from candidates import connected_sum_interval
import os
"""Merge rows + augmentation + second-pass identification; assign composite labels; per-knot BJ-diagram check."""
import json, glob, collections, re
from fractions import Fraction
import database_knotinfo as dk
S = WORK
rows = {}
for f in glob.glob(S + '/rows_*.jsonl'):
    for l in open(f): r = json.loads(l); rows[(r['knot'], r['crossing'])] = r
for f in glob.glob(S + '/aug_*.jsonl'):
    for l in open(f): a = json.loads(l); rows[(a['knot'], a['crossing'])].update(a)
T = json.load(open(U_TABLE))                       # reference table + our bounds (consolidate_results.py)
U = {nm: lo for nm, (lo, hi) in T.items() if lo == hi}; LO = {nm: lo for nm, (lo, hi) in T.items()}
U['0_1'] = 0
ir = json.load(open(S + '/identify_rest.json'))
comp_stats = collections.Counter()
for key, members in ir['classes'].items():
    v = ir['results'][key]['verdict']
    for nm, i in members:
        r = rows[(nm, i)]
        if v is None: comp_stats['unresolved'] += 1; continue
        kind, data = v
        if kind in ('sum-blocks', 'sum-jones'):
            parts = data
            r['result'] = '#'.join(parts); r['result_how'] = kind
            r['u_result'] = None
            r['u_range'] = None
            # A product Jones polynomial is only a proposed identification;
            # it cannot supply an exact label or a Bernhard–Jablan witness.
            interval = connected_sum_interval(parts, T, r.get('sigma_c')) if kind == 'sum-blocks' else None
            if interval is not None:
                lo, up = interval
                if lo == up: r['u_result'] = up; comp_stats[kind + ':exact'] += 1
                else: r['u_result'] = None; r['u_range'] = [lo, up]; comp_stats[kind + ':range'] += 1
            else: comp_stats[kind + ':factor-u-unknown'] += 1
        elif kind == 'prime-isometry':
            r['result'] = data; r['result_how'] = kind; r['u_result'] = U.get(data); comp_stats[kind] += 1
        else: comp_stats[kind] += 1
        uc = r['u_result']; uK = r['u']
        r['label'] = None if uc is None else ('good' if uc == uK - 1 else ('neutral' if uc == uK else 'bad'))
print('composite handling:', dict(comp_stats))
print('labels now:', collections.Counter(r['label'] for r in rows.values()))
json.dump(list(rows.values()), open(os.path.join(CC, 'dataset_v2.json'), 'w'))
# per-knot BJ check on the KnotInfo diagram
byk = collections.defaultdict(list)
for r in rows.values(): byk[r['knot']].append(r)
stat = collections.Counter(); nogood = []; nogood_unknown = []
for nm, rs in byk.items():
    c = collections.Counter(r['label'] for r in rs)
    if c['good'] == 0:
        if c[None] == 0: nogood.append((nm, rs[0]['u'], dict(c)))
        else: nogood_unknown.append((nm, rs[0]['u'], dict(c), sorted({(r['result'], r.get('u_range')) for r in rs if r['label'] is None}, key=str)[:6]))
    stat[(rs[0]['u'], c['good'] > 0, c[None] > 0)] += 1
print('\nper-knot (u, has good crossing, has unknown crossing):')
for k in sorted(stat): print('  ', k, stat[k])
print('\nKnotInfo diagrams with NO good crossing and NO unknown crossing (diagram-level BJ failure):', len(nogood))
for x in nogood[:40]: print('  ', x)
print('\nKnotInfo diagrams with NO good crossing but some unknown results:', len(nogood_unknown))
for x in nogood_unknown[:60]: print('  ', x)
json.dump({'nogood': nogood, 'nogood_unknown': nogood_unknown}, open(os.path.join(CC, 'bj_diagram_check.json'), 'w'), indent=1)
