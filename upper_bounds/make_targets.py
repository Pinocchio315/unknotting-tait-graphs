#!/usr/bin/env python
"""Generate data/targets.json and data/known_u.json (run LOCALLY; the outputs are shipped).

Targets: every knot with <= 13 crossings whose unknotting range, AFTER applying the lower-bound
updates of this project, is [1,2], [2,3] or [3,4] (k = lo = 1, 2, 3 crossing changes to try).  Updates applied on top of KnotInfo 2026.8.1:
  - 222 knots proved u = 3 and 732 proved u = 2 (settled -> excluded),
  - the 93 linking-form range improvements [1,b] -> [2,b] and 13n_689 [1,3] -> [2,3]
    (those landing on [2,3] become NEW targets),
  - 13a_4877 [2,5] -> [3,5] and 13n_1587 [1,3] -> [1,2] (neither lands in scope).

known_u.json: every knot with exact u <= 2 after our updates (u = 1 and u = 2 lists) with
determinant and the canonical codes of its KnotInfo diagram and mirror --- the lookup table for
the "j < k changes reduce to a known small knot" route.

    ~/.pyenv/versions/unknot-venv/bin/python make_targets.py
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd, to_pd  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402

RESULTS = os.path.join(HERE, '..', 'results')


def knotinfo_rows():
    import database_knotinfo
    root = os.path.dirname(database_knotinfo.__file__)
    path = next(os.path.join(dp, f) for dp, dn, fn in os.walk(root) for f in fn
                if f == 'knotinfo_data_complete.csv')
    with open(path, newline='', encoding='utf-8') as h:
        return {r['name']: r for r in csv.DictReader(h, delimiter='|')}


def interval(r):
    lo = r['unknotting_number'].strip()
    if re.fullmatch(r'\d+', lo):
        return int(lo), int(lo)
    m = re.fullmatch(r'\[(\d+),(\d+)\]', lo.replace(' ', ''))
    if m:
        return int(m.group(1)), int(m.group(2))
    raise ValueError((r['name'], lo))


def pd_code(r):
    s = r['pd_notation'].strip()
    return [[int(x) for x in q.split(',')] for q in re.findall(r'\[([\d, ]+)\]', s)]


def mirror_pd_graph(pd):
    g = from_pd([list(q) for q in pd])
    h = g.copy()
    for e in h.signs:
        h.signs[e] = -h.signs[e]
    return h


def main() -> None:
    rows = knotinfo_rows()

    # ------- our updates -------------------------------------------------------------
    u3 = {x['name'] for x in json.load(open(os.path.join(RESULTS, 'owens_u3_2026-08-24.json')))}
    u3.add('13a_1786')
    lick = json.load(open(os.path.join(RESULTS, 'lickorish_cw_obstructed_815_knotinfo2026.8.1.json')))
    u2 = {nm for nm, iv in lick if iv.replace(' ', '') == '[1,2]'}
    u2 |= {'13n_1166', '13n_2504', '13n_30', '13n_45', '13n_80', '13n_2379',
           '13n_2809', '13n_2907', '13n_3033', '13n_3589'}
    improved_lo = {nm: 2 for nm, iv in lick if iv.replace(' ', '') != '[1,2]'}
    improved_lo['13n_689'] = 2
    improved_lo['13a_4877'] = 3
    u3.add('13a_647')          # settled by the Section-4 witness chain (2026-08-28)
    settled = u3 | u2

    # ------- targets ------------------------------------------------------------------
    targets = []
    for nm, r in rows.items():
        try:
            lo, hi = interval(r)
        except ValueError:
            continue
        if nm in settled or lo == hi:
            continue
        lo = max(lo, improved_lo.get(nm, lo))
        if nm == '13n_1587':
            hi = 2
        if (lo, hi) not in ((1, 2), (2, 3), (3, 4)):
            continue
        pd = pd_code(r)
        targets.append({'name': nm, 'range': [lo, hi], 'k': lo, 'crossings': int(r['crossing_number']),
                        'det': abs(int(r['determinant'])), 'alternating': r['alternating'].strip(),
                        'new_23': nm in improved_lo, 'pd': pd})
    targets.sort(key=lambda t: (t['k'], t['crossings'], t['name']))
    n12 = sum(1 for t in targets if t['k'] == 1)
    n23 = sum(1 for t in targets if t['k'] == 2)
    n34 = sum(1 for t in targets if t['k'] == 3)

    # ------- known-u table ------------------------------------------------------------
    known = []
    for nm, r in rows.items():
        try:
            lo, hi = interval(r)
        except ValueError:
            continue
        u = lo if lo == hi else (2 if nm in u2 else None)
        if u is None or u > 2:
            continue
        pd = pd_code(r)
        try:
            g = from_pd([list(q) for q in pd])
            codes = [repr(canonical_code(g)), repr(canonical_code(mirror_pd_graph(pd)))]
        except Exception as exc:
            print('  ! canonical code failed for', nm, exc)
            codes = []
        known.append({'name': nm, 'u': u, 'det': abs(int(r['determinant'])),
                      'crossings': int(r['crossing_number']), 'codes': codes})
    known.sort(key=lambda x: (x['u'], x['crossings'], x['name']))
    n1 = sum(1 for x in known if x['u'] == 1)

    # ------- positive controls: settled minimal-witness knots (k=2) and u=1 knots (k=1)
    controls = []
    for nm in ('13n_30', '13n_45', '13n_80', '13n_2379', '13n_2809',
               '13n_2907', '13n_3033', '13n_3589'):
        r = rows[nm]
        controls.append({'name': nm, 'range': [2, 2], 'k': 2, 'crossings': 13,
                         'det': abs(int(r['determinant'])), 'alternating': r['alternating'].strip(),
                         'new_23': False, 'control': True, 'pd': pd_code(r)})
    n_u1 = 0
    for nm, r in sorted(rows.items()):
        if not r['crossing_number'].strip().isdigit() or r['crossing_number'].strip() != '13':
            continue
        try:
            lo, hi = interval(r)
        except ValueError:
            continue
        if lo == hi == 1 and n_u1 < 5:
            n_u1 += 1
            controls.append({'name': nm, 'range': [1, 1], 'k': 1, 'crossings': 13,
                             'det': abs(int(r['determinant'])),
                             'alternating': r['alternating'].strip(),
                             'new_23': False, 'control': True, 'pd': pd_code(r)})

    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    with open(os.path.join(HERE, 'data', 'targets.json'), 'w') as h:
        json.dump(targets, h)
    with open(os.path.join(HERE, 'data', 'controls.json'), 'w') as h:
        json.dump(controls, h)
    with open(os.path.join(HERE, 'data', 'known_u.json'), 'w') as h:
        json.dump(known, h)
    print(f'targets: {len(targets)}  ([1,2]: {n12}, [2,3]: {n23}, [3,4]: {n34}; '
          f'new [2,3] from our lower bounds: {sum(1 for t in targets if t["new_23"])})')
    print(f'controls: {len(controls)} (8 witness knots k=2, {n_u1} known-u=1 knots k=1)')
    print(f'known u<=2 table: {len(known)}  (u=1: {n1}, u=2: {len(known)-n1})')


if __name__ == '__main__':
    main()
