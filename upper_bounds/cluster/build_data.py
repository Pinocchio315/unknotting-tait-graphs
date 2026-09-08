#!/usr/bin/env python
"""Build the data files for the server (run LOCALLY; needs database_knotinfo).

  data/targets_23_34.json   the [2,3] and [3,4] targets (from ../code/upper_bounds/data/targets.json)
  data/known_u_codes.json   every knot with known u <= 2 (KnotInfo 2026.8.1 + this paper): name, u,
                            det, and the canonical codes of ALL diagrams in the flype orbit of its
                            KnotInfo diagram (= all minimal diagrams for alternating knots)
  data/controls.json        the nine knots settled in the paper, as control targets
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xtait.graph import from_pd  # noqa: E402
from xtait.flypes import flype_orbit  # noqa: E402
from alexander import fingerprint  # noqa: E402

CODE = os.path.join(HERE, '..', 'code')


def knotinfo_rows():
    import database_knotinfo
    root = os.path.dirname(database_knotinfo.__file__)
    path = next(os.path.join(dp, f) for dp, dn, fn in os.walk(root) for f in fn if f == 'knotinfo_data_complete.csv')
    with open(path, newline='', encoding='utf-8') as h:
        return {r['name']: r for r in csv.DictReader(h, delimiter='|')}


def pd_of(r):
    return [[int(x) for x in q.split(',')] for q in re.findall(r'\[([\d, ]+)\]', r['pd_notation'])]


def orbit_job(job):
    name, u, det, crossings, pd = job
    orb = flype_orbit(from_pd([list(q) for q in pd]))
    return {'name': name, 'u': u, 'det': det, 'crossings': crossings, 'codes': sorted(orb), 'alex': list(fingerprint(pd))}


def main() -> None:
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    rows = knotinfo_rows()
    T = json.load(open(os.path.join(CODE, 'upper_bounds', 'data', 'targets.json')))
    T23 = [t for t in T if t['k'] in (2, 3)]
    json.dump(T23, open(os.path.join(HERE, 'data', 'targets_23_34.json'), 'w'))
    print('targets_23_34.json:', len(T23), 'knots')

    known = json.load(open(os.path.join(CODE, 'upper_bounds', 'data', 'known_u.json')))
    jobs = [(x['name'], x['u'], x['det'], x['crossings'], pd_of(rows[x['name']])) for x in known]
    out, t0 = [], time.time()
    from multiprocessing import Pool
    with Pool(max(1, (os.cpu_count() or 2) - 2)) as pool:
        for i, rec in enumerate(pool.imap(orbit_job, jobs, chunksize=16)):
            out.append(rec)
            if (i + 1) % 1000 == 0:
                print(f'  {i + 1}/{len(known)} orbits ({time.time() - t0:.0f}s)', flush=True)
    json.dump({'knotinfo_version': '2026.8.1', 'knots': out}, open(os.path.join(HERE, 'data', 'known_u_codes.json'), 'w'))
    print('known_u_codes.json:', len(out), 'knots,', sum(len(x['codes']) for x in out), 'codes')

    ctrl = []
    for nm, k in [('13n_30', 2), ('13n_45', 2), ('13n_80', 2), ('13n_2379', 2), ('13n_2809', 2), ('13n_2907', 2),
                  ('13n_3033', 2), ('13n_3589', 2), ('13a_647', 3)]:
        r = rows[nm]
        ctrl.append({'name': nm, 'k': k, 'range': [k, k + 1], 'crossings': int(r['crossing_number']),
                     'det': int(r['determinant']), 'alternating': r['alternating'], 'control': True, 'pd': pd_of(r)})
    json.dump(ctrl, open(os.path.join(HERE, 'data', 'controls.json'), 'w'))
    print('controls.json:', [c['name'] for c in ctrl])


if __name__ == '__main__':
    main()
