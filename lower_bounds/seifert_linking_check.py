#!/usr/bin/env python
"""Second route to the linking-pairing obstruction (Section 3, "Linking pairings"): KnotInfo Seifert matrices.

For a Seifert matrix V of K, the symmetric matrix M = V + V^T presents H_1(Sigma_2(K)) and
+-M^{-1} is its linking pairing.  This script re-derives the Lickorish test for every knot whose
recorded range is [1, b] from M instead of from the Goeritz matrix of the Tait graph
(scan_lower_bounds.py), and compares the obstructed set with the deposited 815-knot list.
The test is sign-agnostic (the condition is <+-2/D>), so the overall sign of M is irrelevant.

    python seifert_linking_check.py [--workers 6] [--out seifert_check.json]
"""
from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import os
import re
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))


def knotinfo_rows() -> dict:
    import database_knotinfo
    root = os.path.dirname(database_knotinfo.__file__)
    path = next(os.path.join(dp, f) for dp, dn, fn in os.walk(root) for f in fn
                if f == 'knotinfo_data_complete.csv')
    with open(path, newline='', encoding='utf-8') as h:
        return {r['name']: r for r in csv.DictReader(h, delimiter='|')}


def h1_torsion(M) -> list[int]:
    """Nonunit Smith factors (0 denotes a free summand) of coker(M) = H_1(Sigma_2(K))."""
    from sympy import ZZ
    from sympy.matrices.normalforms import smith_normal_form
    n = M.shape[0]
    snf = smith_normal_form(M, domain=ZZ)
    return [abs(int(snf[i, i])) for i in range(n) if abs(int(snf[i, i])) != 1]


# Shared exact algebra; independence of the Seifert and Goeritz checks is
# independence of their presentations, not of this finite-group routine.
try:
    from .linking_pairing import generator_self_linking
except ImportError:  # direct script invocation
    from linking_pairing import generator_self_linking


def verdict(args):
    """Lickorish test from the Seifert matrix: obstructed iff H_1 is not cyclic or no sign s
    makes the self-linking of a generator equal to s * 2 * k^2 / D."""
    name, seifert, det_recorded = args
    from sympy import Matrix
    V = Matrix(ast.literal_eval(seifert))
    M = V + V.T
    D = abs(int(M.det()))
    rec = {'knot': name, 'det': D, 'det_matches_knotinfo': D == det_recorded}
    if D != det_recorded or D < 1 or D % 2 == 0:
        rec.update(obstructed=None, error='invalid or mismatched knot determinant')
        return rec
    if D == 1:
        rec.update(cyclic=True, allowed_signs=[1, -1], obstructed=False)
        return rec
    tors = h1_torsion(M)
    rec['h1'] = tors
    if tors != [D]:
        rec.update(cyclic=False, allowed_signs=[], obstructed=True)
        return rec
    a = generator_self_linking(M, D)
    if a is None:
        rec.update(cyclic=True, allowed_signs=None, obstructed=None, note='no generator found')
        return rec
    squares = {(k * k) % D for k in range(1, D) if math.gcd(k, D) == 1}
    ia = pow(a, -1, D)
    allowed = [s for s in (1, -1) if (s * 2 * ia) % D in squares]
    rec.update(cyclic=True, self_linking=a, allowed_signs=allowed, obstructed=not allowed)
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--out', default=os.path.join(HERE, 'seifert_check.json'))
    args = ap.parse_args()

    rows = knotinfo_rows()
    targets = [n for n, r in rows.items()
               if re.fullmatch(r'\[1,\d+\]', r['unknotting_number'].replace(' ', ''))
               and r['seifert_matrix'].strip()]
    obs_path = os.path.join(HERE, '..', 'results', 'lickorish_cw_obstructed_815_knotinfo2026.8.1.json')
    deposited = {n for n, _ in json.load(open(obs_path))}
    print(f'{len(targets)} knots with range [1,b]; deposited obstructed list: {len(deposited)}', flush=True)

    jobs = [(n, rows[n]['seifert_matrix'], int(rows[n]['determinant'])) for n in targets]
    t0 = time.time()
    results = []
    with Pool(args.workers) as pool:
        for k, rec in enumerate(pool.imap_unordered(verdict, jobs, chunksize=8), 1):
            results.append(rec)
            if k % 500 == 0:
                print(f'  {k}/{len(jobs)}  {time.time() - t0:.0f}s', flush=True)
    results.sort(key=lambda r: r['knot'])

    obstructed = {r['knot'] for r in results if r['obstructed']}
    undecided = [r['knot'] for r in results if r['obstructed'] is None]
    det_bad = [r['knot'] for r in results if not r['det_matches_knotinfo']]
    print(f'obstructed from Seifert matrices: {len(obstructed)}; undecided: {len(undecided)}; '
          f'det mismatches: {len(det_bad)}')
    print(f'  = deposited list: {obstructed == deposited}'
          f'  (missing {sorted(deposited - obstructed)[:10]}, extra {sorted(obstructed - deposited)[:10]})')
    with open(args.out, 'w') as h:
        json.dump({'targets': len(targets), 'obstructed': sorted(obstructed),
                   'agrees_with_deposited_815': obstructed == deposited,
                   'undecided': undecided, 'det_mismatches': det_bad, 'rows': results}, h, indent=1)
    print('wrote', args.out)


if __name__ == '__main__':
    main()
