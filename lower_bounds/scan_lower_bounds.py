#!/usr/bin/env python
"""Classical u=1 obstructions over all KnotInfo rows with a gap in the unknotting interval.

For each gapped knot: Lickorish linking-form test (from the Tait graph Goeritz matrix) and the
Casson–Walker (Mullins/Boyer–Lines) integrality test (Jones polynomial via Regina from the PD code,
signature from KnotInfo); both signs handled; a knot is "obstructed" when no surgery sign survives.
Section 3, "Linking pairings", reports the Lickorish bound. Casson--Walker
values remain diagnostic auxiliary output. Signs from a checkerboard form
and from a Jones/signature formula cannot be intersected without identifying
their orientations. Only failure of an individual necessary condition is
used below; passing either test never proves u = 1.

    python lower_bounds/scan_lower_bounds.py --out lower_bounds_scan.csv [--workers 10]
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs'))


def worker(name: str) -> dict:
    import tait  # noqa: F401
    from spherogram import Link
    from tait import knotinfo, graph as tg, invariants as inv
    r = knotinfo.row(name)
    lo, hi = knotinfo.unknotting_interval(name)
    rec = {'knot': name, 'crossings': r['crossing_number'], 'interval': f'[{lo},{hi}]', 'det': r['determinant'],
           'signature': r['signature'], 'alternating': r['alternating']}
    try:
        g = tg.from_link(Link(knotinfo.pd_code(name)), name)
        D = inv.determinant(g)
        rec['det_goeritz'] = D
        if D != int(r['determinant']):
            raise ValueError('Goeritz and recorded determinants disagree')
        SL = inv.lickorish_allowed(g)
        J = inv.jones_regina(g)
        SC = inv.casson_walker_allowed(J, int(r['signature']))
        rec['lickorish_allowed'] = sorted(SL) if SL is not None else None
        rec['cw_allowed'] = sorted(SC)
        rec['lickorish_obstructs'] = (SL == set())
        rec['cw_obstructs'] = (SC == set())
        rec['combined_obstructs'] = None  # no orientation identification has been certified
        rec['u1_obstructed'] = bool(rec['lickorish_obstructs'] or rec['cw_obstructs'])
        if lo == 1 and rec['u1_obstructed']:
            rec['new_interval'] = f'[2,{hi}]' if hi > 2 else '2'
        else:
            rec['new_interval'] = rec['interval']
    except Exception as e:
        rec['error'] = f'{type(e).__name__}: {e}'
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lower_bounds_scan.csv'))
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--only-lower1', action='store_true', help='only rows with lower bound 1')
    args = ap.parse_args()
    from tait import knotinfo
    names = knotinfo.gapped_names()
    if args.only_lower1:
        names = [n for n in names if knotinfo.unknotting_interval(n)[0] == 1]
    print(f'{len(names)} gapped rows', flush=True)
    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(worker, n) for n in names]
        for i, f in enumerate(as_completed(futs), 1):
            results.append(f.result())
            if i % 500 == 0:
                print(f'  {i}/{len(names)} ({time.time() - t0:.0f}s)', flush=True)
    results.sort(key=lambda r: (int(r['crossings']) if str(r['crossings']).isdigit() else 99, r['knot']))
    keys = ['knot', 'crossings', 'alternating', 'interval', 'det', 'signature', 'lickorish_allowed', 'cw_allowed',
            'lickorish_obstructs', 'cw_obstructs', 'combined_obstructs', 'u1_obstructed', 'new_interval', 'error']
    with open(args.out, 'w', newline='') as h:
        w = csv.DictWriter(h, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        for r in results:
            w.writerow({k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()})
    n_err = sum(1 for r in results if 'error' in r)
    low1 = [r for r in results if r['interval'].startswith('[1,')]
    obs = [r for r in low1 if r.get('u1_obstructed')]
    print(f'errors {n_err}; rows with lower bound 1: {len(low1)}; u=1 obstructed: {len(obs)} '
          f'(Lickorish {sum(1 for r in low1 if r.get("lickorish_obstructs"))}, CW {sum(1 for r in low1 if r.get("cw_obstructs"))}); '
          f'new exact u=2: {sum(1 for r in obs if r["new_interval"] == "2")}; -> {args.out}; {time.time() - t0:.0f}s')


if __name__ == '__main__':
    main()
