#!/usr/bin/env python
"""Sweep torsion-order lower bounds for the unknotting number over gapped KnotInfo rows.

For every row with unknotting interval [lo, b] (default lo = 1) compute the maximal U-torsion order of
HFK^-(K; F_2) (`tait.invariants.hfk_minus_torsion_order`; Alishahi–Eftekhary: Ord_U <= u).  Rows with
Ord_U > lo get a new lower bound (e.g. [1,2] with Ord 2 -> u = 2 exact).  Results are appended to a
resumable JSONL and summarised. Section 3, "Torsion bounds", uses
this U-torsion invariant. Rational Lee X-torsion and characteristic-two
Bar--Natan h-torsion are reviewed in the Background but are not outputs of
this script. A failed computation is recorded and retried on resumption;
it supplies no lower bound.

    python lower_bounds/scan_torsion_lower_bounds.py --lo 1 --workers 6 --out results_torsion_sweep
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs'))


def one(job):
    nm, prime = job
    import tait  # noqa: F401
    from tait import knotinfo, invariants as inv
    t0 = time.time()
    try:
        pd = knotinfo.pd_code(nm)
        o = inv.hfk_minus_torsion_order(pd, prime=prime)
        return {'name': nm, 'prime': prime, 'ord_u': o, 'seconds': round(time.time() - t0, 2)}
    except Exception as e:
        return {'name': nm, 'prime': prime, 'error': repr(e)}


def latest_records(records):
    """One current record per knot; appended retries replace earlier failures."""
    return {record['name']: record for record in records}


def summarise(records, lower, prime=2):
    """Exclude failures from the order distribution and count each knot once."""
    from collections import Counter
    current = latest_records(records)
    successful = [row for row in current.values() if 'error' not in row and 'ord_u' in row]
    news = sorted((row for row in successful if row['ord_u'] > lower), key=lambda row: row['name'])
    return {'interval_lo': lower, 'prime': prime, 'rows': len(current),
            'successful_rows': len(successful), 'failed_rows': len(current) - len(successful),
            'ord_distribution': dict(Counter(row['ord_u'] for row in successful)),
            'new_lower_bounds': [{'name': row['name'], 'ord_u': row['ord_u']} for row in news]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--lo', type=int, default=1)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--out', default='results_torsion_sweep')
    ap.add_argument('--names-file')
    ap.add_argument('--prime', type=int, default=2)
    args = ap.parse_args()
    import tait  # noqa: F401
    from tait import knotinfo
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs')
    out = os.path.join(root, args.out)
    os.makedirs(out, exist_ok=True)
    res_path = os.path.join(out, f'hfk_torsion_lo{args.lo}' + ('' if args.prime == 2 else f'_p{args.prime}') + '.jsonl')
    done = set()
    if os.path.exists(res_path):
        with open(res_path) as handle:
            current = latest_records(json.loads(line) for line in handle if line.strip())
        done = {name for name, row in current.items() if 'error' not in row and 'ord_u' in row}
    if args.names_file:
        names = [x.strip() for x in open(args.names_file) if x.strip()]
    else:
        names = []
        for nm, r in knotinfo.rows().items():
            m = re.fullmatch(r'\[%d[,;](\d+)\]' % args.lo, r['unknotting_number'].replace(' ', ''))
            if m:
                names.append(nm)
    todo = [nm for nm in names if nm not in done]
    print(f'{len(names)} rows with interval [{args.lo},b], {len(todo)} to compute', flush=True)
    t0 = time.time()
    found = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex, open(res_path, 'a') as h:
        futs = {ex.submit(one, (nm, args.prime)): nm for nm in todo}
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            h.write(json.dumps(r) + '\n')
            h.flush()
            if r.get('ord_u', 0) > args.lo:
                found.append(r)
                print(f'  ** {r["name"]}: Ord_U = {r["ord_u"]} > {args.lo}  -> new lower bound', flush=True)
            if i % 200 == 0:
                print(f'  {i}/{len(todo)} ({time.time()-t0:.0f}s)', flush=True)
    # summary over the full jsonl
    rows = [json.loads(ln) for ln in open(res_path)]
    summary = summarise(rows, args.lo, args.prime)
    news = summary['new_lower_bounds']
    with open(os.path.join(out, f'summary_lo{args.lo}' + ('' if args.prime == 2 else f'_p{args.prime}') + '.json'), 'w') as h:
        json.dump(summary, h, indent=1)
    print(f"done: {summary['rows']} knots, {summary['failed_rows']} failures, "
          f"Ord distribution {summary['ord_distribution']}, new lower bounds: {len(news)}")
    for r in news[:40]:
        print('   ', r['name'], 'Ord_U =', r['ord_u'])
    print(f'-> {out}/summary_lo{args.lo}.json ({time.time()-t0:.0f}s)')


if __name__ == '__main__':
    main()
