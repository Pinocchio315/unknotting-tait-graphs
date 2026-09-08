#!/usr/bin/env python
"""Offline HFK backstop for the j = k route (run LOCALLY, unknot-venv).

The in-run reducer proves the unknot only by reaching 0 crossings; pass moves alone can
stall on a true unknot (ReAPR escapes such stalls by geometric re-embedding).  We instead
close the gap with something stronger: every det-1 candidate whose reduction stalled at a
small diagram is stored in results_*.jsonl (`unknot_near`), and this script decides each one
OUTRIGHT with knot Floer homology --- Seifert genus 0 detects the unknot, no move path
needed.  Any hit settles u(target) = k.

    ~/.pyenv/versions/unknot-venv/bin/python scan_nearmiss_hfk.py <run_dir> [--max-crossings 16]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

from xtait.graph import from_pd  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--max-crossings', type=int, default=16)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    out_path = args.out or os.path.join(args.run_dir, 'nearmiss_hfk.jsonl')

    import tait  # noqa: F401  (sqlite shim)
    from tait import invariants as inv
    from spherogram import Link

    jobs = {}
    for path in sorted(glob.glob(os.path.join(args.run_dir, 'results_*.jsonl'))):
        for line in open(path):
            if not line.strip():
                continue
            rec = json.loads(line)
            for c in rec.get('unknot_near', []):
                if c['n_reduced'] > args.max_crossings:
                    continue
                code = repr(canonical_code(from_pd([list(q) for q in c['reduced_pd']])))
                jobs.setdefault((rec['name'], code), (rec, c))
    print(f'{len(jobs)} distinct stored near-misses to decide by HFK')

    settles = []
    verdict = {}                      # canonical code -> HFK verdict (shared across targets)
    with open(out_path, 'w') as out:
        for n_done, ((target, code), (rec, c)) in enumerate(sorted(jobs.items()), 1):
            if code in verdict:
                is_u = verdict[code]
            else:
                try:
                    is_u = inv.is_unknot(Link([list(q) for q in c['reduced_pd']]))
                except Exception as exc:
                    is_u = None
                    print(f'  ! {target}: HFK failed ({exc})')
                verdict[code] = is_u
            row = {'target': target, 'k': rec['k'], 'rows': c['rows'],
                   'n_reduced': c['n_reduced'], 'hfk_unknot': is_u}
            if is_u:
                row['claim'] = f"u({target}) = {rec['k']}  (det-1 near-miss IS the unknot, HFK-proved)"
                settles.append(row)
                print('  ***', row['claim'])
            out.write(json.dumps(row) + '\n')
            if n_done % 1000 == 0:
                print(f'  {n_done}/{len(jobs)} ({len(verdict)} distinct HFK calls)', flush=True)

    print(f'\nHFK backstop: {len(settles)} settles among {len(jobs)} near-misses '
          f'-> {out_path}')


if __name__ == '__main__':
    main()
