#!/usr/bin/env python
"""Resolve the stored candidates of an enum_search run (run LOCALLY, unknot-venv).

A candidate is a crossing change whose result has the determinant and Alexander fingerprint of a
known knot with small enough u (or of the unknot) but was not identified exactly on the server.
Here the reduced diagram is compared with each claimed partner by SnapPy isometry (mirror-
insensitive), or proved trivial by knot Floer homology, and every confirmed candidate is written
as a SUCCESS_<knot>.json certificate (same format as the server's), verifiable with
verify_enum_certificates.py.

    python resolve_candidates.py <run_dir>
"""
from __future__ import annotations

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'code', 'tait_graphs'))


def isometric(pd_a, pd_b, tries=10):
    from spherogram import Link
    Ea = Link([list(q) for q in pd_a]).exterior()
    Eb = Link([list(q) for q in pd_b]).exterior()
    for _ in range(tries):
        try:
            if Ea.is_isometric_to(Eb):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False


def main() -> None:
    import tait  # noqa: F401
    from tait import knotinfo, invariants as inv
    from spherogram import Link
    run_dir = sys.argv[1]
    known_u = {x['name']: x['u'] for x in json.load(open(os.path.join(HERE, 'data', 'known_u_codes.json')))['knots']}
    n_knots = n_cand = n_ok = 0
    for path in sorted(glob.glob(os.path.join(run_dir, 'results_*.jsonl'))):
        for line in open(path):
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get('success') or not r.get('candidates'):
                continue
            if os.path.exists(os.path.join(run_dir, f"SUCCESS_{r['name']}.json")):
                continue
            n_knots += 1
            for c in r['candidates']:
                n_cand += 1
                partner, pu, good = None, None, False
                if c['kind'] == 'unknot':
                    good = inv.is_unknot(Link([list(q) for q in c['reduced_pd']]))
                else:
                    for nm in c['partners']:
                        if isometric(c['reduced_pd'], knotinfo.pd_code(nm)):
                            partner, pu, good = nm, known_u[nm], True
                            break
                if not good:
                    continue
                n_ok += 1
                u = c['j'] if c['kind'] == 'unknot' else c['j'] + pu
                cert = {'name': r['name'], 'kind': c['kind'], 'j': c['j'], 'u': u, 'partner': partner,
                        'partner_u': pu, 'partner_mirror': None, 'start_pd': None, 'moves': c['moves'],
                        'witness_pd': c['witness_pd'], 'n_witness': len(c['witness_pd']),
                        'change_rows': c['rows'], 'reduced_pd': c['reduced_pd'],
                        'identification': 'resolved locally (' + ('knot Floer homology' if c['kind'] == 'unknot' else 'SnapPy isometry') + ')'}
                targets = {t['name']: t for t in json.load(open(os.path.join(HERE, 'data', 'targets_23_34.json')))}
                controls = {t['name']: t for t in json.load(open(os.path.join(HERE, 'data', 'controls.json')))}
                cert['start_pd'] = (targets.get(r['name']) or controls[r['name']])['pd']
                cp = os.path.join(run_dir, f"SUCCESS_{r['name']}.json")
                json.dump(cert, open(cp, 'w'), indent=1)
                print(f"  *** {r['name']}: u <= {u} via candidate ({c['kind']}, j={c['j']}, partner={partner}) -> {cp}")
                break
    print(f'{n_knots} knots with candidates, {n_cand} candidates checked, {n_ok} settled')


if __name__ == '__main__':
    main()
