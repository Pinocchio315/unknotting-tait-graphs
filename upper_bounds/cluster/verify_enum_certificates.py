#!/usr/bin/env python
"""Verify enumeration certificates for the upper-bound searches (v1.6 Section 4).

Replay authenticates the isotopy path and the exact crossing-row correspondence.
Unknot endpoints use genus detection; named partners use exterior isometry checks
and independent upper bounds from the frozen manuscript snapshot. Numerical
isometry checks are not presented as formal proof certificates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent / 'tait_graphs'))
sys.path.insert(0, str(HERE))
from xtait.graph import from_pd, to_pd
from xtait import moves as mv
from moves_io import replay
from certificate_checks import validate_marks, check_partner_bound
from verify_certificates import start_matches_name, is_unknot_hfk, isometric_to


def verify_enum_upper_bound(cert, known):
    """Return the upper bound after checking all certificate premises."""
    if cert['kind'] not in ('unknot', 'known'):
        raise ValueError('unknown enumeration-certificate kind')
    rows = validate_marks(cert['change_rows'], len(cert['witness_pd']), cert['j'])
    if cert['n_witness'] != len(cert['witness_pd']):
        raise ValueError('recorded crossing count disagrees with the diagram')
    if not start_matches_name(cert['start_pd'], cert['name']):
        raise ValueError('source diagram does not identify as the named knot')
    graph = replay(from_pd([list(q) for q in cert['start_pd']]), cert['moves'])
    # Sorting each quadruple destroys crossing information: over/under switches
    # leave its four labels unchanged. Preserve their cyclic order and rows.
    if to_pd(graph) != [list(q) for q in cert['witness_pd']]:
        raise ValueError('replayed diagram differs from the stored witness')
    edges = sorted(graph.edges)
    for row in rows:
        graph = mv.crossing_change(graph, edges[row])
    changed_pd = to_pd(graph)
    if cert['kind'] == 'unknot':
        if not is_unknot_hfk(changed_pd):
            raise ValueError('changed diagram is not the unknot')
        return len(rows)
    partner = cert['partner']
    bound = check_partner_bound(cert.get('partner_u'), known[partner])
    if not isometric_to(changed_pd, partner):
        raise ValueError('partner identification failed')
    return len(rows) + bound


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_dir', type=Path)
    parser.add_argument('--known-bounds', type=Path,
                        default=HERE.parent.parent / 'results' / 'paper_v1_1_snapshot.json')
    args = parser.parse_args(argv)
    with args.known_bounds.open() as handle:
        known = json.load(handle)
    paths = sorted(args.run_dir.glob('SUCCESS_*.json'))
    if not paths:
        raise SystemExit('No certificates found; no verification performed.')
    failures = 0
    for path in paths:
        try:
            with path.open() as handle:
                cert = json.load(handle)
            bound = verify_enum_upper_bound(cert, known)
            print(f'  ok   {cert["name"]}: u <= {bound}', flush=True)
        except Exception as exc:
            failures += 1
            print(f'  FAIL {path.name}: {exc}', flush=True)
    print(f'{len(paths) - failures}/{len(paths)} certificates verified')
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
