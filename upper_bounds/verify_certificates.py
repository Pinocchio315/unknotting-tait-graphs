#!/usr/bin/env python
"""Independently verify campaign certificates (run LOCALLY, unknot-venv).

Handles both schemas:
  SUCCESS_*.json  (in-run successes)         --- expansion replay + changes, then
      kind 'unknot': knot Floer homology of the changed diagram detects the unknot;
      kind 'known' : SnapPy isometry of the changed diagram to the claimed partner.
  CHAIN_*.json    (re-derived candidate routes) --- additionally replays the recorded
      reduction path (changed diagram -> partner diagram, exact PD match), then
      kind 'known' : SnapPy isometry of the partner diagram to the claimed partner;
                     the bound is j + (KnotInfo upper bound of the partner);
      kind 'chain' : applies the recorded partner_changes and checks the result is the
                     unknot by knot Floer homology; the bound is j + #partner_changes.

The verified conclusion is an upper bound. A certificate's reported lower end is
not a lower-bound proof. Partner upper bounds come from the frozen manuscript
snapshot (or an explicit --known-bounds file), never from partner_u alone.

    ~/.pyenv/versions/unknot-venv/bin/python verify_certificates.py <run_dir> [--out verified]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))  # the tait package

from xtait.graph import from_pd, to_pd  # noqa: E402
from xtait import moves as mv  # noqa: E402
from expand import decode_path  # noqa: E402
from certificate_checks import validate_marks, check_partner_bound  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402


def replay(pd, path):
    """Replay isotopy only: an uncounted crossing change would invalidate the bound."""
    g = from_pd([list(q) for q in pd])
    for kind, params in decode_path(path):
        if kind not in {'r1-', 'r2-p', 'r2-s', 'r3D', 'r3Y', 'flype',
                        'r1+loop', 'r1+pend', 'r2+p', 'r2+s', 'pass'}:
            raise ValueError(f'non-isotopy move in certificate path: {kind}')
        g = mv.apply_move(g, kind, params)
    return g


def apply_rows(g, rows):
    edge_ids = sorted(g.edges)
    for i in validate_marks(rows, len(edge_ids)):
        g = mv.crossing_change(g, edge_ids[i])
    return g


def is_unknot_hfk(pd) -> bool:
    if not pd:
        return True  # connected knot diagrams only; the empty graph denotes the unknot
    import tait  # noqa: F401  (sqlite shim)
    from tait import invariants as inv
    from spherogram import Link
    return inv.is_unknot(Link([list(q) for q in pd]))


def isometric_to(pd, partner_name, tries=12) -> bool:
    import tait  # noqa: F401
    from tait import knotinfo
    from spherogram import Link
    Ea = Link([list(q) for q in pd]).exterior()
    Eb = Link(knotinfo.pd_code(partner_name)).exterior()
    for _ in range(tries):
        try:
            if Ea.is_isometric_to(Eb):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False


def start_matches_name(pd, name) -> bool:
    """Authenticate the source independently of its self-reported certificate name."""
    import tait  # noqa: F401
    from tait import knotinfo
    reference = from_pd(knotinfo.pd_code(name))
    source = from_pd([list(q) for q in pd])
    if canonical_code(source) == canonical_code(reference):
        return True
    return isometric_to(pd, name)


def verify_upper_bound(cert, known_bounds):
    """Replay and independently verify an upper bound; raise on every missing premise."""
    kind = cert['kind']
    if kind not in ('unknot', 'known', 'chain'):
        raise ValueError(f'unknown certificate kind: {kind}')
    rows = validate_marks(cert['change_rows'], len(cert['witness_pd']), cert['j'])
    if cert.get('n_witness', len(cert['witness_pd'])) != len(cert['witness_pd']):
        raise ValueError('recorded witness size disagrees with PD')
    if not start_matches_name(cert['start_pd'], cert['name']):
        raise ValueError('starting diagram was not identified as the named knot')
    g = replay(cert['start_pd'], cert['expansion_path'])
    if to_pd(g) != [list(q) for q in cert['witness_pd']]:
        raise ValueError('expansion replay mismatch')
    changed = apply_rows(g, rows)
    check_pd = to_pd(changed)
    if 'reduction_path' in cert:
        red = replay(check_pd, cert['reduction_path'])
        if to_pd(red) != [list(q) for q in cert['partner_pd']]:
            raise ValueError('reduction replay mismatch')
        check_pd = to_pd(red)
    j = len(rows)
    if kind == 'unknot':
        if not is_unknot_hfk(check_pd):
            raise ValueError('changed diagram is not the unknot')
        return j, 'HFK genus 0 (unknot)'
    if kind == 'chain':
        # Use the actual replayed endpoint. An unrelated partner_pd must never
        # authenticate a chain if the connecting reduction path is absent.
        partner_rows = validate_marks(cert['partner_changes'], len(check_pd))
        final = apply_rows(from_pd(check_pd), partner_rows)
        if not is_unknot_hfk(to_pd(final)):
            raise ValueError('partner changes do not give the unknot')
        return j + len(partner_rows), 'replayed chain ends at HFK genus 0'
    partner = cert['partner']
    pu = check_partner_bound(cert.get('partner_u'), known_bounds[partner])
    if not isometric_to(check_pd, partner):
        raise ValueError('claimed partner was not confirmed by exterior isometry')
    return j + pu, f'isometric to {partner} (independent u <= {pu})'


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--out', default='verified', help='directory for verified witness files')
    ap.add_argument('--known-bounds', default=os.path.join(HERE, '..', 'results', 'paper_v1_1_snapshot.json'),
                    help='independent name -> [lo,hi] or upper-bound JSON')
    ap.add_argument('--paper', action='store_true', help='require exactly the eight main_v1.1 constructions; skip separate archived certificates')
    args = ap.parse_args()
    with open(args.known_bounds) as handle:
        known_bounds = json.load(handle)

    os.makedirs(args.out, exist_ok=True)
    certs = sorted(glob.glob(os.path.join(args.run_dir, 'SUCCESS_*.json'))
                   + glob.glob(os.path.join(args.run_dir, 'CHAIN_*.json')))
    if args.paper:
        from verify_presentation_diagrams import PAPER_KNOTS
        selected, names = [], []
        for path in certs:
            with open(path) as handle:
                name = json.load(handle)['name']
            if name in PAPER_KNOTS:
                selected.append(path)
                names.append(name)
        if sorted(names) != sorted(PAPER_KNOTS):
            raise SystemExit('Paper verification requires exactly one certificate for each of the eight knots.')
        certs = selected
    print(f'{len(certs)} certificates in {args.run_dir}')
    if not certs:
        raise SystemExit('No certificates found; no verification performed.')
    n_ok = 0
    for path in certs:
        cert = json.load(open(path))
        nm, j = cert['name'], cert['j']
        try:
            u, why = verify_upper_bound(cert, known_bounds)
        except (ValueError, KeyError, IndexError, RuntimeError, AssertionError) as exc:
            print(f'  FAIL {nm}: {exc}')
            continue
        n_ok += 1
        print(f'  ok   {nm} (kind {cert["kind"]}, j={j}): {why}; u <= {u}')
        rec = {'knot': nm, 'source': 'upper_bounds_202608', 'claim_u_le': u,
               'kind': cert['kind'], 'j': j, 'partner': cert.get('partner'),
               'best': {'pd_knot': cert['witness_pd'],
                        'change_rows_0based': cert['change_rows'],
                        'crossings': cert['n_witness']}}
        with open(os.path.join(args.out, f'{nm}.json'), 'w') as h:
            json.dump(rec, h, indent=1)
    print(f'{n_ok}/{len(certs)} verified; witness files in {args.out}/')
    if n_ok != len(certs):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
