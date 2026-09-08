#!/usr/bin/env python
"""Re-derive FULL certificates for the targets settled through the j < k candidate routes.

The campaign's per-knot candidate records store only the reduced partner diagram, not the
expanded witness diagram it came from, so the claims cannot be replayed from them alone.
This script re-runs a targeted hunt for each settled target and records the complete chain:

    KnotInfo PD --expansion path--> witness diagram --j changes--> changed diagram
               --reduction path--> partner diagram --[identification | k-j more changes]--> claim

Two claim kinds (matching resolve_unidentified.py):
  'known':  the partner diagram canonically matches a KnotInfo knot P, so
            u(target) <= j + hi(P) (KnotInfo's upper bound for P is certified);
  'chain':  the partner diagram itself is unknotted by <= k-j recorded changes, so
            u(target) <= j + (k-j) = k with no identification at all.

Claims are read from identify_candidates.py / resolve_unidentified.py outputs.  Certificates
are written as CHAIN_<target>.json; verify with verify_certificates.py (handles both kinds).

    ~/.pyenv/versions/unknot-venv/bin/python rehunt_certificates.py \
        xup_runs/candidates_identified.jsonl xup_runs/unidentified_resolved.jsonl \
        --out-dir xup_runs
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

from xtait.graph import from_pd, to_pd  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.reduce import greedy_reduce, search_reduce, ReduceConfig  # noqa: E402
from flipdet import FlipDet, det_pairs  # noqa: E402
from expand import expand_diagrams, jsonable_path  # noqa: E402
from resolve_unidentified import direct_unknotting  # noqa: E402


def reduce_with_path(g, nodes=600, seconds=1.5):
    red, path = greedy_reduce(g)
    if 13 < red.n_crossings() <= 20:      # cheap push into identification range only
        cfg = ReduceConfig(target=13, max_extra=2, max_nodes=nodes, time_limit=seconds,
                           heap_max=20000, log_every=10 ** 9)
        red2, more, _ = search_reduce(red, cfg)
        if red2.n_crossings() < red.n_crossings():
            red, path = red2, path + more
    return red, path


def collect_claims(files):
    """One claim per target: {'target', 'k', 'j', 'det', 'kind', 'partner'?}."""
    claims = {}
    for path in files:
        if not os.path.exists(path):
            continue
        for line in open(path):
            x = json.loads(line)
            if 'claim' not in x:
                continue
            kind = 'chain' if x.get('mechanism') == 'direct' else 'known'
            c = {'target': x['target'], 'k': x['k'], 'j': x['j'], 'det': x['det'],
                 'kind': kind, 'partner': x.get('partner')}
            prev = claims.get(x['target'])
            if prev is None or (prev['kind'] == 'chain' and kind == 'known'):
                claims[x['target']] = c
    return list(claims.values())


def _isometric(pd, partner_exterior, tries=8):
    from spherogram import Link
    Ea = Link([list(q) for q in pd]).exterior()
    for _ in range(tries):
        try:
            if Ea.is_isometric_to(partner_exterior):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        partner_exterior.randomize()
    return False


def hunt_one(t, claim, allcodes, seed=0, max_diagrams=20000, time_cap=1200.0):
    k, j, det_p = t['k'], claim['j'], claim['det']
    pairs = det_pairs([det_p])
    g0 = from_pd([list(q) for q in t['pd']], t['name'])
    partner_codes = partner_ext = None
    if claim['kind'] == 'known':
        partner_codes = set(next(x for x in allcodes if x['name'] == claim['partner'])['codes'])
        import tait  # noqa: F401
        from tait import knotinfo
        from spherogram import Link
        partner_ext = Link(knotinfo.pd_code(claim['partner'])).exterior()
    t0 = time.time()
    tried = 0
    # n_target stays at 40 (NOT the campaign cutoff of 50): the deposited 13a647 chain was
    # re-derived at this size, and its witness diagram has 39 crossings (40 minus the parity
    # correction).  Raising it here would re-derive a different, larger certificate.
    for g, path in itertools.chain([(g0, [])],
                                   expand_diagrams(g0, n_target=40, count=max_diagrams,
                                                   seed=seed)):
        if time.time() - t0 > time_cap:
            return None
        try:
            F = FlipDet(g)
        except ZeroDivisionError:
            continue
        m = len(F.edge_ids)
        for rows in itertools.combinations(range(m), j):
            if F.det_after(rows) not in pairs:
                continue
            h = g
            for i in rows:
                h = mv.crossing_change(h, F.edge_ids[i])
            red, rpath = reduce_with_path(h)
            tried += 1
            if not 0 < red.n_crossings() <= 13:
                continue
            cert = {'name': t['name'], 'range': t['range'], 'k': k, 'j': j,
                    'kind': claim['kind'], 'start_pd': t['pd'],
                    'expansion_path': jsonable_path(path), 'witness_pd': to_pd(g),
                    'change_rows': list(rows), 'n_witness': m, 'det_start': t['det'],
                    'reduction_path': jsonable_path(rpath), 'partner_pd': to_pd(red)}
            if claim['kind'] == 'known':
                if (repr(canonical_code(red)) in partner_codes
                        or _isometric(cert['partner_pd'], partner_ext)):
                    cert['partner'] = claim['partner']
                    return cert
            else:
                changes = direct_unknotting(red, k - j)
                if changes is not None:
                    cert['partner_changes'] = changes
                    cert['u'] = j + len(changes)
                    return cert
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('claim_files', nargs='+')
    ap.add_argument('--out-dir', default='xup_runs')
    ap.add_argument('--names', nargs='*', help='restrict to these targets')
    args = ap.parse_args()

    targets = {t['name']: t for t in json.load(open(os.path.join(HERE, 'data', 'targets.json')))}
    allcodes = json.load(open(os.path.join(HERE, 'data', 'all_knot_codes.json')))
    claims = collect_claims(args.claim_files)
    if args.names:
        wanted = {n.replace('_', '') for n in args.names}
        claims = [c for c in claims if c['target'].replace('_', '') in wanted]
    print(f'{len(claims)} settled targets to certify')

    for claim in claims:
        t = targets[claim['target']]
        found = None
        for seed in range(3):
            found = hunt_one(t, claim, allcodes, seed=seed)
            if found:
                break
        if found:
            if 'u' not in found:
                found['u'] = None      # filled by verify_certificates from KnotInfo hi
            p = os.path.join(args.out_dir, f"CHAIN_{claim['target']}.json")
            with open(p, 'w') as h:
                json.dump(found, h)
            print(f"  ok   {claim['target']} ({claim['kind']}"
                  + (f" -> {claim['partner']}" if claim['partner'] else '') + f') -> {p}',
                  flush=True)
        else:
            print(f"  MISS {claim['target']} --- rerun with more diagrams/seeds", flush=True)


if __name__ == '__main__':
    main()
