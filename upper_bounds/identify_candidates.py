#!/usr/bin/env python
"""Batch-identify the unresolved candidates of a campaign run (run LOCALLY, unknot-venv).

Every results_*.jsonl line may carry up to 20 `candidates`: j < k crossing changes whose
result reduced to <= 13 crossings but whose canonical code matched no KnotInfo diagram.
This script de-duplicates them, identifies each reduced PD with the SnapPy census (isometry
fallback against determinant-matched KnotInfo knots), looks up the partner's unknotting
number (KnotInfo + our updates via data/known_u.json), and reports every case where

    u(target) <= j + u(partner) <= k        --- i.e. the target is SETTLED.

Non-hyperbolic / unidentified results are listed for manual inspection (a composite partner
cannot have u = 1 by Scharlemann, so it can only matter for k = 3 routes).

    ~/.pyenv/versions/unknot-venv/bin/python identify_candidates.py <run_dir_or_jsonl...> \
        [--out candidates_identified.jsonl]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

from xtait.graph import from_pd  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from flipdet import exact_det  # noqa: E402


def census_to_knotinfo(s, known_names):
    """'9_25(0,0)' / 'K11a105(0,0)' / 'K13n13(0,0)' -> KnotInfo-style name, or None.

    CAUTION (Perko): SnapPy's Rolfsen census keeps the pre-Perko numbering, so its
    10-crossing names from 10_162 on are shifted against KnotInfo (census 10_164 is
    KnotInfo 10_163).  Any name-based match is therefore only a HYPOTHESIS; a claim
    must be confirmed by explicit isometry against the KnotInfo diagram before use."""
    s = re.sub(r'\(.*\)$', '', s.strip())
    cands = [s]
    if re.fullmatch(r'K\d+[an]\d+', s):
        t = s[1:]
        cands += [t, re.sub(r'([an])', r'\1_', t, count=1)]
    for c in cands:
        if c in known_names:
            return c
    return None


_EXT_CACHE = {}


def confirmed_partner(pd, name, knotinfo, Link, tries=8):
    """Explicit SnapPy isometry of the reduced diagram against the KnotInfo diagram of
    `name` --- the only accepted proof of partner identity (census names can mislead)."""
    if name not in _EXT_CACHE:
        _EXT_CACHE[name] = Link(knotinfo.pd_code(name)).exterior()
    Eb = _EXT_CACHE[name]
    Ea = Link([list(q) for q in pd]).exterior()
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
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+', help='run directories or results_*.jsonl files')
    ap.add_argument('--out', default='candidates_identified.jsonl')
    args = ap.parse_args()

    import tait  # noqa: F401  (sqlite shim)
    from tait import knotinfo
    from spherogram import Link

    known = json.load(open(os.path.join(HERE, 'data', 'known_u.json')))
    u_of = {x['name']: x['u'] for x in known}          # exact u <= 2 (our updates included)
    all_names = set(knotinfo.rows())

    files = []
    for r in args.runs:
        files += sorted(glob.glob(os.path.join(r, 'results_*.jsonl'))) if os.path.isdir(r) else [r]

    # gather and de-duplicate candidates by (target, canonical code of the reduced diagram)
    jobs = {}
    for path in files:
        for line in open(path):
            try:
                rec = json.loads(line)
            except Exception:
                continue
            for c in rec.get('candidates', []):
                code = repr(canonical_code(from_pd([list(q) for q in c['reduced_pd']])))
                jobs.setdefault((rec['name'], code), (rec, c))
    print(f'{len(jobs)} distinct (target, reduced-knot) pairs from {len(files)} result files')

    settled, open_partner, no_gain, unidentified = [], [], [], []
    id_cache = {}          # canonical code -> (partner_or_None, census ids)
    done = 0
    with open(args.out, 'w') as out:
        for (target, code), (rec, c) in sorted(jobs.items()):
            k, j = rec['k'], c['j']
            pd = [list(q) for q in c['reduced_pd']]
            det = exact_det(from_pd(pd))
            if code in id_cache:
                partner, ids = id_cache[code]
            else:
                partner = None
                try:
                    ids = [str(m) for m in Link(pd).exterior().identify()]
                except Exception:
                    ids = []
                for s in ids:
                    partner = census_to_knotinfo(s, all_names)
                    if partner:
                        break
                id_cache[code] = (partner, ids)
            done += 1
            if done % 500 == 0:
                print(f'  {done}/{len(jobs)} identified ({len(id_cache)} distinct)', flush=True)
            row = {'target': target, 'k': k, 'j': j, 'rows': c['rows'], 'det': det,
                   'n_reduced': c['n_reduced'], 'partner': partner, 'census': ids[:2]}
            if partner is None:
                unidentified.append(row)
            else:
                u = u_of.get(partner)
                if u is None:                          # partner itself open (or u >= 3)
                    lohi = knotinfo.unknotting_interval(partner)
                    row['partner_u'] = list(lohi)
                    open_partner.append(row)
                elif j + u <= k:
                    if not confirmed_partner(pd, partner, knotinfo, Link):
                        row['partner'] = None
                        row['census_hypothesis'] = partner
                        unidentified.append(row)
                        out.write(json.dumps(row) + '\n')
                        continue
                    row['partner_u'] = u
                    row['claim'] = f'u({target}) <= {j + u}  SETTLES {rec["range"]}  (isometry-confirmed)'
                    settled.append(row)
                else:
                    row['partner_u'] = u
                    no_gain.append(row)
            out.write(json.dumps(row) + '\n')

    print(f'\nsettled: {len(settled)}   no-gain (j+u > k): {len(no_gain)}   '
          f'open partner: {len(open_partner)}   unidentified: {len(unidentified)}')
    for row in settled:
        print('  *** ', row['claim'], f"(1 change of row {row['rows']} -> {row['partner']})")
    if unidentified:
        print('\nunidentified (non-hyperbolic / composite / census miss) --- inspect manually:')
        for row in unidentified[:20]:
            print(f"   {row['target']}: j={row['j']} det={row['det']} n={row['n_reduced']}")
    if open_partner:
        pairs = sorted({(r['target'], r['partner']) for r in open_partner})
        print(f'\n{len(pairs)} target/open-partner adjacencies (u differs by <= 1 across a '
              f'crossing change) recorded in {args.out}')


if __name__ == '__main__':
    main()
