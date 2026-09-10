#!/usr/bin/env python
"""Freeze the input data for the correction-term sweep over knots whose double branched cover is an L-space.

The sweep applies the obstructions of Section 3 to non-alternating knots. Reduced Khovanov homology
over F_2 must have rank equal to the determinant, which forces the double branched cover to be an
L-space. Greene's spanning-tree model then supplies finite sets containing the correction terms.
The present implementation also requires cyclic first homology and the signature conditions below.

A target is a knot whose recorded range is still open and whose signature matches the shortest sequence
allowed by that range. For n = 1 and |sigma| <= 2 the test is the half-integral surgery pattern of Ni
and Wu. For n = 2, 3, 4 with |sigma| = 2n it is the enumeration of definite forms of that rank.
Excluding the tested length raises the lower bound to n + 1.

A control is a knot with the same hypotheses whose unknotting number is already known and equals n.  The
obstruction must not fire on a control; an OBSTRUCTED control is an implementation error.

To reproduce the deposited cohort, pass --table with a separately reconstructed v1.2 table, which
includes the two earlier Greene results but precedes the sweep. The default generated/u_table.json
may already contain the sweep's conclusions and is unsuitable for reproducing its input selection.
See README.md for the reconstruction command. Replay workers can instead read the deposited frozen
greene_targets_2026-09-08.json directly, without rebuilding the cohort.

    python build_greene_targets.py                       # writes data/greene_targets.json
    python build_greene_targets.py --max-det 400         # restrict to smaller discriminant groups
"""
from __future__ import annotations
import argparse, ast, hashlib, json, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def thin_over_f2(row):
    """Rank of the reduced mod-2 Khovanov homology and its diagonals, or None when the entry is unusable."""
    raw = row.get('khovanov_reduced_mod2_vector', '')
    if not str(raw).strip():
        return None
    try:
        entries = ast.literal_eval(raw)
    except (SyntaxError, ValueError):
        return None
    if not entries or any(len(e) != 4 or e[0] != 2 or any(type(v) is not int for v in e) or e[1] < 0
                          for e in entries):
        return None
    return {'rank': sum(e[1] for e in entries),
            'diagonals_q_minus_2h': sorted({e[3] - 2 * e[2] for e in entries}),
            'raw': raw}


def cyclic_first_homology(row):
    """Invariant factors of H_1 of the double branched cover from the tabulated torsion numbers.

    Greene's model indexes the Spin^c structures by their first Chern class in a cyclic group, so a
    non-cyclic group is outside the present implementation.  For unknotting number one such a group is
    already excluded by Lickorish's condition, which is applied separately.
    """
    raw = row.get('torsion_numbers', '')
    if not str(raw).strip():
        return None
    try:
        table = dict((int(m), list(factors)) for m, factors in ast.literal_eval(raw))
    except (SyntaxError, ValueError, TypeError):
        return None
    factors = table.get(2)
    return None if factors is None else len(factors) <= 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table', type=Path, default=ROOT / 'generated' / 'u_table.json',
                    help='consolidated ranges; the sweep is measured against these')
    ap.add_argument('--out', type=Path, default=HERE / 'data' / 'greene_targets.json')
    ap.add_argument('--max-det', type=int, default=0, help='skip knots of larger determinant (0 = no limit)')
    args = ap.parse_args()

    import database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    table = json.loads(args.table.read_text())

    targets, controls, skipped = {}, {}, {'thick_or_missing': 0, 'alternating': 0, 'determinant': 0,
                                         'signature_mismatch': 0, 'closed_range': 0,
                                         'noncyclic_homology': 0}
    for name, row in rows.items():
        if name not in table:
            continue
        if str(row.get('alternating', '')).upper().startswith('Y'):
            skipped['alternating'] += 1            # the Goeritz form already gives the correction terms
            continue
        det = int(row['determinant'])
        if args.max_det and det > args.max_det:
            skipped['determinant'] += 1
            continue
        kh = thin_over_f2(row)
        if kh is None or kh['rank'] != det:
            skipped['thick_or_missing'] += 1        # the L-space hypothesis is not certified
            continue
        if cyclic_first_homology(row) is not True:
            skipped['noncyclic_homology'] += 1      # outside the cyclic Spin^c labelling used here
            continue
        lo, hi = table[name]
        sigma = int(row['signature'])
        n = lo if lo >= 1 else 1
        if abs(sigma) != 2 * n and not (n == 1 and abs(sigma) <= 2):
            skipped['signature_mismatch'] += 1      # the shortest allowed sequence is not sign-determined
            continue
        entry = {'name': name, 'pd': [list(map(int, c)) for c in ast.literal_eval(row['pd_notation'])],
                 'determinant': det, 'signature': sigma, 'range': [lo, hi],
                 'crossing_number': int(row['crossing_number']),
                 'test': f'rank{n}' if n >= 2 else 'u1',
                 'reduced_mod2_vector_raw': kh['raw'],
                 'khovanov_rank': kh['rank'], 'diagonals_q_minus_2h': kh['diagonals_q_minus_2h']}
        if lo == hi:
            if lo == n:
                entry['known_unknotting_number'] = lo
                controls[name] = entry                # the obstruction must not fire
            else:
                skipped['closed_range'] += 1
        else:
            targets[name] = entry

    order = sorted(targets, key=lambda n: (targets[n]['crossing_number'], targets[n]['determinant']))
    corder = sorted(controls, key=lambda n: (controls[n]['crossing_number'], controls[n]['determinant']))
    payload = {'built': '2026-09-08',
               'source': {'database': 'KnotInfo', 'package': 'database_knotinfo',
                          'package_version': __import__('importlib.metadata').metadata.version('database_knotinfo'),
                          'ranges': str(args.table.relative_to(ROOT.parent)) if args.table.is_relative_to(ROOT.parent) else str(args.table)},
               'hypothesis': 'reduced Khovanov homology over F_2 of rank equal to the determinant, '
                             'so the double branched cover is an L-space',
               'skipped': skipped,
               'target_order': order, 'control_order': corder,
               'targets': targets, 'controls': controls}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, sort_keys=True, indent=1) + '\n')

    from collections import Counter
    print(f'wrote {args.out} ({args.out.stat().st_size // 1024} KiB)')
    print(f'  targets  {len(targets):5d}  ' + str(dict(sorted(Counter(t["test"] for t in targets.values()).items()))))
    print(f'  controls {len(controls):5d}  ' + str(dict(sorted(Counter(c["test"] for c in controls.values()).items()))))
    print('  skipped   ' + str(skipped))
    print('  sha256   ' + hashlib.sha256(args.out.read_bytes()).hexdigest())


if __name__ == '__main__':
    main()
