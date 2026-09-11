#!/usr/bin/env python
"""Verify the nineteen presentation diagrams of manuscript v1.8.

The v1.3 profile retains the eight earlier diagrams. The v1.8 profile adds
eleven: ten directly give the unknot, and 13n447 goes to 10_91 with u=1.
Named knots require a returned isometry extending over meridians. A failed
isometry search is inconclusive and never used to exclude a knot type.
The fixed required lists prevent an empty or incomplete input from passing.

"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'tait_graphs'))
RESULTS = HERE.parent / 'results'
PAPER_KNOTS = ('13n_30', '13n_45', '13n_80', '13n_2379',
               '13n_2809', '13n_2907', '13n_3033', '13n_3589')
ADDITIONAL_KNOTS = ('13n_221', '13n_436', '13n_447', '13n_636', '13n_1439',
                    '13n_2251', '13n_2639', '13n_3108', '13n_3733', '13n_4025', '13n_4237')
from certificate_checks import validate_marks


def load_paper_diagrams(results=RESULTS):
    """Read exactly the eight required constructions, checking their structural claims."""
    diagrams = []
    for name in PAPER_KNOTS:
        with (Path(results) / 'best' / f'{name}.json').open() as handle:
            data = json.load(handle)
        best = data.get('best', data)
        if data.get('knot', name) != name:
            raise ValueError(f'{name}: file names a different knot')
        pd = best['pd_knot']
        rows = validate_marks(best['change_rows_0based'], len(pd), 2)
        if best.get('crossings', len(pd)) != len(pd):
            raise ValueError(f'{name}: recorded diagram size disagrees with PD')
        diagrams.append((name, pd, rows))
    return diagrams


def load_additional_diagrams(results=RESULTS):
    """Authenticate the new frozen inputs and validate the one-partner exception."""
    sys.path.insert(0, str(HERE.parent))
    from v18_results import inputs
    inputs(results)
    data = json.loads((Path(results)/'extensions_2026-09-10/presentation.json').read_text())
    bounds = json.loads((Path(results)/'paper_v1_1_snapshot.json').read_text())
    if set(data) != set(ADDITIONAL_KNOTS):
        raise ValueError('Missing or unexpected v1.8 presentation diagram')
    diagrams = []
    for name in ADDITIONAL_KNOTS:
        row = data[name]
        best = row['best']
        endpoint = '10_91' if name == '13n_447' else 'unknot'
        count = 1 if endpoint == '10_91' else 2
        if row['knot'] != name or row['endpoint'] != endpoint or row['claim_u_le'] != 2:
            raise ValueError(f'{name}: inconsistent construction claim')
        pd = best['pd_knot']
        if len(pd) != best['crossings']:
            raise ValueError(f'{name}: inconsistent PD size')
        marks = validate_marks(best['change_rows_0based'], len(pd), count)
        if endpoint != 'unknot' and bounds[endpoint][1] != 1:
            raise ValueError('Partner lacks an independent upper bound one')
        diagrams.append((name, pd, marks, endpoint))
    return diagrams


def isometric(pd_a, pd_b, tries=12) -> bool:
    import snappy  # noqa: F401 (enables Link.exterior)
    from spherogram import Link
    Ea = Link([list(q) for q in pd_a]).exterior()
    Eb = Link([list(q) for q in pd_b]).exterior()
    for _ in range(tries):
        try:
            maps = Ea.is_isometric_to(Eb, return_isometries=True)
            if any(f.extends_to_link() for f in maps):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False  # failure to identify is not proof of different knot types


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=RESULTS)
    parser.add_argument('--manuscript', choices=('v1.3', 'v1.8'), default='v1.8')
    args = parser.parse_args(argv)
    import tait  # noqa: F401 (sqlite compatibility)
    from tait import knotinfo, invariants as inv
    from spherogram import Link
    from xtait.graph import pd_crossing_change
    ok_all = True

    def check(label, ok):
        nonlocal ok_all
        ok_all &= bool(ok)
        print(('  ok  ' if ok else '  FAIL') + ' ' + label, flush=True)

    diagrams = [(n, pd, marks, 'unknot') for n, pd, marks in load_paper_diagrams(args.results)]
    if args.manuscript == 'v1.8':
        diagrams += load_additional_diagrams(args.results)
    for name, pd, rows, endpoint in diagrams:
        check(f'{name}: {len(pd)}-crossing diagram identifies as {name} (SnapPy isometry)',
              isometric(pd, knotinfo.pd_code(name)))
        changed = pd
        for row in rows:
            changed = pd_crossing_change(changed, row)
        if endpoint == 'unknot':
            check(f'{name}: changing rows {list(rows)} gives the unknot (HFK genus 0)',
                  inv.is_unknot(Link([list(q) for q in changed])))
        else:
            check(f'{name}: changing rows {list(rows)} gives {endpoint}, with known u=1',
                  isometric(changed, knotinfo.pd_code(endpoint)))
    expected = 19 if args.manuscript == 'v1.8' else 8
    check(f'Paper constructions: {len(diagrams)}/{expected} checked.', len(diagrams) == expected)

    print('ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED')
    if not ok_all:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
