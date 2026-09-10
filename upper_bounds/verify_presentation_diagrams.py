#!/usr/bin/env python
"""Verify the eight marked diagrams in Section 4 and Appendix G of manuscript v1.6.

These are the same eight deposited constructions retained from earlier versions.

Each required file must identify the prescribed knot, have two distinct valid
marked rows, and pass both the exterior-isometry and unknot-recognition checks.
The fixed list prevents an empty or incomplete results directory from passing.
SnapPy's isometry computation is a numerical identification check; the unknot
endpoint is checked by simplification and exact knot Floer genus detection.

An optional --legacy check examines an archived neighboring-knot diagram only.
It is separate from the paper and supplies no unknotting-number conclusion.
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


def isometric(pd_a, pd_b, tries=12) -> bool:
    import snappy  # noqa: F401 (enables Link.exterior)
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
    return False  # failure to identify is not proof of different knot types


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=RESULTS)
    parser.add_argument('--legacy', action='store_true', help='also check a separate archived diagram, without a u claim')
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

    diagrams = load_paper_diagrams(args.results)
    for name, pd, rows in diagrams:
        check(f'{name}: {len(pd)}-crossing diagram identifies as {name} (SnapPy isometry)',
              isometric(pd, knotinfo.pd_code(name)))
        changed = pd
        for row in rows:
            changed = pd_crossing_change(changed, row)
        check(f'{name}: changing rows {list(rows)} gives the unknot (HFK genus 0)',
              inv.is_unknot(Link([list(q) for q in changed])))
    print(f'Paper constructions: {len(diagrams)}/8 checked.', flush=True)

    if args.legacy:
        print('Separate archived diagram check; no unknotting-number inference.', flush=True)
        with (args.results / 'reduced_13a_647.json').open() as handle:
            best = json.load(handle)['best']
        pd = best['pd_knot']
        row, = validate_marks(best['change_rows_0based'], len(pd), 1)
        check('archived source diagram identification', isometric(pd, knotinfo.pd_code('13a_647')))
        check('archived changed diagram identification',
              isometric(pd_crossing_change(pd, row), knotinfo.pd_code('13a_650')))
    print('ALL CHECKS PASSED' if ok_all else 'SOME CHECKS FAILED')
    if not ok_all:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
