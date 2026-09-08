#!/usr/bin/env python
"""Replay the two published upper bounds used in the Bernhard--Jablan result.

Brittenham--Hermiller already proved u(12n491),u(13n3370)<=2. These are
literature inputs, not new upper bounds of this paper. We replay their DT
crossing changes and independently check that the partner has u=1. Source
and partner identifications use numerical SnapPy exterior isometry, as in
the other upper-bound verifiers; they are not interval-certified isometries.
The partner's unknotting crossing is checked by knot Floer genus detection.

Run from any directory with the topology dependencies installed. Nothing is
written unless --out is supplied.
"""

import argparse
import ast
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'tait_graphs'))


def verify():
    import snappy  # noqa: F401: registers Spherogram's exterior method
    from spherogram import Link
    import database_knotinfo as dk
    import tait  # noqa: F401: installs the local KnotInfo adapter
    from tait import invariants
    from xtait.graph import from_pd, to_pd
    from xtait.moves import crossing_change
    from flipdet import exact_det
    from verify_braid_certificate import isometric, require

    source = ROOT / 'results/bernhard_jablan/brittenham_hermiller_upper_bounds.json'
    data = json.loads(source.read_text())
    rows = {r['name']: r for r in dk.link_list()}
    results = {}
    for name, record in data['records'].items():
        dt = record['dt_code']
        require(sorted(abs(x) for x in dt) == list(range(2, 2 * len(dt) + 1, 2)),
                'DT entries must have distinct even absolute values')
        index = record['change_entry_0based']
        require(type(index) is int and 0 <= index < len(dt), 'invalid crossing index')
        changed = dt.copy()
        changed[index] = -changed[index]
        link, child = Link('DT: ' + repr(dt)), Link('DT: ' + repr(changed))
        reference = Link(ast.literal_eval(rows[name]['pd_notation']))
        child_name = record['child']
        child_pd = ast.literal_eval(rows[child_name]['pd_notation'])
        require(isometric(link, reference), f'{name}: source diagram identification failed')
        require(isometric(child, Link(child_pd)), f'{name}: partner identification failed')
        graph = from_pd(child_pd)
        require(abs(exact_det(graph)) > 1, f'{child_name}: nontriviality not established')
        unknotting_row = None
        for row, edge in enumerate(sorted(graph.edges)):
            neighbor = crossing_change(graph, edge)
            if (abs(exact_det(neighbor)) == 1
                    and invariants.is_unknot(Link(to_pd(neighbor)))):
                unknotting_row = row
                break
        require(unknotting_row is not None, f'{child_name}: no verified unknotting crossing')
        require(record['child_upper_bound'] == 1 and record['upper_bound'] == 2,
                'published upper bound differs from the verified crossing count')
        results[name] = {'upper_bound': 2, 'child': child_name,
                         'published_dt_change_entry_0based': index,
                         'child_reference_unknotting_row_0based': unknotting_row,
                         'identification': 'SnapPy exterior isometry (numerical)',
                         'unknot_verification': 'knot Floer genus zero',
                         'attribution': 'Brittenham--Hermiller, arXiv:1705.05985v2'}
        print(f'{name}: published crossing change to {child_name}; '
              f'u({child_name})=1 verified; u({name})<=2', flush=True)
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, help='optional JSON report')
    args = parser.parse_args()
    report = verify()
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + '\n')
