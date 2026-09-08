#!/usr/bin/env python3
"""Bounded geometric spot checks of newly decisive crossing neighbours.

Eight distinct children are chosen deterministically from the first newly
closed parent records. SnapPy matching is numerical and is reported as such;
the exact determinant and signature are also recomputed independently of
the stored crossing-change rows. No search failure establishes a lower bound.
"""
import ast
import gzip
import json
from pathlib import Path
import signal
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'upper_bounds'))
sys.path.insert(0, str(ROOT / 'crossing_changes'))
import database_knotinfo as dk
import snappy
from spherogram import Link
from xtait.graph import from_pd, to_pd
from xtait.moves import crossing_change
from xtait.reduce import greedy_reduce
from flipdet import exact_det
from sig import signature

def main():
    def timeout(*_):
        raise TimeoutError('55-second bounded verification expired')
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(55)
    rows = {(r['knot'], r['crossing']): r for r in map(json.loads,
        gzip.open(ROOT / 'results/open23/rows833.jsonl.gz', 'rt'))}
    audit = json.loads((HERE / 'updated_dichotomy_audit.json').read_text())
    table = {r['name']: r for r in dk.link_list()}
    selected, seen = [], set()
    for parent in audit['new_parent_dependencies']:
        for change in parent['new_exclusions']:
            row = rows[(parent['parent'], change['crossing'])]
            if row['result_how'] == 'isometry-jones' and row['result'] not in seen:
                selected.append(row)
                seen.add(row['result'])
    results = []
    for row in selected[:8]:
        pd = ast.literal_eval(table[row['knot']]['pd_notation'])
        graph = from_pd(pd)
        changed = crossing_change(graph, sorted(graph.edges)[row['crossing']])
        changed_pd = to_pd(changed)
        determinant = exact_det(changed)
        sigma = signature(Link(changed_pd).PD_code())
        assert determinant == row['det_changed']
        assert abs(sigma) == abs(row['sigma_c'])
        reduced, moves = greedy_reduce(changed)
        child_pd = ast.literal_eval(table[row['result']]['pd_notation'])
        matched = Link(to_pd(reduced)).exterior().is_isometric_to(Link(child_pd).exterior())
        assert matched, (row['knot'], row['crossing'], row['result'])
        results.append({'parent': row['knot'], 'crossing': row['crossing'],
                        'child': row['result'], 'determinant': determinant,
                        'signature_abs': abs(sigma), 'isometry_numerical': matched,
                        'source_pd': pd, 'changed_pd': changed_pd,
                        'reduction_moves': moves})
    signal.alarm(0)
    (HERE / 'new_dichotomy_controls.json').write_text(json.dumps(results, indent=2) + '\n')
    print('All eight new crossing-neighbour identifications and invariant checks passed.')

if __name__ == '__main__':
    main()
