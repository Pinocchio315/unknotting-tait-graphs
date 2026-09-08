#!/usr/bin/env python3
"""Replay every crossing change in a named minimal diagram of 14a2539.

The source diagram comes independently from Spherogram's named knot table.
Identifications are reported as numerical SnapPy matches unless exact Regina
diagram signatures match a named reference. No polynomial-only identification
is accepted. Flype completeness is the separate mathematical justification
that these fourteen crossings determine the full Bernhard--Jablan set.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'upper_bounds'))
import snappy  # registers Link.exterior
import regina
from spherogram import Link
from xtait.graph import from_pd
from flipdet import exact_det

def signature(link):
    return regina.Link.fromPD(link.PD_code(min_strand_index=1)).sig()

def main():
    source = Link('14a2539')
    assert source.is_alternating() and len(source.crossings) == 14
    pd = source.PD_code()
    names = ['10_11', '12a178', '12a183', '12a196', '12a684', '8_6', '4_1']
    reference = {name: signature(Link(name)) for name in names}
    report = {'source_pd': pd, 'rows': [], 'reference_signatures': reference}
    for index in range(14):
        changed = [list(q) for q in pd]
        changed[index] = changed[index][1:] + changed[index][:1]
        child = Link(changed)
        child.simplify('basic')
        factors = []
        for factor in child.deconnect_sum():
            sig = signature(factor)
            exact = [name for name, ref in reference.items() if sig == ref]
            numerical = list(map(str, factor.exterior().identify()))
            factors.append({'pd': factor.PD_code(), 'regina_signature': sig,
                            'exact_named_diagram_matches': exact,
                            'snappy_identifications_numerical': numerical})
        report['rows'].append({'crossing_zero_based': index,
                               'determinant': exact_det(from_pd(changed)),
                               'factors': factors})
        print(index, report['rows'][-1]['determinant'],
              [r['snappy_identifications_numerical'] for r in factors], flush=True)
    Path(__file__).with_name('bjset14a2539_replay.json').write_text(
        json.dumps(report, indent=2) + '\n')

if __name__ == '__main__':
    main()
