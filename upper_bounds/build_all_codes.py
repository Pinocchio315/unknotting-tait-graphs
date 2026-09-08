#!/usr/bin/env python
"""Generate data/all_knot_codes.json: canonical codes (diagram + mirror) and unknotting
intervals for EVERY KnotInfo knot <= 13 crossings --- the full identification table used by
resolve_unidentified.py (the campaign's in-run table only covered knots with u <= 2).

    ~/.pyenv/versions/unknot-venv/bin/python build_all_codes.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from make_targets import knotinfo_rows, interval, pd_code, mirror_pd_graph  # noqa: E402


def main() -> None:
    rows = knotinfo_rows()
    out = []
    for nm, r in rows.items():
        if not r['crossing_number'].strip().isdigit():   # header/metadata rows
            continue
        try:
            lo, hi = interval(r)
        except ValueError:
            lo = hi = None
        try:
            pd = pd_code(r)
            g = from_pd([list(q) for q in pd])
            codes = [repr(canonical_code(g)), repr(canonical_code(mirror_pd_graph(pd)))]
        except Exception as exc:
            print('  !', nm, exc)
            continue
        out.append({'name': nm, 'lo': lo, 'hi': hi, 'det': abs(int(r['determinant'])),
                    'crossings': int(r['crossing_number']), 'codes': codes})
    with open(os.path.join(HERE, 'data', 'all_knot_codes.json'), 'w') as h:
        json.dump(out, h)
    print(f'{len(out)} knots -> data/all_knot_codes.json')


if __name__ == '__main__':
    main()
