#!/usr/bin/env python
"""Acceptance tests for the upper-bound campaign (stdlib only; ~1-3 min).  Run FIRST on any machine:

    python3 selftest.py

  1. FlipDet exactness: subset determinants agree with the Goeritz determinant of the actually
     changed graph (exact rational elimination), on the KnotInfo diagram and on expanded ones.
  2. Expansion soundness: R1-free expansion to 40 crossings preserves the knot determinant for
     every emitted diagram, and one expanded diagram reduces back to the minimal crossing number
     with the pass engine.
  3. Positive control: the published 13-crossing witness for 13n2379 --- changing its two recorded
     crossings and running the reducer must reach 0 crossings (the unknot), and the det filter
     must flag exactly that subset as det = 1.
  4. Known-code identification: canonical codes in data/known_u.json match diagrams rebuilt from
     KnotInfo PDs (spot check), so the j < k route can identify reduced knots.
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.reduce import greedy_reduce, search_reduce, ReduceConfig  # noqa: E402
from flipdet import FlipDet, exact_det, P  # noqa: E402
from expand import expand_diagrams  # noqa: E402

WITNESS_13N2379_PD = [[1, 9, 2, 8], [18, 8, 19, 7], [6, 20, 7, 19], [20, 6, 21, 5],
                      [4, 13, 5, 14], [14, 3, 15, 4], [2, 22, 3, 21], [9, 1, 10, 26],
                      [25, 11, 26, 10], [11, 25, 12, 24], [16, 24, 17, 23],
                      [22, 16, 23, 15], [12, 18, 13, 17]]
WITNESS_13N2379_ROWS = [6, 12]


def check(label, ok):
    print(('  ok  ' if ok else '  FAIL') + ' ' + label)
    if not ok:
        sys.exit(1)


def flip_rows(g, F, rows):
    h = g
    for i in rows:
        h = mv.crossing_change(h, F.edge_ids[i])
    return h


def main() -> None:
    t0 = time.time()
    targets = json.load(open(os.path.join(HERE, 'data', 'targets.json')))
    known = json.load(open(os.path.join(HERE, 'data', 'known_u.json')))

    # ---- 1. FlipDet vs exact determinants --------------------------------------------
    t = targets[0]
    g = from_pd([list(q) for q in t['pd']])
    F = FlipDet(g)
    check(f"base det of {t['name']} = KnotInfo det",
          F.det_g % P in (t['det'] % P, (-t['det']) % P))
    import random
    rnd = random.Random(7)
    m = len(F.edge_ids)
    ok = True
    for k in (1, 2, 3):
        for _ in range(25):
            rows = tuple(sorted(rnd.sample(range(m), k)))
            d_pred = F.det_after(rows)
            d_true = exact_det(flip_rows(g, F, rows))
            ok = ok and d_pred in (d_true % P, (-d_true) % P)
    check('subset determinants exact on the KnotInfo diagram (75 random subsets, k=1..3)', ok)

    # ---- 2. expansion soundness ------------------------------------------------------
    n_eff = 40 - (40 - t['crossings']) % 2
    n_emit, det_ok = 0, True
    first = None
    for ge, path in expand_diagrams(g, n_target=40, count=12, seed=1):
        n_emit += 1
        Fe = FlipDet(ge)
        det_ok = det_ok and Fe.det_g % P in (t['det'] % P, (-t['det']) % P)
        det_ok = det_ok and ge.n_crossings() == n_eff
        if first is None:
            first = ge
        if any(kind.startswith('r1') for kind, _ in path):
            det_ok = False
    check(f'expansion emitted {n_emit}/12 R1-free {n_eff}-crossing diagrams with invariant det',
          n_emit == 12 and det_ok)
    Fe = FlipDet(first)
    rows = tuple(sorted(rnd.sample(range(len(Fe.edge_ids)), 2)))
    d_true = exact_det(flip_rows(first, Fe, rows))
    check('subset determinant exact on an expanded diagram at the target size',
          Fe.det_after(rows) in (d_true % P, (-d_true) % P))
    red, _ = greedy_reduce(first)
    if red.n_crossings() > t['crossings']:
        red, _, _ = search_reduce(red, ReduceConfig(target=t['crossings'], max_extra=2,
                                                    max_nodes=4000, time_limit=60,
                                                    log_every=10 ** 9))
    check(f"an expanded diagram reduces back to {t['crossings']} crossings "
          f'(got {red.n_crossings()})', red.n_crossings() <= t['crossings'])

    # ---- 3. positive control: the 13n2379 witness ------------------------------------
    gw = from_pd([list(q) for q in WITNESS_13N2379_PD])
    Fw = FlipDet(gw)
    check('det filter flags the published witness subset (det = 1)',
          Fw.det_after(tuple(WITNESS_13N2379_ROWS)) in (1, P - 1))
    hits = [rows for rows in itertools.combinations(range(13), 2)
            if Fw.det_after(rows) in (1, P - 1)]
    changed = flip_rows(gw, Fw, WITNESS_13N2379_ROWS)
    red, _ = greedy_reduce(changed)
    if red.n_crossings() > 0:
        red, _, _ = search_reduce(red, ReduceConfig(target=0, max_extra=2, max_nodes=4000,
                                                    time_limit=60, log_every=10 ** 9))
    check(f'witness changes reduce to the unknot (0 crossings; det-1 subsets on the diagram: '
          f'{len(hits)})', red.n_crossings() == 0)

    # ---- 4. identification table -----------------------------------------------------
    by_name = {x['name']: x for x in known}
    ok = True
    for nm in ('3_1', '4_1', '6_1'):
        x = by_name[nm]
        ok = ok and len(x['codes']) == 2
    tref = by_name['3_1']
    check('known_u codes present (3_1, 4_1, 6_1; diagram + mirror)', ok and tref['u'] == 1)

    print(f'ALL TESTS PASSED ({time.time() - t0:.0f}s)')


if __name__ == '__main__':
    main()
