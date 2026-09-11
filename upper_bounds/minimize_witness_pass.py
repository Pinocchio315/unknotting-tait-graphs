#!/usr/bin/env python
"""Minimize marked diagrams using Reidemeister moves, flypes, and safe pass moves.

Same best-first search and the same endpoint re-verification (input and result) as
code/upper_bounds/minimize_witness.py; the only change is that children() is replaced by
children_pass.children_with_pass.  The input may be a campaign certificate (witness_pd,
change_rows) or an earlier minimisation record (best.pd_knot, best.change_rows_0based), so a
stalled result can be continued.  The search stops at the crossing number of the knot.

    python upper_bounds/minimize_witness_pass.py \
        results/extensions_2026-09-10/min11/13n_436.json \
        --endpoint unknot --time-limit 1800 --out /tmp/13n_436_reduced.json
"""
from __future__ import annotations

import argparse
import heapq
import itertools
import json
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

import tait  # noqa: F401,E402
from tait import graph as tg  # noqa: E402
from reduce_witnesses import marked_key  # noqa: E402
from minimize_witness import endpoint_ok  # noqa: E402
from children_pass import children_with_pass, changed_det, STATS  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate')
    ap.add_argument('--endpoint', required=True)
    ap.add_argument('--time-limit', type=float, default=1800.0)
    ap.add_argument('--max-extra', type=int, default=2)
    ap.add_argument('--max-states', type=int, default=400000)
    ap.add_argument('--pass-extra', type=int, default=0, help='allow pass routes up to |run| + this')
    ap.add_argument('--pass-routes', type=int, default=4, help='routes tried per run')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    rec = json.load(open(args.certificate))
    if 'best' in rec:
        name, pd, rows = rec['knot'], rec['best']['pd_knot'], rec['best']['change_rows_0based']
    else:
        name, pd, rows = rec['name'], rec['witness_pd'], rec['change_rows']
    floor = int(re.match(r'(\d+)', name).group(1))
    g0 = tg.from_pd([list(q) for q in pd])
    ids = sorted(g0.edges)
    S0 = {ids[i] for i in rows}
    det_target = changed_det(g0, S0)
    print(f'{name}: start {g0.n_crossings()} crossings, marks {sorted(S0)}, endpoint {args.endpoint}, '
          f'floor {floor}, det after change {det_target}', flush=True)
    assert endpoint_ok(g0, S0, args.endpoint), 'endpoint check fails on the input witness!'

    t0 = time.time()
    seen = {marked_key(g0, S0)}
    cnt = itertools.count()
    heap = [(g0.n_crossings(), 0, next(cnt), g0, S0, [])]
    best = (g0.n_crossings(), g0, S0, [])
    expanded = 0
    cap = g0.n_crossings() + args.max_extra
    while heap and time.time() - t0 < args.time_limit and expanded < args.max_states:
        n, depth, _, g, S, path = heapq.heappop(heap)
        expanded += 1
        if n < best[0]:
            best = (n, g, S, path)
            cap = n + args.max_extra
            print(f'  {n} crossings after {expanded} states, {time.time()-t0:.0f}s, '
                  f'{sum(1 for d in path if d[0] == "pass")} pass moves on the path', flush=True)
            if n <= floor:
                break
        for h, S2, desc in children_with_pass(g, S, det_target, args.pass_extra, args.pass_routes):
            m = h.n_crossings()
            if m > cap:
                continue
            key = marked_key(h, S2)
            if key in seen:
                continue
            seen.add(key)
            heapq.heappush(heap, (m, depth + 1, next(cnt), h, S2, path + [desc]))

    n, g, S, path = best
    print(f'best: {n} crossings ({expanded} states expanded, {time.time()-t0:.0f}s)')
    print('pass statistics:', STATS)
    assert endpoint_ok(g, S, args.endpoint), 'endpoint check fails on the reduced witness!'
    print('endpoint re-verified on the reduced witness')
    out_rec = {'knot': name, 'endpoint': args.endpoint, 'start_crossings': g0.n_crossings(),
               'best': {'pd_knot': tg.to_pd(g), 'change_rows_0based': [sorted(g.edges).index(e) for e in sorted(S)],
                        'crossings': n},
               'moves': [[d[0], d[1] if not isinstance(d[1], tuple) else list(d[1])] for d in path],
               'pass_moves_on_path': sum(1 for d in path if d[0] == 'pass'),
               'pass_statistics': dict(STATS)}
    out = args.out or f'reduced_pass_{name}.json'
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, 'w') as h:
        json.dump(out_rec, h, indent=1, default=str)
    print('wrote', out)


if __name__ == '__main__':
    main()
