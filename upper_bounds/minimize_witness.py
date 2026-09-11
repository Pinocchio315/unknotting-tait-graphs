#!/usr/bin/env python
"""Marked-crossing-preserving minimisation of a campaign witness (run LOCALLY, unknot-venv).

Reuses the validated marked-move engine of ../scripts/reduce_witnesses.py (R1-/R2- away from
marks, flypes, mark-respecting R3, temporary R2+), but the endpoint check is generalised:

  --endpoint unknot     changing the marks must give the unknot (HFK genus 0)     [j = k]
  --endpoint <partner>  changing the marks must give the named KnotInfo knot
                        (SnapPy isometry; u bound then comes from the partner)    [j < k]

Input: a certificate JSON (SUCCESS_/CHAIN_ schema: witness_pd + change_rows).

    python upper_bounds/minimize_witness.py \
        results/extensions_2026-09-10/runs190/SUCCESS_13n_221.json \
        --endpoint unknot --time-limit 600 --out /tmp/reduced_13n_221.json
"""
from __future__ import annotations

import argparse
import heapq
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

import tait  # noqa: F401,E402
from tait import graph as tg, moves as mv, invariants as inv  # noqa: E402
from reduce_witnesses import children, marked_key  # noqa: E402  (marked-move engine)


def endpoint_ok(g, marks, endpoint):
    from spherogram import Link
    h = g
    for e in marks:
        h = mv.crossing_change(h, e)
    L = tg.to_link(h)
    if endpoint == 'unknot':
        return inv.is_unknot(L)
    from tait import knotinfo
    Ea = L.exterior()
    Eb = Link(knotinfo.pd_code(endpoint)).exterior()
    for _ in range(10):
        try:
            if any(f.extends_to_link() for f in Ea.is_isometric_to(Eb, return_isometries=True)):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('certificate')
    ap.add_argument('--endpoint', required=True, help="'unknot' or a KnotInfo name like 13a_650")
    ap.add_argument('--time-limit', type=float, default=600.0)
    ap.add_argument('--max-extra', type=int, default=2)
    ap.add_argument('--max-states', type=int, default=400000)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    cert = json.load(open(args.certificate))
    # tait is a superset of the frozen xtait code: from_pd builds identical graphs with
    # identical edge ids, and to_pd rows correspond to sorted edge ids in both.
    g0 = tg.from_pd([list(q) for q in cert['witness_pd']])
    edge_ids = sorted(g0.edges)
    S0 = {edge_ids[i] for i in cert['change_rows']}
    print(f"{cert['name']}: witness {g0.n_crossings()} crossings, marks {sorted(S0)}, "
          f"endpoint {args.endpoint}")
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
            print(f'  {n} crossings after {expanded} states, {time.time()-t0:.0f}s', flush=True)
        for h, S2, desc in children(g, S):
            m = h.n_crossings()
            if m > cap:
                continue
            key = marked_key(h, S2)
            if key in seen:
                continue
            seen.add(key)
            heapq.heappush(heap, (m, depth + 1, next(cnt), h, S2, path + [desc]))

    n, g, S, path = best
    print(f'best: {n} crossings ({expanded} states expanded)')
    assert endpoint_ok(g, S, args.endpoint), 'endpoint check fails on the reduced witness!'
    print('endpoint re-verified on the reduced witness')
    pd = tg.to_pd(g)
    rows = [sorted(g.edges).index(e) for e in sorted(S)]
    rec = {'knot': cert['name'], 'endpoint': args.endpoint, 'claim_u_le': cert.get('u'),
           'best': {'pd_knot': pd, 'change_rows_0based': rows, 'crossings': n},
           'from_crossings': g0.n_crossings(), 'moves': len(path)}
    out = args.out or f"reduced_{cert['name']}.json"
    with open(out, 'w') as h:
        json.dump(rec, h, indent=1)
    print('wrote', out)


if __name__ == '__main__':
    main()
