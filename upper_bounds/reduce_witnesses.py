#!/usr/bin/env python
"""Reduce the size of verified witness diagrams while keeping the marked (to-be-changed) crossings.

Input: results/best/<knot>.json — an independently verified compact diagram D of the knot
with a set S of marked crossings such that changing S gives the unknot (u <= |S|).  These
files are presentation certificates; the forward-search provenance is stored separately in
results/witness_chains/.  We search over diagrams
with marks using only moves that commute with the crossing changes at S:
  * R1-/R2- removals that do not involve a marked crossing,
  * flypes (crossing ids and signs preserved),
  * R3 moves; a marked crossing may take part only if the triangle/star is a valid R3 configuration both in
    D and in D with S changed (then the move is the same combinatorial Y-Δ on both, and the marks are mapped
    to the corresponding new edges),
  * temporary R2+ additions (new crossings unmarked) up to --max-extra above the current best,
best-first on the crossing number (marked-aware canonical code for de-duplication).  The crossing number of
the knot is a hard floor.  Every improved diagram is re-verified (HFK of the changed diagram = unknot, SnapPy
isometry with the census knot) and written to --out (the previous results are left untouched).

    python upper_bounds/reduce_witnesses.py --knots 13n_30 13n_45 13n_80 --time-limit 180 --max-extra 2
"""
from __future__ import annotations

import argparse
import heapq
import itertools
import json
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs'))


def marked_key(g, marks):
    """Canonical code of the diagram with the marked crossings distinguished (sign *2)."""
    from tait.canonical import canonical_code
    h = g.copy()
    for e in marks:
        h.signs[e] = 2 * h.signs[e]
    # colour_invariant=False: no dual (dual() validates signs); all states of one search share the colour class
    return canonical_code(h, colour_invariant=False)


def children(g, marks):
    """Yield (h, marks', description) for all allowed moves."""
    from tait import moves as mv
    from tait.graph import ribbon_faces
    marks = set(marks)
    # R1- on unmarked kinks
    for kind, e in mv.find_r1_removals(g):
        if e not in marks:
            yield mv.r1_remove(g, e), marks, ('r1-', e)
    # R2- not touching marks
    for kind, prm in mv.find_r2_removals(g):
        if kind == 'parallel':
            if not (set(prm) & marks):
                yield mv.r2_remove_parallel(g, *prm), marks, ('r2-p', prm)
        else:
            w = prm[0]
            if not (set(g.rotation[w]) & marks):
                yield mv.r2_remove_series(g, w), marks, ('r2-s', w)
    # flypes (ids preserved)
    for e, pi in mv.find_flypes(g):
        try:
            yield mv.flype(g, e, pi), marks, ('flype', (e, pi))
        except ValueError:
            continue
    # R3 Δ->Y on triangle faces
    faces = ribbon_faces(g)
    for fi in mv.find_r3_triangles(g):
        f = faces[fi]
        tri = [d[1] for d in f]
        verts = [d[0] for d in f]
        if marks & set(tri):
            # validity in the changed diagram: signs (with marks flipped) not all equal
            fl = {e: (-g.signs[e] if e in marks else g.signs[e]) for e in tri}
            if len(set(fl.values())) != 2:
                continue
        try:
            h = mv.r3_delta_to_y(g, fi)
        except ValueError:
            continue
        # new edge base+i at vertex verts[i] corresponds to the triangle edge opposite verts[i]
        new_ids = sorted(set(h.edges) - set(g.edges))
        base = min(new_ids)
        m2 = set(m for m in marks if m not in tri)
        for i, V in enumerate(verts):
            opp = [te for te in tri if V not in g.edges[te]][0]
            if opp in marks:
                m2.add(base + i)
        yield h, m2, ('r3Y', fi)
    # R3 Y->Δ at degree-3 stars
    for t in mv.find_r3_stars(g):
        es = list(g.rotation[t])
        if marks & set(es):
            fl = {e: (-g.signs[e] if e in marks else g.signs[e]) for e in es}
            if len(set(fl.values())) != 2:
                continue
        try:
            h = mv.r3_y_to_delta(g, t)
        except ValueError:
            continue
        new_ids = sorted(set(h.edges) - set(g.edges))
        base = min(new_ids)
        m2 = set(m for m in marks if m not in es)
        for k, e in enumerate(es):
            if e in marks:
                m2.add(base + k)          # new edge base+k has sign -sign(e_k): it is the moved crossing e_k
        yield h, m2, ('r3Δ', t)
    # temporary R2+ (new crossings unmarked)
    for (fi, i, j) in mv.find_r2_parallel_additions(g):
        for s in (1, -1):
            try:
                yield mv.r2_add_parallel(g, fi, i, j, s), marks, ('r2+p', (fi, i, j, s))
            except ValueError:
                continue
    for (v, i, j) in mv.find_r2_series_additions(g):
        for s in (1, -1):
            try:
                yield mv.r2_add_series(g, v, i, j, s), marks, ('r2+s', (v, i, j, s))
            except ValueError:
                continue


def verify(g, marks, knot):
    import snappy
    from tait import graph as tg, moves as mv, invariants as inv
    U = g
    for e in marks:
        U = mv.crossing_change(U, e)
    ok_unknot = bool(inv.is_unknot(U))
    try:
        iso = bool(tg.to_link(g).exterior().is_isometric_to(snappy.Manifold('K' + knot.replace('_', ''))))
    except Exception:
        iso = False
    return ok_unknot, iso


def reduce_one(knot, rec, time_limit, max_extra, max_states, floor):
    from tait import graph as tg, moves as mv
    b = rec['best']
    g0 = tg.from_dict(b['tait_knot'])
    S0 = set(b['marked_crossing_ids'])
    n0 = g0.n_crossings()
    t0 = time.time()
    seen = {marked_key(g0, S0)}
    cnt = itertools.count()
    heap = [(n0, 0, next(cnt), g0, S0, [])]
    best = (n0, g0, S0, [])
    expanded = 0
    cap = n0 + max_extra
    while heap and time.time() - t0 < time_limit and expanded < max_states:
        n, depth, _, g, S, path = heapq.heappop(heap)
        expanded += 1
        if n < best[0]:
            best = (n, g, S, path)
            cap = n + max_extra
            if n <= floor:
                break
        for h, S2, desc in children(g, S):
            m = h.n_crossings()
            if m > cap:
                continue
            key = marked_key(h, S2)
            if key in seen:
                continue
            seen.add(key)
            heapq.heappush(heap, (m, depth + 1, next(cnt), h, S2, path + [desc]))
    n1, g1, S1, path = best
    out = {'knot': knot, 'k': len(S1), 'crossings_before': n0, 'crossings_after': n1, 'crossing_number_floor': floor,
           'expanded': expanded, 'distinct_states': len(seen), 'seconds': round(time.time() - t0, 1),
           'path': [[d[0], list(d[1]) if isinstance(d[1], (tuple, list)) else d[1]] for d in path]}
    if n1 < n0:
        ok_unknot, iso = verify(g1, S1, knot)
        order = sorted(g1.edges)
        out.update({'verified': bool(ok_unknot and iso), 'unknot_after_change': ok_unknot, 'isometric_to_census': iso,
                    'pd_knot': tg.to_oriented_pd(g1), 'change_rows_0based': [order.index(e) for e in S1],
                    'tait_knot': tg.to_dict(g1), 'marked_crossing_ids': sorted(S1)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--knots', nargs='*')
    ap.add_argument('--src', default='../results/best')
    ap.add_argument('--out', default='../results/reduced')
    ap.add_argument('--time-limit', type=float, default=180.0)
    ap.add_argument('--max-extra', type=int, default=2)
    ap.add_argument('--max-states', type=int, default=200000)
    args = ap.parse_args()
    import tait  # noqa: F401
    from tait import knotinfo
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tait_graphs')
    src = os.path.join(root, args.src)
    out = os.path.join(root, args.out)
    os.makedirs(out, exist_ok=True)
    knots = args.knots or sorted(f[:-5] for f in os.listdir(src) if f.endswith('.json'))
    lines = []
    for knot in knots:
        rec = json.load(open(os.path.join(src, f'{knot}.json')))
        floor = int(knotinfo.row(knot)['crossing_number'])
        n0 = rec['best']['crossings_after']
        if n0 <= floor:
            msg = f'{knot}: witness already has {n0} crossings = crossing number; nothing to reduce'
            print(msg, flush=True)
            lines.append(msg)
            continue
        res = reduce_one(knot, rec, args.time_limit, args.max_extra, args.max_states, floor)
        with open(os.path.join(out, f'{knot}.json'), 'w') as h:
            json.dump(res, h, indent=1)
        msg = (f"{knot}: {res['crossings_before']} -> {res['crossings_after']} crossings (floor {floor}); "
               f"expanded {res['expanded']}, states {res['distinct_states']}, {res['seconds']}s"
               + (f"; verified={res['verified']} rows={res['change_rows_0based']} path={len(res['path'])} moves" if 'verified' in res else ''))
        print(msg, flush=True)
        lines.append(msg)
    with open(os.path.join(out, 'summary_reduced.md'), 'a') as h:
        h.write(f'\n## run {time.strftime("%Y-%m-%d %H:%M")} (time limit {args.time_limit}s, max extra {args.max_extra})\n\n')
        h.write('\n'.join('- ' + ln for ln in lines) + '\n')


if __name__ == '__main__':
    main()
