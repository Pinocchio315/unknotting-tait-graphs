"""Serialisable move records for the enumeration search (stdlib only).

A record is a JSON-able pair [kind, params]:
    ["r3D", face_index]            Reidemeister III, triangle -> Y
    ["r3Y", vertex]                Reidemeister III, Y -> triangle
    ["flypeB", {"dual": bool, "e": edge, "edges": [tangle edges]}]
                                   flype of crossing e across the tangle with the given edges,
                                   found on the graph itself or on its dual checkerboard graph
    ["pass", {"run_index": i, "route": [[arc, side], ...], "gain": g}]
                                   pass move (xtait.passmove parameters relative to to_pd(g))
    ["r2P", [face_index, i, j, sign]]   Reidemeister II addition: parallel pair in the Tait graph
    ["r2S", [vertex, i, j, sign]]       Reidemeister II addition: series pair (vertex split)
Replaying the records from the KnotInfo PD reproduces every intermediate graph exactly (all
operations are deterministic), which is what the certificates rely on.
"""
from __future__ import annotations

from xtait import moves as mv
from xtait.graph import dual
from xtait.flypes import find_flypes_all, flype_block


def apply_record(g, rec):
    kind, params = rec
    if kind == 'r3D':
        return mv.r3_delta_to_y(g, int(params))
    if kind == 'r3Y':
        return mv.r3_y_to_delta(g, int(params))
    if kind == 'pass':
        return mv.apply_move(g, 'pass', (params,))
    if kind == 'r2P':
        fi, i, j, sign = (int(x) for x in params)
        return mv.r2_add_parallel(g, fi, i, j, sign)
    if kind == 'r2S':
        v, i, j, sign = (int(x) for x in params)
        return mv.r2_add_series(g, v, i, j, sign)
    if kind == 'flypeB':
        base = dual(g) if params['dual'] else g
        e, edges = int(params['e']), set(int(x) for x in params['edges'])
        for e2, P in find_flypes_all(base):
            if e2 == e and set(P['edges']) == edges:
                return flype_block(base, e2, P)
        raise ValueError('flypeB: recorded tangle not found')
    raise ValueError(f'unknown move record {kind}')


def replay(g0, records):
    g = g0
    for rec in records:
        g = apply_record(g, rec)
    return g
