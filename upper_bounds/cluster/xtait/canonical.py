"""Canonical code of a signed planar map — de-duplication key for diagram searches (stdlib only).

code(g1) == code(g2)  <=>  same diagram on S^2 (same chirality), independent of vertex/edge labels,
rotation-list starting points and of which checkerboard colour class was stored.

Method: for every start dart, BFS-label the map (vertices in discovery order; each vertex's rotation read
from its entry dart; every dart records (opposite vertex label, slot offset at the opposite vertex, edge
sign)); take the lexicographic minimum over all starts, and over the dual (colour invariance).
"""
from __future__ import annotations

from .graph import EmbeddedTaitGraph, _rho_theta, dual


def _code_from(g, d0, theta, dart_lists):
    label = {d0[0]: 0}
    entry = {d0[0]: d0}
    order = [d0[0]]
    queue = [d0[0]]
    qi = 0
    code = []
    while qi < len(queue):
        v = queue[qi]
        qi += 1
        dl = dart_lists[v]
        start = dl.index(entry[v])
        n = len(dl)
        row = [n]
        for i in range(n):
            d = dl[(start + i) % n]
            op = theta[d]
            w = op[0]
            if w not in label:
                label[w] = len(order)
                order.append(w)
                entry[w] = op
                queue.append(w)
            dw = dart_lists[w]
            off = (dw.index(op) - dw.index(entry[w])) % len(dw)
            row.append((label[w], off, g.signs[d[1]]))
        code.append(tuple(row))
    return tuple(code)


def _map_code(g):
    rho, theta = _rho_theta(g)
    dart_lists = {v: g.dart_list(v) for v in g.vertices}
    best = None
    for d0 in theta:
        c = _code_from(g, d0, theta, dart_lists)
        if best is None or c < best:
            best = c
    return best


def canonical_code(g: EmbeddedTaitGraph, include_mirror: bool = False, colour_invariant: bool = True):
    if not g.edges:
        return ('unknot',)
    best = _map_code(g)
    if colour_invariant:
        best = min(best, _map_code(dual(g)))
    if include_mirror:
        best = min(best, canonical_code(g.mirror(), False, colour_invariant))
    return best


def canonical_hash(g: EmbeddedTaitGraph) -> int:
    """64-bit stable hash of the canonical code (visited-set key; low memory)."""
    import hashlib
    b = repr(canonical_code(g)).encode()
    return int.from_bytes(hashlib.blake2b(b, digest_size=8).digest(), 'big')
