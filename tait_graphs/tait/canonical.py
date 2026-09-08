"""Canonical code of a signed planar map (embedded Tait graph) — for search de-duplication.

Why: in a search over diagrams (R-moves, flypes, crossing changes) the same diagram is reached by many
different move sequences, under arbitrary vertex/edge labels.  A *canonical code* is a string/tuple that
depends only on the diagram (not on labels or on where the rotation lists start), so that
    code(g1) == code(g2)  <=>  g1 and g2 are the same signed map (same diagram on S^2, same chirality).
It is used as the key of a visited-set / transposition table, exactly like the canonical (cyclically
reduced, minimal-rotation) words in the AC-conjecture search.

Method: BFS labelling of a rooted map.  For every start dart d0 we traverse the map deterministically:
vertices are numbered in discovery order; at each vertex we walk its rotation starting from the dart we
entered through; for each dart we record (label of the opposite vertex, slot offset of the opposite dart
relative to that vertex's entry dart, sign of the edge).  That description determines the map.  The
canonical code is the lexicographically smallest description over all 2E start darts.  O(E^2), a few ms
for diagrams with <= 60 crossings.
"""
from __future__ import annotations

from .graph import EmbeddedTaitGraph, _rho_theta, Dart


def _code_from(g: EmbeddedTaitGraph, d0: Dart, theta, dart_lists) -> tuple:
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
            # offset of op within w's rotation relative to w's entry dart
            dw = dart_lists[w]
            off = (dw.index(op) - dw.index(entry[w])) % len(dw)
            row.append((label[w], off, g.signs[d[1]]))
        code.append(tuple(row))
    return tuple(code)


def _map_code(g: EmbeddedTaitGraph) -> tuple:
    rho, theta = _rho_theta(g)
    dart_lists = {v: g.dart_list(v) for v in g.vertices}
    best = None
    for d0 in theta:
        c = _code_from(g, d0, theta, dart_lists)
        if best is None or c < best:
            best = c
    return best


def canonical_code(g: EmbeddedTaitGraph, include_mirror: bool = False, colour_invariant: bool = True) -> tuple:
    """Canonical code (hashable tuple) of the diagram represented by g.

    colour_invariant=True (default): the code does not depend on which checkerboard colour class was
    chosen as the vertex set (min over g and its dual), so graphs produced by different routes
    (e.g. from_link on a spherogram diagram vs. our own moves) compare correctly.
    include_mirror=True additionally identifies a diagram with its mirror image."""
    if not g.edges:
        return ('unknot',)
    from .graph import dual
    best = _map_code(g)
    if colour_invariant:
        best = min(best, _map_code(dual(g)))
    if include_mirror:
        m = canonical_code(g.mirror(), include_mirror=False, colour_invariant=colour_invariant)
        best = min(best, m)
    return best


def canonical_string(g: EmbeddedTaitGraph, include_mirror: bool = False, colour_invariant: bool = True) -> str:
    """Compact string form of the canonical code (for file keys / hashing)."""
    return repr(canonical_code(g, include_mirror, colour_invariant))
