"""Complete flype enumeration on embedded Tait graphs, and the flype orbit of a diagram.

A flype moves a crossing c across a 2-tangle T (a disc meeting the diagram in four points) and
rotates T by 180 degrees.  On the checkerboard graph in which the two shaded regions crossed by
the boundary of T are N and S, the crossing c is an edge N--S and T is a union of CONSECUTIVE
"beads" (bridges) of the 2-separation {N, S} of g - c, adjacent to c.  Two things are needed for
completeness that `moves.find_flypes` (used by the random walks of the search) does not do:

  * T may be a union of several consecutive beads, not just one;
  * the tangle boundary may cross two UNshaded regions instead, in which case the same
    configuration is visible only on the dual checkerboard graph.

`flype_orbit` closes a diagram under all flypes on both checkerboard graphs.  For a reduced
alternating diagram of a prime alternating knot the orbit is the set of ALL its minimal diagrams
(Kauffman-Murasugi-Thistlethwaite + Menasco-Thistlethwaite), as oriented maps on the sphere;
`turn_over` gives the diagram seen from the other side of the sphere.
"""
from __future__ import annotations

from . import moves as mv
from .graph import dual
from .canonical import canonical_code


def _merge(pieces):
    P = {'vertices': set(), 'edges': set(), 'N_edges': [], 'S_edges': []}
    for Q in pieces:
        P['vertices'] |= Q['vertices']
        P['edges'] |= Q['edges']
        P['N_edges'] += Q['N_edges']
        P['S_edges'] += Q['S_edges']
    return P


def find_flypes_all(g):
    """Yield (e, P): a crossing e = (N, S) and a tangle P (merged piece dict) adjacent to e."""
    for e, (N, S) in g.edges.items():
        if N == S:
            continue
        pieces = mv._pieces(g, N, S, e)
        rN = g.rotation[N]
        n, iN = len(rN), rN.index(e)
        beads = []
        for P in pieces:
            b = mv._block(rN, P['N_edges'])
            if b is not None:
                beads.append(((b[0] - iN) % n, P))          # position after e, going round N
        beads.sort(key=lambda t: t[0])
        k = len(beads)
        seen = set()
        for m in range(1, k + 1):
            for sel in (beads[:m], beads[k - m:]):
                P = _merge([p for _, p in sel])
                key = frozenset(P['edges'])
                if key in seen or not P['vertices']:
                    continue
                seen.add(key)
                if mv._flype_params(g, e, N, S, P) is not None:
                    yield e, P


def flype_block(g, e, P):
    """The flype of crossing e across the tangle P (a merged piece from find_flypes_all)."""
    N, S = g.edges[e]
    params = mv._flype_params(g, e, N, S, P)
    if params is None:
        raise ValueError('flype_block: tangle is not adjacent to the crossing')
    side, (sN, mN), (sS, mS) = params
    h = g.copy()
    rN, rS = g.rotation[N], g.rotation[S]
    blockN = [rN[(sN + k) % len(rN)] for k in range(mN)]
    blockS = [rS[(sS + k) % len(rS)] for k in range(mS)]

    def rebuild(rot, start, m, newblock, e_first):
        n = len(rot)
        first = (start - 1) % n if e_first else start % n
        order = [rot[(first + k) % n] for k in range(n)]
        return newblock + order[m + 1:]

    if side == 'N_before':
        newN = rebuild(rN, sN, mN, list(reversed(blockS)) + [e], e_first=True)
        newS = rebuild(rS, sS, mS, [e] + list(reversed(blockN)), e_first=False)
    else:
        newN = rebuild(rN, sN, mN, [e] + list(reversed(blockS)), e_first=False)
        newS = rebuild(rS, sS, mS, list(reversed(blockN)) + [e], e_first=True)
    h.rotation[N] = newN
    h.rotation[S] = newS
    for f_ in P['edges']:
        a, b = h.edges[f_]
        a2 = S if a == N else (N if a == S else a)
        b2 = S if b == N else (N if b == S else b)
        h.edges[f_] = (a2, b2)
    for v in P['vertices']:
        h.rotation[v] = list(reversed(h.rotation[v]))
    return mv._check(h, 'flype')


def turn_over(g):
    """The same diagram seen from the other side of the sphere (rotations reversed, signs kept)."""
    h = g.copy()
    for v in h.vertices:
        h.rotation[v] = list(reversed(h.rotation[v]))
    return h


def flype_orbit(g0) -> dict:
    """All diagrams reachable by flypes (both checkerboard graphs): canonical code -> graph."""
    code = lambda h: repr(canonical_code(h))  # noqa: E731
    orbit = {code(g0): g0}
    frontier = [g0]
    while frontier:
        g = frontier.pop()
        for h0 in (g, dual(g)):
            for e, P in find_flypes_all(h0):
                h = flype_block(h0, e, P)
                c = code(h)
                if c not in orbit:
                    orbit[c] = h
                    frontier.append(h)
    return orbit


def orbit_up_to_turning_over(orbit: dict) -> list:
    """Group an orbit into classes {D, turn_over(D)}; returns a list of lists of codes."""
    classes = {}
    for c, h in orbit.items():
        ct = repr(canonical_code(turn_over(h)))
        classes.setdefault(min(c, ct), []).append(c)
    return list(classes.values())
