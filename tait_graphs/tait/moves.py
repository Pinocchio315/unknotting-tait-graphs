"""Reidemeister moves, flypes and crossing changes on embedded Tait graphs.

All functions return a NEW EmbeddedTaitGraph (inputs are never modified) and raise ValueError when the
requested move is not applicable.  `find_*` functions enumerate the applicable move parameters.

Dictionary (white graph = one checkerboard colour class; vertices = white faces, edges = crossings):

  crossing change   : flip the sign of one edge                              (changes the knot)
  R1  remove        : delete a loop edge whose two ends are adjacent in the rotation (empty kink),
                      or delete a pendant edge together with its degree-1 vertex (white monogon)
  R1  add           : insert a loop at a slot of a vertex, or attach a new pendant vertex at a slot
  R2  remove        : delete two parallel edges of opposite sign bounding a bigon face,
                      or contract the two opposite-sign edges of a degree-2 vertex (series pair)
  R2  add           : two parallel opposite-sign edges between two vertices of one face,
                      or split a vertex into two joined by a 2-path through a new degree-2 vertex
  R3                : Y–Δ / Δ–Y with all three signs negated (triangle edge opposite leaf L gets
                      -sign(edge t–L)); applicable iff the three signs are not all equal
  flype             : for an edge e=(N,S) and a tangle piece P of the 2-separation {N,S} adjacent to
                      e: reflect P (reverse its rotations, swap its N/S attachments) and move e to the
                      other side of P; signs unchanged

The sign rules were derived from the spherogram/Gordon–Litherland convention used in graph.py and are
checked in tests/test_tait.py against knot invariants (HFK, Jones) and against spherogram's own R3.
"""
from __future__ import annotations

import itertools
from collections import Counter

from .graph import EmbeddedTaitGraph, ribbon_faces, Dart


# ----------------------------------------------------------------------------- helpers
def _positions(rot: list[int], e: int) -> list[int]:
    return [i for i, x in enumerate(rot) if x == e]


def _cyclic_adjacent(i: int, j: int, n: int) -> bool:
    return n >= 2 and ((j - i) % n == 1 or (i - j) % n == 1)


def _replace_slot(rot: list[int], idx: int, new: list[int]) -> list[int]:
    return rot[:idx] + list(new) + rot[idx + 1:]


def _check(h: EmbeddedTaitGraph, what: str) -> EmbeddedTaitGraph:
    try:
        h.validate()
    except ValueError as exc:
        raise ValueError(f'{what}: result invalid ({exc})') from exc
    return h


# ----------------------------------------------------------------------------- crossing change
def crossing_change(g: EmbeddedTaitGraph, e: int) -> EmbeddedTaitGraph:
    h = g.copy()
    h.signs[e] = -h.signs[e]
    return h


# ----------------------------------------------------------------------------- R1
def find_r1_removals(g: EmbeddedTaitGraph) -> list[tuple[str, int]]:
    """[('loop', e), ('pendant', e)] removable by R1."""
    out = []
    for e, (u, w) in g.edges.items():
        if u == w:
            pos = _positions(g.rotation[u], e)
            if len(pos) == 2 and _cyclic_adjacent(pos[0], pos[1], len(g.rotation[u])):
                out.append(('loop', e))
        else:
            if g.degree(u) == 1 or g.degree(w) == 1:
                out.append(('pendant', e))
    return out


def r1_remove(g: EmbeddedTaitGraph, e: int) -> EmbeddedTaitGraph:
    h = g.copy()
    u, w = h.edges[e]
    if u == w:
        pos = _positions(h.rotation[u], e)
        if len(pos) != 2 or not _cyclic_adjacent(pos[0], pos[1], len(h.rotation[u])):
            raise ValueError('R1: loop is not an empty kink')
        h.rotation[u] = [x for x in h.rotation[u] if x != e]
    else:
        if h.degree(u) == 1:
            pend, other = u, w
        elif h.degree(w) == 1:
            pend, other = w, u
        else:
            raise ValueError('R1: edge is neither a kink loop nor pendant')
        h.rotation[other] = [x for x in h.rotation[other] if x != e]
        del h.rotation[pend]
        h.vertices.remove(pend)
    del h.edges[e]
    del h.signs[e]
    if len(h.vertices) == 0:  # cannot happen (pendant removal keeps `other`)
        raise ValueError('R1: empty graph')
    return _check(h, 'R1 remove')


def r1_add_loop(g: EmbeddedTaitGraph, v: int, pos: int, sign: int) -> EmbeddedTaitGraph:
    """Insert an empty kink (loop) at vertex v before rotation index pos (0 <= pos <= deg)."""
    h = g.copy()
    e = h.new_edge_id()
    r = h.rotation[v]
    pos %= max(1, len(r) + 1) if r else 1
    h.rotation[v] = r[:pos] + [e, e] + r[pos:]
    h.edges[e] = (v, v)
    h.signs[e] = sign
    return _check(h, 'R1 add loop')


def r1_add_pendant(g: EmbeddedTaitGraph, v: int, pos: int, sign: int) -> EmbeddedTaitGraph:
    """Attach a new degree-1 vertex to v before rotation index pos (white monogon kink)."""
    h = g.copy()
    e = h.new_edge_id()
    w = h.new_vertex_id()
    r = h.rotation[v]
    pos %= max(1, len(r) + 1) if r else 1
    h.rotation[v] = r[:pos] + [e] + r[pos:]
    h.vertices.append(w)
    h.rotation[w] = [e]
    h.edges[e] = (v, w)
    h.signs[e] = sign
    return _check(h, 'R1 add pendant')


# ----------------------------------------------------------------------------- R2
def find_r2_removals(g: EmbeddedTaitGraph) -> list[tuple[str, tuple]]:
    """[('parallel', (e1, e2)), ('series', (w,))] removable by R2."""
    out = []
    for f in ribbon_faces(g):
        if len(f) == 2:
            (v1, e1, _), (v2, e2, _) = f
            if e1 != e2 and not g.is_loop(e1) and not g.is_loop(e2) and g.signs[e1] == -g.signs[e2] \
                    and set(g.edges[e1]) == set(g.edges[e2]) and g.edges[e1][0] != g.edges[e1][1]:
                out.append(('parallel', tuple(sorted((e1, e2)))))
    for w in g.vertices:
        r = g.rotation[w]
        if len(r) == 2 and r[0] != r[1] and g.signs[r[0]] == -g.signs[r[1]] \
                and not g.is_loop(r[0]) and not g.is_loop(r[1]):
            out.append(('series', (w,)))
    return out


def r2_remove_parallel(g: EmbeddedTaitGraph, e1: int, e2: int) -> EmbeddedTaitGraph:
    if ('parallel', tuple(sorted((e1, e2)))) not in find_r2_removals(g):
        raise ValueError('R2: edges do not form an opposite-sign bigon')
    h = g.copy()
    u, w = h.edges[e1]
    for v in (u, w):
        h.rotation[v] = [x for x in h.rotation[v] if x not in (e1, e2)]
    for e in (e1, e2):
        del h.edges[e]
        del h.signs[e]
    return _check(h, 'R2 remove parallel')


def r2_remove_series(g: EmbeddedTaitGraph, w: int) -> EmbeddedTaitGraph:
    if ('series', (w,)) not in find_r2_removals(g):
        raise ValueError('R2: vertex is not an opposite-sign series pair')
    h = g.copy()
    e1, e2 = h.rotation[w]
    u, v = h.other_end(e1, w), h.other_end(e2, w)
    if u != v:
        ru, rv = h.rotation[u], h.rotation[v]
        iu = ru.index(e1)
        iv = rv.index(e2)
        spliced = rv[iv + 1:] + rv[:iv]          # v's edges starting after e2, cyclically
        h.rotation[u] = ru[:iu] + spliced + ru[iu + 1:]
        for e in spliced:
            a, b = h.edges[e]
            h.edges[e] = (u if a == v else a, u if b == v else b)
        del h.rotation[v]
        h.vertices.remove(v)
    else:
        h.rotation[u] = [x for x in h.rotation[u] if x not in (e1, e2)]
    del h.rotation[w]
    h.vertices.remove(w)
    for e in (e1, e2):
        del h.edges[e]
        del h.signs[e]
    return _check(h, 'R2 remove series')


def find_r2_parallel_additions(g: EmbeddedTaitGraph) -> list[tuple[int, int, int]]:
    """(face index, i, j): darts i<j of that face at distinct vertices -> two new parallel edges."""
    out = []
    for fi, f in enumerate(ribbon_faces(g)):
        for i, j in itertools.combinations(range(len(f)), 2):
            if f[i][0] != f[j][0]:
                out.append((fi, i, j))
    return out


def r2_add_parallel(g: EmbeddedTaitGraph, face_index: int, i: int, j: int, sign: int = 1) -> EmbeddedTaitGraph:
    """Push the arc of face F at dart i across the arc at dart j: two new parallel edges (opposite signs)
    between the white vertices of the two darts, inserted inside face F (before the dart's edge)."""
    faces = ribbon_faces(g)
    f = faces[face_index]
    (u, eu, ku), (v, ev, kv) = f[i], f[j]
    if u == v:
        raise ValueError('R2 add parallel: darts must be at distinct vertices')
    h = g.copy()
    f1, f2 = h.new_edge_id(), h.new_edge_id() + 1
    ru = h.rotation[u]
    pu = [p for p in _positions(ru, eu)][ku]
    h.rotation[u] = ru[:pu] + [f1, f2] + ru[pu:]
    rv = h.rotation[v]
    pv = [p for p in _positions(rv, ev)][kv]
    h.rotation[v] = rv[:pv] + [f2, f1] + rv[pv:]
    h.edges[f1] = (u, v)
    h.edges[f2] = (u, v)
    h.signs[f1] = sign
    h.signs[f2] = -sign
    return _check(h, 'R2 add parallel')


def find_r2_series_additions(g: EmbeddedTaitGraph) -> list[tuple[int, int, int]]:
    """(v, i, j): split rotation[v] into arcs [i:j) and [j:i) (both non-empty)."""
    out = []
    for v in g.vertices:
        d = g.degree(v)
        if d >= 2:
            for i in range(d):
                for j in range(d):
                    if i != j:
                        out.append((v, i, j))
    return out


def r2_add_series(g: EmbeddedTaitGraph, v: int, i: int, j: int, sign: int = 1) -> EmbeddedTaitGraph:
    """Push a strand across white region v: v splits into v (keeping rotation[i:j)) and v2 (the rest),
    joined by a 2-path v –f1– w –f2– v2 through a new degree-2 vertex w; signs (sign, -sign)."""
    r = g.rotation[v]
    d = len(r)
    if d < 2 or i == j:
        raise ValueError('R2 add series: need two distinct split positions')
    arcA = [r[(i + k) % d] for k in range((j - i) % d)]
    arcB = [r[(j + k) % d] for k in range((i - j) % d)]
    h = g.copy()
    f1, f2 = h.new_edge_id(), h.new_edge_id() + 1
    w, v2 = h.new_vertex_id(), h.new_vertex_id() + 1
    h.vertices += [w, v2]
    h.rotation[v] = arcA + [f1]
    h.rotation[w] = [f1, f2]
    h.rotation[v2] = arcB + [f2]
    for e in arcB:
        a, b = h.edges[e]
        if a == v and b == v:          # loop at v split across the two arcs? both ends must be in arcB
            if arcB.count(e) == 2:
                h.edges[e] = (v2, v2)
            else:
                raise ValueError('R2 add series: split would cut a loop')
        else:
            h.edges[e] = (v2 if a == v else a, v2 if b == v else b)
    h.edges[f1] = (v, w)
    h.edges[f2] = (w, v2)
    h.signs[f1] = sign
    h.signs[f2] = -sign
    return _check(h, 'R2 add series')


# ----------------------------------------------------------------------------- R3
def find_r3_triangles(g: EmbeddedTaitGraph) -> list[int]:
    """Indices of faces of length 3 (distinct vertices and edges) whose signs are not all equal."""
    out = []
    for fi, f in enumerate(ribbon_faces(g)):
        if len(f) == 3:
            vs = [d[0] for d in f]
            es = [d[1] for d in f]
            if len(set(vs)) == 3 and len(set(es)) == 3 and len({g.signs[e] for e in es}) == 2:
                out.append(fi)
    return out


def find_r3_stars(g: EmbeddedTaitGraph) -> list[int]:
    """Degree-3 vertices with three distinct non-loop edges to three distinct leaves, signs not all equal."""
    out = []
    for t in g.vertices:
        r = g.rotation[t]
        if len(r) == 3 and len(set(r)) == 3 and not any(g.is_loop(e) for e in r):
            leaves = [g.other_end(e, t) for e in r]
            if len(set(leaves)) == 3 and len({g.signs[e] for e in r}) == 2:
                out.append(t)
    return out


def r3_y_to_delta(g: EmbeddedTaitGraph, t: int) -> EmbeddedTaitGraph:
    """Star at t (rotation [e0,e1,e2], leaves L0,L1,L2) -> triangle: at leaf Li replace ei by
    [edge(Li,L_{i+1}), edge(Li,L_{i-1})]; the edge opposite leaf Lk gets sign -sign(ek)."""
    if t not in find_r3_stars(g):
        raise ValueError('R3: not an applicable star')
    h = g.copy()
    es = list(h.rotation[t])
    leaves = [h.other_end(e, t) for e in es]
    base = h.new_edge_id()
    # new edge between Li and Lj (i<j) opposite leaf Lk
    new = {}
    for k in range(3):
        i, j = [x for x in range(3) if x != k]
        new[frozenset((i, j))] = base + k
        h.edges[base + k] = (leaves[i], leaves[j])
        h.signs[base + k] = -h.signs[es[k]]
    for i in range(3):
        Li = leaves[i]
        f_next = new[frozenset((i, (i + 1) % 3))]
        f_prev = new[frozenset((i, (i - 1) % 3))]
        r = h.rotation[Li]
        idx = r.index(es[i])
        h.rotation[Li] = _replace_slot(r, idx, [f_next, f_prev])
    for e in es:
        del h.edges[e]
        del h.signs[e]
    del h.rotation[t]
    h.vertices.remove(t)
    return _check(h, 'R3 Y->Δ')


def r3_delta_to_y(g: EmbeddedTaitGraph, face_index: int) -> EmbeddedTaitGraph:
    """Triangle face [(V0,f0),(V1,f1),(V2,f2)] -> star: at Vi the pair [prev_Vi(fi), fi] is replaced by
    a new edge ei to the new centre t; t's rotation is the reverse of the face order; ei gets
    -sign(triangle edge opposite Vi)."""
    faces = ribbon_faces(g)
    f = faces[face_index]
    if face_index not in find_r3_triangles(g):
        raise ValueError('R3: not an applicable triangle')
    h = g.copy()
    t = h.new_vertex_id()
    base = h.new_edge_id()
    tri_edges = [d[1] for d in f]
    verts = [d[0] for d in f]
    new_e = {}
    for i, (V, fe, _) in enumerate(f):
        r = h.rotation[V]
        idx = r.index(fe)
        prev_e = r[(idx - 1) % len(r)]
        if prev_e not in tri_edges or prev_e == fe:
            raise ValueError('R3: triangle edges are not consecutive at a vertex')
        e_new = base + i
        new_e[V] = e_new
        opposite = [te for te in tri_edges if V not in h.edges[te]]
        if len(opposite) != 1:
            raise ValueError('R3: could not find opposite edge')
        h.edges[e_new] = (V, t)
        h.signs[e_new] = -h.signs[opposite[0]]
        # replace the consecutive pair [prev_e, fe] by e_new (handle wrap-around)
        pidx = (idx - 1) % len(r)
        if pidx < idx:
            h.rotation[V] = r[:pidx] + [e_new] + r[idx + 1:]
        else:  # pidx == len(r)-1, idx == 0
            h.rotation[V] = [e_new] + r[1:pidx]
    h.vertices.append(t)
    h.rotation[t] = [new_e[V] for V in reversed(verts)]
    for te in tri_edges:
        del h.edges[te]
        del h.signs[te]
    return _check(h, 'R3 Δ->Y')


# ----------------------------------------------------------------------------- flype
def _pieces(g: EmbeddedTaitGraph, N: int, S: int, e: int):
    """Connected pieces of G - e - {N,S}: list of dicts with 'vertices' (set), 'edges' (set, incl. attachments),
    'N_edges', 'S_edges' (edges attached at N / S)."""
    others = [v for v in g.vertices if v not in (N, S)]
    parent = {v: v for v in others}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for f_, (a, b) in g.edges.items():
        if f_ == e or a in (N, S) or b in (N, S):
            continue
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    comps: dict[int, dict] = {}
    for v in others:
        comps.setdefault(find(v), {'vertices': set(), 'edges': set(), 'N_edges': [], 'S_edges': []})['vertices'].add(v)
    for f_, (a, b) in g.edges.items():
        if f_ == e:
            continue
        if a in (N, S) and b in (N, S):
            # trivial piece: an N-S edge (or loop at N or S) — treat separately
            comps.setdefault(('triv', f_), {'vertices': set(), 'edges': {f_}, 'N_edges': [], 'S_edges': []})
            c = comps[('triv', f_)]
            if a == N or b == N:
                c['N_edges'].append(f_)
            if a == S or b == S:
                c['S_edges'].append(f_)
            continue
        x = a if a not in (N, S) else b
        c = comps[find(x)]
        c['edges'].add(f_)
        if N in (a, b):
            c['N_edges'].append(f_)
        if S in (a, b):
            c['S_edges'].append(f_)
    return list(comps.values())


def _block(rot: list[int], edges: list[int]):
    """If `edges` occupy a cyclically contiguous block of rot, return (start, length) else None."""
    idx = sorted(rot.index(x) for x in edges)  # edges at N/S are not loops here, each once
    n, m = len(rot), len(idx)
    if m == 0:
        return None
    for s in idx:
        if all(((s + k) % n) in idx for k in range(m)):
            return s, m
    return None


def find_flypes(g: EmbeddedTaitGraph) -> list[tuple[int, int]]:
    """(e, piece_index) pairs where the edge e=(N,S) can be flyped across the tangle piece."""
    out = []
    for e, (N, S) in g.edges.items():
        if N == S:
            continue
        pieces = _pieces(g, N, S, e)
        for pi, P in enumerate(pieces):
            if not P['vertices'] or not P['N_edges'] or not P['S_edges']:
                continue  # trivial piece or one-sided lobe
            if _flype_params(g, e, N, S, P) is not None:
                out.append((e, pi))
    return out


def _flype_params(g, e, N, S, P):
    rN, rS = g.rotation[N], g.rotation[S]
    bN = _block(rN, P['N_edges'])
    bS = _block(rS, P['S_edges'])
    if bN is None or bS is None:
        return None
    sN, mN = bN
    sS, mS = bS
    iN, iS = rN.index(e), rS.index(e)
    nN, nS = len(rN), len(rS)
    before_N = (iN == (sN - 1) % nN)
    after_N = (iN == (sN + mN) % nN)
    before_S = (iS == (sS - 1) % nS)
    after_S = (iS == (sS + mS) % nS)
    if before_N and after_S and not (after_N and before_S):
        side = 'N_before'
    elif after_N and before_S and not (before_N and after_S):
        side = 'N_after'
    elif (before_N and after_S) and (after_N and before_S):
        side = 'N_before'   # degenerate (only e and P around N,S): both readings coincide
    else:
        return None
    return side, (sN, mN), (sS, mS)


def flype(g: EmbeddedTaitGraph, e: int, piece_index: int) -> EmbeddedTaitGraph:
    N, S = g.edges[e]
    pieces = _pieces(g, N, S, e)
    P = pieces[piece_index]
    params = _flype_params(g, e, N, S, P)
    if params is None:
        raise ValueError('flype: edge/piece not adjacent beads of the 2-separation')
    side, (sN, mN), (sS, mS) = params
    h = g.copy()
    rN, rS = g.rotation[N], g.rotation[S]
    blockN = [rN[(sN + k) % len(rN)] for k in range(mN)]
    blockS = [rS[(sS + k) % len(rS)] for k in range(mS)]
    # rotation at N: remove e and blockN, insert new block (reversed blockS) with e on the other side
    restN = [x for x in rN if x != e and x not in blockN]
    # position: keep the cyclic context: find the element before the removed region
    def rebuild(rot, removed_block, e_idx, start, m, newblock, e_first):
        n = len(rot)
        region = set(removed_block) | {e}
        # walk the cycle starting right after the region to collect the remaining elements in order
        # the region occupies [start-1 .. start+m-1] (e before) or [start .. start+m] (e after)
        if e_first:
            first = (start - 1) % n
        else:
            first = start % n
        order = [rot[(first + k) % n] for k in range(n)]
        # order begins with the region (m+1 elements), then the rest
        rest = order[m + 1:]
        assert set(order[:m + 1]) == region, 'flype: block bookkeeping failed'
        return newblock + rest
    if side == 'N_before':     # N: [e, blockN]  S: [blockS, e]   ->  N: [rev(blockS), e]   S: [e, rev(blockN)]
        newN = rebuild(rN, blockN, None, sN, mN, list(reversed(blockS)) + [e], e_first=True)
        newS = rebuild(rS, blockS, None, sS, mS, [e] + list(reversed(blockN)), e_first=False)
    else:                      # N: [blockN, e]  S: [e, blockS]   ->  N: [e, rev(blockS)]   S: [rev(blockN), e]
        newN = rebuild(rN, blockN, None, sN, mN, [e] + list(reversed(blockS)), e_first=False)
        newS = rebuild(rS, blockS, None, sS, mS, list(reversed(blockN)) + [e], e_first=True)
    h.rotation[N] = newN
    h.rotation[S] = newS
    # re-attach P's edges: N <-> S
    for f_ in P['edges']:
        a, b = h.edges[f_]
        a2 = S if a == N else (N if a == S else a)
        b2 = S if b == N else (N if b == S else b)
        h.edges[f_] = (a2, b2)
    # reflect the tangle: reverse rotations of P's internal vertices
    for v in P['vertices']:
        h.rotation[v] = list(reversed(h.rotation[v]))
    return _check(h, 'flype')


# ----------------------------------------------------------------------------- enumeration
def all_isotopy_moves(g: EmbeddedTaitGraph, include_increasing: bool = True,
                      include_flypes: bool = True) -> list[tuple[str, tuple]]:
    """List of (kind, params) for every applicable knot-type-preserving move.
    kinds: r1-, r2-, r3, flype, (r1+, r2+ when include_increasing)."""
    out: list[tuple[str, tuple]] = []
    for kind, e in find_r1_removals(g):
        out.append(('r1-', (e,)))
    for kind, params in find_r2_removals(g):
        out.append(('r2-' + kind[0], params))
    for fi in find_r3_triangles(g):
        out.append(('r3Δ', (fi,)))
    for t in find_r3_stars(g):
        out.append(('r3Y', (t,)))
    if include_flypes:
        for e, pi in find_flypes(g):
            out.append(('flype', (e, pi)))
    if include_increasing:
        for v in g.vertices:
            for pos in range(max(1, g.degree(v))):
                for s in (1, -1):
                    out.append(('r1+loop', (v, pos, s)))
                    out.append(('r1+pend', (v, pos, s)))
        for fi, i, j in find_r2_parallel_additions(g):
            for s in (1, -1):
                out.append(('r2+p', (fi, i, j, s)))
        for v, i, j in find_r2_series_additions(g):
            for s in (1, -1):
                out.append(('r2+s', (v, i, j, s)))
    return out


def twist(g: EmbeddedTaitGraph, r1kind: str, v: int, pos: int, s1: int, k: int, s2: int) -> EmbeddedTaitGraph:
    """Composite R1+R2 "twist": add a kink (loop/pendant, sign s1) at (v, pos) and immediately push its loop
    over a neighbouring arc (the k-th R2-parallel addition whose face pair touches the new edge, sign s2).
    Net: +3 crossings, writhe changes by s1, and the result has no R1-removable crossing (else ValueError).
    This reaches the writhe-changed diagrams that otherwise need kinky intermediate states."""
    h1 = r1_add_loop(g, v, pos, s1) if r1kind == 'loop' else r1_add_pendant(g, v, pos, s1)
    e_new = max(h1.edges)                                  # id assigned by the r1 move
    faces = ribbon_faces(h1)
    cands = []
    for (fi, i, j) in find_r2_parallel_additions(h1):
        f = faces[fi]
        if f[i][1] == e_new or f[j][1] == e_new:
            cands.append((fi, i, j))
    cands.sort()
    if k >= len(cands):
        raise ValueError('twist: no such R2 candidate')
    fi, i, j = cands[k]
    h2 = r2_add_parallel(h1, fi, i, j, s2)
    if find_r1_removals(h2):
        raise ValueError('twist: result still has a kink')
    return h2


def apply_move(g: EmbeddedTaitGraph, kind: str, params: tuple) -> EmbeddedTaitGraph:
    if kind == 'twist':
        return twist(g, *params)
    if kind == 'r1-':
        return r1_remove(g, *params)
    if kind == 'r2-p':
        return r2_remove_parallel(g, *params)
    if kind == 'r2-s':
        return r2_remove_series(g, *params)
    if kind == 'r3Δ':
        return r3_delta_to_y(g, *params)
    if kind == 'r3Y':
        return r3_y_to_delta(g, *params)
    if kind == 'flype':
        return flype(g, *params)
    if kind == 'r1+loop':
        return r1_add_loop(g, *params)
    if kind == 'r1+pend':
        return r1_add_pendant(g, *params)
    if kind == 'r2+p':
        return r2_add_parallel(g, *params)
    if kind == 'r2+s':
        return r2_add_series(g, *params)
    if kind == 'cc':
        return crossing_change(g, *params)
    raise ValueError(f'unknown move kind {kind}')
