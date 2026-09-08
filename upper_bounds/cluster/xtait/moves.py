"""Exact isotopy moves on embedded Tait graphs (stdlib only; frozen copy of the verified project code).

  crossing change    flip the sign of one edge                                    (changes the knot)
  R1-  r1_remove     delete an empty-kink loop, or a pendant edge with its vertex
  R1+  r1_add_loop / r1_add_pendant
  R2-  r2_remove_parallel (opposite-sign bigon) / r2_remove_series (degree-2 vertex, opposite signs)
  R2+  r2_add_parallel (two opposite-sign edges inside a face) / r2_add_series (vertex split by a 2-path)
  R3   r3_delta_to_y / r3_y_to_delta  (Y–Δ with all three signs negated; edge opposite leaf L <- -sign(t–L))
  flype              reflect a tangle piece of a 2-separation and move the adjacent crossing across

All functions return a NEW graph and raise ValueError when not applicable; results are validated
(planarity via Euler).  Exactness (knot type preserved) is checked in selftest.py with the Kauffman
bracket as an independent oracle.
"""
from __future__ import annotations

import itertools

from .graph import EmbeddedTaitGraph, ribbon_faces


def _positions(rot, e):
    return [i for i, x in enumerate(rot) if x == e]


def _cyclic_adjacent(i, j, n):
    return n >= 2 and ((j - i) % n == 1 or (i - j) % n == 1)


def _replace_slot(rot, idx, new):
    return rot[:idx] + list(new) + rot[idx + 1:]


def _check(h, what):
    try:
        h.validate()
    except ValueError as exc:
        raise ValueError(f'{what}: result invalid ({exc})') from exc
    return h


# ----------------------------------------------------------------------------- crossing change
def crossing_change(g, e):
    h = g.copy()
    h.signs[e] = -h.signs[e]
    return h


# ----------------------------------------------------------------------------- R1
def find_r1_removals(g):
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


def r1_remove(g, e):
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
    return _check(h, 'R1 remove')


def r1_add_loop(g, v, pos, sign):
    h = g.copy()
    e = h.new_edge_id()
    r = h.rotation[v]
    pos = pos % (len(r) + 1) if r else 0
    h.rotation[v] = r[:pos] + [e, e] + r[pos:]
    h.edges[e] = (v, v)
    h.signs[e] = sign
    return _check(h, 'R1 add loop')


def r1_add_pendant(g, v, pos, sign):
    h = g.copy()
    e = h.new_edge_id()
    w = h.new_vertex_id()
    r = h.rotation[v]
    pos = pos % (len(r) + 1) if r else 0
    h.rotation[v] = r[:pos] + [e] + r[pos:]
    h.vertices.append(w)
    h.rotation[w] = [e]
    h.edges[e] = (v, w)
    h.signs[e] = sign
    return _check(h, 'R1 add pendant')


# ----------------------------------------------------------------------------- R2
def find_r2_removals(g):
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


def r2_remove_parallel(g, e1, e2):
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


def r2_remove_series(g, w):
    if ('series', (w,)) not in find_r2_removals(g):
        raise ValueError('R2: vertex is not an opposite-sign series pair')
    h = g.copy()
    e1, e2 = h.rotation[w]
    u, v = h.other_end(e1, w), h.other_end(e2, w)
    if u != v:
        ru, rv = h.rotation[u], h.rotation[v]
        iu = ru.index(e1)
        iv = rv.index(e2)
        spliced = rv[iv + 1:] + rv[:iv]
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


def find_r2_parallel_additions(g):
    out = []
    for fi, f in enumerate(ribbon_faces(g)):
        for i, j in itertools.combinations(range(len(f)), 2):
            if f[i][0] != f[j][0]:
                out.append((fi, i, j))
    return out


def r2_add_parallel(g, face_index, i, j, sign=1):
    faces = ribbon_faces(g)
    f = faces[face_index]
    (u, eu, ku), (v, ev, kv) = f[i], f[j]
    if u == v:
        raise ValueError('R2 add parallel: darts must be at distinct vertices')
    h = g.copy()
    f1 = h.new_edge_id()
    f2 = f1 + 1
    ru = h.rotation[u]
    pu = _positions(ru, eu)[ku]
    h.rotation[u] = ru[:pu] + [f1, f2] + ru[pu:]
    rv = h.rotation[v]
    pv = _positions(rv, ev)[kv]
    h.rotation[v] = rv[:pv] + [f2, f1] + rv[pv:]
    h.edges[f1] = (u, v)
    h.edges[f2] = (u, v)
    h.signs[f1] = sign
    h.signs[f2] = -sign
    return _check(h, 'R2 add parallel')


def find_r2_series_additions(g):
    out = []
    for v in g.vertices:
        d = g.degree(v)
        if d >= 2:
            for i in range(d):
                for j in range(d):
                    if i != j:
                        out.append((v, i, j))
    return out


def r2_add_series(g, v, i, j, sign=1):
    r = g.rotation[v]
    d = len(r)
    if d < 2 or i == j:
        raise ValueError('R2 add series: need two distinct split positions')
    arcA = [r[(i + k) % d] for k in range((j - i) % d)]
    arcB = [r[(j + k) % d] for k in range((i - j) % d)]
    h = g.copy()
    f1 = h.new_edge_id()
    f2 = f1 + 1
    w = h.new_vertex_id()
    v2 = w + 1
    h.vertices += [w, v2]
    h.rotation[v] = arcA + [f1]
    h.rotation[w] = [f1, f2]
    h.rotation[v2] = arcB + [f2]
    for e in set(arcB):
        a, b = h.edges[e]
        if a == v and b == v:
            if arcB.count(e) == 2:
                h.edges[e] = (v2, v2)
            elif arcB.count(e) == 1:
                raise ValueError('R2 add series: split would cut a loop')
        else:
            h.edges[e] = (v2 if a == v else a, v2 if b == v else b)
    h.edges[f1] = (v, w)
    h.edges[f2] = (w, v2)
    h.signs[f1] = sign
    h.signs[f2] = -sign
    return _check(h, 'R2 add series')


# ----------------------------------------------------------------------------- R3
def find_r3_triangles(g):
    out = []
    for fi, f in enumerate(ribbon_faces(g)):
        if len(f) == 3:
            vs = [d[0] for d in f]
            es = [d[1] for d in f]
            if len(set(vs)) == 3 and len(set(es)) == 3 and len({g.signs[e] for e in es}) == 2:
                out.append(fi)
    return out


def find_r3_stars(g):
    out = []
    for t in g.vertices:
        r = g.rotation[t]
        if len(r) == 3 and len(set(r)) == 3 and not any(g.is_loop(e) for e in r):
            leaves = [g.other_end(e, t) for e in r]
            if len(set(leaves)) == 3 and len({g.signs[e] for e in r}) == 2:
                out.append(t)
    return out


def r3_y_to_delta(g, t):
    if t not in find_r3_stars(g):
        raise ValueError('R3: not an applicable star')
    h = g.copy()
    es = list(h.rotation[t])
    leaves = [h.other_end(e, t) for e in es]
    base = h.new_edge_id()
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


def r3_delta_to_y(g, face_index):
    faces = ribbon_faces(g)
    if face_index not in find_r3_triangles(g):
        raise ValueError('R3: not an applicable triangle')
    f = faces[face_index]
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
        pidx = (idx - 1) % len(r)
        if pidx < idx:
            h.rotation[V] = r[:pidx] + [e_new] + r[idx + 1:]
        else:
            h.rotation[V] = [e_new] + r[1:pidx]
    h.vertices.append(t)
    h.rotation[t] = [new_e[V] for V in reversed(verts)]
    for te in tri_edges:
        del h.edges[te]
        del h.signs[te]
    return _check(h, 'R3 Δ->Y')


# ----------------------------------------------------------------------------- flype
def _pieces(g, N, S, e):
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
    comps = {}
    for v in others:
        comps.setdefault(find(v), {'vertices': set(), 'edges': set(), 'N_edges': [], 'S_edges': []})['vertices'].add(v)
    for f_, (a, b) in g.edges.items():
        if f_ == e:
            continue
        if a in (N, S) and b in (N, S):
            comps.setdefault(('triv', f_), {'vertices': set(), 'edges': {f_}, 'N_edges': [], 'S_edges': []})
            c = comps[('triv', f_)]
            if N in (a, b):
                c['N_edges'].append(f_)
            if S in (a, b):
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


def _block(rot, edges):
    idx = sorted(rot.index(x) for x in edges)
    n, m = len(rot), len(idx)
    if m == 0:
        return None
    idxset = set(idx)
    for s in idx:
        if all(((s + k) % n) in idxset for k in range(m)):
            return s, m
    return None


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
    if before_N and after_S:
        side = 'N_before'
    elif after_N and before_S:
        side = 'N_after'
    else:
        return None
    return side, (sN, mN), (sS, mS)


def find_flypes(g):
    out = []
    for e, (N, S) in g.edges.items():
        if N == S:
            continue
        for pi, P in enumerate(_pieces(g, N, S, e)):
            if not P['vertices'] or not P['N_edges'] or not P['S_edges']:
                continue
            if _flype_params(g, e, N, S, P) is not None:
                out.append((e, pi))
    return out


def flype(g, e, piece_index):
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
    return _check(h, 'flype')


# ----------------------------------------------------------------------------- dispatch
def apply_move(g, kind, params):
    if kind == 'r1-':
        return r1_remove(g, *params)
    if kind == 'r2-p':
        return r2_remove_parallel(g, *params)
    if kind == 'r2-s':
        return r2_remove_series(g, *params)
    if kind == 'r3D':
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
    if kind == 'pass':
        from .passmove import apply_pass_params
        return apply_pass_params(g, params)
    raise ValueError(f'unknown move kind {kind}')
