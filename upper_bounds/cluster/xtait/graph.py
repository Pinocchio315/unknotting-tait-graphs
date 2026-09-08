"""Embedded Tait graphs (signed planar maps) and PD codes — stdlib only.

Diagram model
-------------
A knot diagram on S^2 is stored as one checkerboard colour class ("white" faces):

    vertices : list[int]              white faces
    edges    : dict[int, (int,int)]   crossing id -> (u, v)  (loops allowed)
    signs    : dict[int, int]         +1 / -1 per crossing (Gordon–Litherland style: +1 iff the white
                                      corners of the crossing are PD corners 0 and 2)
    rotation : dict[int, list[int]]   cyclic order of incident crossings around each vertex
                                      (clockwise; a loop appears twice)

PD codes are lists of quadruples [a,b,c,d] of arc labels, each label occurring exactly twice; the strand
through positions 0->2 is the UNDER strand, positions are listed counterclockwise around the crossing.

`to_pd` and `from_pd` are mutually consistent (see selftest: round trip preserves the canonical code and
the Jones polynomial, including chirality).
"""
from __future__ import annotations

import json
import random as _random
from dataclasses import dataclass, field

Dart = tuple  # (vertex, edge, slot)


# =============================================================================
# The graph object
# =============================================================================
@dataclass
class EmbeddedTaitGraph:
    vertices: list
    edges: dict
    signs: dict
    rotation: dict
    knot_name: str | None = None
    meta: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- basics
    def copy(self) -> "EmbeddedTaitGraph":
        return EmbeddedTaitGraph(list(self.vertices), dict(self.edges), dict(self.signs),
                                 {v: list(r) for v, r in self.rotation.items()}, self.knot_name, dict(self.meta))

    def n_crossings(self) -> int:
        return len(self.edges)

    def degree(self, v) -> int:
        return len(self.rotation[v])

    def other_end(self, e, v):
        u, w = self.edges[e]
        return w if v == u else u

    def is_loop(self, e) -> bool:
        u, w = self.edges[e]
        return u == w

    def new_vertex_id(self) -> int:
        return (max(self.vertices) + 1) if self.vertices else 0

    def new_edge_id(self) -> int:
        return (max(self.edges) + 1) if self.edges else 0

    def dart_list(self, v) -> list:
        seen = {}
        out = []
        for e in self.rotation[v]:
            k = seen.get(e, 0)
            seen[e] = k + 1
            out.append((v, e, k))
        return out

    def euler_characteristic(self) -> int:
        return len(self.vertices) - len(self.edges) + len(ribbon_faces(self))

    def is_connected(self) -> bool:
        if not self.vertices:
            return True
        seen = {self.vertices[0]}
        stack = [self.vertices[0]]
        while stack:
            v = stack.pop()
            for e in self.rotation[v]:
                w = self.other_end(e, v)
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return len(seen) == len(self.vertices)

    def validate(self) -> None:
        if not self.edges:
            if len(self.vertices) != 1 or self.rotation.get(self.vertices[0]) not in ([],):
                raise ValueError('edgeless diagram must be a single empty vertex (the 0-crossing unknot)')
            return
        for v in self.vertices:
            if v not in self.rotation:
                raise ValueError(f'vertex {v} has no rotation')
        for e, (u, w) in self.edges.items():
            if u not in self.rotation or w not in self.rotation:
                raise ValueError(f'edge {e} has unknown endpoint')
            if u == w:
                if self.rotation[u].count(e) != 2:
                    raise ValueError(f'loop {e} must appear twice at {u}')
            else:
                if self.rotation[u].count(e) != 1 or self.rotation[w].count(e) != 1:
                    raise ValueError(f'edge {e} must appear once at each endpoint')
            if self.signs.get(e) not in (1, -1):
                raise ValueError(f'edge {e} has no valid sign')
        for v, r in self.rotation.items():
            for e in r:
                if e not in self.edges or v not in self.edges[e]:
                    raise ValueError(f'rotation of {v} mentions non-incident edge {e}')
        if not self.is_connected():
            raise ValueError('graph is disconnected (split diagram)')
        if self.euler_characteristic() != 2:
            raise ValueError(f'rotation system is not spherical (chi={self.euler_characteristic()})')

    def mirror(self) -> "EmbeddedTaitGraph":
        h = self.copy()
        for v in h.vertices:
            h.rotation[v] = list(reversed(h.rotation[v]))
        for e in h.signs:
            h.signs[e] = -h.signs[e]
        return h


# =============================================================================
# Ribbon faces (rotation system -> face cycles of darts)
# =============================================================================
def _rho_theta(g: EmbeddedTaitGraph):
    rho = {}
    for v in g.vertices:
        dv = g.dart_list(v)
        for i, d in enumerate(dv):
            rho[d] = dv[(i + 1) % len(dv)]
    theta = {}
    for e, (u, w) in g.edges.items():
        if u != w:
            theta[(u, e, 0)] = (w, e, 0)
            theta[(w, e, 0)] = (u, e, 0)
        else:
            theta[(u, e, 0)] = (u, e, 1)
            theta[(u, e, 1)] = (u, e, 0)
    return rho, theta


def ribbon_faces(g: EmbeddedTaitGraph) -> list:
    rho, theta = _rho_theta(g)
    phi = {d: rho[theta[d]] for d in rho}
    seen, faces = set(), []
    for d in rho:
        if d in seen:
            continue
        cyc, cur = [], d
        while cur not in seen:
            seen.add(cur)
            cyc.append(cur)
            cur = phi[cur]
        faces.append(cyc)
    return faces


def dual(g: EmbeddedTaitGraph) -> EmbeddedTaitGraph:
    """The other checkerboard colour class of the same diagram (reversed face cycles, negated signs)."""
    faces = ribbon_faces(g)
    edge_faces = {e: [] for e in g.edges}
    rotation = {}
    for fi, f in enumerate(faces):
        seq = [d[1] for d in reversed(f)]
        rotation[fi] = seq
        for e in seq:
            edge_faces[e].append(fi)
    edges = {}
    for e, fl in edge_faces.items():
        if len(fl) != 2:
            raise ValueError(f'edge {e} incident to {len(fl)} face-slots')
        edges[e] = (fl[0], fl[1])
    h = EmbeddedTaitGraph(list(range(len(faces))), edges, {e: -s for e, s in g.signs.items()},
                          rotation, g.knot_name, dict(g.meta))
    h.validate()
    return h


# =============================================================================
# Graph -> PD  (medial construction; see the main project for the derivation)
# =============================================================================
def to_pd(g: EmbeddedTaitGraph) -> list:
    """PD code of the diagram.  Quadruple at edge e=(u,v): [arc_u_in, arc_u_out, arc_v_in, arc_v_out]
    listed counterclockwise; positions 0,2 = under strand.  sign +1 keeps the quadruple (white faces at
    corners 0,2); sign -1 rotates it by one position."""
    arc, prev, label = {}, {}, 1
    for v in g.vertices:
        dv = g.dart_list(v)
        for i, d in enumerate(dv):
            arc[d] = label
            label += 1
            prev[d] = dv[(i - 1) % len(dv)]
    pd = []
    for e in sorted(g.edges):
        u, w = g.edges[e]
        du = (u, e, 0)
        dw = (w, e, 1) if u == w else (w, e, 0)
        q = [arc[prev[du]], arc[du], arc[prev[dw]], arc[dw]]
        if g.signs[e] == -1:
            q = [q[1], q[2], q[3], q[0]]
        pd.append(q)
    return pd


# =============================================================================
# PD -> graph
# =============================================================================
def pd_faces(pd) -> tuple:
    """Faces of a PD as cycles of darts (crossing_index, position); face_of[(c,i)] = face index.

    The next corner of (c, i) is (c2, i2+1 mod 4) where (c2, i2) is the other occurrence of the arc
    label pd[c][i].  Euler's formula (faces == crossings + 2) is enforced."""
    occ = {}
    for c, q in enumerate(pd):
        if len(q) != 4:
            raise ValueError('PD quadruples must have 4 labels')
        for i, l in enumerate(q):
            occ.setdefault(int(l), []).append((c, i))
    for l, ds in occ.items():
        if len(ds) != 2:
            raise ValueError(f'label {l} occurs {len(ds)} times')
    partner = {}
    for l, (d1, d2) in occ.items():
        partner[d1] = d2
        partner[d2] = d1
    face_of, faces = {}, []
    for start in partner:
        if start in face_of:
            continue
        fi = len(faces)
        cyc, cur = [], start
        while cur not in face_of:
            face_of[cur] = fi
            cyc.append(cur)
            c2, i2 = partner[cur]
            cur = (c2, (i2 + 1) % 4)
        if cur != start:
            raise ValueError('face walk did not close at its start')
        faces.append(cyc)
    if len(faces) != len(pd) + 2:
        raise ValueError(f'PD is not a planar knot diagram (faces {len(faces)} != crossings+2)')
    return faces, face_of


def from_pd(pd, knot_name: str | None = None) -> EmbeddedTaitGraph:
    """Embedded Tait graph of a PD code (one checkerboard colour class).

    Corner (c,i) of the face structure is the corner between strand-positions i-1 and i; the checkerboard
    classes are {faces of corners 1,3} vs {faces of corners 0,2} of each crossing.  With the face-walk
    convention of `pd_faces`, taking the colour class of the corners (1,3) with sign +1 and rotation lists
    equal to the face cycles reproduces exactly the convention of `to_pd` — verified by the round-trip and
    Jones tests in selftest.py."""
    pd = [[int(x) for x in q] for q in pd]
    if not pd:
        return EmbeddedTaitGraph([0], {}, {}, {0: []}, knot_name)   # 0-crossing unknot
    faces, face_of = pd_faces(pd)
    n = len(pd)
    # 2-colour the faces: corners (c,1) and (c,3) same colour; (c,0),(c,2) same; (c,0) vs (c,1) opposite.
    colour = {}
    adj_same = {fi: set() for fi in range(len(faces))}
    adj_opp = {fi: set() for fi in range(len(faces))}
    for c in range(n):
        f0, f1, f2, f3 = (face_of[(c, i)] for i in range(4))
        adj_same[f0].add(f2); adj_same[f2].add(f0)
        adj_same[f1].add(f3); adj_same[f3].add(f1)
        for a, b in ((f0, f1), (f1, f2), (f2, f3), (f3, f0)):
            adj_opp[a].add(b); adj_opp[b].add(a)
    stack = [0]
    colour[0] = 0
    while stack:
        f = stack.pop()
        for gset, rel in ((adj_same[f], 0), (adj_opp[f], 1)):
            for h in gset:
                want = colour[f] ^ rel
                if h in colour:
                    if colour[h] != want:
                        raise ValueError('faces are not checkerboard colourable')
                else:
                    colour[h] = want
                    stack.append(h)
    # the class containing the (c,1)/(c,3) corners of crossing 0 is our "white" class
    white_col = colour[face_of[(0, 1)]]
    white = sorted(fi for fi in range(len(faces)) if colour[fi] == white_col)
    edges, signs = {}, {}
    for c in range(n):
        f0, f1, f2, f3 = (face_of[(c, i)] for i in range(4))
        if colour[f1] == white_col:
            edges[c] = (f1, f3)
            signs[c] = 1
        else:
            edges[c] = (f0, f2)
            signs[c] = -1
    rotation = {}
    for fi in white:
        rotation[fi] = [d[0] for d in faces[fi]]
    g = EmbeddedTaitGraph(white, edges, signs, rotation, knot_name)
    g.validate()
    return g


# =============================================================================
# PD utilities: relabel, splice (connected sum), crossing change
# =============================================================================
def pd_labels(pd) -> set:
    return {int(x) for q in pd for x in q}


def pd_shift(pd, s: int) -> list:
    return [[int(x) + s for x in q] for q in pd]


def pd_normalize_labels(pd) -> list:
    """Relabel arcs to 1..2n (order-preserving)."""
    labs = sorted(pd_labels(pd))
    m = {l: i + 1 for i, l in enumerate(labs)}
    return [[m[int(x)] for x in q] for q in pd]


def pd_connected_sum(pd1, arc1: int, pd2, arc2: int, variant: int = 0) -> list:
    """Splice two knot PDs along the arcs arc1 (of pd1) and arc2 (of pd2): cut both arcs and reconnect.
    For knots the result is a diagram of the connected sum with n1+n2 crossings.  The two variants are the
    two ways of reconnecting (they differ by the relative orientation of the summands)."""
    s = max(pd_labels(pd1)) + 1
    q2 = pd_shift(pd2, s)
    a, b = int(arc1), int(arc2) + s
    occ1 = [(c, i) for c, q in enumerate(pd1) for i, l in enumerate(q) if l == a]
    occ2 = [(c, i) for c, q in enumerate(q2) for i, l in enumerate(q) if l == b]
    if len(occ1) != 2 or len(occ2) != 2:
        raise ValueError('arcs must occur exactly twice')
    out1 = [list(q) for q in pd1]
    out2 = [list(q) for q in q2]
    fresh = max(pd_labels(pd1) | pd_labels(q2)) + 1
    (c1, i1), (c1b, i1b) = occ1
    (c2, i2), (c2b, i2b) = occ2
    if variant % 2 == 0:
        out2[c2][i2] = a            # a runs pd1-side -> pd2-side
        out2[c2b][i2b] = fresh
        out1[c1b][i1b] = fresh      # fresh runs pd2-side -> pd1-side
    else:
        out2[c2b][i2b] = a
        out2[c2][i2] = fresh
        out1[c1b][i1b] = fresh
    return pd_normalize_labels(out1 + out2)


def pd_crossing_change(pd, c: int) -> list:
    out = [list(q) for q in pd]
    a, b, cc, d = out[c]
    out[c] = [b, cc, d, a]
    return out


# =============================================================================
# Strand structure of a PD (shared by jones / passmove)
# =============================================================================
def pd_passages(pd) -> list:
    """Trace the knot: cyclic list of passages (crossing, in_position); in_position odd => OVER passage.
    Each crossing appears exactly twice (once under: in_position 0 or 2; once over: 1 or 3)."""
    occ = {}
    for c, q in enumerate(pd):
        for i, l in enumerate(q):
            occ.setdefault(int(l), []).append((c, i))
    partner = {}
    for l, (d1, d2) in occ.items():
        partner[d1] = d2
        partner[d2] = d1
    # walk: enter crossing c at position i, leave at (i+2)%4, follow that arc to its other occurrence
    start = (0, 0)
    walk, cur = [], start
    while True:
        walk.append(cur)
        c, i = cur
        out = (c, (i + 2) % 4)
        cur = partner[out]
        if cur == start:
            break
    if len(walk) != 2 * len(pd):
        raise ValueError('PD is not a single-component knot diagram')
    return walk


# =============================================================================
# IO
# =============================================================================
def to_dict(g: EmbeddedTaitGraph) -> dict:
    return {'knot_name': g.knot_name,
            'vertices': list(g.vertices),
            'edges': {str(e): list(uv) for e, uv in g.edges.items()},
            'signs': {str(e): s for e, s in g.signs.items()},
            'rotation': {str(v): list(r) for v, r in g.rotation.items()},
            'meta': g.meta}


def from_dict(d: dict) -> EmbeddedTaitGraph:
    return EmbeddedTaitGraph([int(v) for v in d['vertices']],
                             {int(e): tuple(uv) for e, uv in d['edges'].items()},
                             {int(e): int(s) for e, s in d['signs'].items()},
                             {int(v): [int(x) for x in r] for v, r in d['rotation'].items()},
                             d.get('knot_name'), d.get('meta', {}))


def save_json(g: EmbeddedTaitGraph, path: str) -> None:
    with open(path, 'w') as h:
        json.dump(to_dict(g), h, indent=1)


def load_json(path: str) -> EmbeddedTaitGraph:
    with open(path) as h:
        return from_dict(json.load(h))


def relabel(g: EmbeddedTaitGraph, seed=None) -> EmbeddedTaitGraph:
    rnd = _random.Random(seed)
    vs = list(g.vertices)
    vp = vs[:]
    rnd.shuffle(vp)
    vmap = dict(zip(vs, vp))
    es = sorted(g.edges)
    ep = es[:]
    rnd.shuffle(ep)
    emap = dict(zip(es, ep))
    h = EmbeddedTaitGraph(sorted(vmap[v] for v in g.vertices),
                          {emap[e]: (vmap[u], vmap[w]) for e, (u, w) in g.edges.items()},
                          {emap[e]: s for e, s in g.signs.items()}, {}, g.knot_name, dict(g.meta))
    for v in g.vertices:
        r = [emap[e] for e in g.rotation[v]]
        if r:
            k = rnd.randrange(len(r))
            r = r[k:] + r[:k]
        h.rotation[vmap[v]] = r
    return h
