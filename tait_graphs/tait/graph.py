"""Embedding-preserving Tait graphs (signed planar maps) for knot diagrams.

An abstract signed multigraph determines a knot only up to mutation, and for non-3-connected
graphs different planar embeddings give different knots.  The *rotation system* (cyclic order of
edge-ends around each vertex) together with the crossing signs determines the diagram on S^2 exactly.

Data model
----------
EmbeddedTaitGraph
    vertices : list[int]                  faces of one checkerboard colour class ("white" faces)
    edges    : dict[int, (int, int)]      crossing id -> (u, v)   (loops u == v allowed)
    signs    : dict[int, int]             crossing id -> +1 / -1   (Gordon–Litherland convention of
                                          spherogram.white_graph: +1 iff the white corners are PD
                                          corners 0 and 2, i.e. are swept going counterclockwise from
                                          the under-strand to the over-strand)
    rotation : dict[int, list[int]]       vertex -> cyclic list of incident crossing ids in the order
                                          spherogram.faces() returns corners (clockwise); a loop
                                          appears twice.

All rotation-system conventions inside this package are self-consistent; the only external contract
is: from_link(L) followed by to_link() reproduces the same diagram (verified on all KnotInfo knots
<= 12 crossings with Regina's canonical diagram signature, including chirality).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import networkx as nx
from spherogram import Link
from spherogram.links.links_base import CrossingStrand

Dart = tuple[int, int, int]  # (vertex, edge, slot)  slot = 0, or 1 for the second end of a loop


@dataclass
class EmbeddedTaitGraph:
    vertices: list[int]
    edges: dict[int, tuple[int, int]]
    signs: dict[int, int]
    rotation: dict[int, list[int]]
    knot_name: str | None = None
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ basics
    def copy(self) -> "EmbeddedTaitGraph":
        return EmbeddedTaitGraph(list(self.vertices), dict(self.edges), dict(self.signs),
                                 {v: list(r) for v, r in self.rotation.items()}, self.knot_name, dict(self.meta))

    def n_crossings(self) -> int:
        return len(self.edges)

    def degree(self, v: int) -> int:
        return len(self.rotation[v])

    def other_end(self, e: int, v: int) -> int:
        u, w = self.edges[e]
        return w if v == u else u

    def is_loop(self, e: int) -> bool:
        u, w = self.edges[e]
        return u == w

    def new_vertex_id(self) -> int:
        return (max(self.vertices) + 1) if self.vertices else 0

    def new_edge_id(self) -> int:
        return (max(self.edges) + 1) if self.edges else 0

    def incident(self, v: int) -> list[int]:
        return list(self.rotation[v])

    def darts(self) -> list[Dart]:
        out = []
        for v in self.vertices:
            seen: dict[int, int] = {}
            for e in self.rotation[v]:
                k = seen.get(e, 0)
                seen[e] = k + 1
                out.append((v, e, k))
        return out

    def dart_list(self, v: int) -> list[Dart]:
        seen: dict[int, int] = {}
        out = []
        for e in self.rotation[v]:
            k = seen.get(e, 0)
            seen[e] = k + 1
            out.append((v, e, k))
        return out

    def euler_characteristic(self) -> int:
        return len(self.vertices) - len(self.edges) + len(ribbon_faces(self))

    def is_planar(self) -> bool:
        """The rotation system describes a map on the sphere iff V - E + F == 2 (connected graph)."""
        return self.is_connected() and self.euler_characteristic() == 2

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
        """Raise ValueError if the structure is inconsistent (incidences, planarity)."""
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
            if e not in self.signs or self.signs[e] not in (1, -1):
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
        """Mirror image of the diagram: reverse all rotations and flip all signs."""
        h = self.copy()
        for v in h.vertices:
            h.rotation[v] = list(reversed(h.rotation[v]))
        for e in h.signs:
            h.signs[e] = -h.signs[e]
        return h

    def writhe_like_sum(self) -> int:
        return sum(self.signs.values())


# ----------------------------------------------------------------------------
# Ribbon-graph faces (rotation system -> faces as dart cycles)
# ----------------------------------------------------------------------------
def _rho_theta(g: EmbeddedTaitGraph):
    rho: dict[Dart, Dart] = {}
    for v in g.vertices:
        dv = g.dart_list(v)
        for i, d in enumerate(dv):
            rho[d] = dv[(i + 1) % len(dv)]
    theta: dict[Dart, Dart] = {}
    for e, (u, w) in g.edges.items():
        if u != w:
            theta[(u, e, 0)] = (w, e, 0)
            theta[(w, e, 0)] = (u, e, 0)
        else:
            theta[(u, e, 0)] = (u, e, 1)
            theta[(u, e, 1)] = (u, e, 0)
    return rho, theta


def ribbon_faces(g: EmbeddedTaitGraph) -> list[list[Dart]]:
    """Faces of the embedded graph as cycles of darts; phi = rho o theta."""
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


def face_of_dart(g: EmbeddedTaitGraph) -> dict[Dart, int]:
    out = {}
    for i, f in enumerate(ribbon_faces(g)):
        for d in f:
            out[d] = i
    return out


# ----------------------------------------------------------------------------
# Link -> embedded Tait graph
# ----------------------------------------------------------------------------
def from_link(link: Link, knot_name: str | None = None) -> EmbeddedTaitGraph:
    """Checkerboard graph of one colour class with the rotation system read off link.faces().
    Vertex ids = indices into link.faces(); crossing ids = indices into link.crossings."""
    faces = link.faces()
    face_of = {corner: n for n, face in enumerate(faces) for corner in face}
    cidx = {c: i for i, c in enumerate(link.crossings)}

    G = nx.MultiGraph()
    for c in link.crossings:
        G.add_edge(face_of[CrossingStrand(c, 0)], face_of[CrossingStrand(c, 2)], crossing=cidx[c], sign=1)
        G.add_edge(face_of[CrossingStrand(c, 1)], face_of[CrossingStrand(c, 3)], crossing=cidx[c], sign=-1)
    comps = sorted(nx.connected_components(G), key=lambda s: sorted(s))
    if len(comps) > 2:
        raise ValueError('split diagram')
    white = set(comps[1]) if len(comps) == 2 else set(comps[0])

    edges, signs = {}, {}
    for u, v, d in G.edges(data=True):
        if u in white and v in white:
            edges[d['crossing']] = (u, v)
            signs[d['crossing']] = d['sign']
    if len(edges) != len(link.crossings):
        raise ValueError('every crossing must give exactly one white edge')
    rotation = {fidx: [cidx[cs.crossing] for cs in faces[fidx]] for fidx in sorted(white)}
    g = EmbeddedTaitGraph(sorted(white), edges, signs, rotation, knot_name)
    g.validate()
    return g


def from_pd(pd, knot_name: str | None = None) -> EmbeddedTaitGraph:
    return from_link(Link([list(q) for q in pd]), knot_name)


# ----------------------------------------------------------------------------
# Embedded Tait graph -> PD code (medial graph) -> Link
# ----------------------------------------------------------------------------
def to_pd(g: EmbeddedTaitGraph) -> list[list[int]]:
    """PD code realising the signed map.

    The medial quadruple at edge e=(u,v) is [arc_u_in, arc_u_out, arc_v_in, arc_v_out], listed
    counterclockwise (rotation lists are clockwise around faces).  In PD convention strands 0->2 are
    the under-strand and corner j lies between strands j, j+1, so corners 0 and 2 are exactly the
    two faces u, v of our colour class.  The white-graph sign is +1 when the white edge joins
    corners 0 and 2, hence: sign +1 -> keep the quadruple, sign -1 -> rotate it by one position."""
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


def to_link(g: EmbeddedTaitGraph) -> Link:
    if not g.edges:
        # The closure of sigma_1 is a one-crossing unknot. The cancelling
        # two-letter word [1,-1] instead closes to the TWO-component unlink.
        return Link(braid_closure=[1])
    return Link(to_pd(g))


# ----------------------------------------------------------------------------
# Export / import
# ----------------------------------------------------------------------------
def to_networkx(g: EmbeddedTaitGraph) -> nx.MultiGraph:
    """MultiGraph with edge attrs crossing/sign and node attr 'rotation' (JSON list of crossing ids)."""
    G = nx.MultiGraph(knot_name=g.knot_name or '')
    for v in g.vertices:
        G.add_node(v, rotation=json.dumps(g.rotation[v]))
    for e, (u, w) in sorted(g.edges.items()):
        G.add_edge(u, w, key=e, crossing=e, sign=g.signs[e])
    return G


def to_dict(g: EmbeddedTaitGraph) -> dict:
    return {
        'knot_name': g.knot_name,
        'vertices': g.vertices,
        'edges': {str(e): list(uv) for e, uv in g.edges.items()},
        'signs': {str(e): s for e, s in g.signs.items()},
        'rotation': {str(v): r for v, r in g.rotation.items()},
        'meta': g.meta,
    }


def from_dict(d: dict) -> EmbeddedTaitGraph:
    return EmbeddedTaitGraph(
        vertices=[int(v) for v in d['vertices']],
        edges={int(e): tuple(uv) for e, uv in d['edges'].items()},
        signs={int(e): int(s) for e, s in d['signs'].items()},
        rotation={int(v): [int(x) for x in r] for v, r in d['rotation'].items()},
        knot_name=d.get('knot_name'),
        meta=d.get('meta', {}),
    )


def save_json(g: EmbeddedTaitGraph, path: str) -> None:
    with open(path, 'w') as h:
        json.dump(to_dict(g), h, indent=1)


def load_json(path: str) -> EmbeddedTaitGraph:
    with open(path) as h:
        return from_dict(json.load(h))


def relabel(g: EmbeddedTaitGraph, vmap: dict[int, int] | None = None, emap: dict[int, int] | None = None,
            rotate_rotations: bool = False, seed: int | None = None) -> EmbeddedTaitGraph:
    """Relabel vertices/edges (and optionally cyclically shift rotation lists) — same diagram."""
    import random
    rnd = random.Random(seed)
    if vmap is None:
        vs = list(g.vertices)
        perm = vs[:]
        rnd.shuffle(perm)
        vmap = dict(zip(vs, perm))
    if emap is None:
        es = sorted(g.edges)
        perm = es[:]
        rnd.shuffle(perm)
        emap = dict(zip(es, perm))
    h = EmbeddedTaitGraph(sorted(vmap[v] for v in g.vertices),
                          {emap[e]: (vmap[u], vmap[w]) for e, (u, w) in g.edges.items()},
                          {emap[e]: s for e, s in g.signs.items()},
                          {}, g.knot_name, dict(g.meta))
    for v in g.vertices:
        r = [emap[e] for e in g.rotation[v]]
        if rotate_rotations and r:
            k = rnd.randrange(len(r))
            r = r[k:] + r[:k]
        h.rotation[vmap[v]] = r
    return h


def dual(g: EmbeddedTaitGraph, reverse: bool = True) -> EmbeddedTaitGraph:
    """The other checkerboard colour class of the same diagram: vertices = faces of g (ribbon faces),
    edges = the same crossings with negated signs, rotation at a face = the crossings around the face in
    the *reversed* face-cycle order (reverse=True, verified against spherogram's other colour class)."""
    faces = ribbon_faces(g)
    vertices = list(range(len(faces)))
    edge_faces: dict[int, list[int]] = {e: [] for e in g.edges}
    rotation = {}
    for fi, f in enumerate(faces):
        seq = [d[1] for d in f]
        if reverse:
            seq = list(reversed(seq))
        rotation[fi] = seq
        for e in seq:
            edge_faces[e].append(fi)
    edges = {}
    for e, fl in edge_faces.items():
        if len(fl) != 2:
            raise ValueError(f'edge {e} is incident to {len(fl)} face-slots (expected 2)')
        edges[e] = (fl[0], fl[1])
    h = EmbeddedTaitGraph(vertices, edges, {e: -s for e, s in g.signs.items()}, rotation, g.knot_name, dict(g.meta))
    h.validate()
    return h


def to_oriented_pd(g: EmbeddedTaitGraph) -> list[list[int]]:
    """PD code in the KnotTheory convention (strand labels increase along the orientation, position 0 of
    every crossing is the incoming under-strand), directly from the map — no spherogram needed.  Row i of
    the result is the crossing sorted(g.edges)[i] (same order as to_pd)."""
    pd = to_pd(g)
    n = len(pd)
    if n == 0:
        return []
    # arc label -> list of (row, position)
    where: dict[int, list[tuple[int, int]]] = {}
    for r, q in enumerate(pd):
        for pos, a in enumerate(q):
            where.setdefault(a, []).append((r, pos))
    # traverse the (single) knot component: start at row 0 position 0, leave through position 2, ...
    new_label: dict[tuple[int, int], int] = {}   # (row, pos) of the *head* end of an arc -> new label
    head_of_arc: dict[int, tuple[int, int]] = {}
    label = 1
    r, pos = 0, 0
    # we enter row 0 through position 0 (so that arc at pos 0 is incoming): choose the arc pd[0][0] as arc #1
    # whose head is (0,0); its tail is the other occurrence.
    start = (0, 0)
    cur_head = start
    visited_heads = set()
    order: list[tuple[int, int, int]] = []   # (row, entry_pos, new label of the arc entering at (row, entry_pos))
    while cur_head not in visited_heads:
        visited_heads.add(cur_head)
        rr, pp = cur_head
        a = pd[rr][pp]
        head_of_arc[a] = cur_head
        new_label[a] = label
        label += 1
        # exit the crossing through the opposite position, follow that arc to its other end
        exit_pos = (pp + 2) % 4
        b = pd[rr][exit_pos]
        occ = [x for x in where[b] if x != (rr, exit_pos)]
        if not occ:   # arc b both ends at the same crossing: the other occurrence is the other position
            occ = [x for x in where[b] if x != (rr, exit_pos)]
        nxt = occ[0] if occ else None
        if nxt is None:
            raise ValueError('could not follow strand')
        cur_head = nxt
    if len(new_label) != 2 * n:
        raise ValueError('not a single-component knot diagram')
    out = []
    for r, q in enumerate(pd):
        a, b, c, d = q
        # under strand is positions 0 -> 2 or 2 -> 0; position 0 must be the incoming (head) end
        if head_of_arc[a] == (r, 0):
            row = [new_label[a], new_label[b], new_label[c], new_label[d]]
        elif head_of_arc[c] == (r, 2):
            row = [new_label[c], new_label[d], new_label[a], new_label[b]]
        else:
            raise ValueError('under-strand orientation inconsistent')
        out.append(row)
    return out
