"""Early forward-search prototype over diagrams of a knot.

This module is retained as an experimental utility but is not used for any result reported
in the paper.  The production change-reduce-identify pipeline is implemented in
``code/upper_bounds``.

Goal: for a knot K (start diagram g0 from KnotInfo) and a target k, find a diagram D isotopic to g0 and a
set S of crossings with |S| = k such that
   (A) flipping S gives the unknot  — proved by knot Floer homology (Seifert genus 0), or
   (B) |S| = 1 and flipping gives a knot K' with KnotInfo u(K') = k-1 — candidate by (det, Jones), proved by
       SnapPy isometry with the census knot (HFK detection for 3_1 / 4_1).
Either way u(K) <= k.  The diagram D is reached from g0 by exact isotopy moves (R1±, R2±, R3, flype), so the
witness is certified by (path of moves, crossing set, unknot/identification proof) and can be re-verified
with `verify_witness`.

Search: best-first with a priority = crossings + depth_weight*depth + jitter, a crossing budget
n0 + max_extra, canonical codes as visited set (colour-class invariant), and a per-crossing-count cap to
keep diversity (QD-like).  A learned value function can be plugged in via `config.value_fn`.
"""
from __future__ import annotations

import heapq
import itertools
import json
import math
import random
import time
from dataclasses import dataclass, field, asdict
from typing import Callable, Optional

from . import knotinfo
from .graph import EmbeddedTaitGraph, to_link, to_pd, to_dict, from_dict, to_oriented_pd
from . import moves as mv
from .canonical import canonical_code
from .invariants import goeritz_matrix, jones_regina, hfk, is_unknot


# ----------------------------------------------------------------------------- fast determinant
def _det_bareiss(M: list[list[int]]) -> int:
    n = len(M)
    if n == 0:
        return 1
    A = [row[:] for row in M]
    sign, prev = 1, 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        piv = A[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * piv - A[i][k] * A[k][j]) // prev
        prev = piv
    return sign * A[-1][-1]


def goeritz_list(g: EmbeddedTaitGraph, flipped: set[int] = frozenset()) -> list[list[int]]:
    verts = g.vertices
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    M = [[0] * n for _ in range(n)]
    for e, (u, w) in g.edges.items():
        s = -g.signs[e] if e in flipped else g.signs[e]
        i, j = idx[u], idx[w]
        M[i][i] += s
        M[j][j] += s
        M[i][j] -= s
        M[j][i] -= s
    return [row[:-1] for row in M[:-1]]


def det_after_flips(g: EmbeddedTaitGraph, flipped: set[int]) -> int:
    return abs(_det_bareiss(goeritz_list(g, flipped)))


_P61 = (1 << 61) - 1          # Mersenne prime; |det| of our diagrams is far below p/2, so residues are exact


def _inverse_mod(M: list[list[int]], p: int) -> list[list[int]] | None:
    """Inverse of an integer matrix modulo the prime p (Gauss-Jordan); None if singular mod p."""
    n = len(M)
    A = [[x % p for x in row] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c]), None)
        if piv is None:
            return None
        A[c], A[piv] = A[piv], A[c]
        inv = pow(A[c][c], p - 2, p)
        A[c] = [(x * inv) % p for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]:
                f = A[r][c]
                A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


class FlipDet:
    """Determinants of the diagram after flipping any small set of crossings, in O(1) per query.

    Flipping crossing e = (u, w) with sign s changes the (reduced) Goeritz matrix by c·x·xᵀ with
    c = -2s and x = e_u - e_w (the deleted vertex contributes the zero vector, loops give x = 0).  By the
    matrix determinant lemma det(G + X C Xᵀ) = det(G)·det(I + C·XᵀG⁻¹X), and XᵀG⁻¹X needs only four
    entries of G⁻¹ per pair of flips.  G⁻¹ is computed once modulo the prime 2^61-1; since the true
    determinants are tiny compared with p, the centred residue is the exact signed determinant.
    `det(S)` returns |det| of the diagram with the crossings S flipped (S empty: the knot determinant)."""

    def __init__(self, g: EmbeddedTaitGraph):
        verts = g.vertices
        idx = {v: i for i, v in enumerate(verts)}
        self.m = len(verts) - 1
        M = goeritz_list(g)
        self.d0 = _det_bareiss(M)
        self.p = _P61
        self.inv = _inverse_mod(M, self.p) if (self.m > 0 and self.d0 % self.p) else None
        self.x = {}
        for e, (u, w) in g.edges.items():
            i, j = idx[u], idx[w]
            if i == j:
                continue                                  # loop: flipping a nugatory crossing changes nothing
            self.x[e] = ((i if i < self.m else None), (j if j < self.m else None), (-2 * g.signs[e]) % self.p)
        self._g = g

    def _q(self, a: int, b: int) -> int:
        """x_aᵀ G⁻¹ x_b mod p."""
        ia, ja, _ = self.x[a]
        ib, jb, _ = self.x[b]
        inv = self.inv
        t = 0
        for r, sr in ((ia, 1), (ja, -1)):
            if r is None:
                continue
            row = inv[r]
            for c, sc in ((ib, 1), (jb, -1)):
                if c is None:
                    continue
                t += sr * sc * row[c]
        return t % self.p

    def det(self, S) -> int:
        S = [e for e in S if e in self.x]               # loops do not change the determinant
        if not S or self.m == 0:
            return abs(self.d0)
        if self.inv is None:                              # singular mod p (never for knots): exact fallback
            return det_after_flips(self._g, set(S))
        p = self.p
        k = len(S)
        # A = I + C·Q,  Q_ab = x_aᵀ G⁻¹ x_b
        A = [[((self.x[S[a]][2] * self._q(S[a], S[b])) + (1 if a == b else 0)) % p for b in range(k)] for a in range(k)]
        if k == 1:
            dA = A[0][0]
        elif k == 2:
            dA = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % p
        else:
            dA = 1
            for c in range(k):
                piv = next((r for r in range(c, k) if A[r][c]), None)
                if piv is None:
                    dA = 0
                    break
                if piv != c:
                    A[c], A[piv] = A[piv], A[c]
                    dA = -dA
                dA = (dA * A[c][c]) % p
                inv = pow(A[c][c], p - 2, p)
                for r in range(c + 1, k):
                    if A[r][c]:
                        f = (A[r][c] * inv) % p
                        A[r] = [(x - f * y) % p for x, y in zip(A[r], A[c])]
        d = (self.d0 * dA) % p
        if d > p // 2:
            d -= p
        return abs(d)


# ----------------------------------------------------------------------------- adjacency table
_ADJ_CACHE: dict[int, dict] = {}


def _jones_key(d: dict[int, int]) -> tuple:
    return tuple(sorted(d.items()))


def _mirror_key(key: tuple) -> tuple:
    return tuple(sorted((-e, c) for e, c in key))


def adjacency_table(u_value: int) -> dict[tuple, list[str]]:
    """(det, Jones key) -> KnotInfo names with exact unknotting number u_value (mirror keys included)."""
    if u_value in _ADJ_CACHE:
        return _ADJ_CACHE[u_value]
    table: dict[tuple, list[str]] = {}
    for name, r in knotinfo.rows().items():
        if r['unknotting_number'].strip() != str(u_value) or not r['jones_polynomial']:
            continue
        try:
            D = int(r['determinant'])
        except ValueError:
            continue
        key = _jones_key(knotinfo.parse_jones_string(r['jones_polynomial']))
        try:
            vol = float(r.get('volume') or 0.0)
        except ValueError:
            vol = 0.0
        for kk in (key, _mirror_key(key)):
            table.setdefault((D, kk), [])
            if name not in [n for n, _v in table[(D, kk)]]:
                table[(D, kk)].append((name, vol))
    table['_dets'] = {D for (D, _k) in table if isinstance(D, int)}
    _ADJ_CACHE[u_value] = table
    return table


def _snappy_name_to_knotinfo(name: str) -> str | None:
    import re
    m = re.fullmatch(r'K(\d+)([an])(\d+)', name)
    if m:
        return f'{m.group(1)}{m.group(2)}_{m.group(3)}'
    m = re.fullmatch(r'(\d+)_(\d+)', name)
    if m:
        return name
    return None


def identify_knot(link) -> list[str]:
    """KnotInfo names of census identifications of the knot (SnapPy); [] if none."""
    try:
        ids = link.exterior().identify()
    except Exception:
        return []
    out = []
    for M in ids:
        nm = _snappy_name_to_knotinfo(M.name().split('(')[0])
        if nm and nm in knotinfo.rows():
            out.append(nm)
    return out


def confirm_identification(link, candidates: list[str]) -> str | None:
    """Return the candidate name proved by SnapPy isometry (or HFK detection for 3_1/4_1), else None."""
    import snappy
    H = None
    for nm in candidates:
        if nm in ('3_1', '4_1'):   # HFK detects the trefoil and figure-eight
            H = H or hfk(link)
            ref = hfk(to_link_from_name(nm))
            if H['ranks'] == ref['ranks'] and H['tau'] == ref['tau'] and H['seifert_genus'] == ref['seifert_genus']:
                return nm
            continue
    names = identify_knot(link)
    for nm in candidates:
        if nm in names:
            return nm
    # fall back: explicit isometry test with the census manifold
    try:
        M = link.exterior()
        for nm in candidates:
            r = knotinfo.row(nm)
            cen = 'K' + nm.replace('_', '') if int(r['crossing_number']) >= 11 else nm
            try:
                N = snappy.Manifold(cen)
                if M.is_isometric_to(N):
                    return nm
            except Exception:
                continue
    except Exception:
        pass
    return None


def to_link_from_name(nm: str):
    from spherogram import Link
    return Link(knotinfo.pd_code(nm))


# ----------------------------------------------------------------------------- witnesses
@dataclass
class Witness:
    target: str
    k: int
    kind: str                       # 'unknot' or 'adjacent'
    diagram: dict                   # to_dict(EmbeddedTaitGraph)
    crossings: list[int]
    path: list                      # [(kind, params), ...] from the start diagram
    adjacent_knot: str | None = None
    adjacent_u: int | None = None
    checks: dict = field(default_factory=dict)

    def save(self, path: str) -> None:
        with open(path, 'w') as h:
            json.dump(asdict(self), h, indent=1)


def verify_witness(g0: EmbeddedTaitGraph, w: Witness) -> bool:
    """Replay the move path from g0, compare canonical codes, and re-check the unknot/identification."""
    cur = g0
    for kind, params in w.path:
        cur = mv.apply_move(cur, kind, tuple(params))
    D = from_dict(w.diagram)
    if canonical_code(cur) != canonical_code(D):
        return False
    h = D
    for e in w.crossings:
        h = mv.crossing_change(h, e)
    if w.kind == 'unknot':
        return is_unknot(h)
    L = to_link(h)
    nm = confirm_identification(L, [w.adjacent_knot])
    return nm == w.adjacent_knot and knotinfo.unknotting_interval(nm) == (w.k - 1, w.k - 1)


# ----------------------------------------------------------------------------- config / search
@dataclass
class SearchConfig:
    k: int = 2                          # target upper bound
    max_extra: int = 4                  # crossing budget above the start diagram
    max_nodes: int = 5000               # expanded nodes
    time_limit: float = 600.0           # seconds
    depth_weight: float = 0.3
    jitter: float = 0.5
    per_level_cap: int = 2000           # max expanded nodes per crossing count
    increasing_samples: int = 12        # sampled R1+/R2+ moves per node
    use_unknot_goal: bool = True        # k-subset flips -> unknot (HFK)
    use_adjacency_goal: bool = True     # single flip -> knot with u = k-1
    seed: int = 0
    value_fn: Optional[Callable[[EmbeddedTaitGraph], float]] = None   # lower is better (added to priority)
    log_every: int = 500
    dedup: str = 'sig'                   # 'sig' (Regina canonical signature) or 'canonical' (tait.canonical)
    r1_add: Optional[bool] = None       # R1+ moves; None = only when max_extra >= 3 (R2/R3/flypes preserve the
                                        # writhe, so a kink is only useful as an intermediate of an R1+R2 twist)
    skip_kinky: bool = True             # no goal tests on diagrams with an R1-removable crossing (their reduced
                                        # diagram has the same flip behaviour and is expanded first)
    kink_penalty: float = 2.0           # priority penalty per R1-removable crossing: a kink is only a stepping stone
                                        # to an R1+R2 twist two crossings higher, so it is ordered like that level
    heap_max: int = 2_000_000           # prune the lazy heap to half this size when it grows beyond it (long runs
                                        # would otherwise exhaust memory: ~30 pushes per expansion, ~250 B each)


def _goal_unknot(g: EmbeddedTaitGraph, k: int, fd: "FlipDet | None" = None) -> list[int] | None:
    """Any k-subset of crossings whose flip gives the unknot.  Prefilter: determinant after the flips must be 1
    (FlipDet, O(1) per subset); survivors: Jones polynomial == 1 in Regina (crossing changes applied to a single
    oriented-PD link), then HFK (Seifert genus 0) as the proof."""
    import regina
    edges = sorted(g.edges)
    if k > len(edges) or not edges:
        return None
    fd = fd or FlipDet(g)
    R = None
    for S in itertools.combinations(range(len(edges)), k):
        if fd.det([edges[i] for i in S]) != 1:
            continue
        if R is None:
            R = _regina_link(g)
        R2 = regina.Link(R)
        for i in S:
            R2.change(R2.crossing(i))
        J = R2.jones()
        if J.minExp() == 0 and J.maxExp() == 0 and str(J[0]) == '1':
            h = g
            for i in S:
                h = mv.crossing_change(h, edges[i])
            if is_unknot(h):
                return [edges[i] for i in S]
    return None


def _regina_link(g: EmbeddedTaitGraph):
    import regina
    return regina.Link.fromPD(to_oriented_pd(g))


def _jones_dict_regina(J) -> dict[int, int]:
    out = {}
    for e in range(J.minExp(), J.maxExp() + 1):
        c = int(str(J[e]))
        if c:
            out[e // 2] = c
    return out


def _goal_adjacent(g: EmbeddedTaitGraph, k: int, table: dict, fd: "FlipDet | None" = None) -> tuple[list[int], str] | None:
    """Single flip giving a knot with u = k-1: determinant prefilter (FlipDet), Jones via Regina (crossing change
    done in Regina on the oriented PD), table lookup by (det, Jones), then SnapPy/HFK confirmation."""
    import regina
    edges = sorted(g.edges)
    dets = table.get('_dets', set())
    fd = fd or FlipDet(g)
    R = None
    for i, e in enumerate(edges):
        D = fd.det([e])
        if D not in dets:
            continue
        if R is None:
            R = _regina_link(g)
        R2 = regina.Link(R)
        R2.change(R2.crossing(i))
        key = (D, _jones_key(_jones_dict_regina(R2.jones())))
        cands = table.get(key)
        if not cands:
            continue
        # Jones collisions are common (e.g. 10_34 vs 12n_468): filter by hyperbolic volume before SnapPy identify
        h = mv.crossing_change(g, e)
        L = to_link(h)
        try:
            vol = float(L.exterior().volume())
        except Exception:
            vol = 0.0
        names = [nm for nm, v in cands if (v < 0.5 and vol < 0.5) or abs(v - vol) < 1e-5]
        if not names:
            continue
        nm = confirm_identification(L, names)
        if nm is not None:
            return [e], nm
    return None


def diagram_key(g: EmbeddedTaitGraph, method: str = 'sig'):
    """De-duplication key of a diagram: Regina's canonical signature ('sig', fast, identifies a diagram with
    its 180-degree rotation and orientation reversal, not with its mirror) or our canonical_code."""
    if method == 'sig':
        if not g.edges:
            return 'unknot'
        return _regina_link(g).sig(False, True)
    return canonical_code(g)


def _expansions(g: EmbeddedTaitGraph, budget: int, rnd: random.Random, n_inc: int,
                allow_r1: bool = True) -> list[tuple[str, tuple]]:
    out: list[tuple[str, tuple]] = []
    for kind, e in mv.find_r1_removals(g):
        out.append(('r1-', (e,)))
    for kind, p in mv.find_r2_removals(g):
        out.append(('r2-' + kind[0], p))
    for fi in mv.find_r3_triangles(g):
        out.append(('r3Δ', (fi,)))
    for t in mv.find_r3_stars(g):
        out.append(('r3Y', (t,)))
    for e, pi in mv.find_flypes(g):
        out.append(('flype', (e, pi)))
    n = g.n_crossings()
    if n + 2 <= budget:
        par = mv.find_r2_parallel_additions(g)
        if not allow_r1:
            from .graph import ribbon_faces
            faces = ribbon_faces(g)
            par = [(fi, i, j) for (fi, i, j) in par if faces[fi][i][0] != faces[fi][j][0]]   # same-vertex pair = kink pair
        ser = mv.find_r2_series_additions(g)
        if n_inc <= 0:                                    # exhaustive: every R2+ addition, both signs
            for (fi, i, j) in par:
                out.append(('r2+p', (fi, i, j, 1)))
                out.append(('r2+p', (fi, i, j, -1)))
            for (v, i, j) in ser:
                out.append(('r2+s', (v, i, j, 1)))
                out.append(('r2+s', (v, i, j, -1)))
        else:
            for _ in range(n_inc):
                if par and rnd.random() < 0.6:
                    fi, i, j = rnd.choice(par)
                    out.append(('r2+p', (fi, i, j, rnd.choice((1, -1)))))
                elif ser:
                    v, i, j = rnd.choice(ser)
                    out.append(('r2+s', (v, i, j, rnd.choice((1, -1)))))
    if n + 3 <= budget:
        # composite R1+R2 twists: writhe +-1 without kinky intermediate states (see mv.twist); the k-th
        # qualifying R2-over-the-kink candidate is resolved lazily inside apply_move.
        placements = [(v, pos) for v in g.vertices for pos in range(max(1, g.degree(v)))]
        if n_inc > 0 and len(placements) > n_inc:
            placements = rnd.sample(placements, n_inc)
        for (v, pos) in placements:
            for r1kind in ('loop', 'pend'):
                for s1 in (1, -1):
                    for k in (0, 1):
                        for s2 in (1, -1):
                            out.append(('twist', (r1kind, v, pos, s1, k, s2)))
    if allow_r1 and n + 1 <= budget:
        for _ in range(max(1, (n_inc if n_inc > 0 else 12) // 4)):
            v = rnd.choice(g.vertices)
            pos = rnd.randrange(max(1, g.degree(v)))
            out.append((rnd.choice(('r1+loop', 'r1+pend')), (v, pos, rnd.choice((1, -1)))))
    return out


_MOVE_DELTA = {'r1-': -1, 'r2-p': -2, 'r2-s': -2, 'r3Δ': 0, 'r3Y': 0, 'flype': 0, 'r2+p': 2, 'r2+s': 2,
               'r1+loop': 1, 'r1+pend': 1, 'twist': 3}
_KINK_DELTA = {'r1-': -1, 'r1+loop': 1, 'r1+pend': 1, 'twist': 0}     # estimated change of the number of R1-removable crossings


def _flatten_path(node) -> list:
    out = []
    while node is not None:
        node, step = node
        out.append(step)
    out.reverse()
    return out


def search_upper_bound(g0: EmbeddedTaitGraph, target: str, config: SearchConfig,
                       log: Callable[[str], None] | None = print) -> tuple[Witness | None, dict]:
    """Best-first search for a u(K) <= k witness starting from diagram g0 of knot `target`.

    Children are pushed lazily as (parent, move): the child diagram is only built, de-duplicated (canonical
    key) and goal-tested when it is popped — most generated children are never expanded, so this saves the
    dominant move-application + key cost and keeps memory small (the heap holds moves, not diagrams).  With a
    value function the children must be built eagerly (value_fn needs the diagram)."""
    rnd = random.Random(config.seed)
    t0 = time.time()
    n0 = g0.n_crossings()
    budget = n0 + config.max_extra
    table = adjacency_table(config.k - 1) if (config.use_adjacency_goal and config.k >= 1) else {}
    allow_r1 = config.r1_add if config.r1_add is not None else (config.max_extra >= 3)
    lazy = config.value_fn is None
    stats = {'expanded': 0, 'generated': 0, 'dedup_hits': 0, 'by_crossings': {}, 'goal_tests': 0, 'skipped_kinky': 0,
             'heap_prunes': 0}
    visited = {diagram_key(g0, config.dedup)}
    counter = itertools.count()
    heap: list = []

    def priority(n, depth, kinks, g=None):
        p = n + config.kink_penalty * max(0, kinks) + config.depth_weight * depth + config.jitter * rnd.random()
        if config.value_fn is not None and g is not None:
            p += config.value_fn(g)
        return p

    k0 = len(mv.find_r1_removals(g0))
    # heap entries: (priority, counter, parent_graph_or_child, kind, params, depth, path_node, est_kinks)
    heapq.heappush(heap, (priority(n0, 0, k0, g0), next(counter), g0, None, None, 0, None, k0))
    level_count: dict[int, int] = {}
    while heap:
        if stats['expanded'] >= config.max_nodes or time.time() - t0 > config.time_limit:
            break
        _, _, gp, kind, params, depth, pnode, _ek = heapq.heappop(heap)
        if kind is None:
            g = gp                                        # root, or eagerly built child
        else:
            try:
                g = mv.apply_move(gp, kind, params)
            except ValueError:
                continue
            code = diagram_key(g, config.dedup)
            if code in visited:
                stats['dedup_hits'] += 1
                continue
            visited.add(code)
        n = g.n_crossings()
        if level_count.get(n, 0) >= config.per_level_cap:
            continue
        level_count[n] = level_count.get(n, 0) + 1
        stats['expanded'] += 1
        stats['by_crossings'][n] = stats['by_crossings'].get(n, 0) + 1
        # ---- goal tests (skipped on diagrams with an R1-removable crossing: same flip behaviour as the reduced
        # diagram, which has fewer crossings and is expanded before this one)
        n_kinks = len(mv.find_r1_removals(g))
        kinky = config.skip_kinky and n_kinks > 0
        if kinky:
            stats['skipped_kinky'] += 1
        else:
            stats['goal_tests'] += 1
        fd = None if kinky else FlipDet(g)
        if config.use_unknot_goal and not kinky:
            S = _goal_unknot(g, config.k, fd)
            if S is not None:
                w = Witness(target, config.k, 'unknot', to_dict(g), S, _flatten_path(pnode),
                            checks={'hfk_genus_0': True, 'crossings_in_diagram': n, 'depth': depth,
                                    'expanded': stats['expanded'], 'seconds': round(time.time() - t0, 1)})
                stats['visited'] = len(visited)
                return w, stats
        if config.use_adjacency_goal and table and not kinky:
            r = _goal_adjacent(g, config.k, table, fd)
            if r is not None:
                S, nm = r
                w = Witness(target, config.k, 'adjacent', to_dict(g), S, _flatten_path(pnode), adjacent_knot=nm,
                            adjacent_u=config.k - 1,
                            checks={'identified_by': 'snappy_isometry_or_hfk', 'crossings_in_diagram': n,
                                    'depth': depth, 'expanded': stats['expanded'], 'seconds': round(time.time() - t0, 1)})
                stats['visited'] = len(visited)
                return w, stats
        # ---- expand
        for kind2, params2 in _expansions(g, budget, rnd, config.increasing_samples, allow_r1):
            stats['generated'] += 1
            step = (kind2, list(params2))
            ek = n_kinks + _KINK_DELTA.get(kind2, 0)
            if lazy:
                heapq.heappush(heap, (priority(n + _MOVE_DELTA[kind2], depth + 1, ek), next(counter), g, kind2, params2,
                                      depth + 1, (pnode, step), ek))
            else:
                try:
                    h = mv.apply_move(g, kind2, params2)
                except ValueError:
                    continue
                code = diagram_key(h, config.dedup)
                if code in visited:
                    stats['dedup_hits'] += 1
                    continue
                visited.add(code)
                kh = len(mv.find_r1_removals(h))
                heapq.heappush(heap, (priority(h.n_crossings(), depth + 1, kh, h), next(counter), h, None, None,
                                      depth + 1, (pnode, step), kh))
        if len(heap) > config.heap_max:
            heap = sorted(heap)[:config.heap_max // 2]      # a sorted list is a valid heap
            stats['heap_prunes'] += 1
        if log and stats['expanded'] % config.log_every == 0:
            log(f'  [{target}] expanded {stats["expanded"]} visited {len(visited)} heap {len(heap)} '
                f'levels {dict(sorted(stats["by_crossings"].items()))} {time.time() - t0:.0f}s')
    stats['visited'] = len(visited)
    stats['seconds'] = round(time.time() - t0, 1)
    return None, stats


# ----------------------------------------------------------------------------- random inflate/simplify walks
def _random_step(g: EmbeddedTaitGraph, rnd: random.Random, budget: int, p_increase: float, allow_r1: bool = True):
    """One random isotopy move; with probability p_increase prefer crossing-increasing moves (if budget allows)."""
    n = g.n_crossings()
    inc, other = [], []
    if n + 2 <= budget:
        par = mv.find_r2_parallel_additions(g)
        ser = mv.find_r2_series_additions(g)
        for _ in range(3):
            if par:
                fi, i, j = rnd.choice(par)
                inc.append(('r2+p', (fi, i, j, rnd.choice((1, -1)))))
            if ser:
                v, i, j = rnd.choice(ser)
                inc.append(('r2+s', (v, i, j, rnd.choice((1, -1)))))
    if allow_r1 and n + 1 <= budget and rnd.random() < 0.25:    # kinks rarely help: keep them rare in walks
        v = rnd.choice(g.vertices)
        inc.append((rnd.choice(('r1+loop', 'r1+pend')), (v, rnd.randrange(max(1, g.degree(v))), rnd.choice((1, -1)))))
    for fi in mv.find_r3_triangles(g):
        other.append(('r3Δ', (fi,)))
    for t in mv.find_r3_stars(g):
        other.append(('r3Y', (t,)))
    for e, pi in mv.find_flypes(g):
        other.append(('flype', (e, pi)))
    for kind, prm in mv.find_r2_removals(g):
        other.append(('r2-' + kind[0], prm))
    for kind, e in mv.find_r1_removals(g):
        other.append(('r1-', (e,)))
    pool = inc if (inc and rnd.random() < p_increase) else (other or inc)
    if not pool:
        return g
    kind, params = rnd.choice(pool)
    try:
        return mv.apply_move(g, kind, params)
    except ValueError:
        return g


def search_random_walks(g0: EmbeddedTaitGraph, target: str, config: SearchConfig, n_walks: int = 200,
                        walk_len: int = 30, p_increase: float = 0.6,
                        log: Callable[[str], None] | None = None) -> tuple[Witness | None, dict]:
    """Sampling strategy (Brittenham–Hermiller / DKT style): from g0 do random walks that mostly inflate the
    diagram up to the crossing budget, then R3/flype shuffles and simplifications; every distinct diagram on the
    way is goal-tested.  Complements the best-first search, which stays near minimal diagrams."""
    rnd = random.Random(config.seed + 12345)
    t0 = time.time()
    budget = g0.n_crossings() + config.max_extra
    table = adjacency_table(config.k - 1) if config.use_adjacency_goal else {}
    allow_r1 = config.r1_add if config.r1_add is not None else (config.max_extra >= 3)
    seen = set()
    stats = {'walks': 0, 'tested': 0, 'by_crossings': {}, 'skipped_kinky': 0}
    for wk in range(n_walks):
        if time.time() - t0 > config.time_limit:
            break
        stats['walks'] += 1
        g = g0
        for step in range(walk_len):
            p_inc = p_increase if step < walk_len * 0.6 else 0.15   # inflate first, then shuffle / simplify
            h = _random_step(g, rnd, budget, p_inc, allow_r1)
            if h is g:
                continue
            g = h
            key = diagram_key(g, config.dedup)
            if key in seen:
                continue
            seen.add(key)
            if config.skip_kinky and mv.find_r1_removals(g):
                stats['skipped_kinky'] += 1
                continue
            stats['tested'] += 1
            n = g.n_crossings()
            stats['by_crossings'][n] = stats['by_crossings'].get(n, 0) + 1
            fd = FlipDet(g)
            S = _goal_unknot(g, config.k, fd) if config.use_unknot_goal else None
            if S is not None:
                return Witness(target, config.k, 'unknot', to_dict(g), S, [],
                               checks={'hfk_genus_0': True, 'crossings_in_diagram': n, 'walk': wk, 'step': step,
                                       'path_recorded': False, 'seconds': round(time.time() - t0, 1)}), stats
            if table:
                r = _goal_adjacent(g, config.k, table, fd)
                if r is not None:
                    S, nm = r
                    return Witness(target, config.k, 'adjacent', to_dict(g), S, [], adjacent_knot=nm, adjacent_u=config.k - 1,
                                   checks={'identified_by': 'snappy_isometry_or_hfk', 'crossings_in_diagram': n, 'walk': wk,
                                           'step': step, 'path_recorded': False, 'seconds': round(time.time() - t0, 1)}), stats
        if log and wk % 50 == 0:
            log(f'  [{target}] walk {wk} tested {stats["tested"]} {time.time() - t0:.0f}s')
    stats['seconds'] = round(time.time() - t0, 1)
    stats['visited'] = len(seen)
    return None, stats


def verify_witness_light(g0: EmbeddedTaitGraph, w: Witness) -> bool:
    """Verification when the move path was not recorded: the witness diagram must be a diagram of the same knot
    (SnapPy isometry of exteriors; fallback HFK+Jones+det agreement), then the flip proof is re-checked."""
    D = from_dict(w.diagram)
    same = False
    try:
        same = bool(to_link(g0).exterior().is_isometric_to(to_link(D).exterior()))
    except Exception:
        same = False
    if not same:
        from .invariants import invariant_fingerprint
        same = invariant_fingerprint(g0) == invariant_fingerprint(D)
    if not same:
        return False
    h = D
    for e in w.crossings:
        h = mv.crossing_change(h, e)
    if w.kind == 'unknot':
        return is_unknot(h)
    nm = confirm_identification(to_link(h), [w.adjacent_knot])
    return nm == w.adjacent_knot and knotinfo.unknotting_interval(nm) == (w.k - 1, w.k - 1)
