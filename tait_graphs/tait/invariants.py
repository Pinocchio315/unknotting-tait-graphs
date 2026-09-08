"""Invariants and obstructions computed from embedded Tait graphs / PD codes (no Sage required).

  goeritz_matrix(g)          reduced Goeritz matrix (signed Laplacian of the Tait graph, one row/col deleted)
  determinant(g)
  h1_torsion(g)              invariant factors of H1(Sigma_2 K) = coker(Goeritz)
  linking_form_class(g)      a with lambda(gen,gen) = a/D (mod 1) when H1 is cyclic
  lickorish_allowed(g)       Lickorish 1985: u=1 => linking form ≅ <±2/D>; returns allowed surgery signs
  jones_regina(g|pd|link)    Jones polynomial {exp: coeff} via Regina (exact)
  casson_walker_allowed(...) Mullins + Boyer–Lines integrality test for u=1
  hfk(link) / is_unknot(link)  knot Floer homology (Ozsváth–Szabó program); genus 0 <=> unknot
  sigma2_snappy(link)        double branched cover via the order-two meridional orbifold filling
  identify(link)             SnapPy identification (census names) of the knot exterior
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction

import snappy  # noqa: F401  (must be imported before spherogram.Link.exterior works)
from spherogram import Link
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from .graph import EmbeddedTaitGraph, to_pd, to_link, from_link, to_oriented_pd


# ----------------------------------------------------------------------------- Goeritz / H1 / linking form
def goeritz_matrix(g: EmbeddedTaitGraph) -> Matrix:
    verts = list(g.vertices)
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    M = [[0] * n for _ in range(n)]
    for e, (u, w) in g.edges.items():
        s = g.signs[e]
        i, j = idx[u], idx[w]
        M[i][i] += s
        M[j][j] += s
        M[i][j] -= s
        M[j][i] -= s
    return Matrix([row[:-1] for row in M[:-1]]) if n > 1 else Matrix([])


def determinant(g: EmbeddedTaitGraph) -> int:
    G = goeritz_matrix(g)
    return abs(int(G.det())) if G.shape[0] else 1


def h1_torsion(g: EmbeddedTaitGraph) -> list[int]:
    G = goeritz_matrix(g)
    n = G.shape[0]
    if n == 0:
        return []
    snf = smith_normal_form(G, domain=ZZ)
    return [abs(int(snf[i, i])) for i in range(n) if abs(int(snf[i, i])) != 1]


def linking_form_class(g: EmbeddedTaitGraph) -> int | None:
    """A generator's self-pairing numerator for cyclic H1, using +G^-1.

    Use the shared exhaustive primary-component construction rather than
    searching a few short vectors. The presentation still comes independently
    from this graph's Goeritz form. The global orientation sign is unspecified.
    """
    G = goeritz_matrix(g)
    D = determinant(g)
    if D == 1 or h1_torsion(g) != [D]:
        return None
    import sys
    from pathlib import Path
    root = str(Path(__file__).resolve().parents[2])
    if root not in sys.path:
        sys.path.insert(0, root)
    from lower_bounds.linking_pairing import generator_self_linking
    return generator_self_linking(G, D)


def lickorish_allowed(g: EmbeddedTaitGraph) -> set[int] | None:
    """Set of signs s with linking form ≅ <s*2/D>. Empty set => u(K) != 1. {±1} when D = 1 (no info)."""
    D = determinant(g)
    if D == 1:
        return {1, -1}
    if h1_torsion(g) != [D]:
        return set()
    a = linking_form_class(g)
    if a is None:
        return None
    squares = {(k * k) % D for k in range(1, D) if math.gcd(k, D) == 1}
    ia = pow(a, -1, D)
    return {s for s in (1, -1) if (s * 2 * ia) % D in squares}


# ----------------------------------------------------------------------------- Jones / Casson–Walker
def jones_regina(obj) -> dict[int, int]:
    """Jones polynomial {t-exponent: coeff} of an EmbeddedTaitGraph, PD code or spherogram Link (Regina)."""
    import regina
    # Regina's fromPD expects a consistently *oriented* PD code (KnotTheory convention).  Our to_pd()
    # and arbitrary user PDs are not oriented, so always pass the PD through spherogram first, which
    # re-orients the strands without changing the diagram.
    if isinstance(obj, EmbeddedTaitGraph):
        if not obj.edges:
            return {0: 1}
        pd = to_oriented_pd(obj)
    elif isinstance(obj, Link):
        pd = obj.PD_code(min_strand_index=1)
    else:
        pd0 = [list(q) for q in obj]
        if not pd0:
            return {0: 1}
        pd = Link(pd0).PD_code(min_strand_index=1)
    J = regina.Link.fromPD(pd).jones()
    out = {}
    for e in range(J.minExp(), J.maxExp() + 1):
        c = int(str(J[e]))
        if c:
            if e % 2:
                raise ValueError('odd sqrt(t) exponent: not a knot')
            out[e // 2] = c
    return out


def casson_walker_allowed(jones: dict[int, int], signature: int) -> set[int]:
    """u=1 => D/2 * (lambda_CW(Sigma_2 K) -/+ lambda_CW(L(D,2))) is an integer for the surgery sign +/-
    (Mullins' formula lambda = sigma/8 - V'(-1)/(12 V(-1)); Boyer–Lines).

    These are the signs in the Jones/signature orientation convention. Do not
    intersect them with signs from an arbitrary checkerboard Goeritz matrix
    without identifying the corresponding orientation of the branched cover.
    """
    V = sum(c * (-1) ** (e % 2) for e, c in jones.items())
    dV = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in jones.items())
    D = abs(V)
    lam = Fraction(signature, 8) - Fraction(dV, 12 * V)
    lens = -Fraction((D - 1) * (D - 5), 48 * D)
    plus = Fraction(D, 2) * (lam - lens)
    minus = Fraction(D, 2) * (lam + lens)
    return {s for s, val in ((1, plus), (-1, minus)) if val.denominator == 1}


# ----------------------------------------------------------------------------- HFK / unknot / identification
def hfk(link_or_graph):
    """Knot Floer homology (dict from knot_floer_homology). The HFK program rejects diagrams with a
    reducible R1 kink, so the diagram is first simplified by R1/R2 moves (knot type unchanged)."""
    L = to_link(link_or_graph) if isinstance(link_or_graph, EmbeddedTaitGraph) else link_or_graph
    L2 = L.copy()
    L2.simplify('basic')
    return L2.knot_floer_homology()


def is_unknot(link_or_graph) -> bool:
    """HFK detects the unknot: Seifert genus 0 <=> unknot (rigorous)."""
    if isinstance(link_or_graph, EmbeddedTaitGraph) and not link_or_graph.edges:
        return True
    return hfk(link_or_graph)['seifert_genus'] == 0


def knot_sig(link_or_graph, use_reflection: bool = False, use_reversal: bool = True) -> str:
    """Regina's canonical signature of the *diagram* (diagram-level invariant, for tests)."""
    import regina
    L = to_link(link_or_graph) if isinstance(link_or_graph, EmbeddedTaitGraph) else link_or_graph
    R = regina.Link.fromPD(L.PD_code(min_strand_index=1))
    try:
        return R.sig(use_reflection, use_reversal)
    except Exception:
        return R.knotSig(use_reflection, use_reversal)


def identify(link_or_graph):
    L = to_link(link_or_graph) if isinstance(link_or_graph, EmbeddedTaitGraph) else link_or_graph
    return L.exterior().identify()


def sigma2_snappy(link_or_graph):
    """Double branched cover, preserving the knot's meridian through the cover.

    The exterior from Link has the knot meridian as (1,0). Filling with (2,0)
    first gives the order-two knot orbifold. Its connected double cyclic cover
    is a manifold: the lifted filling is primitive and removes the singular
    locus. SnapPy carries this filling through the covering construction.
    Homology is only a consistency check; its order does not identify a slope.
    """
    L = to_link(link_or_graph) if isinstance(link_or_graph, EmbeddedTaitGraph) else link_or_graph
    if len(L.link_components) + L.unlinked_unknot_components != 1:
        raise ValueError('expected a one-component knot')
    g = from_link(L)
    D = determinant(g)
    orbifold = L.exterior()
    orbifold.dehn_fill((2, 0))
    covers = orbifold.covers(2, cover_type='cyclic')
    if len(covers) != 1:
        raise ValueError('expected a unique 2-fold cover')
    N = covers[0]
    for a, b in N.cusp_info('filling'):
        if a != int(a) or b != int(b) or math.gcd(int(a), int(b)) != 1:
            raise ValueError('the orbifold cover did not lift to a primitive manifold filling')
    H = N.homology()
    if H.betti_number() != 0 or H.order() != D:
        raise ValueError('branched-cover homology disagrees with the Goeritz determinant')
    return N


def invariant_fingerprint(link_or_graph) -> tuple:
    """Cheap knot-type fingerprint for tests/search: (det, Jones, HFK ranks, tau)."""
    L = to_link(link_or_graph) if isinstance(link_or_graph, EmbeddedTaitGraph) else link_or_graph
    g = link_or_graph if isinstance(link_or_graph, EmbeddedTaitGraph) else from_link(L)
    H = hfk(L)
    return (determinant(g), tuple(sorted(jones_regina(L).items())), tuple(sorted(H['ranks'].items())), H['tau'])


# ----------------------------------------------------------------------------- knot Floer torsion order
def hfk_minus_torsion_order(pd, prime: int = 2, nmax: int = 8) -> int:
    """Maximal order of U-torsion in HFK^-(K) (F[U]-module), a lower bound for the unknotting number
    (Alishahi–Eftekhary 2020; see also Juhász–Miller–Zemke 2020: Ord_U(K) <= bridge index - 1).

    HFK^- is the homology of the n_z = 0 (pure U-power) part of the "UV = 0" complex returned by
    knot_floer_homology (arrow a -> b with U^i: i = A(b) - A(a) >= 0, M(b) = M(a) - 1 + 2i).  With
    f(N) = dim_F H_*(C / U^N) = r*N + 2*sum_k m_k*min(k, N) (r = 1 for a knot), the largest N with
    f(N) - f(N-1) > r is the maximal torsion order.  Validated against the list of <= 10-crossing knots with
    order 2 (8_19, 10_124, 10_128, 10_139, 10_152, 10_154, 10_161)."""
    import numpy as np
    from knot_floer_homology import pd_to_hfk
    data = pd_to_hfk(pd, prime=prime, complex=True)
    gens, diff = data['generators'], data['differentials']
    names = list(gens)
    idx = {g: i for i, g in enumerate(names)}
    arrows = []
    for (a, b), c in diff.items():
        Aa, Ma = gens[a]
        Ab, Mb = gens[b]
        i = Ab - Aa
        if i >= 0 and Mb == Ma - 1 + 2 * i:
            arrows.append((idx[a], idx[b], i, c % prime))
    n = len(names)

    def rank_mod_p(A):
        A = A % prime
        r = 0
        rows, cols = A.shape
        for cc in range(cols):
            piv = next((rr for rr in range(r, rows) if A[rr, cc]), None)
            if piv is None:
                continue
            A[[r, piv]] = A[[piv, r]]
            A[r] = (A[r] * pow(int(A[r, cc]), prime - 2, prime)) % prime
            nz = np.nonzero(A[:, cc])[0]
            for rr in nz:
                if rr != r:
                    A[rr] = (A[rr] - A[rr, cc] * A[r]) % prime
            r += 1
            if r == rows:
                break
        return r

    def f(N):
        dim = n * N
        M = np.zeros((dim, dim), dtype=np.int64)
        for (a, b, i, c) in arrows:
            for k in range(N - i):
                M[b * N + k + i, a * N + k] = (M[b * N + k + i, a * N + k] + c) % prime
        return dim - 2 * rank_mod_p(M)

    prev = None
    best = 0
    for N in range(1, nmax + 1):
        cur = f(N)
        inc = cur if prev is None else cur - prev
        if inc > 1:
            best = N
        prev = cur
    return best
