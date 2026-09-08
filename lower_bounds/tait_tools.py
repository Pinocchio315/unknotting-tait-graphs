"""Sage-free Tait-graph / Goeritz / linking-pairing utilities for the unknot-venv
(spherogram without Sage cannot build white_graph/goeritz_matrix/jones; regina supplies Jones).

Usage (inside ~/.pyenv/versions/unknot-venv):
    from tait_tools import *
    L = link_from_knotinfo('12n_491')          # spherogram Link from KnotInfo PD
    G, W = goeritz_matrix(L)                     # signed Tait (white) graph + Goeritz matrix
    det_K = abs(int(G.det()))
    lickorish_u1_allowed(L)                      # set of allowed surgery signs; empty set => u != 1
    jones_regina(L)                              # Jones polynomial dict {exponent: coeff} in t (regina)
    sigma2_homology(L)                           # H1 of double branched cover via SnapPy (order check)
"""
from __future__ import annotations

import ast
import csv
import math
import os
from fractions import Fraction

import snappy  # must be imported before spherogram.Link.exterior() works
import networkx as nx
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
from spherogram import Link
from spherogram.links.links_base import CrossingStrand

# ---------------------------------------------------------------------------
# KnotInfo access
# ---------------------------------------------------------------------------
_KNOTINFO_ROWS = None


def knotinfo_rows() -> dict:
    """name -> row dict from the installed database_knotinfo CSV ('|'-delimited)."""
    global _KNOTINFO_ROWS
    if _KNOTINFO_ROWS is None:
        import database_knotinfo
        root = os.path.dirname(database_knotinfo.__file__)
        path = next(os.path.join(dp, f) for dp, dn, fn in os.walk(root) for f in fn
                    if f == 'knotinfo_data_complete.csv')
        with open(path, newline='', encoding='utf-8') as h:
            _KNOTINFO_ROWS = {r['name']: r for r in csv.DictReader(h, delimiter='|')}
    return _KNOTINFO_ROWS


def link_from_knotinfo(name: str) -> Link:
    """spherogram Link built from the KnotInfo PD code (KnotInfo names like '12n_491', '7_4')."""
    return Link(ast.literal_eval(knotinfo_rows()[name]['pd_notation']))


# ---------------------------------------------------------------------------
# Tait (checkerboard) graph and Goeritz matrix — port of spherogram.white_graph
# ---------------------------------------------------------------------------
def white_graph(link: Link) -> nx.MultiGraph:
    """Signed checkerboard multigraph (one colour class) of a non-split diagram.
    Vertices = face indices (as in link.faces()), edges = crossings with 'sign' in {+1,-1}
    following the Gordon–Litherland convention used by spherogram."""
    face_of = {corner: n for n, face in enumerate(link.faces()) for corner in face}
    G = nx.MultiGraph()
    for c in link.crossings:
        G.add_edge(face_of[CrossingStrand(c, 0)], face_of[CrossingStrand(c, 2)], crossing=c, sign=1)
        G.add_edge(face_of[CrossingStrand(c, 1)], face_of[CrossingStrand(c, 3)], crossing=c, sign=-1)
    comps = sorted(nx.connected_components(G), key=lambda s: sorted(s))
    if len(comps) > 2:
        raise ValueError('split diagram')
    # spherogram keeps the *second* component (sorted by vertex set); either colour class works
    return G.subgraph(comps[1]).copy() if len(comps) == 2 else G


def goeritz_matrix(link: Link):
    """Return (G, W): reduced Goeritz matrix (sympy Matrix) of the white graph W."""
    W = white_graph(link)
    verts = sorted(W.nodes())
    idx = {v: i for i, v in enumerate(verts)}
    n = len(verts)
    M = [[0] * n for _ in range(n)]
    for u, v, d in W.edges(data=True):
        s = d['sign']
        i, j = idx[u], idx[v]
        M[i][i] += s
        M[j][j] += s
        M[i][j] -= s
        M[j][i] -= s
    return Matrix([row[:-1] for row in M[:-1]]), W


def determinant(link: Link) -> int:
    G, _ = goeritz_matrix(link)
    return abs(int(G.det()))


# ---------------------------------------------------------------------------
# H1 of the double branched cover and Lickorish's unknotting-number-one test
# ---------------------------------------------------------------------------
def h1_torsion(G: Matrix) -> list[int]:
    """Nonunit Smith factors (0 denotes a free summand) of coker(G) = H1(Sigma_2(K))."""
    n = G.shape[0]
    if n == 0:
        return []
    snf = smith_normal_form(G, domain=ZZ)
    return [abs(int(snf[i, i])) for i in range(n) if abs(int(snf[i, i])) != 1]


# Shared exact algebra; independence of the Seifert and Goeritz checks is
# independence of their presentations, not of this finite-group routine.
try:
    from .linking_pairing import generator_self_linking
except ImportError:  # direct script invocation
    from linking_pairing import generator_self_linking


def lickorish_u1_allowed(link: Link) -> set[int] | None:
    """Lickorish (1985): u(K)=1 => H1(Sigma_2 K) = Z/D and the linking pairing is <±2/D>.
    Returns the set of signs s in {+1,-1} with linking pairing ≅ <s*2/D>; empty set => u(K) != 1;
    {+1,-1} is returned for D=1 (no information); None if a generator was not found."""
    G, _ = goeritz_matrix(link)
    D = abs(int(G.det()))
    if D == 1:
        return {1, -1}
    if h1_torsion(G) != [D]:
        return set()
    a = generator_self_linking(G, D)
    if a is None:
        return None
    squares = {(k * k) % D for k in range(1, D) if math.gcd(k, D) == 1}
    ia = pow(a, -1, D)
    return {s for s in (1, -1) if (s * 2 * ia) % D in squares}


# ---------------------------------------------------------------------------
# Jones polynomial via regina (exact, fast), Casson–Walker test
# ---------------------------------------------------------------------------
def jones_regina(link: Link) -> dict[int, int]:
    """Jones polynomial as {exponent_of_t: coeff}. regina returns a Laurent polynomial in sqrt(t);
    exponents are halved (all even for knots)."""
    import regina
    R = regina.Link.fromPD(link.PD_code(min_strand_index=1))
    J = R.jones()
    d = {}
    for e in range(J.minExp(), J.maxExp() + 1):
        c = int(str(J[e]))
        if c:
            if e % 2:
                raise ValueError('odd sqrt(t) exponent — not a knot?')
            d[e // 2] = c
    return d


def casson_walker_u1_allowed(jones: dict[int, int], signature: int) -> set[int]:
    """Mullins + Boyer–Lines: u=1 => D/2 * (lambda_CW(Sigma) -/+ lambda_CW(L(D,2))) is an integer
    for the surgery sign +/-.  Returns the set of signs for which the integrality holds."""
    V = sum(c * (-1) ** (e % 2) for e, c in jones.items())
    dV = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in jones.items())
    D = abs(V)
    lam = Fraction(signature, 8) - Fraction(dV, 12 * V)
    lens = -Fraction((D - 1) * (D - 5), 48 * D)
    plus = Fraction(D, 2) * (lam - lens)
    minus = Fraction(D, 2) * (lam + lens)
    S = set()
    if plus.denominator == 1:
        S.add(1)
    if minus.denominator == 1:
        S.add(-1)
    return S


# ---------------------------------------------------------------------------
# Double branched cover via SnapPy (independent check of H1)
# ---------------------------------------------------------------------------
def sigma2_snappy(link: Link):
    """The branched double cover, with the filling lifted by SnapPy.

    Fill the original knot exterior along twice its *known* meridian. This
    is the orbifold with cone angle pi along K. Its connected cyclic double
    cover has primitive lifted filling and is Sigma_2(K). Filling a cover
    along an arbitrary slope with the right |H1| does not identify a manifold.
    """
    if len(link.link_components) != 1:
        raise ValueError('the branched-cover helper requires a knot')
    exterior = link.exterior()
    exterior.dehn_fill((2, 0))
    covers = exterior.covers(2, cover_type='cyclic')
    if len(covers) != 1:
        raise ValueError('expected one connected double cover of the knot orbifold')
    cover = covers[0]
    for cusp in cover.cusp_info():
        meridian, longitude = cusp['filling']
        if not float(meridian).is_integer() or not float(longitude).is_integer() or math.gcd(int(meridian), int(longitude)) != 1:
            raise ValueError('the lifted filling is not primitive')
    homology = cover.homology()
    if homology.betti_number() or homology.order() != determinant(link):
        raise ArithmeticError('branched-cover homology does not match the Goeritz determinant')
    return cover


def sigma2_homology(link: Link):
    return sigma2_snappy(link).homology()
