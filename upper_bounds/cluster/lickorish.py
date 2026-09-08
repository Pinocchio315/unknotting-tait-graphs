"""Exact Lickorish pre-filter for partner routes (stdlib only).

A knot with unknotting number one has cyclic H_1 of its double branched cover and a generator of
self-linking +-2/D (Lickorish).  For a diagram whose Goeritz matrix is G, changing the crossings
of a subset S gives G' = G + sum_{e in S} (-2 eta_e) b_e b_e^T; its exact rational inverse follows
from G^{-1} by Sherman-Morrison updates in O(n^2), so the two conditions can be tested for every
determinant hit before any diagram reduction is attempted.  Hits that fail cannot have a u = 1
partner (no loss).
"""
from __future__ import annotations

import itertools
import math
from fractions import Fraction

_SQUARES = {}


def _squares(D):
    if D not in _SQUARES:
        _SQUARES[D] = {(k * k) % D for k in range(1, D) if math.gcd(k, D) == 1}
    return _SQUARES[D]


class LickorishFilter:
    def __init__(self, g):
        self.edge_ids = sorted(g.edges)
        verts = sorted(g.vertices, key=str)
        drop = verts[-1]
        idx = {v: i for i, v in enumerate(verts) if v != drop}
        self.n = n = len(verts) - 1
        self.eta = [g.signs[e] for e in self.edge_ids]
        self.b = []
        for e in self.edge_ids:
            u, v = g.edges[e]
            vec = [0] * n
            if u != v:
                if u in idx:
                    vec[idx[u]] += 1
                if v in idx:
                    vec[idx[v]] -= 1
            self.b.append(vec)
        G = [[Fraction(0)] * n for _ in range(n)]
        for eta, vec in zip(self.eta, self.b):
            nz = [(i, c) for i, c in enumerate(vec) if c]
            for i, ci in nz:
                for j, cj in nz:
                    G[i][j] += eta * ci * cj
        self.Ginv = self._inverse(G)

    @staticmethod
    def _inverse(G):
        n = len(G)
        A = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(G)]
        for c in range(n):
            piv = next(r for r in range(c, n) if A[r][c] != 0)
            A[c], A[piv] = A[piv], A[c]
            f = A[c][c]
            A[c] = [x / f for x in A[c]]
            for r in range(n):
                if r != c and A[r][c] != 0:
                    m = A[r][c]
                    A[r] = [x - m * y for x, y in zip(A[r], A[c])]
        return [row[n:] for row in A]

    def inverse_after(self, subset):
        """exact inverse of G' after flipping the given edge indices (Sherman-Morrison per edge)"""
        n = self.n
        Ginv = [row[:] for row in self.Ginv]
        for a in subset:
            c = -2 * self.eta[a]
            bvec = self.b[a]
            w = [sum(Ginv[i][j] * bvec[j] for j in range(n) if bvec[j]) for i in range(n)]
            denom = 1 + c * sum(bvec[i] * w[i] for i in range(n) if bvec[i])
            if denom == 0:
                return None
            f = Fraction(c) / denom
            for i in range(n):
                if w[i] == 0:
                    continue
                wi = f * w[i]
                for j in range(n):
                    if w[j] != 0:
                        Ginv[i][j] -= wi * w[j]
        return Ginv

    def passes(self, subset, D):
        """Conservative necessary test for a partner with unknotting number one.

        True includes an inconclusive generator search; it never certifies u=1.
        D=1 gives a trivial linking pairing and cannot obstruct u=1: a nontrivial
        knot with determinant one need not be the unknot. D must be the exact
        determinant, independently of the modular determinant screening stage.
        """
        if D == 1:
            return True
        Ginv = self.inverse_after(subset)
        if Ginv is None:
            return False
        n = self.n
        # cyclic iff gcd of the entries of the adjugate D * G'^{-1} is 1
        gcd = 0
        for row in Ginv:
            for x in row:
                y = x * D
                if y.denominator != 1:
                    return False                          # inconsistent D (should not happen)
                gcd = math.gcd(gcd, abs(int(y)))
                if gcd == 1:
                    break
            if gcd == 1:
                break
        if gcd != 1:
            return False
        squares = _squares(D)

        def order_and_selflink(x):
            col = [sum(Ginv[i][j] * x[j] for j in range(n) if x[j]) for i in range(n)]
            den = 1
            for e in col:
                den = den * e.denominator // math.gcd(den, e.denominator)
            if den != D:
                return None
            lam = sum(x[i] * col[i] for i in range(n) if x[i])
            return (lam.numerator * (D // lam.denominator)) % D
        for k in (1, 2, 3):
            for S in itertools.combinations(range(n), k):
                x = [1 if j in S else 0 for j in range(n)]
                a = order_and_selflink(x)
                if a is None:
                    continue
                if math.gcd(a, D) != 1:
                    return False
                ia = pow(a, -1, D)
                return any((s * 2 * ia) % D in squares for s in (1, -1))
        return True                                       # generator not found among small vectors: keep
