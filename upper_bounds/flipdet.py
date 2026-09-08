"""FlipDet — determinant of the knot after flipping any subset of edge signs, O(1) per subset.

The Goeritz matrix of an embedded Tait graph is the signed graph Laplacian with one vertex
deleted: G = sum_e eta(e) b_e b_e^T with b_e = e_u - e_v (coordinates of the deleted vertex
dropped; loops contribute 0), and det K = |det G|.  Changing the crossing at edge e flips
eta(e), i.e. G' = G - 2 eta(e) b_e b_e^T --- a rank-one update per changed crossing.  By the
matrix determinant lemma,

    det G' = det G * det( I_k - 2 M Lambda ),   M[a,b] = b_{e_a}^T G^{-1} b_{e_b},
                                                Lambda = diag(eta(e_a)),

so after one modular inversion of G per diagram, the determinant of the knot obtained by
changing ANY k crossings costs a k x k determinant (k <= 3 here: a handful of multiplications).
All arithmetic is mod p = 2^61 - 1 (stdlib ints); determinants are compared as the pair
{d mod p, -d mod p} since the Goeritz determinant is defined up to sign.
An unequal residue excludes a target determinant; an equal residue is only a
candidate and must not be used as an exact determinant or unknot certificate.

Edge order: sorted(g.edges) --- the SAME order as xtait.graph.to_pd rows, so subset indices
reported by this module are PD row indices of to_pd(g).
"""
from __future__ import annotations

P = (1 << 61) - 1


def _inv_mod(a: int) -> int:
    return pow(a, P - 2, P)


class FlipDet:
    def __init__(self, g):
        self.edge_ids = sorted(g.edges)
        verts = sorted(g.vertices, key=str)
        drop = verts[-1]
        idx = {v: i for i, v in enumerate(verts) if v != drop}
        n = len(verts) - 1
        self.eta = [g.signs[e] for e in self.edge_ids]
        # b-vectors as sparse (index, coeff) lists
        self.b = []
        for e in self.edge_ids:
            u, v = g.edges[e]
            sp = []
            if u != v:
                if u in idx:
                    sp.append((idx[u], 1))
                if v in idx:
                    sp.append((idx[v], -1))
            self.b.append(sp)
        # Goeritz matrix mod P
        G = [[0] * n for _ in range(n)]
        for eta, sp in zip(self.eta, self.b):
            for i, ci in sp:
                for j, cj in sp:
                    G[i][j] = (G[i][j] + eta * ci * cj) % P
        self.det_g = self._inv_and_det(G, n)
        m = len(self.edge_ids)
        # M[a][b] = b_a^T G^{-1} b_b from the inverse's entries
        Gi = self.Gi
        self.M = [[0] * m for _ in range(m)]
        for a in range(m):
            for bb in range(a, m):
                s = 0
                for i, ci in self.b[a]:
                    for j, cj in self.b[bb]:
                        s += ci * cj * Gi[i][j]
                s %= P
                self.M[a][bb] = s
                self.M[bb][a] = s

    def _inv_and_det(self, G, n):
        """Gauss-Jordan mod P; sets self.Gi, returns det(G) mod P."""
        A = [row[:] + [1 if i == j else 0 for j in range(n)] for i, row in enumerate(G)]
        det = 1
        for c in range(n):
            piv = next((r for r in range(c, n) if A[r][c] % P), None)
            if piv is None:
                raise ZeroDivisionError('Goeritz matrix singular modulo P; determinant may be a nonzero multiple of P')
            if piv != c:
                A[c], A[piv] = A[piv], A[c]
                det = -det
            det = det * A[c][c] % P
            inv = _inv_mod(A[c][c])
            A[c] = [x * inv % P for x in A[c]]
            for r in range(n):
                if r != c and A[r][c]:
                    f = A[r][c]
                    A[r] = [(x - f * y) % P for x, y in zip(A[r], A[c])]
        self.Gi = [row[n:] for row in A]
        return det % P

    # ------------------------------------------------------------------ subset factors
    def det_after(self, subset) -> int:
        """det of the knot after flipping the edges at the given sorted PD-row indices, mod P."""
        k = len(subset)
        M, eta = self.M, self.eta
        if k == 1:
            a, = subset
            f = (1 - 2 * eta[a] * M[a][a]) % P
        elif k == 2:
            a, b = subset
            x = (1 - 2 * eta[a] * M[a][a]) % P
            y = (1 - 2 * eta[b] * M[b][b]) % P
            f = (x * y - 4 * eta[a] * eta[b] * M[a][b] * M[a][b]) % P
        elif k == 3:
            a, b, c = subset
            ea, eb, ec = eta[a], eta[b], eta[c]
            m11 = (1 - 2 * ea * M[a][a]) % P
            m22 = (1 - 2 * eb * M[b][b]) % P
            m33 = (1 - 2 * ec * M[c][c]) % P
            m12 = -2 * eb * M[a][b]
            m13 = -2 * ec * M[a][c]
            m21 = -2 * ea * M[b][a]
            m23 = -2 * ec * M[b][c]
            m31 = -2 * ea * M[c][a]
            m32 = -2 * eb * M[c][b]
            f = (m11 * (m22 * m33 - m23 * m32)
                 - m12 * (m21 * m33 - m23 * m31)
                 + m13 * (m21 * m32 - m22 * m31)) % P
        else:
            raise ValueError('k <= 3 only')
        return self.det_g * f % P

    def matches(self, subset, det_pairs) -> bool:
        return self.det_after(subset) in det_pairs


def det_pair(d: int):
    """The two residues mod P that |det| = d can take (sign of det G is not controlled)."""
    return frozenset((d % P, (-d) % P))


def det_pairs(values):
    s = set()
    for d in values:
        s.add(d % P)
        s.add((-d) % P)
    return frozenset(s)


def exact_det(g) -> int:
    """|det(Goeritz)| by exact rational Gaussian elimination, for verification."""
    verts = sorted(g.vertices, key=str)
    idx = {v: i for i, v in enumerate(verts[:-1])}
    n = len(verts) - 1
    if n == 0:
        return 1
    A = [[0] * n for _ in range(n)]
    for e, (u, v) in g.edges.items():
        s = g.signs[e]
        if u == v:
            continue
        iu, iv = idx.get(u), idx.get(v)
        if iu is not None:
            A[iu][iu] += s
        if iv is not None:
            A[iv][iv] += s
        if iu is not None and iv is not None:
            A[iu][iv] -= s
            A[iv][iu] -= s
    from fractions import Fraction
    A = [[Fraction(x) for x in row] for row in A]
    det = Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    assert det.denominator == 1
    return abs(int(det))
