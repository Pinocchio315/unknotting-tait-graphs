"""Exact correction-term obstruction for alternating knots (paper §§2.4, 3).

Owens, *Unknotting information from Heegaard Floer homology* (2008),
Theorems 1 and 5, assumes a sequence with n = sigma(K)/2 negative
crossing changes after choosing the representative with positive signature.
For |sigma| = 4 and a sequence of length two, the signature forces this
hypothesis. For |sigma| = 2, the pattern (p,n) = (0,2) is not covered:
``obstruct_u2`` must first exclude it by Traczyk's criterion before testing
(1,1). Failure for just one permitted sign pattern is not an obstruction to u=2.

The surgery form Qt must admit a group isomorphism phi with
    m_Qt(g) >= m_G(phi(g)),   m_Qt(g) - m_G(phi(g)) in 2 Z.
The sharp positive Goeritz form computes the boundary correction terms.
For odd determinant we label characteristic classes by [xi] in coker(Q).
Conjugation sends [xi] to [-xi], so the unique spin structure is labelled 0;
there is no subtraction of diag(Q). This convention must be shared by both
forms. The spin correction term -sigma/4 selects the appropriate cover
orientation (Manolescu–Owens).

All quadratic forms, characteristic-vector bounds, and correction-term
comparisons use exact integers or Fractions. A completed obstruction is a
necessary condition: PASS only means that this test leaves a surgery form
possible. A bounded isomorphism search returns UNDECIDED, never OBSTRUCTED.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction

import numpy as np


# ----------------------------------------------------------------------------- Goeritz forms from a PD code
def definite_goeritz_candidates(pd):
    """The two positive-definite checkerboard Goeritz forms of an alternating knot diagram: the reduced graph
    Laplacians of the two checkerboard (Tait) graphs.  One belongs to the diagram, the other to its mirror;
    the caller picks the correct side via d(Sigma_2, spin) = -sigma/4."""
    import _compat  # noqa: F401
    from spherogram import Link
    import taitgraph
    g = taitgraph.from_link(Link([list(q) for q in pd]))
    signs = set(g.signs.values())
    if len(signs) != 1:
        raise ValueError('diagram is not alternating (mixed checkerboard signs)')

    def lap_reduced(h):
        verts = sorted(h.vertices)
        idx = {v: i for i, v in enumerate(verts)}
        n = len(verts)
        M = np.zeros((n, n), dtype=np.int64)
        for e, (u, w) in h.edges.items():
            if u == w:
                continue
            i, j = idx[u], idx[w]
            M[i, i] += 1
            M[j, j] += 1
            M[i, j] -= 1
            M[j, i] -= 1
        return M[:-1, :-1].copy()

    return [lap_reduced(g), lap_reduced(taitgraph.dual(g))]


# ----------------------------------------------------------------------------- exact linear algebra helpers
def det_int(M) -> int:
    """Exact determinant of an integer matrix (Bareiss)."""
    A = [[int(x) for x in row] for row in M]
    n = len(A)
    if n == 0:
        return 1
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


def adjugate_int(M):
    """Exact adjugate of an integer matrix via Fraction inverse."""
    n = len(M)
    D = det_int(M)
    A = [[Fraction(int(M[i][j])) for j in range(n)] for i in range(n)]
    I = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        I[c], I[piv] = I[piv], I[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        I[c] = [x / pv for x in I[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
                I[r] = [x - f * y for x, y in zip(I[r], I[c])]
    adj = [[I[i][j] * D for j in range(n)] for i in range(n)]
    out = np.zeros((n, n), dtype=object)
    for i in range(n):
        for j in range(n):
            assert adj[i][j].denominator == 1
            out[i, j] = int(adj[i][j])
    return out


def snf_U(M):
    """Diagonalize over Z, retaining the left unimodular change of basis.

    Return U,d with U M V diagonal and diagonal entries of absolute value d.
    The entries need not satisfy the divisibility convention of Smith normal
    form. They still present the cokernel as a product of cyclic groups;
    ClassMap splits them into prime powers to obtain canonical group factors.
    """
    A = [[int(x) for x in row] for row in M]
    n = len(A)
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]
        U[i], U[j] = U[j], U[i]

    def add_row(i, j, c):
        A[i] = [x + c * y for x, y in zip(A[i], A[j])]
        U[i] = [x + c * y for x, y in zip(U[i], U[j])]

    def swap_cols(i, j):
        for r in range(n):
            A[r][i], A[r][j] = A[r][j], A[r][i]

    def add_col(i, j, c):
        for r in range(n):
            A[r][i] += c * A[r][j]

    for t in range(n):
        while True:
            # find pivot: smallest nonzero |entry| in submatrix
            best = None
            for i in range(t, n):
                for j in range(t, n):
                    if A[i][j] != 0 and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                        best = (i, j)
            if best is None:
                break
            i0, j0 = best
            swap_rows(t, i0)
            swap_cols(t, j0)
            done = True
            for i in range(t + 1, n):
                if A[i][t] != 0:
                    add_row(i, t, -(A[i][t] // A[t][t]))
                    if A[i][t] != 0:
                        done = False
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    add_col(j, t, -(A[t][j] // A[t][t]))
                    if A[t][j] != 0:
                        done = False
            # divisibility of the rest is not needed for coordinates; ensure zeros in row/col t
            if done:
                break
        # nothing
    d = [abs(A[i][i]) for i in range(n)]
    return U, d


def _prime_powers(n):
    out = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            q = 1
            while n % p == 0:
                n //= p
                q *= p
            out.append(q)
        p += 1
    if n > 1:
        out.append(n)
    return out


class ClassMap:
    """Coordinates of Z^r/(Q Z^r) in the canonical primary decomposition prod Z/p^k (via the SNF left
    transform and CRT), so that two isomorphic groups always get the same modulus list."""

    def __init__(self, Q):
        self.Q = _integral_symmetric_matrix(Q)
        if det_int(self.Q) == 0:
            raise ValueError('a discriminant group requires a nonsingular form')
        self.r = len(Q)
        U, d = snf_U(Q)
        rows, mods = [], []
        for i in range(self.r):
            if d[i] == 1:
                continue
            for q in _prime_powers(d[i]):
                rows.append(U[i])
                mods.append(q)
        order_ = sorted(range(len(mods)), key=lambda t: mods[t])
        self.U = [rows[t] for t in order_]
        self.d = [mods[t] for t in order_]
        self.order = 1
        for x in self.d:
            self.order *= x

    def coords(self, xi):
        return tuple(sum(u * x for u, x in zip(row, xi)) % dd for row, dd in zip(self.U, self.d))

    def all_elements(self):
        return itertools.product(*[range(dd) for dd in self.d])


# ----------------------------------------------------------------------------- exact quadratic-form enumeration
class IncompleteCorrectionTerms(RuntimeError):
    """The configured ellipsoid search did not reach every characteristic class."""


class IsomorphismSearchLimit(RuntimeError):
    """The finite group-isomorphism search exceeds its configured work limit."""


def _integral_symmetric_matrix(Q):
    """Validate a form without silently truncating rational matrix entries."""
    rows = [list(row) for row in Q]
    n = len(rows)
    if any(len(row) != n for row in rows):
        raise ValueError('the form must be square')
    if any(value != int(value) for row in rows for value in row):
        raise ValueError('the form must have integral entries')
    matrix = np.array([[int(value) for value in row] for row in rows], dtype=object).reshape(n, n)
    if any(matrix[i, j] != matrix[j, i] for i in range(n) for j in range(i)):
        raise ValueError('the form must be symmetric')
    return matrix


def positive_definite(Q):
    """Sylvester's criterion with exact determinants, not an eigenvalue tolerance."""
    Q = _integral_symmetric_matrix(Q)
    return all(det_int(Q[:i, :i]) > 0 for i in range(1, len(Q) + 1))


def _ldl_upper(A):
    """Return exact U,d with A=U.T diag(d) U and U upper unit triangular.

    Descending coordinate enumeration then makes each squared term depend
    only on the current coordinate and coordinates already chosen.
    """
    n = len(A)
    L = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    d = []
    for i in range(n):
        pivot = A[i][i] - sum(L[i][k] ** 2 * d[k] for k in range(i))
        if pivot <= 0:
            raise ValueError('the form must be positive definite')
        d.append(pivot)
        for j in range(i + 1, n):
            L[j][i] = (A[j][i] - sum(L[j][k] * L[i][k] * d[k]
                                      for k in range(i))) / pivot
    return [[L[j][i] for j in range(n)] for i in range(n)], d


def _enumerate_characteristic(U, weights, parity, bound, cm, best):
    """Enumerate an exact ellipsoid, including its boundary.

    At coordinate i the condition is (q*x+p)^2 <= q^2*remaining/d_i,
    where p/q is the contribution of the later coordinates. Integer square
    roots give rigorous endpoints; no floating-point pruning can discard a
    shortest representative. This matters because an overestimated m_G can
    otherwise produce a false obstruction.
    """
    rank = len(parity)
    xi = [0] * rank

    def visit(i, partial):
        if i < 0:
            label = cm.coords(xi)
            value = (partial - rank) / 4
            if label not in best or value < best[label]:
                best[label] = value
            return
        remaining = bound - partial
        if remaining < 0:
            return
        offset = sum((U[i][j] * xi[j] for j in range(i + 1, rank)), Fraction(0))
        p, q = offset.numerator, offset.denominator
        radius_squared = remaining * q * q / weights[i]
        radius = math.isqrt(radius_squared.numerator // radius_squared.denominator)
        lower = -((radius + p) // q)  # ceil((-radius-p)/q)
        upper = (radius - p) // q
        lower += (parity[i] - lower) % 2
        for x in range(lower, upper + 1, 2):
            xi[i] = x
            visit(i - 1, partial + weights[i] * (x + offset) ** 2)

    visit(rank - 1, Fraction(0))


def m_Q(Q, need_all_classes=True, b0=None, max_grow=24):
    """Compute m_Q on the odd-order discriminant group with exact minima.

    Characteristic covectors satisfy xi_i=Q_ii mod 2. Their images [xi] in
    coker(Q) label every spin-c structure because det(Q) is odd. An exhaustive
    ellipsoid containing a representative of each class determines every
    minimum: any vector outside it has strictly larger norm than the recorded
    representatives. Even-order forms require a different affine labelling
    and are deliberately rejected by this implementation.

    ``need_all_classes=False`` returns the partial table inside the initial
    ellipsoid and is useful for diagnostics only. Exhausting ``max_grow`` in
    complete mode raises IncompleteCorrectionTerms; a partial table must never
    be interpreted as a completed obstruction.
    """
    Q = _integral_symmetric_matrix(Q)
    rank = len(Q)
    determinant = det_int(Q)
    if determinant <= 0 or determinant % 2 == 0:
        raise ValueError('m_Q requires a positive-definite form of odd determinant')
    if max_grow < 1:
        raise ValueError('max_grow must be positive')
    adj = adjugate_int(Q)
    inverse = [[Fraction(int(adj[i, j]), determinant) for j in range(rank)]
               for i in range(rank)]
    U, weights = _ldl_upper(inverse)
    cm = ClassMap(Q)
    parity = [int(Q[i, i]) % 2 for i in range(rank)]
    bound = Fraction(rank + 8) if b0 is None else Fraction(str(b0))
    if bound < 0:
        raise ValueError('the norm bound must be nonnegative')
    best = {}
    for _ in range(max_grow):
        _enumerate_characteristic(U, weights, parity, bound, cm, best)
        if not need_all_classes or len(best) == cm.order:
            return best, cm
        # Ensure a zero initial bound still grows to a positive ellipsoid.
        bound = max(bound * Fraction(8, 5), Fraction(1))
    raise IncompleteCorrectionTerms(
        f'found {len(best)} of {cm.order} characteristic classes after {max_grow} bounds')


# ----------------------------------------------------------------------------- candidates and matching
def qtilde(m1, m2, a):
    return np.array([[m1, 1, a, 0], [1, 2, 0, 0], [a, 0, m2, 1], [0, 0, 1, 2]], dtype=np.int64)


def candidates(D, n_even):
    """All (m1, m2, a): 0 <= a < m1 <= m2, (2m1-1)(2m2-1) - 4a^2 = D, exactly n_even of {m1,m2} even."""
    out = []
    m1 = 1
    while (2 * m1 - 1) ** 2 <= D + 4 * (m1 - 1) ** 2:
        for a in range(0, m1):
            num = D + 4 * a * a
            if num % (2 * m1 - 1) == 0:
                t = num // (2 * m1 - 1)
                if t % 2 == 1:
                    m2 = (t + 1) // 2
                    if m2 >= m1:
                        ne = (m1 % 2 == 0) + (m2 % 2 == 0)
                        if ne == n_even:
                            out.append((m1, m2, a))
            # also allow a >= m1? theorem says 0 <= a < m1 <= m2 (reduced form)
        m1 += 1
    return out


def group_isos(cmA: ClassMap, cmB: ClassMap, cap=2_000_000):
    """Yield every isomorphism A -> B from possible images of generators.

    Each column obeys the order-divisibility relations in the target group;
    bijectivity is then checked on the entire finite group. Different primary
    decompositions yield no maps. IsomorphismSearchLimit means the prescribed
    cap was exceeded, not that no isomorphism exists.
    """
    if sorted(cmA.d) != sorted(cmB.d):
        return
    dA, dB = list(cmA.d), list(cmB.d)
    k = len(dA)
    if k == 0:
        yield lambda g: ()
        return
    # match up factors of equal order: enumerate hom matrices M (k x k): column j = image of e_j^A in B,
    # entries M[i][j] in Z/dB[i] with dA[j] * M[i][j] == 0 mod dB[i]  (i.e. M[i][j] multiple of dB[i]/gcd)
    choices = []
    total = 1
    for j in range(k):
        col = []
        for i in range(k):
            step = dB[i] // math.gcd(dB[i], dA[j])
            col.append(list(range(0, dB[i], step)))
        sz = 1
        for c in col:
            sz *= len(c)
        total *= sz
        choices.append(col)
    if total > cap:
        raise IsomorphismSearchLimit('group-isomorphism search cap exceeded')
    ele = list(itertools.product(*[range(d) for d in dA]))

    def is_iso(M):
        seen = set()
        for g in ele:
            img = tuple(sum(M[i][j] * g[j] for j in range(k)) % dB[i] for i in range(k))
            if img in seen:
                return False
            seen.add(img)
        return True

    for cols in itertools.product(*[itertools.product(*choices[j]) for j in range(k)]):
        M = [[cols[j][i] for j in range(k)] for i in range(k)]
        # quick: injectivity on generators' orders
        ok = True
        for j in range(k):
            img = tuple(M[i][j] for i in range(k))
            o = 1
            for i in range(k):
                if img[i]:
                    o = o * (dB[i] // math.gcd(img[i], dB[i])) // math.gcd(o, dB[i] // math.gcd(img[i], dB[i]))
            if o != dA[j]:
                ok = False
                break
        if not ok:
            continue
        if is_iso(M):
            def phi(g, M=M):
                return tuple(sum(M[i][j] * g[j] for j in range(k)) % dB[i] for i in range(k))
            yield phi


def matching_exists(mQt: dict, cmQt: ClassMap, mG: dict, cmG: ClassMap) -> bool | None:
    """Test all group isomorphisms in the shared spin-origin convention.

    True means at least one candidate satisfies the correction-term necessary
    conditions, not that a surgery or an unknotting sequence exists. False
    requires a completed search. None means the search cap or incomplete
    input tables leave the comparison undecided. Congruence mod 2 means an
    even *integer* difference of the exact rational correction terms.
    """
    if cmQt.d != cmG.d:
        return False
    elems = list(cmQt.all_elements())
    if any(g not in mQt for g in elems) or any(g not in mG for g in cmG.all_elements()):
        return None
    try:
        isos = group_isos(cmQt, cmG)
        for phi in isos:
            ok = True
            for g in elems:
                a = mQt[g]
                b = mG[phi(g)]
                dif = a - b
                if dif < 0 or dif.denominator != 1 or dif.numerator % 2 != 0:
                    ok = False
                    break
            if ok:
                return True
        return False
    except IsomorphismSearchLimit:
        return None  # An unfinished search never excludes a surgery form.


# ----------------------------------------------------------------------------- top level per knot
def obstruct_u2(pd, sigma, det_k, jones=None, verbose=False):
    """Obstruction to u = 2 for an alternating knot.
    |sigma| = 4: Owens Theorem 1 alone (sign pattern (0,2) forced).
    |sigma| = 2: Owens Theorem 1 with n = 1 rules out the (1,1) pattern; Traczyk's e^{i pi/3} criterion must
    also rule out (0,2) (needs jones = {exp: coeff} of the KnotInfo diagram; requires dim3 = 2, eps = -1).
    Returns dict with verdict: 'OBSTRUCTED' (u >= 3), 'PASS', 'NOT_APPLICABLE', 'UNDECIDED', or 'ERROR'."""
    sigma_signed = int(sigma)
    sigma = abs(sigma_signed)
    if sigma == 2:
        if jones is None:
            return {'verdict': 'NOT_APPLICABLE', 'reason': 'sigma 2 needs the Jones polynomial (Traczyk)'}
        from traczyk import mirror_jones, traczyk_excludes_00_2
        Jpos = jones if sigma_signed > 0 else mirror_jones(jones)
        try:
            excl, eps, d3 = traczyk_excludes_00_2(Jpos)
        except ValueError as e:
            return {'verdict': 'ERROR', 'reason': f'traczyk: {e}'}
        if d3 > 2:
            return {'verdict': 'OBSTRUCTED', 'by': 'wendt', 'eps': eps, 'dim3': d3}
        if not excl:
            return {'verdict': 'PASS', 'by': 'traczyk allows (0,2)', 'eps': eps, 'dim3': d3}
        res = _owens_theorem1(pd, sigma, det_k, n=1)
        res.update({'eps': eps, 'dim3': d3, 'by': 'owens(1,1) + traczyk(0,2)'})
        return res
    if sigma != 4:
        return {'verdict': 'NOT_APPLICABLE', 'reason': f'|sigma| = {sigma} not in (2, 4)'}
    return _owens_theorem1(pd, sigma, det_k, n=2)


def _owens_theorem1(pd, sigma, det_k, n, verbose=False):
    # Goeritz side: pick the definite form with m_G(0) == -sigma/4  (d(Sigma_2, spin) = -sigma/4)
    forms = definite_goeritz_candidates(pd)
    target0 = Fraction(-sigma, 4)
    chosen = None
    info = []
    for G in forms:
        Dg = det_int(G)
        if Dg != det_k:
            info.append(('det mismatch', Dg))
            continue
        mG, cmG = m_Q(G)
        # spin class = identity of Z^r/Q Z^r (odd order group: the unique self-conjugate class is g = 0)
        g_spin = cmG.coords([0] * len(G))
        info.append(('d_spin', str(mG.get(g_spin))))
        if mG.get(g_spin) == target0:
            chosen = (G, mG, cmG)
            break
    if chosen is None:
        return {'verdict': 'ERROR', 'reason': f'no Goeritz side with d(spin) = -sigma/4; info {info}'}
    G, mG, cmG = chosen
    cands = candidates(det_k, n)
    if not cands:
        return {'verdict': 'OBSTRUCTED', 'candidates': 0}
    results = []
    for (m1, m2, a) in cands:
        Qt = qtilde(m1, m2, a)
        assert det_int(Qt) == det_k
        mQt, cmQt = m_Q(Qt)
        r = matching_exists(mQt, cmQt, mG, cmG)
        results.append(((m1, m2, a), r))
        if r is True:
            return {'verdict': 'PASS', 'candidates': len(cands), 'witness_form': (m1, m2, a)}
    if any(r is None for _, r in results):
        return {'verdict': 'UNDECIDED', 'candidates': len(cands)}
    return {'verdict': 'OBSTRUCTED', 'candidates': len(cands)}
