"""Owens-type obstruction to u(K) = 4 for alternating knots with |sigma| = 8 (four same-sign changes).

Same scheme as owens_u3.py one rank up: Sigma_2(K) must bound the positive-definite 8x8 plumbing of four
(m_i, 2) chains linked by a_ij (all m_i even), det = det K, and its d-invariants must match the sharp
Goeritz form.  Quaternary forms Q' = 2M - I of determinant D are enumerated through Minkowski-reduced
representatives (a11 a22 a33 a44 <= 4 D) and the GL(4, F_2) lifts that carry the mod-2 form to the
identity, which reaches every Gamma(2)-class with the plumbing parity pattern.
"""
from __future__ import annotations
import itertools, math, sys, collections
from fractions import Fraction
import numpy as np
from owens_obstruction import definite_goeritz_candidates, det_int, m_Q, matching_exists, positive_definite, adjugate_int
from linkform import linking_invariants, compatible

N = 4


def reduced_quaternary_forms(D):
    """Symmetric integer 4x4 forms with a<=b<=c<=d, |2x_ij|<=a_ii (i<j), a b c d <= 4D, det = D, positive
    definite.  Every GL(4,Z)-class of positive quaternary forms of determinant D has a Minkowski-reduced
    representative in this set.  The last row is solved from the 3x3 block: d = (D + v^T adj(F3) v)/det F3."""
    out = []
    lim = 4 * D
    for a in range(1, int(lim ** 0.25) + 2):
        for b in range(a, int((lim / a) ** (1 / 3)) + 2):
            cmax = math.isqrt(lim // (a * b)) + 1
            for c in range(b, cmax + 1):
                dmax = lim // (a * b * c)
                if dmax < c:
                    break
                r12 = range(-(a // 2), a // 2 + 1); r23 = range(-(b // 2), b // 2 + 1)
                # v-grid for the last column
                v1 = np.arange(-(a // 2), a // 2 + 1); v2 = np.arange(-(b // 2), b // 2 + 1); v3 = np.arange(-(c // 2), c // 2 + 1)
                V = np.array(np.meshgrid(v1, v2, v3, indexing='ij')).reshape(3, -1).T.astype(np.int64)   # (N,3)
                for x12 in r12:
                    if a * b - x12 * x12 <= 0:
                        continue
                    for x13 in r12:
                        for x23 in r23:
                            F3 = np.array([[a, x12, x13], [x12, b, x23], [x13, x23, c]], dtype=np.int64)
                            d3 = det_int(F3)
                            if d3 <= 0:
                                continue
                            adj3 = np.array(adjugate_int(F3), dtype=np.int64)
                            num = D + np.einsum('ni,ij,nj->n', V, adj3, V)      # d * d3 = D + v^T adj v
                            ok = (num % d3 == 0)
                            dd = num[ok] // d3
                            sel = (dd >= c) & (dd <= dmax)
                            for v, d in zip(V[ok][sel], dd[sel]):
                                F = np.array([[a, x12, x13, v[0]], [x12, b, x23, v[1]], [x13, x23, c, v[2]], [v[0], v[1], v[2], int(d)]], dtype=np.int64)
                                out.append(F)
    return out


def gl4f2_lifts():
    """integer lift (product of elementary matrices, det 1) for every element of GL(4, F_2)"""
    gens = []
    for i in range(N):
        for j in range(N):
            if i != j:
                E = np.eye(N, dtype=np.int64); E[i, j] = 1; gens.append(E)
    start = np.eye(N, dtype=np.int64)
    key = lambda P: tuple((P % 2).flatten())
    lifts = {key(start): start}
    frontier = [start]
    while frontier:
        nxt = []
        for P in frontier:
            for E in gens:
                Q = (P @ E)
                k = key(Q)
                if k not in lifts:
                    lifts[k] = Q; nxt.append(Q)
        frontier = nxt
    assert len(lifts) == 20160, len(lifts)
    return lifts



_SYM = None


def _canonical_key(m, A):
    """lexicographically smallest (m, a12, a13, a14, a23, a24, a34) over permutations and sign changes"""
    global _SYM
    if _SYM is None:
        _SYM = [(perm, signs) for perm in itertools.permutations(range(N)) for signs in itertools.product((1, -1), repeat=N)]
    best = None
    for perm, signs in _SYM:
        mm = (m[perm[0]], m[perm[1]], m[perm[2]], m[perm[3]])
        if best is not None and mm > best[0]:
            continue
        aa = tuple(signs[i] * signs[j] * A[perm[i]][perm[j]] for i in range(N) for j in range(i + 1, N))
        kk = (mm, aa)
        if best is None or kk < best:
            best = kk
    return best


def plumbing(m, A):
    Q = np.zeros((2 * N, 2 * N), dtype=np.int64)
    for i in range(N):
        Q[2 * i, 2 * i] = m[i]; Q[2 * i, 2 * i + 1] = Q[2 * i + 1, 2 * i] = 1; Q[2 * i + 1, 2 * i + 1] = 2
        for j in range(N):
            if i != j:
                Q[2 * i, 2 * j] = A[i][j]
    return Q


_LIFTS = None
_ORTHO = {}


def candidate_plumbings(D, n_even=4, log=None):
    """Enumerate surgery candidates, preserving the number of even m_i.

    The reduction search covers GL(4,Z)-classes. Lifts of GL(4,F_2) recover
    the admissible bases modulo Gamma(2), where diagonal residues modulo four
    are fixed. The public obstruction uses n_even=4; the parameter is also
    useful for testing the parity enumeration and must not be silently ignored.
    """
    if D <= 0 or D % 2 == 0 or not 0 <= n_even <= N:
        raise ValueError('require positive odd determinant and 0 <= n_even <= 4')
    global _LIFTS
    if _LIFTS is None:
        _LIFTS = gl4f2_lifts()
    keys_all = list(_LIFTS.keys())
    seen = set(); out = []
    forms = reduced_quaternary_forms(D)
    if log: log(f'  {len(forms)} reduced quaternary forms of determinant {D}')
    I2 = tuple(np.eye(N, dtype=np.int64).flatten())
    for F in forms:
        Fb = tuple((F % 2).flatten())
        if Fb not in _ORTHO:
            Fm = np.array(Fb, dtype=np.int64).reshape(N, N)
            sel = []
            for k in keys_all:
                Pm = np.array(k, dtype=np.int64).reshape(N, N)
                if tuple(((Pm.T @ Fm @ Pm) % 2).flatten()) == I2:
                    sel.append(k)
            _ORTHO[Fb] = sel
        for k in _ORTHO[Fb]:
            P = _LIFTS[k]
            Qp = P.T @ F @ P
            m = [(int(Qp[i, i]) + 1) // 2 for i in range(N)]
            if sum(value % 2 == 0 for value in m) != n_even:
                continue
            A = [[int(Qp[i, j]) // 2 if i != j else 0 for j in range(N)] for i in range(N)]
            best = _canonical_key(m, A)
            if best in seen:
                continue
            seen.add(best)
            Q = plumbing(m, A)
            if det_int(Q) != D or not positive_definite(Q):
                continue
            out.append((best, Q))
    return out


def obstruct_u4(pd, sigma, det_k, verbose=False):
    sig = abs(int(sigma))
    if sig != 8:
        return {'verdict': 'NOT_APPLICABLE', 'reason': f'|sigma| = {sig} != 8'}
    forms = definite_goeritz_candidates(pd)
    target0 = Fraction(-sig, 4)
    chosen = None
    for G in forms:
        if det_int(G) != det_k:
            continue
        mG, cmG = m_Q(G)
        if mG.get(cmG.coords([0] * len(G))) == target0:
            chosen = (G, mG, cmG); break
    if chosen is None:
        return {'verdict': 'ERROR', 'reason': 'no Goeritz side with d(spin) = -sigma/4'}
    G, mG, cmG = chosen
    log = (lambda s: print(s, file=sys.stderr, flush=True)) if verbose else None
    cands = candidate_plumbings(det_k, log=log)
    if log: log(f'  {len(cands)} candidate plumbings')
    if not cands:
        return {'verdict': 'OBSTRUCTED', 'candidates': 0}
    # A surviving form is only a failure to obstruct. To certify OBSTRUCTED,
    # every admissible form must fail a completed exact comparison.
    undecided = 0
    invG = linking_invariants(G)
    cands = [(key, Q) for key, Q in cands if compatible(linking_invariants(Q), invG)]   # H_1 and linking form must agree
    if log: log(f'  {len(cands)} candidates after the linking-form filter')
    for key, Q in cands:
        mQ, cmQ = m_Q(Q)
        r = matching_exists(mQ, cmQ, mG, cmG)
        if r is True:
            return {'verdict': 'PASS', 'candidates': len(cands), 'witness_form': key}
        if r is None:
            undecided += 1
    if undecided:
        return {'verdict': 'UNDECIDED', 'candidates': len(cands), 'undecided': undecided}
    return {'verdict': 'OBSTRUCTED', 'candidates': len(cands)}


if __name__ == '__main__':
    import ast, json, time, database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for nm in sys.argv[1:]:
        r = rows[nm]; t0 = time.time()
        res = obstruct_u4(ast.literal_eval(r['pd_notation']), int(r['signature']), int(r['determinant']), verbose=True)
        print(json.dumps({'name': nm, 'u': str(r['unknotting_number']).strip(), 'sigma': int(r['signature']), 'det': int(r['determinant']), **res, 'seconds': round(time.time() - t0, 1)}), flush=True)
