"""Owens-type obstruction to u(K) = 3 for alternating knots with |sigma| = 6 (signature-sharp case).

If u(K) = 3 and |sigma(K)| = 6, all three crossing changes have the same sign, so the Montesinos trick gives
Sigma_2(K) as half-integer surgery on a 3-component link and Sigma_2(K) bounds a positive-definite plumbing
X with intersection form
    Qt = plumbing of three chains (m_i, 2), i = 1,2,3, with linking numbers a_ij between the m_i-vertices,
    m_i even, det Qt = det [[2m1-1, 2a12, 2a13], [2a12, 2m2-1, 2a23], [2a13, 2a23, 2m3-1]] = det K.
Exactly as in Owens' Theorem 1 for two changes, the d-invariants of Sigma_2(K) (sharp Goeritz form for
alternating K) must satisfy m_Qt(g) >= m_G(phi g), == mod 2, for some group isomorphism phi.
The rank-3 forms Q' = 2M - I are enumerated up to the equivalences that preserve the plumbing structure
(GL(3,Z) matrices congruent to I mod 2, permutations, sign changes): every GL(3,Z)-class of positive-definite
ternary forms of determinant D has a Minkowski-reduced representative with a d f <= 2 D, and the
Gamma(2)-classes inside a GL(3,Z)-class are reached through lifts of the 168 elements of GL(3, F_2).
If no form admits a matching, u(K) != 3; with u(K) >= 3 from the signature this gives u(K) >= 4.
"""
from __future__ import annotations
import itertools, math, sys
from fractions import Fraction
import numpy as np
from owens_obstruction import definite_goeritz_candidates, det_int, m_Q, matching_exists, positive_definite


def reduced_ternary_forms(D):
    """symmetric integer 3x3 forms with 1<=a<=d<=f, |2b|<=a, |2c|<=a, |2e|<=d, a*d*f<=2D, det=D, pos.def."""
    out = []
    for a in range(1, int(round((2 * D) ** (1 / 3))) + 2):
        for d in range(a, int(2 * D // a) + 1):
            if a * d * d > 2 * D:
                break
            fmax = (2 * D) // (a * d)
            for f in range(d, fmax + 1):
                for b in range(-(a // 2), a // 2 + 1):
                    for c in range(-(a // 2), a // 2 + 1):
                        for e in range(-(d // 2), d // 2 + 1):
                            F = np.array([[a, b, c], [b, d, e], [c, e, f]], dtype=np.int64)
                            if det_int(F) != D:
                                continue
                            if a * d - b * b <= 0:
                                continue
                            out.append(F)
    return out


def gl3f2_lifts():
    """one GL(3,Z) lift (entries in {-1,0,1}, det +-1) for each element of GL(3, F_2)"""
    lifts = {}
    for entries in itertools.product((-1, 0, 1), repeat=9):
        P = np.array(entries, dtype=np.int64).reshape(3, 3)
        dt = det_int(P)
        if dt not in (1, -1):
            continue
        key = tuple((P % 2).flatten())
        if key not in lifts:
            lifts[key] = P
    assert len(lifts) == 168, len(lifts)
    return list(lifts.values())


def plumbing(m, A):
    """6x6 plumbing form for chains (m_i, 2) with linking numbers A[i][j]"""
    Q = np.zeros((6, 6), dtype=np.int64)
    for i in range(3):
        Q[2 * i, 2 * i] = m[i]; Q[2 * i, 2 * i + 1] = Q[2 * i + 1, 2 * i] = 1; Q[2 * i + 1, 2 * i + 1] = 2
        for j in range(3):
            if i != j:
                Q[2 * i, 2 * j] = A[i][j]
    return Q


def candidate_plumbings(D, n_even=3):
    """Enumerate every admissible parity class of the half-integral surgery form.

    Reduction is under GL(3,Z), but the prescribed meridians require bases
    modulo Gamma(2). Lifting all GL(3,F_2) elements accounts for those bases.
    The parity of each diagonal modulo four is unchanged by Gamma(2), so it
    suffices to test one integral lift per residue class. Signed permutations
    only identify equivalent plumbing components; they remove no obstruction.
    """
    if D <= 0 or D % 2 == 0 or not 0 <= n_even <= 3:
        raise ValueError('require positive odd determinant and 0 <= n_even <= 3')
    lifts = gl3f2_lifts()
    seen = set(); out = []
    for F in reduced_ternary_forms(D):
        for P in lifts:
            Qp = P.T @ F @ P
            if any(Qp[i, i] % 2 == 0 for i in range(3)) or any(Qp[i, j] % 2 for i in range(3) for j in range(3) if i != j):
                continue
            m = [(int(Qp[i, i]) + 1) // 2 for i in range(3)]
            if sum(1 for x in m if x % 2 == 0) != n_even:
                continue
            A = [[int(Qp[i, j]) // 2 if i != j else 0 for j in range(3)] for i in range(3)]
            # canonical key up to permutations and sign changes of the components
            keys = []
            for perm in itertools.permutations(range(3)):
                for signs in itertools.product((1, -1), repeat=3):
                    mm = tuple(m[p] for p in perm)
                    aa = tuple(signs[i] * signs[j] * A[perm[i]][perm[j]] for i in range(3) for j in range(i + 1, 3))
                    keys.append((mm, aa))
            key = min(keys)
            if key in seen:
                continue
            seen.add(key)
            Q = plumbing(m, A)
            if det_int(Q) != D:
                continue
            if not positive_definite(Q):
                continue
            out.append((key, Q))
    return out


def obstruct_u3(pd, sigma, det_k, verbose=False):
    sig = abs(int(sigma))
    if sig != 6:
        return {'verdict': 'NOT_APPLICABLE', 'reason': f'|sigma| = {sig} != 6'}
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
    cands = candidate_plumbings(det_k, n_even=3)
    if verbose:
        print(f'  {len(cands)} candidate plumbing forms', file=sys.stderr)
    if not cands:
        return {'verdict': 'OBSTRUCTED', 'candidates': 0}
    # A surviving form is only a failure to obstruct. To certify OBSTRUCTED,
    # every admissible form must fail a completed exact comparison.
    undecided = 0
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
    import ast, database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    names = sys.argv[1:]
    for nm in names:
        r = rows[nm]
        res = obstruct_u3(ast.literal_eval(r['pd_notation']), int(r['signature']), int(r['determinant']), verbose=True)
        print(nm, 'u =', r['unknotting_number'], 'sigma', r['signature'], 'det', r['determinant'], '->', res, flush=True)
