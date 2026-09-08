"""Donaldson lattice-embedding refinement of Owens' u = 2 obstruction (alternating K, |sigma| = 4).

If u(K) = 2 with both changes of the same sign, the positive-definite Goeritz lattice L_G (the side used in
owens_obstruction, d(spin) = -sigma/4) and the plumbing lattice Q~ = Q~(m1, m2, a) (m1, m2 even, det = det K)
embed orthogonally into Z^N, N = rank G + 4, as mutually orthogonal sublattices (Donaldson's theorem applied
to the closed positive-definite manifold obtained by gluing the two fillings).  Hence:
    some embedding L_G -> Z^N has orthogonal complement isometric to a plumbing Q~(m1, m2, a) with m1, m2 even,
    and that Q~ passes the d-invariant matching.
The script enumerates all embeddings of L_G up to signed permutations (backtracking, most-constrained vertex
first), computes each orthogonal complement, tests isometry with Owens' candidate list, and reports
PASS / OBSTRUCTED / UNDECIDED (time or count cap).  Complements of non-primitive embeddings are reported
separately and never used to obstruct.
"""
from __future__ import annotations
import itertools, math, sys, time
from fractions import Fraction
import numpy as np
from owens_obstruction import definite_goeritz_candidates, det_int, m_Q, matching_exists, candidates, qtilde, adjugate_int


def hermite_kernel(V):
    """integer basis (columns) of {x in Z^N : V x = 0} for an integer r x N matrix V of rank r"""
    V = [list(map(int, row)) for row in V]
    r, N = len(V), len(V[0])
    U = [[1 if i == j else 0 for j in range(N)] for i in range(N)]     # column operations accumulated
    A = [row[:] for row in V]
    col = 0
    for i in range(r):
        # bring the gcd of A[i][col:] to position (i, col) by column operations
        while True:
            nz = [c for c in range(col, N) if A[i][c] != 0]
            if not nz:
                break
            p = min(nz, key=lambda c: abs(A[i][c]))
            if p != col:
                for row in A: row[col], row[p] = row[p], row[col]
                for row in U: row[col], row[p] = row[p], row[col]
            done = True
            for c in range(col + 1, N):
                if A[i][c]:
                    q = A[i][c] // A[i][col]
                    for row in A: row[c] -= q * row[col]
                    for row in U: row[c] -= q * row[col]
                    if A[i][c]: done = False
            if done:
                break
        if any(A[i][c] for c in range(col, N)):
            col += 1
    K = [[U[n][c] for c in range(col, N)] for n in range(N)]     # N x (N - col)
    return np.array(K, dtype=np.int64)


def embeddings(G, N, time_limit, max_embeddings):
    """all r-tuples of vectors in Z^N with Gram matrix G, up to signed permutations of Z^N"""
    G = np.array(G, dtype=np.int64); r = len(G)
    # order: greedy most-constrained (largest number of edges to already placed vertices)
    order = [int(np.argmax(np.diag(G)))]
    while len(order) < r:
        best = max((i for i in range(r) if i not in order), key=lambda i: (sum(abs(int(G[i, j])) for j in order), int(G[i, i])))
        order.append(best)
    vecs = []           # placed vectors (in Z^N) in 'order'
    out = []
    t0 = time.time()
    state = {'stop': False}

    def rec(k, used):
        if state['stop']:
            return
        if k == r:
            out.append([v.copy() for v in vecs])
            if len(out) >= max_embeddings: state['stop'] = True
            return
        if time.time() - t0 > time_limit:
            state['stop'] = True; return
        i = order[k]; norm = int(G[i, i])
        targets = [int(G[i, order[j]]) for j in range(k)]
        ucoords = sorted(used); fcoords = [c for c in range(N) if c not in used]
        x = np.zeros(N, dtype=np.int64)
        # enumerate the part on used coordinates with incremental pruning
        def rec_used(idx, rem_norm, sums):
            if state['stop']:
                return
            if idx == len(ucoords):
                if any(sums[j] != targets[j] for j in range(k)):
                    return
                # free part: multiset of positive integers with squares summing to rem_norm, non-increasing
                def parts(rem, maxv, acc):
                    if rem == 0:
                        yield acc; return
                    for a in range(min(maxv, math.isqrt(rem)), 0, -1):
                        yield from parts(rem - a * a, a, acc + [a])
                for p in parts(rem_norm, math.isqrt(rem_norm) if rem_norm else 0, []):
                    if len(p) > len(fcoords):
                        continue
                    y = x.copy()
                    for a, c in zip(p, fcoords): y[c] = a
                    vecs.append(y)
                    rec(k + 1, used | {c for c in range(N) if y[c] != 0})
                    vecs.pop()
                return
            c = ucoords[idx]
            # remaining capacity bounds for each linear constraint
            rest = ucoords[idx + 1:]
            B = math.isqrt(rem_norm)
            for val in range(-B, B + 1):
                new = [sums[j] + val * int(vecs[j][c]) for j in range(k)]
                ok = True
                for j in range(k):
                    cap = sum(abs(int(vecs[j][cc])) for cc in rest) * math.isqrt(max(rem_norm - val * val, 0))
                    if abs(targets[j] - new[j]) > cap:
                        ok = False; break
                if not ok:
                    continue
                x[c] = val
                rec_used(idx + 1, rem_norm - val * val, new)
                x[c] = 0
        rec_used(0, norm, [0] * k)

    rec(0, set())
    return out, state['stop'], order


def short_vectors(C, bound):
    """all vectors (as coefficient rows) of a positive-definite integer form C with norm <= bound, up to sign"""
    C = np.array(C, dtype=float); n = len(C)
    L = np.linalg.cholesky(C)          # C = L L^T ; norm(x) = |L^T x|^2
    R = np.linalg.qr(L.T, mode='r')
    for i in range(n):
        if R[i, i] < 0: R[i, :] = -R[i, :]
    out = []
    x = np.zeros(n)
    def rec(i, partial):
        if i < 0:
            v = x.astype(np.int64)
            if np.any(v):
                # up to sign: first nonzero positive
                fz = next(k for k in range(n) if v[k] != 0)
                if v[fz] > 0: out.append(v.copy())
            return
        # solve for x[i]: contribution (R[i,i] x_i + sum_{j>i} R[i,j] x_j)^2 <= bound - partial
        s = sum(R[i, j] * x[j] for j in range(i + 1, n))
        rad = math.sqrt(max(bound - partial, 0) + 1e-9)
        lo = math.ceil((-rad - s) / R[i, i] - 1e-9); hi = math.floor((rad - s) / R[i, i] + 1e-9)
        for v in range(lo, hi + 1):
            x[i] = v
            t = (R[i, i] * v + s) ** 2
            if partial + t <= bound + 1e-9:
                rec(i - 1, partial + t)
        x[i] = 0
    rec(n - 1, 0.0)
    return out


def is_plumbing(C, cands):
    """C (4x4 Gram) isometric to Q~(m1,m2,a) for some (m1,m2,a) in cands?  Returns the triple or None."""
    C = np.array(C, dtype=np.int64)
    ip = lambda u, v: int(u @ C @ v)
    two = [v for v in short_vectors(C, 2) if ip(v, v) == 2]
    if len(two) < 2:
        return None
    maxm = max(max(m1, m2) for m1, m2, a in cands)
    longv = {}
    def vectors_of_norm(m):
        if m not in longv:
            longv[m] = [v for v in short_vectors(C, m) if ip(v, v) == m] if m <= maxm else []
        return longv[m]
    for (m1, m2, a) in cands:
        for f1, f2 in itertools.permutations(two, 2):
            if ip(f1, f2) != 0:
                continue
            for e1 in vectors_of_norm(m1):
                for s1 in (1, -1):
                    E1 = s1 * e1
                    if ip(E1, f1) != 1 or ip(E1, f2) != 0:
                        continue
                    for e2 in vectors_of_norm(m2):
                        for s2 in (1, -1):
                            E2 = s2 * e2
                            if ip(E2, f2) != 1 or ip(E2, f1) != 0 or ip(E1, E2) != a:
                                continue
                            B = np.array([E1, f1, E2, f2])
                            if abs(int(round(np.linalg.det(B.astype(float))))) == 1:
                                return (m1, m2, a)
    return None


def donaldson_u2(pd, sigma, det_k, time_limit=300.0, max_embeddings=200000, verbose=False):
    sig = abs(int(sigma))
    if sig != 4:
        return {'verdict': 'NOT_APPLICABLE'}
    forms = definite_goeritz_candidates(pd); target0 = Fraction(-sig, 4); chosen = None
    for G in forms:
        if det_int(G) != det_k: continue
        mG, cmG = m_Q(G)
        if mG.get(cmG.coords([0] * len(G))) == target0:
            chosen = (G, mG, cmG); break
    if chosen is None:
        return {'verdict': 'ERROR', 'reason': 'no Goeritz side'}
    G, mG, cmG = chosen
    # Donaldson gluing: the plumbing X (positive-definite filling of Sigma_2) is glued to the OTHER checkerboard
    # Goeritz filling, the positive-definite filling of -Sigma_2; its form G2 is the other definite candidate.
    others = [g for g in forms if det_int(g) == det_k and not (g.shape == G.shape and np.array_equal(g, G))]
    if not others:
        return {'verdict': 'ERROR', 'reason': 'other Goeritz side not found'}
    G2 = others[0]
    cands = candidates(det_k, 2)
    dmatch = {}
    for c in cands:
        mQ, cmQ = m_Q(qtilde(*c)); dmatch[c] = matching_exists(mQ, cmQ, mG, cmG)
    good = [c for c in cands if dmatch[c] is True]
    if not good:
        return {'verdict': 'OBSTRUCTED', 'stage': 'd-invariants', 'candidates': len(cands)}
    r = len(G2); N = r + 4
    t0 = time.time()
    embs, capped, order = embeddings(G2, N, time_limit, max_embeddings)
    n_prim = n_nonprim = 0; found = None
    for vecs in embs:
        V = np.array(vecs, dtype=np.int64)
        K = hermite_kernel(V)
        if K.shape[1] != 4:
            continue
        C = K.T @ K
        dC = det_int(C)
        if dC != det_k:
            n_nonprim += 1; continue
        n_prim += 1
        hit = is_plumbing(C, good)
        if hit is not None:
            found = hit; break
    res = {'rank_G': r, 'N': N, 'embeddings': len(embs), 'primitive': n_prim, 'nonprimitive': n_nonprim,
           'd_pass_candidates': len(good), 'seconds': round(time.time() - t0, 1)}
    if found is not None:
        res.update({'verdict': 'PASS', 'witness_form': found})
    elif capped:
        res.update({'verdict': 'UNDECIDED', 'reason': 'embedding enumeration capped'})
    elif n_nonprim and not n_prim:
        res.update({'verdict': 'UNDECIDED', 'reason': 'only non-primitive embeddings'})
    else:
        res.update({'verdict': 'OBSTRUCTED', 'stage': 'Donaldson'})
    return res


if __name__ == '__main__':
    import ast, json, database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for nm in sys.argv[1:]:
        r = rows[nm]
        res = donaldson_u2(ast.literal_eval(r['pd_notation']), int(r['signature']), int(r['determinant']))
        print(json.dumps({'name': nm, 'u': str(r['unknotting_number']).strip(), 'det': int(r['determinant']), **res}), flush=True)
