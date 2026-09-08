"""Owens' u = 2 obstruction strengthened by the Donaldson embedding WITH the gluing-induced identification.

Alternating K, |sigma| = 4.  Let Gamma be one checkerboard Tait graph of the diagram and Gamma* the other;
the cut lattice of Gamma and the cut lattice of Gamma* (= cycle lattice of Gamma) are orthogonal complements
in Z^c (c = crossings), which encodes the double cover of B^4 over the two checkerboard surfaces and gives the
geometric identification psi : disc(G2) -> disc(G1) between the two Goeritz discriminant groups (= H_1(Sigma_2)).
If u(K) = 2 (same-sign changes) the plumbing X (form Q~(m1,m2,a), m_i even) glues to the Goeritz filling of
-Sigma_2, so Q~ and G2 embed orthogonally in Z^N, N = rank G2 + 4, and the gluing gives phi2 : disc(Q~) -> disc(G2).
Owens' d-invariant inequalities must then hold for the SPECIFIC identification phi = psi o phi2:
    m_Q~(g) >= m_G1(phi(g)),  m_Q~(g) - m_G1(phi(g)) in 2Z   for all g.
Verdict PASS if some embedding satisfies this; OBSTRUCTED if none does (u >= 3); UNDECIDED if capped.
"""
from __future__ import annotations
import itertools, math, sys, time, collections
from fractions import Fraction
import numpy as np
sys.path.insert(0, '/Users/pinocchio/Documents/05_unknotting_number/0001_writing_with_pavel/server_code')
from xtait.graph import from_pd, ribbon_faces
from owens_obstruction import det_int, m_Q, ClassMap, candidates, qtilde
from donaldson_u2 import hermite_kernel, embeddings, short_vectors


def tait_data(pd):
    """cut vectors (vertices of Gamma) and face vectors (vertices of Gamma*) in Z^c, both Laplacians"""
    g = from_pd([list(q) for q in pd])
    edges = sorted(g.edges); eidx = {e: i for i, e in enumerate(edges)}; c = len(edges)
    verts = list(g.vertices)
    cut = {}
    for w in verts:
        v = np.zeros(c, dtype=np.int64)
        for e in edges:
            u, x = g.edges[e]
            if u == w: v[eidx[e]] += 1
            if x == w: v[eidx[e]] -= 1
        cut[w] = v
    faces = ribbon_faces(g)
    def face_vec(f, conv):
        v = np.zeros(c, dtype=np.int64)
        for d in f:
            vtx, e = d[0], d[1]
            u, x = g.edges[e]
            s = 1 if vtx == u else -1
            if u == x: s = 0
            v[eidx[e]] += s * conv
        return v
    for conv in (1, -1):
        fv = [face_vec(f, conv) for f in faces]
        ok = all(int(cut[w] @ f) == 0 for w in verts for f in fv)
        if ok:
            break
    if not ok:
        raise RuntimeError('cut/face orthogonality failed')
    C1 = np.array([cut[w] for w in verts[:-1]], dtype=np.int64)      # reduced: delete last vertex
    F2 = np.array(fv[:-1], dtype=np.int64)                            # reduced: delete last face
    return C1, F2


def glue_map(A_rows, B_rows, GA, GB):
    """A_rows (rA x M), B_rows (rB x M): orthogonal integer lattices in Z^M with Gram GA, GB.
    Returns the isomorphism disc(A) -> disc(B) induced by Z^M (dict coords->coords), or None."""
    cmA, cmB = ClassMap(GA), ClassMap(GB)
    if cmA.order != cmB.order:
        return None
    M = A_rows.shape[1]
    gens = []
    for n in range(M):
        gens.append((cmA.coords(A_rows[:, n]), cmB.coords(B_rows[:, n])))
    def add(x, y, d): return tuple((a + b) % m for a, b, m in zip(x, y, d))
    zero_a = tuple(0 for _ in cmA.d); zero_b = tuple(0 for _ in cmB.d)
    phi = {zero_a: zero_b}; queue = collections.deque([zero_a])
    while queue:
        a = queue.popleft(); b = phi[a]
        for pa, pb in gens:
            a2 = add(a, pa, cmA.d); b2 = add(b, pb, cmB.d)
            if a2 in phi:
                if phi[a2] != b2:
                    return None
            else:
                phi[a2] = b2; queue.append(a2)
    if len(phi) != cmA.order or len(set(phi.values())) != cmB.order:
        return None
    return phi, cmA, cmB


def plumbing_bases(C, cands):
    """all (triple, basis rows B (4x4 in the C-basis)) with B C B^T = qtilde(triple)"""
    C = np.array(C, dtype=np.int64)
    ip = lambda u, v: int(u @ C @ v)
    two = [v for v in short_vectors(C, 2) if ip(v, v) == 2]
    out = []
    if len(two) < 2:
        return out
    cache = {}
    def vecs(m):
        if m not in cache:
            cache[m] = [v for v in short_vectors(C, m) if ip(v, v) == m]
        return cache[m]
    for (m1, m2, a) in cands:
        for f1, f2 in itertools.permutations(two, 2):
            if ip(f1, f2) != 0: continue
            for sf1 in (1, -1):
                F1 = sf1 * f1
                for sf2 in (1, -1):
                    F2 = sf2 * f2
                    for e1 in vecs(m1):
                        for s1 in (1, -1):
                            E1 = s1 * e1
                            if ip(E1, F1) != 1 or ip(E1, F2) != 0: continue
                            for e2 in vecs(m2):
                                for s2 in (1, -1):
                                    E2 = s2 * e2
                                    if ip(E2, F2) != 1 or ip(E2, F1) != 0 or ip(E1, E2) != a: continue
                                    B = np.array([E1, F1, E2, F2], dtype=np.int64)
                                    if abs(int(round(np.linalg.det(B.astype(float))))) == 1:
                                        out.append(((m1, m2, a), B))
    return out


def obstruct(pd, sigma, det_k, time_limit=300.0, max_embeddings=200000):
    sig = abs(int(sigma))
    if sig != 4:
        return {'verdict': 'NOT_APPLICABLE'}
    C1, F2 = tait_data(pd)
    L1 = C1 @ C1.T; L2 = F2 @ F2.T
    if det_int(L1) != det_k or det_int(L2) != det_k:
        return {'verdict': 'ERROR', 'reason': f'Laplacian dets {det_int(L1)}, {det_int(L2)} != {det_k}'}
    # sharp side: m(0) = -sigma/4
    sides = []
    for rows_, L in ((C1, L1), (F2, L2)):
        mG, cmG = m_Q(L); sides.append((rows_, L, mG, cmG, mG.get(cmG.coords([0] * len(L)))))
    s1 = [s for s in sides if s[4] == Fraction(-sig, 4)]
    if not s1:
        return {'verdict': 'ERROR', 'reason': 'no sharp side'}
    rows1, G1, mG1, cmG1, _ = s1[0]
    rows2, G2, mG2, cmG2, _ = [s for s in sides if s is not s1[0]][0]
    # psi: disc(G2) -> disc(G1) from the Z^c gluing
    gm = glue_map(rows2, rows1, G2, G1)
    if gm is None:
        return {'verdict': 'ERROR', 'reason': 'Z^c gluing is not a graph'}
    psi, cm2, cm1 = gm
    cands = candidates(det_k, 2)
    mQ = {c: m_Q(qtilde(*c)) for c in cands}
    r = len(G2); N = r + 4
    t0 = time.time()
    embs, capped, order = embeddings(G2, N, time_limit, max_embeddings)
    n_prim = n_nonprim = 0; n_plumb = 0; found = None; elems = None
    for vecs in embs:
        V = np.array(vecs, dtype=np.int64)          # rows: images of G2's basis in the 'order' permutation
        # G2 basis order used in embeddings(): 'order'; reorder V to the original basis order
        Vo = np.zeros_like(V); 
        for k, i in enumerate(order): Vo[i] = V[k]
        K = hermite_kernel(Vo)
        if K.shape[1] != 4: continue
        C = K.T @ K
        if det_int(C) != det_k:
            n_nonprim += 1; continue
        n_prim += 1
        for triple, B in plumbing_bases(C, cands):
            n_plumb += 1
            P = (K @ B.T).T                          # 4 x N: plumbing basis vectors in Z^N, Gram = qtilde(triple)
            Qt = qtilde(*triple)
            gm2 = glue_map(P, Vo, Qt, G2)            # phi2: disc(Q~) -> disc(G2)
            if gm2 is None: continue
            phi2, cmQ, _ = gm2
            mQt, cmQt = mQ[triple]
            ok = True
            for g, b in phi2.items():
                a = mQt.get(g); c = mG1.get(psi[b])
                if a is None or c is None: ok = False; break
                dif = a - c
                if dif < 0 or dif.denominator != 1 or dif.numerator % 2: ok = False; break
            if ok:
                found = triple; break
        if found: break
    res = {'rank_G2': r, 'N': N, 'embeddings': len(embs), 'primitive': n_prim, 'nonprimitive': n_nonprim,
           'plumbing_bases_tested': n_plumb, 'seconds': round(time.time() - t0, 1)}
    if found is not None: res.update({'verdict': 'PASS', 'witness_form': found})
    elif capped: res.update({'verdict': 'UNDECIDED', 'reason': 'capped'})
    elif n_prim == 0 and n_nonprim: res.update({'verdict': 'UNDECIDED', 'reason': 'only non-primitive embeddings'})
    else: res.update({'verdict': 'OBSTRUCTED'})
    return res


if __name__ == '__main__':
    import ast, json, database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for nm in sys.argv[1:]:
        r = rows[nm]
        try:
            res = obstruct(ast.literal_eval(r['pd_notation']), int(r['signature']), int(r['determinant']))
        except Exception as e:
            res = {'verdict': 'ERROR', 'reason': repr(e)[:120]}
        print(json.dumps({'name': nm, 'u': str(r['unknotting_number']).strip(), 'det': int(r['determinant']), **res}), flush=True)
