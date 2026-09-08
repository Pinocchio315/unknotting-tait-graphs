#!/usr/bin/env python
"""Independent cross-check of the Owens u=2 obstruction (|sigma| = 4) for specific knots.

Deliberately shares NO algorithmic code with owens_obstruction.py:
  * checkerboard graphs built directly from spherogram's face list (no Tait-graph machinery),
  * m_Q by full box enumeration over characteristic vectors (-Q_ii <= xi_i <= Q_ii - 2)  [OS bound],
  * homology classes as adjugate labels in (Z/D)^r (no Smith normal form),
  * isomorphisms enumerated from generator images with brute-force bijectivity check.
Checks, per knot: det(G) = det K for both checkerboard forms, d-multisets of the two forms are exact
negatives, side selection d(spin) = -sigma/4, and the final verdict.

    python crosscheck.py 11a_63 13a_16 9_10 8_5
"""
import sys
from fractions import Fraction
from itertools import product

sys.path.insert(0, '.')
import _compat  # noqa
import kinfo


def checkerboard_laplacians(pd):
    from spherogram import Link
    L = Link([list(q) for q in pd])
    faces = L.faces()
    corner_face = {}
    for fi, f in enumerate(faces):
        for cs in f:
            corner_face[(cs.crossing, cs.strand_index)] = fi
    # two faces sharing an arc are adjacent; colour by BFS bipartition of the arc-adjacency graph
    import collections
    adj = collections.defaultdict(set)
    for c in L.crossings:
        for i in range(4):
            # the arc leaving corner (c, i) also appears as the opposite corner on the neighbour
            d, j = c.adjacent[i]
            f1 = corner_face[(c, i)]             # one side of the arc from (c, i) to (d, j)
            f2 = corner_face[(d, j)]             # the other side (the face walked from the far end)
            if f1 != f2:
                adj[f1].add(f2)
                adj[f2].add(f1)
    colour = {0: 0}
    queue = [0]
    while queue:
        x = queue.pop()
        for y in adj[x]:
            if y not in colour:
                colour[y] = 1 - colour[x]
                queue.append(y)
    assert len(colour) == len(faces), 'face adjacency not connected'
    # crossings connect the two same-colour faces diagonally: corners (c,0)&(c,2) / (c,1)&(c,3)
    forms = []
    for col in (0, 1):
        verts = sorted(f for f in range(len(faces)) if colour[f] == col)
        idx = {v: i for i, v in enumerate(verts)}
        n = len(verts)
        M = [[0] * n for _ in range(n)]
        for c in L.crossings:
            pair0 = (corner_face[(c, 0)], corner_face[(c, 2)])
            pair1 = (corner_face[(c, 1)], corner_face[(c, 3)])
            pair = pair0 if colour[pair0[0]] == col else pair1
            assert colour[pair[0]] == col and colour[pair[1]] == col
            u, w = idx[pair[0]], idx[pair[1]]
            if u == w:
                continue
            M[u][u] += 1
            M[w][w] += 1
            M[u][w] -= 1
            M[w][u] -= 1
        forms.append([row[:-1] for row in M[:-1]])
    return forms


def det_frac(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    d = Fraction(1)
    for c in range(n):
        p = next((r for r in range(c, n) if A[r][c] != 0), None)
        if p is None:
            return 0
        if p != c:
            A[c], A[p] = A[p], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return int(d)


def inv_frac(M):
    n = len(M)
    A = [[Fraction(x) for x in row] + [Fraction(1 if i == j else 0) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [x / piv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return [row[n:] for row in A]


def m_and_labels(Q):
    """Full box enumeration: returns dict label -> min (xi Q^-1 xi - r)/4, labels = tuple(adj @ xi mod D)."""
    n = len(Q)
    D = abs(det_frac(Q))
    Qi = inv_frac(Q)
    adj = [[int(Qi[i][j] * D) for j in range(n)] for i in range(n)]
    ranges = []
    for i in range(n):
        lo, hi = -Q[i][i], Q[i][i] - 2
        ranges.append([v for v in range(lo, hi + 1) if (v - Q[i][i]) % 2 == 0])
    best = {}
    for xi in product(*ranges):
        f = Fraction(0)
        for i in range(n):
            s = sum(adj[i][j] * xi[j] for j in range(n))
            f += Fraction(xi[i] * s, D)
        lbl = tuple(sum(adj[i][j] * xi[j] for j in range(n)) % D for i in range(n))
        val = (f - n) / 4
        if lbl not in best or val < best[lbl]:
            best[lbl] = val
    return best, D


def group_elems(best):
    return list(best.keys())


def add(x, y, D):
    return tuple((a + b) % D for a, b in zip(x, y))


def scalar(k, x, D):
    return tuple((k * a) % D for a in x)


def span(gens, D, order_target):
    seen = {tuple([0] * len(gens[0]))}
    frontier = [tuple([0] * len(gens[0]))]
    while frontier:
        x = frontier.pop()
        for g in gens:
            y = add(x, g, D)
            if y not in seen:
                seen.add(y)
                frontier.append(y)
                if len(seen) > order_target:
                    return seen
    return seen


def isomorphisms(A_elems, DA, B_elems, DB):
    """Brute force: generators of A = up to 2 elements spanning A; map to candidate images in B."""
    n = len(A_elems)
    zeroA = tuple([0] * len(A_elems[0]))
    # find a generating pair of A
    gens = None
    for g1 in A_elems:
        if len(span([g1], DA, n)) == n:
            gens = [g1]
            break
    if gens is None:
        for g1 in A_elems:
            for g2 in A_elems:
                if len(span([g1, g2], DA, n)) == n:
                    gens = [g1, g2]
                    break
            if gens:
                break
    assert gens, 'need <= 2 generators'
    def order(x, D, elems_len):
        k, y = 1, x
        zero = tuple([0]*len(x))
        while y != zero:
            y = add(y, x, D)
            k += 1
        return k
    for imgs in product(B_elems, repeat=len(gens)):
        if any(order(g, DA, n) != order(h, DB, n) for g, h in zip(gens, imgs)):
            continue
        # build the map by spanning
        m = {zeroA: tuple([0] * len(B_elems[0]))}
        frontier = [zeroA]
        ok = True
        while frontier and ok:
            x = frontier.pop()
            for g, h in zip(gens, imgs):
                y = add(x, g, DA)
                fy = add(m[x], h, DB)
                if y in m:
                    if m[y] != fy:
                        ok = False
                        break
                else:
                    m[y] = fy
                    frontier.append(y)
        if ok and len(m) == n and len(set(m.values())) == n:
            yield m


def check(nm):
    r = kinfo.row(nm)
    sigma = abs(int(r['signature']))
    detk = int(r['determinant'])
    assert sigma == 4, 'crosscheck covers the |sigma| = 4 branch'
    pd = kinfo.pd_code(nm)
    F1, F2 = checkerboard_laplacians(pd)
    d1, d2 = abs(det_frac(F1)), abs(det_frac(F2))
    assert d1 == detk and d2 == detk, ('det mismatch', d1, d2, detk)
    m1, D = m_and_labels(F1)
    m2, _ = m_and_labels(F2)
    ms1 = sorted(m1.values())
    ms2 = sorted((-v for v in m2.values()))
    assert ms1 == ms2, 'the two checkerboard d-multisets are not exact negatives'
    z1 = tuple([0] * len(F1))
    z2 = tuple([0] * len(F2))
    target = Fraction(-sigma, 4)
    if m1[z1] == target:
        mg, DG = m1, D
    elif m2[z2] == target:
        mg, DG = m2, D
    else:
        raise AssertionError('no side with d(spin) = -sigma/4')
    # candidates by brute force
    cands = []
    m1v = 1
    while (2 * m1v - 1) ** 2 <= detk + 4 * (m1v - 1) ** 2:
        for a in range(m1v):
            num = detk + 4 * a * a
            if num % (2 * m1v - 1) == 0 and (num // (2 * m1v - 1)) % 2 == 1:
                m2v = (num // (2 * m1v - 1) + 1) // 2
                if m2v >= m1v and m1v % 2 == 0 and m2v % 2 == 0:
                    cands.append((m1v, m2v, a))
        m1v += 1
    verdict = 'OBSTRUCTED'
    for (a1, a2, aa) in cands:
        Qt = [[a1, 1, aa, 0], [1, 2, 0, 0], [aa, 0, a2, 1], [0, 0, 1, 2]]
        assert abs(det_frac(Qt)) == detk
        mq, DQ = m_and_labels(Qt)
        A = group_elems(mq)
        B = group_elems(mg)
        found = False
        for iso in isomorphisms(A, DQ, B, DG):
            if all((mq[g] - mg[iso[g]]) >= 0 and (mq[g] - mg[iso[g]]).denominator == 1
                   and (mq[g] - mg[iso[g]]).numerator % 2 == 0 for g in A):
                found = True
                break
        if found:
            verdict = 'PASS'
            break
    print(f'{nm}: det {detk}, sigma -{sigma}, candidates {len(cands)} -> {verdict}   '
          f'[checks: det x2 OK, +-d-multiset OK, spin side OK]')
    return verdict


if __name__ == '__main__':
    for nm in (sys.argv[1:] or ['9_10', '8_5', '11a_63', '13a_16']):
        check(nm)
