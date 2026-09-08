"""Alexander-polynomial fingerprint of a knot diagram from its PD code (stdlib only).

Wirtinger presentation with one generator per over-arc; at a crossing with under-arcs a (in),
c (out), over generator o and sign eps the relation is x_c = x_o^eps x_a x_o^-eps, whose Fox
derivatives (abelianised) give the row (1-t, t, -1) for eps = +1 and (t-1, 1, -t) for eps = -1
(the latter multiplied by t).  Any (n-1)-minor is the Alexander polynomial up to +-t^m.  The
fingerprint removes that ambiguity: F(t) = Delta(t) * Delta(1/t) evaluated modulo the prime
P = 2^61 - 1 at t = 2, 3, 5, together with |Delta(-1)| = det.  Equal fingerprints are a
necessary condition for equal knot types (mirror- and orientation-insensitive); the certificates
that rely on it are confirmed locally by SnapPy isometry.
"""
from __future__ import annotations

P = (1 << 61) - 1


def _det_mod(M, n):
    A = [row[:] for row in M]
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] % P), None)
        if piv is None:
            return 0
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = -det
        det = det * A[c][c] % P
        inv = pow(A[c][c], P - 2, P)
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c] * inv % P
                A[r] = [(x - f * y) % P for x, y in zip(A[r], A[c])]
    return det % P


def alexander_rows(pd):
    """(over generator, under-in arc generator, under-out arc generator, sign) per crossing, with
    arcs merged into over-arc generators"""
    pd = [[int(x) for x in q] for q in pd]
    n = len(pd)
    m = 2 * n
    labels = sorted({x for q in pd for x in q})
    # arcs b and d of every crossing belong to the same over-arc: union-find over labels
    parent = {l: l for l in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b, c, d in pd:
        rb, rd = find(b), find(d)
        if rb != rd:
            parent[rd] = rb
    gens = sorted({find(l) for l in labels})
    gid = {g: i for i, g in enumerate(gens)}
    # orientation by tracing the knot: the under strand enters at position 0 and leaves at 2;
    # the over passage enters at position 1 (b -> d) or 3 (d -> b).  Labels need not be sequential.
    occ = {}
    for c, q in enumerate(pd):
        for i, l in enumerate(q):
            occ.setdefault(l, []).append((c, i))
    partner = {}
    for l, ds in occ.items():
        (c1, i1), (c2, i2) = ds
        partner[(c1, i1)] = (c2, i2)
        partner[(c2, i2)] = (c1, i1)
    over_in = {}
    cur = (0, 0)
    for _ in range(m):
        c, i = cur
        if i % 2 == 1:
            over_in[c] = i
        cur = partner[(c, (i + 2) % 4)]
    rows = []
    for c, (a, b, cc, d) in enumerate(pd):
        pos = over_in[c] == 3                 # over strand runs d -> b: positive crossing
        rows.append((gid[find(b)], gid[find(a)], gid[find(cc)], 1 if pos else -1))
    return rows, len(gens)


def alexander_eval(pd, t):
    """Delta(t) mod P up to +-t^m (a unit); returns (value, sign/power-free? no) -> raw minor"""
    rows, ng = alexander_rows(pd)
    n = len(rows)
    M = [[0] * ng for _ in range(n)]
    for r, (o, a, c, eps) in enumerate(rows):
        if eps == 1:
            M[r][o] = (M[r][o] + 1 - t) % P
            M[r][a] = (M[r][a] + t) % P
            M[r][c] = (M[r][c] - 1) % P
        else:
            M[r][o] = (M[r][o] + t - 1) % P
            M[r][a] = (M[r][a] + 1) % P
            M[r][c] = (M[r][c] - t) % P
    # delete the last row and the last column
    minor = [row[:-1] for row in M[:-1]]
    return _det_mod(minor, n - 1)


def fingerprint(pd):
    out = []
    for t in (2, 3, 5):
        tinv = pow(t, P - 2, P)
        out.append(alexander_eval(pd, t) * alexander_eval(pd, tinv) % P)
    return tuple(out)


def det_from_alexander(pd):
    v = alexander_eval(pd, P - 1)             # t = -1
    return v if v < P // 2 else P - v


UNKNOT_FP = (1, 1, 1)
