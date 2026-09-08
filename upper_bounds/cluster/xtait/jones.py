"""Exact Jones polynomial of a knot PD via the Kauffman bracket state sum (stdlib only).

O(2^n * n) — intended for verification (selftest oracles, jackpot sanity checks), not for search.
Default hard cap n <= 24 (~40 s at 20 crossings in CPython).

Conventions: V(unknot) = 1, variable t, V(right-handed trefoil per the standard PD) matching the KnotInfo
table (calibrated in selftest against the shipped summand data)."""
from __future__ import annotations

from .graph import pd_passages


class _DSU:
    __slots__ = ('p',)

    def __init__(self, items):
        self.p = {x: x for x in items}

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def bracket(pd) -> dict:
    """Kauffman bracket <D> as {A-exponent: coefficient}."""
    labels = {int(x) for q in pd for x in q}
    n = len(pd)
    if n > 24:
        raise ValueError(f'bracket state sum limited to 24 crossings (got {n})')
    total = {}
    for mask in range(1 << n):
        dsu = _DSU(labels)
        a_cnt = 0
        for i, (a, b, c, d) in enumerate(pd):
            if (mask >> i) & 1 == 0:
                a_cnt += 1
                dsu.union(a, b)
                dsu.union(c, d)
            else:
                dsu.union(b, c)
                dsu.union(d, a)
        loops = len({dsu.find(x) for x in labels})
        exp = 2 * a_cnt - n            # A^{a} B^{n-a}, B = A^{-1}
        coeff = 1
        # (-A^2 - A^-2)^(loops-1): expand
        # accumulate into total: coeff * sum_k C(loops-1,k) (-1)^(loops-1) A^{2(loops-1) - 4k}
        m = loops - 1
        sign = -1 if m % 2 else 1
        ck = 1
        for k in range(m + 1):
            E = exp + 2 * m - 4 * k
            total[E] = total.get(E, 0) + sign * ck
            ck = ck * (m - k) // (k + 1)
    return {e: c for e, c in total.items() if c}


def writhe(pd) -> int:
    """Writhe of the (single-component) PD.  At each crossing the under passage enters at position 0 or 2
    and the over passage at 1 or 3; the crossing is positive iff over_in == (under_in + 3) mod 4.
    (Calibrated so that the normalised bracket equals the KnotInfo Jones polynomial; see selftest.)"""
    under_in = {}
    over_in = {}
    for (c, i) in pd_passages(pd):
        if i % 2 == 0:
            under_in[c] = i
        else:
            over_in[c] = i
    w = 0
    for c in range(len(pd)):
        w += 1 if over_in[c] == (under_in[c] + 3) % 4 else -1
    return w


def jones(pd) -> dict:
    """Jones polynomial {t-exponent: coefficient} (exact)."""
    if not pd:
        return {0: 1}
    br = bracket(pd)
    w = writhe(pd)
    sign = -1 if (3 * w) % 2 else 1
    poly = {}
    for aexp, c in br.items():
        E = aexp - 3 * w             # multiply by (-A)^{-3w} = (-1)^{3w} A^{-3w}
        if E % 4:
            raise ValueError('normalised bracket exponent not divisible by 4')
        t = -E // 4                   # t = A^{-4}
        poly[t] = poly.get(t, 0) + sign * c
    poly = {e: c for e, c in poly.items() if c}
    if sum(poly.values()) not in (1, -1):
        # V(1) = (-2)^{#components-1}; for a knot it must be ±1 (sign checks the normalisation)
        raise ValueError('normalisation failed: V(1) != ±1')
    return poly


def jones_mirror(poly: dict) -> dict:
    return {-e: c for e, c in poly.items()}


def jones_mul(p1: dict, p2: dict) -> dict:
    out = {}
    for e1, c1 in p1.items():
        for e2, c2 in p2.items():
            out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
    return {e: c for e, c in out.items() if c}


def jones_span(poly: dict) -> int:
    return max(poly) - min(poly) if poly else 0


def determinant_from_jones(poly: dict) -> int:
    return abs(sum(c * (-1) ** (e % 2) for e, c in poly.items()))
