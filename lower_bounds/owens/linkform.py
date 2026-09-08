"""Cheap necessary condition for the Owens-type matching: the linking forms of the two fillings agree.

For a positive-definite integral form Q the boundary linking form on H = Z^n / Q Z^n is
b(x, y) = x^T adj(Q) y / det(Q)  (mod 1).  For every odd prime p whose p-part of H is cyclic, Z/p^k,
the form restricted to that part is <u / p^k> and its isomorphism class is the Legendre symbol (u/p).
Two forms that must be isomorphic (same oriented 3-manifold, both fillings positive definite) therefore
have equal invariant factors and equal Legendre symbols at every such prime.  Non-cyclic p-parts are not
compared (no pruning), so the filter is always safe.
"""
from __future__ import annotations
import math
import numpy as np
from owens_obstruction import det_int, adjugate_int, ClassMap, _prime_powers, _integral_symmetric_matrix


def _order(cm, x):
    c = cm.coords(x)
    o = 1
    for v, d in zip(c, cm.d):
        if v:
            k = d // math.gcd(v, d)
            o = o * k // math.gcd(o, k)
    return o


def linking_invariants(Q):
    """dict: p -> (k, legendre) for the cyclic p-parts of the discriminant group of Q; plus 'factors'."""
    Q = _integral_symmetric_matrix(Q)
    n = len(Q)
    D = det_int(Q)
    adj = adjugate_int(Q)
    cm = ClassMap(Q)
    factors = sorted(cm.d)
    out = {'factors': tuple(factors)}
    # p-parts: cyclic iff exactly one invariant factor is divisible by p
    primes = sorted({p for d in cm.d for p in _prime_powers(d) for p in [_smallest_prime(p)]})
    for p in primes:
        if p == 2:
            continue
        parts = [d for d in cm.d if d % p == 0]
        if len(parts) != 1:
            continue
        pk = 1
        while parts[0] % (pk * p) == 0:
            pk *= p
        # generator of the p-part: (D/pk) * x for a vector x whose class has order divisible by pk
        # The coordinate classes generate H. If this p-part is cyclic, at
        # least one coordinate projects to an element of full p-power order;
        # otherwise all would lie in the proper subgroup pH. Thus a scan of
        # the standard basis is complete and needs no random trial budget.
        for i in range(n):
            y = np.zeros(n, dtype=object)
            y[i] = D // pk
            if _order(cm, y) == pk:
                val = int(y @ adj @ y)          # b(y,y) = val / D
                # reduce val/D to u/pk
                num, den = val % D, D
                g = math.gcd(num, den) if num else den
                num, den = num // g, den // g
                if den != pk:
                    continue
                u = num % p
                out[p] = (pk, 1 if pow(u, (p - 1) // 2, p) == 1 else -1)
                break
    return out


def _smallest_prime(q):
    for p in range(2, q + 1):
        if q % p == 0:
            return p
    return q


def compatible(inv_a, inv_b):
    """Return False only when recorded necessary linking invariants conflict.

    Noncyclic p-primary pairings and the 2-primary part are not classified
    here. Their omission only weakens this pre-filter; the subsequent exact
    correction-term comparison still decides whether a candidate survives.
    """
    if inv_a['factors'] != inv_b['factors']:
        return False
    for p, v in inv_a.items():
        if p == 'factors':
            continue
        if p in inv_b and inv_b[p] != v:
            return False
    return True
