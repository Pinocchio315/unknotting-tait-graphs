#!/usr/bin/env python
"""Reduced Heegaard Floer homology of the boundary of a negative-definite plumbing with at most one bad vertex,
per Spin^c structure, from lattice cohomology (Nemethi; Ozsvath-Szabo).

For a characteristic element k of the lattice (L, (.,.)) with negative-definite form Q_neg, the weight function
chi_k(x) = -((x,x) + k(x))/2 on x in Z^n has bounded sublevel sets; the graded root is the tree of connected
components (under unit steps) of the sublevel sets {chi_k <= N}, and HF^+(-Y, [k]) is its homology (Nemethi 2005,
Ozsvath-Szabo 2003 for plumbings with at most one bad vertex).  The tower contributes one generator per level from
the minimum upwards, so rank HF_red(Y, [k]) = sum_N (#components(N) - 1); the ranks of HF_red are the same for Y and
-Y.  Classes [k] are labelled as in owens_obstruction.ClassMap (k mod Q Z^n), the labelling used by montesinos_u1.
"""
import itertools, math, sys, os
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, 'lower_bounds', 'owens'))
from owens_obstruction import ClassMap

class TooManyPoints(Exception):
    """A resource bound was reached; no Floer obstruction may use this run."""


class _ExactEllipsoid:
    """Rational LDL^T enumeration for the weight (x^T A x - k.x)/2.

    No floating-point boundary test is used: a missed lattice point could
    change connectivity and produce an invalid unknotting obstruction.
    """
    def __init__(self, matrix, characteristic):
        from sympy import Matrix
        self.A = Matrix(matrix)
        self.k = Matrix(characteristic)
        self.n = self.A.rows
        if self.A.rows != self.A.cols or self.A != self.A.T or self.k.shape != (self.n, 1):
            raise ValueError('a symmetric square form and a characteristic vector are required')
        if any(x.q != 1 for x in self.A) or any(x.q != 1 for x in self.k):
            raise ValueError('the form and characteristic vector must be integral')
        self.diagonal = [int(self.A[i, i]) for i in range(self.n)]
        if any((int(self.k[i]) - self.diagonal[i]) % 2 for i in range(self.n)):
            raise ValueError('k must be characteristic')
        lower, diagonal = self.A.LDLdecomposition(hermitian=False)
        self.L = [[Fraction(lower[i, j]) for j in range(self.n)] for i in range(self.n)]
        self.d = [Fraction(diagonal[i, i]) for i in range(self.n)]
        if any(d <= 0 for d in self.d):
            raise ValueError('the negative plumbing form must be negative definite')
        self.inverse = self.A.inv()
        self.c = [Fraction(x) / 2 for x in self.inverse * self.k]
        self.offset = Fraction((self.k.T * self.inverse * self.k)[0]) / 4

    def chi(self, vector):
        # Python integers avoid overflow in intermediate quadratic products.
        quadratic = sum(int(self.A[i, j]) * vector[i] * vector[j]
                        for i in range(self.n) for j in range(self.n))
        return (quadratic - sum(int(t) * x for t, x in zip(self.k, vector))) // 2

    def points(self, bound, max_points):
        radius = Fraction(bound) + self.offset
        if radius < 0:
            return []
        x = [0] * self.n
        output = []
        def visit(i, remaining):
            if i < 0:
                output.append(tuple(x))
                if len(output) > max_points:
                    raise TooManyPoints
                return
            center = self.c[i] - sum(self.L[j][i] * (x[j] - self.c[j])
                                    for j in range(i + 1, self.n))
            squared_radius = remaining / self.d[i]
            # This deliberately enlarged integer interval contains the exact
            # rational ellipsoid section; the inequality removes excess points.
            width = math.isqrt(squared_radius.numerator // squared_radius.denominator) + 1
            for value in range(math.floor(center) - width, math.ceil(center) + width + 1):
                used = self.d[i] * (value - center) ** 2
                if used <= remaining:
                    x[i] = value
                    visit(i - 1, remaining - used)
        visit(self.n - 1, radius)
        return output

    def local_minimum_levels(self, max_points):
        # A weak local minimum under every unit step satisfies
        # |(2Ax-k)_i| <= A_ii. Characteristic parity reduces this finite box.
        # Above its largest weight no new connected component can appear.
        ranges = [range(-a, a + 1, 2) for a in self.diagonal]
        if math.prod(len(values) for values in ranges) > max_points:
            raise TooManyPoints
        from sympy import Matrix
        weak_minima = {}
        for z in itertools.product(*ranges):
            vector = self.inverse * (Matrix(z) + self.k) / 2
            if all(value.q == 1 for value in vector):
                point = tuple(int(value) for value in vector)
                weak_minima[point] = self.chi(point)
        # Flat plateaux need not create components: a weak minimum can move
        # at equal weight to a point with a descending step. Collapse each
        # equal-weight component of weak minima and retain only those with
        # no such exit. Every birth of a sublevel component occurs here.
        seen = set()
        births = []
        for point, level in weak_minima.items():
            if point in seen:
                continue
            pending = [point]
            seen.add(point)
            escapes = False
            while pending:
                current = pending.pop()
                for j in range(self.n):
                    for sign in (-1, 1):
                        adjacent = list(current)
                        adjacent[j] += sign
                        adjacent = tuple(adjacent)
                        if self.chi(adjacent) != level:
                            continue
                        if adjacent not in weak_minima:
                            escapes = True
                        elif adjacent not in seen:
                            seen.add(adjacent)
                            pending.append(adjacent)
            if not escapes:
                births.append(level)
        return births


def lattice_points(Qpos, k, bound, max_points=400_000):
    """All integer x with x^T Qpos x - k.x <= bound, using exact arithmetic."""
    return _ExactEllipsoid(Qpos, k).points(bound, max_points)


def hf_red_rank(Qneg, k, max_depth=40, max_points=400_000):
    """Return (dim_F2 HF_red, min chi), or (None, min chi) if uncertified.

    Section 3, ``Montesinos knots'': the graded root of the negative-definite
    plumbing gives HF^+(-Y) for a tree with at most one bad vertex
    (Ozsvath--Szabo, On the Floer homology of plumbed three-manifolds,
    Theorem 1.2; Nemethi, On the Ozsvath--Szabo invariant..., Section 11).
    HF_red has the same rank after orientation reversal. Each level adds
    number_of_components - 1 finite generators.

    We enumerate all weak local minima and their equal-weight plateaux.
    A plateau with an exit to a lower weight does not create a new component.
    We continue until the sublevel set is connected and contains all possible
    births. Above that level every point descends through unit steps of
    nonincreasing weight, with eventual strict decrease; positive definiteness
    makes descent terminate. Every later sublevel set is therefore connected.
    Resource limits cause abstention, never an asserted rank from an unproved stable-looking tail.
    """
    try:
        from .montesinos_u1 import validate_plumbing
    except ImportError:
        from montesinos_u1 import validate_plumbing
    validate_plumbing(Qneg)
    ellipsoid = _ExactEllipsoid(-np.array(Qneg, dtype=object), k)
    try:
        minima = ellipsoid.local_minimum_levels(max_points)
    except TooManyPoints:
        return None, None
    minimum, certificate_level = min(minima), max(minima)
    if certificate_level - minimum > max_depth:
        return None, minimum
    top = certificate_level
    while True:
        try:
            points = ellipsoid.points(2 * top, max_points)
        except TooManyPoints:
            return None, minimum
        levels = {}
        for vector in points:
            levels.setdefault(ellipsoid.chi(vector), []).append(vector)
        parent = {}
        components = rank = 0
        def find(vector):
            while parent[vector] != vector:
                parent[vector] = parent[parent[vector]]
                vector = parent[vector]
            return vector
        for level in range(minimum, top + 1):
            for vector in levels.get(level, []):
                parent[vector] = vector
                components += 1
                for j in range(ellipsoid.n):
                    for sign in (-1, 1):
                        adjacent = list(vector)
                        adjacent[j] += sign
                        adjacent = tuple(adjacent)
                        if adjacent in parent:
                            left, right = find(vector), find(adjacent)
                            if left != right:
                                parent[left] = right
                                components -= 1
            rank += components - 1
        if components == 1:
            return rank, minimum
        if top - minimum >= max_depth:
            return None, minimum
        top = min(top + 1, minimum + max_depth)

def hf_red_all(Qneg):
    """rank HF_red(Y, t) for every Spin^c structure t, keyed by ClassMap coordinates (as in m_Q)."""
    Qpos = -np.array(Qneg, dtype=object); n = len(Qpos); cm = ClassMap(Qpos); D = cm.order
    diag = np.array([Qpos[i, i] % 2 for i in range(n)], dtype=np.int64)
    reps = {}
    # representatives: characteristic vectors k = diag + 2y with small y, until every class is hit
    for rad in range(0, 6):
        for y in itertools.product(range(-rad, rad + 1), repeat=n):
            k = diag + 2 * np.array(y, dtype=np.int64); g = cm.coords(k)
            if g not in reps: reps[g] = k
        if len(reps) == D: break
    if len(reps) != D:
        return None, cm  # incomplete Spin^c coverage is not a Floer computation
    out = {}
    for g, k in reps.items():
        r, m = hf_red_rank(Qneg, k)
        if r is None: return None, cm            # too large to enumerate within the memory bound: no Floer test
        out[g] = r
    return out, cm

if __name__ == '__main__':
    import ast, database_knotinfo as dk, time
    from montesinos_u1 import plumbing_for, correction_terms
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for nm in sys.argv[1:]:
        r = rows[nm]; t0 = time.time(); Q, fl = plumbing_for(r['montesinos_notation'].strip(), int(r['determinant']))
        d, cm = correction_terms(Q); ranks, _ = hf_red_all(Q)
        print(nm, r['montesinos_notation'], 'det', r['determinant'], 'plumbing rank', len(Q), 'flipped', fl)
        for g in sorted(d): print('   class', g, 'd(Y)=', d[g], 'rank HF_red =', ranks[g])
        print('   time', round(time.time() - t0, 1), 's')
