"""Exact, conservative Ni--Wu obstruction to unknotting number one.

If u(K)=1 and D=det(K), the Montesinos trick identifies one orientation of
the double branched cover with D/2 surgery on a knot in S^3.  Proposition
1.6 of Ni--Wu (arXiv:1009.4720), with H_j=V_{-j}, therefore prescribes

    d(L(D,2), i) - d(Y, i) = 2 V_min(floor(i/2), floor((D+1-i)/2)).

The V_j are nonnegative integers, with successive differences zero or one.
This module tests EVERY affine bijection i -> a*i+b of Z/D, for both
orientations. It deliberately imposes no linking-pairing or spin-origin
restriction. This larger search includes all identifications induced by a
surgery, independently of conventions for first-Chern-class coordinates.
Failure of the larger search is sufficient for an obstruction. A surviving
fit is only inconclusive: these necessary conditions do not construct surgery
or an unknotting crossing. The supplied d-invariants require their own proof.
"""

from fractions import Fraction
from math import gcd


def lens_half_integral_d(order, index):
    """d(L(D,2),i), for odd D>1, in the positive-surgery convention.

    This is the lens-space recursion expanded through L(2,1). Keeping the
    expression here independent of the plumbing implementation provides a
    second arithmetic check of its lens-space formula.
    """
    if type(order) is not int or order < 3 or order % 2 == 0:
        raise ValueError('D must be an odd integer greater than one')
    if type(index) is not int or not 0 <= index < order:
        raise ValueError('lens-space label must lie in 0,...,D-1')
    return (Fraction((2 * index - order - 1) ** 2, 8 * order)
            - (Fraction(1, 2) if index % 2 == 0 else 0))


def surgery_test(d, order):
    """Return all surviving fits and exhaustive failure counts, using rationals.

    The labels 0,...,D-1 can be any affine cyclic coordinates on Spin^c(Y).
    In particular, odd order makes first-Chern-class coordinates admissible.
    Each reported fit contains the orientation, affine map, and finite V
    prefix. A nonnegative prefix with drops at most one extends to a sequence
    eventually equal to zero; no unwarranted genus bound is imposed.
    """
    # Validate before iteration: D <= 0 would otherwise skip every test and
    # yield an empty fit list that could be mistaken for an obstruction.
    if type(order) is not int or order < 3 or order % 2 == 0:
        raise ValueError('D must be an odd integer greater than one')
    lens = [lens_half_integral_d(order, i) for i in range(order)]
    if (set(d) != set(range(order))
            or any(type(k) is not int for k in d)):
        raise ValueError('exactly D correction terms with integer cyclic labels are required')
    if any(isinstance(v, float) for v in d.values()):
        raise ValueError('correction terms must be exact rational numbers, not floats')
    values = [Fraction(d[k]) for k in range(order)]
    indices = [min(i // 2, (order + 1 - i) // 2) for i in range(order)]
    units = [a for a in range(1, order) if gcd(a, order) == 1]
    fits, orientations = [], {}
    for eps in (1, -1):
        counts = dict(tested=0, invalid_gap=0, inconsistent_V=0,
                      invalid_successive_difference=0, surviving=0)
        for a in units:
            for b in range(order):
                counts['tested'] += 1
                prefix, failure = {}, None
                for i, j in enumerate(indices):
                    v = (lens[i] - eps * values[(a * i + b) % order]) / 2
                    if v < 0 or v.denominator != 1:
                        failure = 'invalid_gap'
                        break
                    if j in prefix and prefix[j] != v:
                        failure = 'inconsistent_V'
                        break
                    prefix[j] = int(v)
                if failure is None:
                    seq = [prefix[j] for j in range(max(prefix) + 1)]
                    if any(seq[j] - seq[j + 1] not in (0, 1)
                           for j in range(len(seq) - 1)):
                        failure = 'invalid_successive_difference'
                if failure is not None:
                    counts[failure] += 1
                else:
                    fits.append((eps, a, b, seq))
                    counts['surviving'] += 1
        orientations[str(eps)] = counts
    return {'det': order, 'identifications': 'all_affine_bijections',
            'linking_filter': False, 'orientations': orientations, 'fits': fits}


def u1_admissible(d, lam, D):
    """Compatibility entry point; ``lam`` is intentionally not used as a filter.

    The pairing is still needed upstream to align the diagram's candidate
    gradings. Here, enumerating every affine map is a stronger verification
    of a negative result than restricting the search by that pairing.
    """
    return surgery_test(d, D)['fits']
