"""Gordon–Litherland signature from a PD code via the xtait checkerboard graph.  Convention calibrated
against KnotInfo (see calibrate())."""
from paths import *
import sys, ast
import numpy as np
from fractions import Fraction
from xtait.graph import from_pd
VARIANT = {'type2_flip': True, 'global': 1}   # calibrated: 12965/12965 KnotInfo signatures

def symmetric_signature(matrix):
    """Exact inertia difference of a rational symmetric form, including null directions.

    Congruence preserves inertia (Sylvester's law). A nonzero diagonal entry
    splits off a one-dimensional form; if all diagonals vanish, a nonzero
    off-diagonal entry splits off a hyperbolic two-dimensional form. Schur
    complements therefore avoid an eigenvalue tolerance in lower-bound tests.
    """
    rows = matrix.tolist() if hasattr(matrix, 'tolist') else matrix
    A = [[Fraction(x) for x in row] for row in rows]
    n = len(A)
    if any(len(row) != n for row in A) or any(A[i][j] != A[j][i] for i in range(n) for j in range(n)):
        raise ValueError('expected a square symmetric matrix')
    signature = 0
    while A:
        n = len(A)
        pivot = next((i for i in range(n) if A[i][i]), None)
        if pivot is not None:
            order = [pivot] + [i for i in range(n) if i != pivot]
            A = [[A[i][j] for j in order] for i in order]
            d = A[0][0]
            signature += 1 if d > 0 else -1
            A = [[A[i][j] - A[i][0] * A[0][j] / d
                  for j in range(1, n)] for i in range(1, n)]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i + 1, n) if A[i][j]), None)
            if pair is None:
                break  # the remaining zero form contributes only to nullity
            order = list(pair) + [i for i in range(n) if i not in pair]
            A = [[A[i][j] for j in order] for i in order]
            b = A[0][1]  # [[0,b],[b,0]] has one positive and one negative eigenvalue
            A = [[A[i][j] - (A[i][0] * A[1][j] + A[i][1] * A[0][j]) / b
                  for j in range(2, n)] for i in range(2, n)]
    return signature

def goeritz_data(pd):
    """Return integral Goeritz form and type-II correction for consecutive oriented PD labels."""
    pd = [[int(x) for x in q] for q in pd]
    g = from_pd(pd); n2 = 2 * len(pd)
    V = sorted(g.vertices); idx = {v: i for i, v in enumerate(V)}
    G = np.zeros((len(V), len(V)), dtype=object)
    mu = 0
    for c, q in enumerate(pd):
        a, b, cc, d = q; eta = g.signs[c]; u, w = g.edges[c]
        G[idx[u], idx[w]] -= eta; G[idx[w], idx[u]] -= eta; G[idx[u], idx[u]] += eta; G[idx[w], idx[w]] += eta
        if (d - b) % n2 == 1: over_b_to_d = True
        elif (b - d) % n2 == 1: over_b_to_d = False
        else: raise ValueError('labels not consecutive')
        type2 = (over_b_to_d if eta == 1 else not over_b_to_d)
        if VARIANT['type2_flip']: type2 = not type2
        if type2: mu += eta
    return G[1:, 1:], mu
def signature(pd):
    G, mu = goeritz_data(pd)
    s = symmetric_signature(G)
    return VARIANT['global'] * (s - mu)
def calibrate():
    import database_knotinfo as dk
    rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit() and r.get('signature', '').strip()]
    best = None
    for flip in (False, True):
        for gl in (1, -1):
            VARIANT['type2_flip'], VARIANT['global'] = flip, gl
            ok = bad = 0; ex = []
            for r in rows:
                try:
                    s = signature(ast.literal_eval(r['pd_notation']))
                except Exception as e:
                    bad += 1; ex.append((r['name'], repr(e))); continue
                if s == int(r['signature']): ok += 1
                else:
                    bad += 1
                    if len(ex) < 3: ex.append((r['name'], s, r['signature']))
            print(f'flip={flip} global={gl}: ok={ok} bad={bad} {ex[:3]}')
            if best is None or ok > best[0]: best = (ok, flip, gl)
    # Calibration is diagnostic; do not leave the last, possibly wrong, convention active.
    VARIANT['type2_flip'], VARIANT['global'] = best[1], best[2]
    return best
if __name__ == '__main__':
    print(calibrate())
