"""Knot diagrams from closed space polygons (PL knots): generic projection to PD code.

This general conversion utility is not used for the results reported in the paper.  A random
rotation is applied first, so repeated calls can give different projections of the same knot.

PD convention: KnotTheory X[i, j, k, l] with i the incoming under-arc, then counterclockwise; arcs are labelled
0..2c-1 in traversal order (arc m runs from crossing visit m to visit m+1), which spherogram.Link accepts.
"""
from __future__ import annotations

import numpy as np


def random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.standard_normal((3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def read_polygon(path: str) -> np.ndarray:
    """Whitespace separated x y z per line (comments '#'); the last point is not repeated."""
    pts = []
    for ln in open(path):
        ln = ln.strip()
        if not ln or ln.startswith('#'):
            continue
        pts.append([float(t) for t in ln.replace(',', ' ').split()[:3]])
    P = np.asarray(pts, float)
    if len(P) > 1 and np.allclose(P[0], P[-1]):
        P = P[:-1]
    return P


def _seg_intersections(X: np.ndarray, Y: np.ndarray, eps: float):
    """All proper intersections between non-adjacent edges of the closed 2-D polygon (X, Y).
    Returns list of (i, ti, j, tj) with i < j, and raises ValueError on a degenerate (non-generic) projection."""
    n = len(X)
    i_idx, j_idx = np.triu_indices(n, 2)
    keep = ~((i_idx == 0) & (j_idx == n - 1))            # edge n-1 is adjacent to edge 0
    i_idx, j_idx = i_idx[keep], j_idx[keep]
    p = np.stack([X, Y], 1)
    d = np.roll(p, -1, 0) - p                             # edge vectors
    P, D = p[i_idx], d[i_idx]
    Q, E = p[j_idx], d[j_idx]
    den = D[:, 0] * E[:, 1] - D[:, 1] * E[:, 0]
    W = Q - P
    par = np.abs(den) < 1e-14
    with np.errstate(divide='ignore', invalid='ignore'):
        t = (W[:, 0] * E[:, 1] - W[:, 1] * E[:, 0]) / den
        s = (W[:, 0] * D[:, 1] - W[:, 1] * D[:, 0]) / den
    hit = (~par) & (t > -eps) & (t < 1 + eps) & (s > -eps) & (s < 1 + eps)
    # degeneracy: an intersection at / extremely near an endpoint, or (anti)parallel overlapping edges
    near_end = hit & ((np.abs(t) < eps) | (np.abs(t - 1) < eps) | (np.abs(s) < eps) | (np.abs(s - 1) < eps))
    if near_end.any():
        raise ValueError('non-generic projection (intersection at a vertex)')
    if par.any():
        # parallel edges: reject if collinear and overlapping (rare); check quickly via cross product of W and D
        cr = np.abs(W[par, 0] * D[par, 1] - W[par, 1] * D[par, 0])
        if (cr < 1e-12).any():
            raise ValueError('non-generic projection (collinear edges)')
    idx = np.nonzero(hit)[0]
    return [(int(i_idx[m]), float(t[m]), int(j_idx[m]), float(s[m])) for m in idx]


def polygon_to_pd(points: np.ndarray, rng: np.random.Generator | None = None, rotation: np.ndarray | None = None,
                  max_tries: int = 20, eps: float = 1e-9):
    """Project the closed polygon (after a random rotation) and return (pd_code, info).
    info = {'crossings': c, 'rotation': R}.  pd_code is a list of 4-lists (0-based arc labels), [] for a
    crossing-free projection."""
    rng = rng or np.random.default_rng()
    P = np.asarray(points, float)
    n = len(P)
    if n < 3:
        raise ValueError('need at least 3 points')
    for _ in range(max_tries):
        R = rotation if rotation is not None else random_rotation(rng)
        Q = P @ R.T
        X, Y, Z = Q[:, 0], Q[:, 1], Q[:, 2]
        try:
            inter = _seg_intersections(X, Y, eps)
        except ValueError:
            if rotation is not None:
                raise
            continue
        # heights along edges
        dz = np.roll(Z, -1) - Z
        visits = []                                       # (edge, t, crossing id, is_under)
        cross_dirs = []
        for cid, (i, ti, j, tj) in enumerate(inter):
            zi = Z[i] + ti * dz[i]
            zj = Z[j] + tj * dz[j]
            if abs(zi - zj) < 1e-12:
                break
            under_i = zi < zj
            visits.append((i, ti, cid, under_i))
            visits.append((j, tj, cid, not under_i))
            di = np.array([X[(i + 1) % n] - X[i], Y[(i + 1) % n] - Y[i]])
            dj = np.array([X[(j + 1) % n] - X[j], Y[(j + 1) % n] - Y[j]])
            cross_dirs.append((di, dj, under_i))
        else:
            c = len(inter)
            if c == 0:
                return [], {'crossings': 0, 'rotation': R}
            visits.sort(key=lambda v: (v[0], v[1]))
            m = len(visits)                               # = 2c
            # arc m runs from visit m to visit m+1 (mod 2c)
            pos = {}                                      # (cid, is_under) -> visit index
            for k, (i, t, cid, under) in enumerate(visits):
                pos[(cid, under)] = k
            pd = []
            for cid in range(c):
                ku, ko = pos[(cid, True)], pos[(cid, False)]
                in_u, out_u = (ku - 1) % m, ku
                in_o, out_o = (ko - 1) % m, ko
                di, dj, under_i = cross_dirs[cid]
                du, do = (di, dj) if under_i else (dj, di)
                cr = du[0] * do[1] - du[1] * do[0]
                if cr > 0:
                    pd.append([in_u, in_o, out_u, out_o])
                else:
                    pd.append([in_u, out_o, out_u, in_o])
            return pd, {'crossings': c, 'rotation': R}
        if rotation is not None:
            raise ValueError('degenerate heights')
    raise ValueError('no generic projection found')


def torus_knot_polygon(p: int, q: int, n: int = 200, R: float = 2.0, r: float = 0.7) -> np.ndarray:
    """Sampled (p, q) torus knot, for tests."""
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = (R + r * np.cos(q * t)) * np.cos(p * t)
    y = (R + r * np.cos(q * t)) * np.sin(p * t)
    z = r * np.sin(q * t)
    return np.stack([x, y, z], 1)
