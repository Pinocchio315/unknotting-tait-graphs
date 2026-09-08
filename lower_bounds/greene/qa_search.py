#!/usr/bin/env python
"""Quasi-alternating certificate search with re-simplification of every resolved link.

A link L is quasi-alternating (QA) if it is the unknot, or some diagram has a crossing whose two resolutions L_0,
L_1 are QA with det(L) = det(L_0) + det(L_1).  Non-split alternating links are QA.  The search resolves crossings at
the level of the PD code, simplifies every resolved diagram with spherogram (global moves), reads the determinant
from the signed Tait graph (det = |sum over spanning trees of the product of the edge signs|), and recurses.
A certificate proves QA; failure proves nothing.
    python qa_search.py 12n_491 8_20 9_43
"""
import ast, sys, os, time, sympy, spherogram
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', '..', 'upper_bounds'))
from xtait.graph import from_pd

def tau_of_pd(pd):
    if not pd: return 1
    g = from_pd([list(q) for q in pd]); V = list(g.vertices); idx = {v: i for i, v in enumerate(V)}; n = len(V)
    if n == 1: return 1
    L = sympy.zeros(n, n)
    for e, (u, v) in g.edges.items():
        if u == v: continue
        s = g.signs[e]; i, j = idx[u], idx[v]; L[i, i] += s; L[j, j] += s; L[i, j] -= s; L[j, i] -= s
    return int(L[1:, 1:].det())

def smoothings(pd, k):
    """the two resolutions of crossing k of the PD code (labels merged); None if a free circle is produced"""
    a, b, c, d = pd[k]; rest = [list(q) for i, q in enumerate(pd) if i != k]; out = []
    for pairs in (((a, b), (c, d)), ((a, d), (b, c))):
        if any(x == y for x, y in pairs): out.append(None); continue        # the smoothing closes a free circle
        sub = {}
        for x, y in pairs: sub[y] = x
        new = [[sub.get(v, v) for v in q] for q in rest]
        # a label pair may have merged two arcs of the same circle: harmless; keep
        out.append(new)
    return out

def simplify(pd):
    if not pd: return [], 1
    L = spherogram.Link([tuple(int(v) for v in q) for q in pd]); L.simplify('global')
    return [list(q) for q in L.PD_code()], len(L.link_components)

def relabel(pd):
    labels = sorted({v for q in pd for v in q}); m = {l: i + 1 for i, l in enumerate(labels)}
    return [[m[v] for v in q] for q in pd]

memo = {}
def qa(pd, stats, depth=0):
    stats['nodes'] += 1
    pd, ncomp = simplify(pd)
    if not pd: return 'unknot' if ncomp == 1 else None
    pd = relabel(pd); key = str(sorted(map(tuple, pd)))
    if key in memo: return memo[key]
    t = tau_of_pd(pd)
    if t == 0: memo[key] = None; return None
    L = spherogram.Link([tuple(q) for q in pd])
    if L.is_alternating(): memo[key] = f'alternating({len(pd)})'; return memo[key]
    for k in range(len(pd)):
        s0, s1 = smoothings(pd, k)
        if s0 is None or s1 is None: continue
        p0, n0 = simplify(s0); p1, n1 = simplify(s1)
        if (not p0 and n0 > 1) or (not p1 and n1 > 1): continue
        d0, d1 = abs(tau_of_pd(relabel(p0)) if p0 else 1), abs(tau_of_pd(relabel(p1)) if p1 else 1)
        if d0 == 0 or d1 == 0 or d0 + d1 != abs(t): continue
        r0 = qa(p0, stats, depth + 1)
        if r0 is None: continue
        r1 = qa(p1, stats, depth + 1)
        if r1 is None: continue
        memo[key] = {'crossing': k, 'det': (abs(t), d0, d1), 'res0': r0, 'res1': r1}; return memo[key]
    memo[key] = None; return None

if __name__ == '__main__':
    import database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for nm in sys.argv[1:]:
        pd = [list(q) for q in ast.literal_eval(rows[nm]['pd_notation'])]; stats = {'nodes': 0}; t0 = time.time(); r = qa(pd, stats)
        print(f'{nm}: det {rows[nm]["determinant"]}, KnotInfo QA={rows[nm]["quasi_alternating"]!r}: {"QUASI-ALTERNATING" if r else "no certificate"} ({stats["nodes"]} nodes, {time.time()-t0:.1f}s)', flush=True)
        if r and nm == sys.argv[-1]: print('  ', r)
