"""Add signature of K_c (spherogram), effective resistance of the crossing's edge in the Tait graph, and
Tait sign to each row of the crossing dataset."""
from paths import *
import sys, json, ast, glob
from fractions import Fraction
import sympy, spherogram, database_knotinfo as dk
from sig import signature as gl_signature
from xtait.graph import from_pd, to_pd, dual
from xtait import moves as mv
S = WORK
si, sn = (int(x) for x in sys.argv[1].split('/'))
rows = {r['name']: r for r in dk.link_list()}
knots = sorted({json.loads(l)['knot'] for f in glob.glob(S + '/rows_*.jsonl') for l in open(f)})[si::sn]
def reff(g):
    V = sorted(g.vertices); idx = {v: i for i, v in enumerate(V)}; n = len(V)
    L = sympy.zeros(n, n)
    for e in g.edges:
        u, w = g.edges[e]; L[idx[u], idx[u]] += 1; L[idx[w], idx[w]] += 1; L[idx[u], idx[w]] -= 1; L[idx[w], idx[u]] -= 1
    Lr = L[1:, 1:].inv()
    out = {}
    for e in g.edges:
        u, w = g.edges[e]; x = sympy.zeros(n - 1, 1)
        if idx[u] > 0: x[idx[u] - 1] += 1
        if idx[w] > 0: x[idx[w] - 1] -= 1
        out[e] = sympy.Rational((x.T * Lr * x)[0, 0])
    return out
with open(S + f'/aug_{si}.jsonl', 'a') as h:
    for nm in knots:
        pd = ast.literal_eval(rows[nm]['pd_notation']); g = from_pd([list(q) for q in pd]); gd = dual(g)
        edges = sorted(g.edges); Rb = reff(g); Rw = reff(gd)
        sigK = int(rows[nm]['signature'])
        for i, e in enumerate(edges):
            hgr = mv.crossing_change(g, e)
            try: sc = int(gl_signature(spherogram.Link([list(q) for q in to_pd(hgr)]).PD_code()))
            except Exception: sc = None
            h.write(json.dumps({'knot': nm, 'crossing': i, 'sigma': sigK, 'sigma_c': sc, 'R_black': str(Rb[e]), 'R_white': str(Rw[e]),
                                'tait_sign': int(g.signs[e]), 'n_black': len(g.vertices), 'n_white': len(gd.vertices)}) + '\n')
        h.flush()
