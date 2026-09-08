"""Per-crossing dataset: for alternating knots with known u <= 4, change each crossing of the KnotInfo
minimal diagram, identify the result rigorously and record u(K_c) together with local features."""
from paths import *
import sys, json, re, ast, time, collections, math
import numpy as np, spherogram, snappy, database_knotinfo as dk
from xtait.graph import from_pd, to_pd, dual
from xtait import moves as mv
from xtait.canonical import canonical_code
from xtait.reduce import greedy_reduce
from flipdet import FlipDet, P
from alexander import det_from_alexander, alexander_eval
out, shard = sys.argv[1], sys.argv[2]; si, sn = (int(x) for x in shard.split('/'))
rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
U = {nm: lo for nm, (lo, hi) in json.load(open(U_TABLE)).items() if lo == hi}   # reference table + our bounds (consolidate_results.py)
U['0_1'] = 0
allk = json.load(open(ALL_CODES)); code2name = {c: x['name'] for x in allk for c in x['codes']}
alt = lambda r: str(r.get('alternating', '')).upper().startswith('Y')
knots = [nm for nm, r in rows.items() if alt(r) and nm in U and 1 <= U[nm] <= 4]
knots.sort(key=lambda n: (int(rows[n]['crossing_number']), n)); knots = knots[si::sn]
def turn(g):
    h = g.copy()
    for v in h.vertices: h.rotation[v] = list(reversed(h.rotation[v]))
    return h
def mirror(g):
    h = g.copy()
    for e in h.signs: h.signs[e] = -h.signs[e]
    return h
def ki_name(census):
    n = census.split('(')[0]
    if n.startswith('K'):
        n = n[1:]
        for s in ('a', 'n'):
            if s in n: return n.replace(s, s + '_', 1)
    return n
extcache = {}
def ext(nm):
    if nm not in extcache: extcache[nm] = spherogram.Link(ast.literal_eval(rows[nm]['pd_notation'])).exterior()
    return extcache[nm]
idcache = {}
def identify(h):
    red, _ = greedy_reduce(h)
    key = repr(canonical_code(red))
    if key in idcache: return idcache[key]
    name = None; how = None
    if red.n_crossings() == 0: name, how = '0_1', 'trivial'
    if name is None:
        for x in (red, turn(red), mirror(red), mirror(turn(red))):
            c = repr(canonical_code(x))
            if c in code2name: name, how = code2name[c], 'code'; break
    if name is None:
        try:
            M = spherogram.Link(to_pd(red)).exterior()
            for cand in M.identify():
                n = ki_name(str(cand))
                if n in rows:
                    try:
                        if M.is_isometric_to(ext(n)): name, how = n, 'isometry'; break
                    except Exception: pass
        except Exception: pass
    idcache[key] = (name, how, red.n_crossings())
    return idcache[key]
with open(out, 'a') as h:
    for nm in knots:
        r = rows[nm]; pd = ast.literal_eval(r['pd_notation']); uK = U[nm]; D = int(r['determinant'])
        g = from_pd([list(q) for q in pd]); gd = dual(g); edges = sorted(g.edges)
        L = spherogram.Link(to_pd(g)); signs = [c.sign for c in L.crossings]
        fd = FlipDet(g)
        deg = {v: g.degree(v) for v in g.vertices}
        # multiplicities: parallel edges between the same endpoints in g and in the dual
        par = collections.Counter(tuple(sorted(g.edges[e])) for e in edges)
        pard = collections.Counter(tuple(sorted(gd.edges[e])) for e in gd.edges)
        degd = {v: gd.degree(v) for v in gd.vertices}
        twist_max = max(max(par.values()), max(pard.values()))
        for i, e in enumerate(edges):
            u, w = g.edges[e]; ud, wd = gd.edges[e]
            hgraph = mv.crossing_change(g, e)
            t0 = time.time()
            name, how, nred = identify(hgraph)
            uc = U.get(name) if name else None
            detc = fd.det_after([i]); detc = detc if detc < P // 2 else P - detc
            pdc = to_pd(hgraph)
            rec = {'knot': nm, 'u': uK, 'det': D, 'crossings': len(edges), 'crossing': i,
                   'sign': int(signs[i]), 'deg_min': min(deg[u], deg[w]), 'deg_max': max(deg[u], deg[w]),
                   'face_min': min(degd[ud], degd[wd]), 'face_max': max(degd[ud], degd[wd]),
                   'twist_par': par[tuple(sorted((u, w)))], 'twist_ser': pard[tuple(sorted((ud, wd)))], 'twist_max_diagram': twist_max,
                   'det_changed': int(detc), 'reduced_crossings': nred, 'result': name, 'result_how': how,
                   'u_result': uc, 'label': (None if uc is None else ('good' if uc == uK - 1 else ('neutral' if uc == uK else 'bad'))),
                   'seconds': round(time.time() - t0, 2)}
            h.write(json.dumps(rec) + '\n'); h.flush()
