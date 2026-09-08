"""Per-crossing feature rows for the open alternating [2,3] knots (hypothesis u = 2), with rigorous
identification of the knot K_c obtained by each crossing change and its known u-range.
    python apply833.py i/n        (shard i of n; writes results/open23/work/rows833_i.jsonl)"""
import sys, os, json, re, ast, collections
from paths import *
from candidates import connected_sum_interval
S = os.path.join(OPEN23, 'work'); os.makedirs(S, exist_ok=True)
import sympy, spherogram, snappy, database_knotinfo as dk, networkx as nx
from xtait.jones import jones
from flipdet import exact_det
def jmul(A, B):
    out = {}
    for i, a in A.items():
        for j, b in B.items(): out[i + j] = out.get(i + j, 0) + a * b
    return {k: v for k, v in out.items() if v}
def jmirror(A): return {-k: v for k, v in A.items()}
jcache = {}
def jones_of(nm):
    if nm not in jcache: jcache[nm] = {k: v for k, v in jones(ast.literal_eval(rows[nm]['pd_notation'])).items() if v}
    return jcache[nm]

from xtait.graph import from_pd, to_pd, dual, EmbeddedTaitGraph
from xtait import moves as mv
from xtait.canonical import canonical_code
from xtait.reduce import greedy_reduce
from flipdet import FlipDet, P
from sig import signature as gl_signature
si, sn = (int(x) for x in sys.argv[1].split('/'))
rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
RANGE = {nm: list(v) for nm, v in json.load(open(U_TABLE)).items()}   # reference table + our bounds
RANGE['0_1'] = [0, 0]
allk = json.load(open(ALL_CODES)); code2name = {c: x['name'] for x in allk for c in x['codes']}
bydet = collections.defaultdict(list)
for _x in allk:
    if _x['crossings'] >= 3: bydet[_x['det']].append(_x['name'])
# targets: every alternating knot whose range in the reference table (with our bounds) is [2,3]
targets = sorted(nm for nm, v in RANGE.items() if v == [2, 3] and nm in rows and str(rows[nm].get('alternating', '')).upper().startswith('Y'))
import os as _os
P3 = {x['name']: x['P_u3'] for x in json.load(open(RANKING))} if _os.path.exists(RANKING) else {}   # optional ML ranking (experiments/ml)
P3 = {nm: P3.get(nm, 0.0) for nm in targets}
targets = targets[si::sn]
def turn(g):
    h = g.copy()
    for v in h.vertices: h.rotation[v] = list(reversed(h.rotation[v]))
    return h
def mirror(g):
    h = g.copy()
    for e in h.signs: h.signs[e] = -h.signs[e]
    return h
def lookup(red):
    for x in (red, turn(red), mirror(red), mirror(turn(red))):
        c = repr(canonical_code(x))
        if c in code2name: return code2name[c]
    return None
extcache = {}
def ext(nm):
    if nm not in extcache: extcache[nm] = spherogram.Link(ast.literal_eval(rows[nm]['pd_notation'])).exterior()
    return extcache[nm]
def ki_name(census):
    n = census.split('(')[0]
    if n.startswith('K'):
        n = n[1:]
        for s in ('a', 'n'):
            if s in n: return n.replace(s, s + '_', 1)
    return n
def split_sum(red):
    for gg in (red, dual(red)):
        G = nx.Graph()
        for x in gg.edges: G.add_edge(*gg.edges[x])
        comps = [set(c) for c in nx.biconnected_components(G)]
        if len(comps) < 2: continue
        names = []; ok = True
        for vs in comps:
            eids = [x for x in gg.edges if set(gg.edges[x]) <= vs and gg.edges[x][0] != gg.edges[x][1]]
            if not eids: continue
            sub = EmbeddedTaitGraph(sorted(vs), {x: gg.edges[x] for x in eids}, {x: gg.signs[x] for x in eids}, {v: [x for x in gg.rotation[v] if x in set(eids)] for v in vs}, None)
            try:
                hr, _ = greedy_reduce(sub); nmx = lookup(hr) if hr.n_crossings() else '0_1'
            except Exception: nmx = None
            if nmx is None: ok = False; break
            if nmx != '0_1': names.append(nmx)
        if ok and names: return names
    return None
idcache = {}
def identify(h):
    red, _ = greedy_reduce(h); key = repr(canonical_code(red))
    if key in idcache: return idcache[key]
    name = how = None; parts = None
    if red.n_crossings() == 0: name, how = '0_1', 'trivial'
    if name is None:
        name = lookup(red)
        if name: how = 'code'
    hyp = None
    if name is None:
        try:
            pdr = to_pd(red); M = spherogram.Link([list(q) for q in pdr]).exterior(); hyp = False
            for attempt in range(4):
                if M.solution_type().startswith('all tetrahedra positively'): hyp = True; break
                M.randomize()
            if hyp:
                # rigorous: isometry against every KnotInfo knot with the same determinant and Jones polynomial (mirror allowed)
                D = exact_det(red); Jp = {k: v for k, v in jones(pdr).items() if v}
                for n in bydet.get(D, []):
                    jn = jones_of(n)
                    if jn == Jp or jmirror(jn) == Jp:
                        try:
                            if M.is_isometric_to(ext(n)): name, how = n, 'isometry-jones'; break
                        except Exception: pass
                if name is None:
                    for cand in M.identify():
                        n = ki_name(str(cand))
                        if n in rows:
                            try:
                                if M.is_isometric_to(ext(n)): name, how = n, 'isometry'; break
                            except Exception: pass
        except Exception: pass
    if name is None:
        parts = split_sum(red)
        if parts: name, how = '#'.join(parts), 'sum-blocks'
    if name is None and hyp is not True:
        pdr = to_pd(red); D = exact_det(red); Jp = {k: v for k, v in jones(pdr).items() if v}; n = red.n_crossings(); cands = set()
        for A in allk:
            if A['det'] == 0 or D % A['det']: continue
            for B in allk:
                if A['name'] > B['name'] or A['det'] * B['det'] != D or A['crossings'] + B['crossings'] > n + 2: continue
                for ma in (0, 1):
                    for mb in (0, 1):
                        JA = jones_of(A['name']); JB = jones_of(B['name'])
                        if ma: JA = jmirror(JA)
                        if mb: JB = jmirror(JB)
                        if jmul(JA, JB) == Jp: cands.add((A['name'], B['name']))
        if len(cands) == 1:
            parts = list(cands)[0]; name, how = '#'.join(parts), 'sum-jones'
        elif cands: name, how = 'composite?', 'sum-ambiguous'; parts = None
    # u-range of the result
    if name is None: ur = None
    elif parts:
        # A decomposition with only one nontrivial block is not composite.
        # Jones-only matches retain their heuristic tag in candidate_status.
        ur = connected_sum_interval(parts, RANGE)
    else: ur = RANGE.get(name)
    idcache[key] = (name, how, red.n_crossings(), ur)
    return idcache[key]
def reff(g):
    V = sorted(g.vertices); idx = {v: i for i, v in enumerate(V)}; n = len(V)
    L = sympy.zeros(n, n)
    for e in g.edges:
        u, w = g.edges[e]; L[idx[u], idx[u]] += 1; L[idx[w], idx[w]] += 1; L[idx[u], idx[w]] -= 1; L[idx[w], idx[u]] -= 1
    Lr = L[1:, 1:].inv(); out = {}
    for e in g.edges:
        u, w = g.edges[e]; x = sympy.zeros(n - 1, 1)
        if idx[u] > 0: x[idx[u] - 1] += 1
        if idx[w] > 0: x[idx[w] - 1] -= 1
        out[e] = sympy.Rational((x.T * Lr * x)[0, 0])
    return out
with open(S + f'/rows833_{si}.jsonl', 'a') as h:
    for nm in targets:
        r = rows[nm]; pd = ast.literal_eval(r['pd_notation']); D = int(r['determinant']); sigK = int(r['signature'])
        g = from_pd([list(q) for q in pd]); gd = dual(g); edges = sorted(g.edges)
        L = spherogram.Link(to_pd(g)); signs = [c.sign for c in L.crossings]
        fd = FlipDet(g); deg = {v: g.degree(v) for v in g.vertices}; degd = {v: gd.degree(v) for v in gd.vertices}
        par = collections.Counter(tuple(sorted(g.edges[e])) for e in edges); pard = collections.Counter(tuple(sorted(gd.edges[e])) for e in gd.edges)
        twist_max = max(max(par.values()), max(pard.values())); Rb = reff(g); Rw = reff(gd)
        for i, e in enumerate(edges):
            u, w = g.edges[e]; ud, wd = gd.edges[e]
            hgraph = mv.crossing_change(g, e)
            name, how, nred, ur = identify(hgraph)
            detc = fd.det_after([i]); detc = detc if detc < P // 2 else P - detc
            try: sc = int(gl_signature(spherogram.Link([list(q) for q in to_pd(hgraph)]).PD_code()))
            except Exception: sc = None
            rec = {'knot': nm, 'P_u3': P3[nm], 'u': 2, 'det': D, 'crossings': len(edges), 'crossing': i, 'sign': int(signs[i]),
                   'deg_min': min(deg[u], deg[w]), 'deg_max': max(deg[u], deg[w]), 'face_min': min(degd[ud], degd[wd]), 'face_max': max(degd[ud], degd[wd]),
                   'twist_par': par[tuple(sorted((u, w)))], 'twist_ser': pard[tuple(sorted((ud, wd)))], 'twist_max_diagram': twist_max,
                   'det_changed': int(detc), 'reduced_crossings': nred, 'result': name, 'result_how': how, 'u_range': ur,
                   'sigma': sigK, 'sigma_c': sc, 'R_black': str(Rb[e]), 'R_white': str(Rw[e]), 'tait_sign': int(g.signs[e]), 'n_black': len(g.vertices), 'n_white': len(gd.vertices)}
            h.write(json.dumps(rec) + '\n'); h.flush()
