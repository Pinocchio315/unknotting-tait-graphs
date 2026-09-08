"""Second pass over unidentified change results: reduced alternating composites via block decomposition
(rigorous), SnapPy retries, and a Jones/det/Alexander composite hypothesis (flagged non-rigorous)."""
from paths import *
import sys, json, ast, glob, re, itertools, collections
import spherogram, snappy, database_knotinfo as dk, networkx as nx
from xtait.graph import from_pd, to_pd, dual, EmbeddedTaitGraph
from xtait import moves as mv
from xtait.canonical import canonical_code
from xtait.reduce import greedy_reduce
from xtait.jones import jones
def jmul(A, B):
    out = {}
    for i, a in A.items():
        for j, b in B.items(): out[i + j] = out.get(i + j, 0) + a * b
    return {k: v for k, v in out.items() if v}
def jmirror(A): return {-k: v for k, v in A.items()}
from flipdet import exact_det
S = WORK
rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
allk = json.load(open(ALL_CODES)); code2name = {c: x['name'] for x in allk for c in x['codes']}
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
def blocks_split(g):
    """Split a signed plane graph at cut vertices into its blocks (as sub-diagrams)."""
    G = nx.MultiGraph()
    for e in g.edges: u, w = g.edges[e]; G.add_edge(u, w, key=e)
    bl = list(nx.biconnected_component_edges(nx.Graph(G))) if G.number_of_nodes() > 1 else []
    return G, bl
def sub_pd(g, edge_set):
    """PD of the sub-diagram consisting of the given edges: build from the original PD by taking the
    connected-sum factor -> we reconstruct via the signed graph restricted to the block."""
    h = g.copy()
    for e in list(h.edges):
        if e not in edge_set: h.delete_edge(e)   # may not exist; fallback below
    return h
def jones_key(pd): return {k: v for k, v in jones(pd).items() if v}
jcache = {}
def jones_of(nm):
    if nm not in jcache: jcache[nm] = jones_key(ast.literal_eval(rows[nm]['pd_notation']))
    return jcache[nm]
def mirror_pd(pd): return [[q[1], q[2], q[3], q[0]] if False else q for q in pd]
todo = collections.OrderedDict()
for f in sorted(glob.glob(S + '/rows_*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        if r['result'] is None: todo.setdefault((r['knot'], r['crossing']), r)
print('unidentified rows', len(todo)); sys.stdout.flush()
res = {}
classes = {}
for (nm, i), r in todo.items():
    pd = ast.literal_eval(rows[nm]['pd_notation']); g = from_pd([list(q) for q in pd]); e = sorted(g.edges)[i]
    red, _ = greedy_reduce(mv.crossing_change(g, e)); key = repr(canonical_code(red))
    classes.setdefault(key, []).append((nm, i))
    if key in res: continue
    verdict = None
    altern = len({red.signs[x] for x in red.edges}) == 1
    if altern:
        # reduced alternating, not prime-table => must be composite; split at cut vertices of the graph or its dual
        for gg in (red, dual(red)):
            G = nx.Graph()
            for x in gg.edges: G.add_edge(*gg.edges[x])
            comps = [set(c) for c in nx.biconnected_components(G)]
            if len(comps) < 2: continue
            names = []; ok = True
            for vs in comps:
                eids = [x for x in gg.edges if set(gg.edges[x]) <= vs and gg.edges[x][0] != gg.edges[x][1]]
                if not eids: continue
                sub = EmbeddedTaitGraph(sorted(vs), {x: gg.edges[x] for x in eids}, {x: gg.signs[x] for x in eids},
                                        {v: [x for x in gg.rotation[v] if x in set(eids)] for v in vs}, None)
                try:
                    hr, _ = greedy_reduce(sub); nmx = lookup(hr) if hr.n_crossings() else '0_1'
                except Exception as ex:
                    nmx = None
                if nmx is None: ok = False; break
                if nmx != '0_1': names.append(nmx)
            if ok and names:
                verdict = ('sum-blocks', sorted(names)); break
    if verdict is None:
        try:
            M = spherogram.Link(to_pd(red)).exterior()
            hyp = M.solution_type().startswith('all tetrahedra positively')
            for attempt in range(6):
                ids = M.identify()
                for cand in ids:
                    n = cand.name().split('(')[0]
                    if n.startswith('K'):
                        n = n[1:]
                        for s in ('a', 'n'):
                            if s in n: n = n.replace(s, s + '_', 1); break
                    if n in rows and M.is_isometric_to(spherogram.Link(ast.literal_eval(rows[n]['pd_notation'])).exterior()):
                        verdict = ('prime-isometry', n); break
                if verdict: break
                M.randomize()
        except Exception as ex:
            hyp = None
    if verdict is None:
        # Jones/det composite hypothesis (non-rigorous): K = A # B with det(A)det(B) = det, J(A)J(B) = J
        pdr = to_pd(red); D = exact_det(red); Jp = jones_key(pdr); n = red.n_crossings()
        cands = []
        for A in allk:
            if A['det'] == 0 or D % A['det'] or A['crossings'] + 3 > n + 2: continue
            for B in allk:
                if A['name'] > B['name'] or A['det'] * B['det'] != D or A['crossings'] + B['crossings'] > n + 2: continue
                for ma in (0, 1):
                    for mb in (0, 1):
                        JA = jones_of(A['name']); JB = jones_of(B['name'])
                        if ma: JA = jmirror(JA)
                        if mb: JB = jmirror(JB)
                        if jmul(JA, JB) == Jp: cands.append((A['name'], ma, B['name'], mb))
        cands = sorted(set(cands))
        if len({(a, b) for a, _, b, _ in cands}) == 1:
            verdict = ('sum-jones', sorted([cands[0][0], cands[0][2]]))
        elif cands: verdict = ('sum-ambiguous', cands)
    res[key] = {'verdict': verdict, 'alternating_reduced': altern, 'hyperbolic': hyp if 'hyp' in dir() else None, 'n': red.n_crossings()}
    print(key[:40], res[key]); sys.stdout.flush()
json.dump({'classes': {k: v for k, v in classes.items()}, 'results': res}, open(S + '/identify_rest.json', 'w'))
