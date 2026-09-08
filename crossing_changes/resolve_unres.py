"""Retry unidentified crossing neighbors, retaining unresolved and heuristic cases.

Failure of a hyperbolic solver does not prove nonhyperbolicity or compositeness;
prime satellite knots also exist. Only an explicit identification can certify
the partner. Jones-product guesses remain in the separate heuristic tier.
"""
from paths import *
import os
import sys, json, ast, collections
S = os.path.join(OPEN23, 'work')
src = open(os.path.join(HERE, 'apply833.py')).read(); head = src.split("with open(S + f'/rows833_{si}.jsonl', 'a') as h:")[0].replace("si, sn = (int(x) for x in sys.argv[1].split('/'))", "si, sn = 0, 1")
exec(head)
from xtait.flypes import flype_orbit
OUT = OPEN23
fin = json.load(open(OUT + '/priority_u23_final.json'))
unres = fin['unresolved_candidates']
tierA = set(fin['knots_with_candidates'])
TORUS = {'3_1', '5_1', '7_1', '9_1', '11a_367', '13a_4878', '8_19', '10_124'}
res = []; cache = {}
orbits = {}
for nm, c in unres:
    g0 = from_pd([list(q) for q in ast.literal_eval(rows[nm]['pd_notation'])])
    if 'diagram' in c:
        if nm not in orbits: orbits[nm] = list(flype_orbit(g0).values())
        g = orbits[nm][c['diagram']]; e = c['crossing']
    else:
        g = g0; e = sorted(g.edges)[c['crossing']]
    hg = mv.crossing_change(g, e); red, _ = greedy_reduce(hg); key = repr(canonical_code(red))
    if key not in cache:
        pdr = to_pd(red); M = spherogram.Link([list(q) for q in pdr]).exterior(); hyp = None; name = None
        for attempt in range(8):
            st = M.solution_type()
            if st.startswith('all tetrahedra positively'):
                hyp = True
                for cand in M.identify():
                    n = ki_name(str(cand))
                    if n in rows:
                        try:
                            if M.is_isometric_to(ext(n)): name = n; break
                        except Exception: pass
                if name: break
            M.randomize()
        if name: verdict = ('prime', name, RANGE.get(name))
        elif c['status'].startswith('candidate:composite?'): verdict = ('composite(jones)', None, [2, 99])   # Jones = product, det multiplicative: composite => u >= 2 (Scharlemann); heuristic
        elif hyp is None: verdict = ('non-hyperbolic-unidentified', None, None)
        else: verdict = ('hyperbolic-unidentified', None, None)
        cache[key] = verdict
    res.append({'knot': nm, 'cand': c, 'verdict': cache[key]})
cnt = collections.Counter(r['verdict'][0] for r in res); print('resolution:', dict(cnt))
print('prime identifications:', collections.Counter((r['verdict'][1], str(r['verdict'][2])) for r in res if r['verdict'][0] == 'prime').most_common(20))
json.dump(res, open(OUT + '/unresolved_candidates_resolution.json', 'w'), indent=1)
# recompute tiers: a knot in tier A/A2 whose only candidates were unresolved and are now all composite/u>=2 moves to tier B (heuristic)
still = collections.defaultdict(list)
for r in res:
    v = r['verdict']
    if v[0] == 'prime' and v[2] and v[2][0] <= 1: still[r['knot']].append((v[1], v[2]))
    elif v[0] in ('non-hyperbolic-unidentified', 'hyperbolic-unidentified'): still[r['knot']].append(('unresolved:' + v[0], None))
fl = {}
for l in open(OUT + '/flype_orbit_check.jsonl'): r = json.loads(l); fl[r['knot']] = r
pri = json.load(open(OUT + '/priority_u23.json')); knots = {k['knot']: k for k in pri['knots']}
def real_cands(nm):
    out = []
    if nm in knots: out += [c for c in knots[nm]['candidates'] if c['status'] == 'candidate']
    if nm in fl: out += [c for c in fl[nm]['candidates'] if c['status'] == 'candidate']
    return out
moved = [nm for nm in set(fin['knots_with_candidates']) if not real_cands(nm) and not still.get(nm)]
fin['unresolved_after_resolution'] = {nm: v for nm, v in still.items() if all(str(x[0]).startswith('unresolved') for x in v)}
newprime = {nm: v for nm, v in still.items() if any(x[0] != 'unresolved' for x in v)}
print(f'knots moved to tier B after resolution (all their candidates were composites): {len(moved)}')
print('knots with newly identified prime [1,2] children:', {k: v for k, v in newprime.items()})
print('knots keeping unresolved candidates:', {nm: v for nm, v in still.items() if all(str(x[0]).startswith('unresolved') for x in v)})
fin['moved_to_tierB_after_resolution'] = sorted(moved); fin['new_children_after_resolution'] = {k: v for k, v in newprime.items()}
json.dump(fin, open(OUT + '/priority_u23_final.json', 'w'), indent=1)
