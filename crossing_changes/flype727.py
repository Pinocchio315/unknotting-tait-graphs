"""For the target knots whose KnotInfo diagram has no u=2 candidate crossing: check every minimal diagram
(complete flype orbit).  A candidate = crossing whose change gives |sigma| <= 2 and a knot with u possibly 1."""
from paths import *
from candidates import candidate_status as status_of
import os
import sys, json, ast, time
S = os.path.join(OPEN23, 'work')
src = open(os.path.join(HERE, 'apply833.py')).read()
head = src.split("with open(S + f'/rows833_{si}.jsonl', 'a') as h:")[0]
exec(head)
from xtait.flypes import flype_orbit
OUT = OPEN23
pri = json.load(open(OUT + '/priority_u23.json'))
zero = sorted(k['knot'] for k in pri['knots'] if k['n_candidates'] == 0)[si::sn]
with open(S + f'/flype727_{si}.jsonl', 'a') as out:
    for nm in zero:
        t0 = time.time(); r = rows[nm]; g0 = from_pd([list(q) for q in ast.literal_eval(r['pd_notation'])])
        orbit = flype_orbit(g0); cands = []; n_cross = 0; heur = 0; stat = {}
        for di, (code, g) in enumerate(orbit.items()):
            for e in sorted(g.edges):
                n_cross += 1; hg = mv.crossing_change(g, e)
                try: sc = int(gl_signature(spherogram.Link([list(q) for q in to_pd(hg)]).PD_code()))
                except Exception: sc = None
                if sc is None: name, how, nred, ur = identify(hg); st = status_of(name, how, ur, 2)   # unknown signature: treat as sigma-dropping
                elif abs(sc) > 2: st = 'excluded:sigma'; name = how = ur = None
                else: name, how, nred, ur = identify(hg); st = status_of(name, how, ur, sc)
                stat[st] = stat.get(st, 0) + 1
                if st.startswith('candidate') or st == 'WITNESS': cands.append({'diagram': di, 'crossing': e, 'result': name, 'how': how, 'u_range': ur, 'status': st})
                if st == 'excluded:u>=2(jones)': heur += 1
        out.write(json.dumps({'knot': nm, 'P_u3': P3[nm], 'diagrams': len(orbit), 'crossings_total': n_cross, 'status': stat, 'heuristic_exclusions': heur, 'candidates': cands, 'seconds': round(time.time() - t0, 1)}) + '\n'); out.flush()
