#!/usr/bin/env python
"""Which knots of unknotting number one show an unknotting crossing in their tabulated minimal diagram?

McCoy's theorem answers this for alternating knots: every alternating diagram of such a knot contains an
unknotting crossing.  Nothing of the kind is known without that hypothesis, and a knot with u(K) = 1 whose
minimal diagrams contain no unknotting crossing satisfies u(K) = 1 < 2 <= u^s_BJ(K) and therefore fails the
strong Bernhard-Jablan equality.  By the descent lemma such a knot sits at the bottom of every chain of
crossing neighbours that a failure of the weak equality produces, so this is the class to look for.

For each knot of unknotting number one, this script changes every crossing of the tabulated diagram.  A
change gives the unknot only if the changed diagram has determinant one; those cases are decided by
simplification and, failing that, by knot Floer homology.  The result for one diagram is not the whole
answer for a non-alternating knot, whose other minimal diagrams are not tested here.

    python u1_minimal_diagram_scan.py                 # every knot with u = 1 in the consolidated table
    python u1_minimal_diagram_scan.py 13n_1587 ...    # named knots
"""
import ast, json, os, sys
from paths import *
import database_knotinfo as dk
from xtait.graph import from_pd, to_pd
from xtait import moves as mv
from xtait.reduce import greedy_reduce
from flipdet import exact_det

rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
table = json.load(open(U_TABLE))
names = [a for a in sys.argv[1:] if not a.startswith('-')]
if not names:
    names = sorted(n for n, v in table.items() if v == [1, 1])
print(f'{len(names)} knots with unknotting number one', flush=True)

out = {}
for nm in names:
    row = rows[nm]
    graph = from_pd([list(q) for q in ast.literal_eval(row['pd_notation'])])
    unknotting, undecided = [], []
    for edge in sorted(graph.edges):
        reduced, _ = greedy_reduce(mv.crossing_change(graph, edge))
        if reduced.n_crossings() == 0:
            unknotting.append(edge); continue
        if exact_det(reduced) != 1:
            continue
        try:
            import spherogram
            from knot_floer_homology import pd_to_hfk
            link = spherogram.Link([tuple(int(x) for x in q) for q in to_pd(reduced)])
            link.simplify('global')
            if not link.crossings or pd_to_hfk(link.PD_code())['seifert_genus'] == 0:
                unknotting.append(edge)
        except Exception as exc:
            undecided.append(edge)
            print(f'  {nm}: crossing {edge} undecided ({exc!r})'[:160], flush=True)
    out[nm] = {'alternating': row['alternating'].upper().startswith('Y'),
               'crossing_number': int(row['crossing_number']), 'determinant': int(row['determinant']),
               'unknotting_crossings': unknotting, 'undecided': undecided,
               'verdict': ('unknotting crossing' if unknotting
                           else 'undecided' if undecided else 'none in this diagram')}
    if not unknotting:
        print(f'  {nm}: {out[nm]["verdict"]}', flush=True)
path = os.path.join(CC, 'u1_minimal_diagram_scan_2026-09-08.json')
json.dump(out, open(path, 'w'), indent=1)
from collections import Counter
print('\n' + str(dict(Counter((v['alternating'], v['verdict']) for v in out.values()))))
print('wrote', path)
