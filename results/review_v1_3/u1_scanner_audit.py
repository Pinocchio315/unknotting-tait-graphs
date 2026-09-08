#!/usr/bin/env python3
"""Audit the scanner's frozen cohorts and certificates, with a bounded HFK replay.

The determinant calculation uses a separately constructed Goeritz matrix
and its exact adjugate, without xtait or greedy reduction. It checks all
7,052 crossing changes of the 546 open-range chosen diagrams. Direct PD
crossing changes are independently simplified with Spherogram; determinant
one cases receive HFK genus computations within an explicit wall-time budget.
"""
import ast, hashlib, json, sys, time
from pathlib import Path
import sympy
import database_knotinfo as dk
import spherogram
from knot_floer_homology import pd_to_hfk
ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'lower_bounds/greene'))
from greene_dinv import Diagram


def changed(pd, crossing):
    out = [list(c) for c in pd]
    out[crossing] = out[crossing][1:] + out[crossing][:1]
    return out


def identify(pd):
    knot = spherogram.Link(pd)
    assert len(knot.link_components) == 1
    knot.simplify('global')
    if not knot.crossings:
        assert knot.unlinked_unknot_components == 1
        return {'method': 'Spherogram', 'seifert_genus': 0}
    hfk = pd_to_hfk(knot.PD_code(), prime=2)
    assert sum(hfk['ranks'].values()) == hfk['total_rank']
    genus = max(a for a, m in hfk['ranks'])
    assert genus == hfk['seifert_genus']
    return {'method': 'HFK', 'seifert_genus': genus, 'total_rank': hfk['total_rank'],
            'hfk_ranks': [[a, m, r] for (a, m), r in sorted(hfk['ranks'].items())]}


def main():
    start = time.monotonic()
    rows = {r['name']: r for r in dk.link_list()}
    table = json.loads((ROOT / 'generated/u_table.json').read_text())
    paths = [ROOT / 'results/crossing_changes' / f for f in
             ('u1_minimal_diagram_scan_2026-09-08.json', 'u1_open_diagram_scan_2026-09-08.json')]
    known, opened = [json.loads(p.read_text()) for p in paths]
    assert set(known) == {n for n, v in table.items() if v == [1, 1]}
    assert set(opened) == {n for n, v in table.items() if v[0] == 1 and v[1] > 1}
    assert not set(known) & set(opened)
    for n, r in {**known, **opened}.items():
        pd = ast.literal_eval(rows[n]['pd_notation'])
        assert len(pd) == r['crossing_number'] == int(rows[n]['crossing_number'])
        assert r['determinant'] == int(rows[n]['determinant'])
        assert all(type(e) is int and 0 <= e < len(pd) for e in r['unknotting_crossings'])
        assert len(set(r['unknotting_crossings'])) == len(r['unknotting_crossings'])
        assert not r['undecided']
    assert all(r['unknotting_crossings'] and r['verdict'] == 'unknotting crossing' for r in known.values())
    assert all(not r['unknotting_crossings'] and r['verdict'] == 'none in this diagram' for r in opened.values())
    determinants, det1 = {}, []
    for name in opened:
        pd = ast.literal_eval(rows[name]['pd_notation'])
        graph = Diagram(pd, mark=pd[0][0])
        values = []
        for crossing in range(graph.n):
            a, b = graph.wcorners[crossing]
            v = [0] * graph.m
            for sign, region in ((1, graph.cf[crossing, a]), (-1, graph.cf[crossing, b])):
                if region in graph.widx:
                    v[graph.widx[region]] += sign
            change = sum(v[j] * int(graph.adj[j, k]) * v[k]
                         for j in range(graph.m) for k in range(graph.m))
            determinant = abs(graph.det + 2 * graph.mu[crossing] * change)
            values.append(determinant)
            if determinant == 1:
                det1.append((name, crossing))
        determinants[name] = values
    # One explicit witness suffices for each positive existential claim. Replay
    # it directly from the PD, without the scanner's Tait-graph reducer.
    positives = {}
    for name in known:
        edge = known[name]['unknotting_crossings'][0]
        result = identify(changed(ast.literal_eval(rows[name]['pd_notation']), edge))
        assert result['seifert_genus'] == 0, name
        positives[name] = {'crossing': edge, **result}
    negatives = []
    deadline = time.monotonic() + 60
    for name, edge in det1:
        if time.monotonic() > deadline:
            break
        result = identify(changed(ast.literal_eval(rows[name]['pd_notation']), edge))
        assert result['seifert_genus'] > 0, (name, edge, result)
        negatives.append({'name': name, 'crossing': edge, **result})
    report = {'scope': 'one tabulated minimal diagram per knot; no enumeration of other minimal diagrams',
              'known_cohort': len(known), 'open_cohort': len(opened),
              'known_alternating': sum(r['alternating'] for r in known.values()),
              'known_nonalternating': sum(not r['alternating'] for r in known.values()),
              'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'crossing_determinants': determinants,
              'open_crossings': sum(map(len, determinants.values())),
              'open_determinant_one_crossings': len(det1), 'independent_positives': positives,
              'independent_hfk_nontrivial': negatives,
              'hfk_remaining': det1[len(negatives):], 'seconds': round(time.monotonic() - start, 1)}
    (HERE / 'u1_scanner_audit.json').write_text(json.dumps(report, indent=1) + '\n')
    print('known', len(known), 'open', len(opened), 'crossings', report['open_crossings'],
          'det1', len(det1), 'HFK checked', len(negatives), 'remaining', len(report['hfk_remaining']),
          'seconds', report['seconds'])


if __name__ == '__main__':
    main()
