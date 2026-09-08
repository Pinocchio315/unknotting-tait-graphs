#!/usr/bin/env python
"""Test every crossing of one tabulated minimal diagram per selected knot.

A positive result gives an unknotting crossing in that diagram. A negative
result says nothing about other minimal diagrams of a non-alternating knot.
Determinant one is only a filter: remaining diagrams are decided by explicit
simplification or knot Floer homology, which detects Seifert genus.

With no --out argument results are printed and no archive is overwritten.
Use --cohort open for the unresolved ranges with lower endpoint one.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

from paths import U_TABLE
from xtait.graph import from_pd, to_pd
from xtait import moves as mv
from xtait.reduce import greedy_reduce
from flipdet import exact_det


def hfk_genus(result):
    """Reject incomplete or internally inconsistent HFK output before use.

    Ozsvath--Szabo's genus-detection theorem identifies the maximal Alexander
    grading supporting HFK with Seifert genus (math/0311496). In S^3, genus
    zero is equivalent to the unknot. No analogous inference follows merely
    from tau=0, signature=0 or determinant one.
    """
    genus, ranks = result['seifert_genus'], result['ranks']
    if type(genus) is not int or genus < 0 or result['modulus'] != 2 or not ranks:
        raise ValueError('HFK did not return a nonnegative integer genus over F_2')
    if any(not isinstance(grade, tuple) or len(grade) != 2
           or any(type(k) is not int for k in grade)
           or type(rank) is not int or rank <= 0 for grade, rank in ranks.items()):
        raise ValueError('HFK ranks or gradings are malformed')
    if (sum(ranks.values()) != result['total_rank']
            or max(a for a, _ in ranks) != genus):
        raise ValueError('HFK genus and rank table disagree')
    if genus == 0 and result['total_rank'] != 1:
        raise ValueError('Genus-zero HFK output does not have the rank of the unknot')
    return genus


def scan_row(row):
    """One specified PD, with an evidence record for every crossing.

    Edge ids of from_pd are the original zero-based PD-row indices. The
    number of PD rows must equal the tabulated crossing number; minimality
    here is inherited from that tabulated invariant, not proved by the scan.
    """
    import spherogram
    pd = [list(q) for q in ast.literal_eval(row['pd_notation'])]
    if len(pd) != int(row['crossing_number']):
        raise ValueError('The supplied PD does not have the tabulated minimal crossing number')
    if len(spherogram.Link(pd).link_components) != 1:
        raise ValueError('A one-component knot diagram is required')
    graph = from_pd(pd)
    if exact_det(graph) != int(row['determinant']):
        raise ValueError('The supplied PD and tabulated determinant disagree')
    unknotting, undecided, evidence = [], [], {}
    for edge in sorted(graph.edges):
        try:
            changed = mv.crossing_change(graph, edge)
            determinant = exact_det(changed)
            record = {'determinant': determinant}
            evidence[edge] = record
            if determinant != 1:
                record['method'] = 'determinant excludes unknot'
                continue
            reduced, path = greedy_reduce(changed)
            record['reduction_path'] = path
            if reduced.n_crossings() == 0:
                record['method'] = 'explicit Tait-graph reduction to zero crossings'
                unknotting.append(edge)
                continue
            link = spherogram.Link([tuple(int(x) for x in q) for q in to_pd(reduced)])
            if len(link.link_components) != 1:
                raise ValueError('The reduction did not preserve one component')
            link.simplify('global')
            if not link.crossings:
                if link.unlinked_unknot_components != 1:
                    raise ValueError('The reduced link is not a single crossing-free component')
                record['method'] = 'Spherogram simplification to zero crossings'
                unknotting.append(edge)
                continue
            from knot_floer_homology import pd_to_hfk
            computed = pd_to_hfk(link.PD_code(), prime=2)
            genus = hfk_genus(computed)
            record.update(method='HFK Seifert genus', seifert_genus=genus,
                          hfk_input_pd=link.PD_code(), total_rank=computed['total_rank'],
                          hfk_ranks=[[a, m, rank] for (a, m), rank in sorted(computed['ranks'].items())])
            if genus == 0:
                unknotting.append(edge)
        except Exception as exc:
            undecided.append(edge)
            evidence.setdefault(edge, {})['error'] = repr(exc)
    return {'alternating': row['alternating'].upper().startswith('Y'),
            'crossing_number': len(pd), 'determinant': int(row['determinant']),
            'pd': pd, 'scope': 'one tabulated minimal diagram',
            'unknotting_crossings': unknotting, 'undecided': undecided,
            'crossing_evidence': evidence,
            'verdict': ('unknotting crossing' if unknotting
                        else 'undecided' if undecided else 'none in this diagram')}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('knots', nargs='*')
    ap.add_argument('--cohort', choices=('known', 'open'), default='known')
    ap.add_argument('--table', type=Path, default=Path(U_TABLE))
    ap.add_argument('--out', type=Path, help='Write this explicit output file; no file is changed by default')
    args = ap.parse_args()
    import database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list()
            if str(r.get('crossing_number', '')).strip().isdigit()}
    table = json.loads(args.table.read_text())
    names = sorted(set(args.knots)) if args.knots else sorted(
        n for n, interval in table.items()
        if (interval == [1, 1] if args.cohort == 'known'
            else interval[0] == 1 and interval[1] > 1))
    unknown = set(names) - rows.keys()
    if unknown:
        ap.error('Unknown knot names: ' + ' '.join(sorted(unknown)))
    out = {name: scan_row(rows[name]) for name in names}
    payload = {'schema_version': 2, 'scope': 'one tabulated minimal diagram per knot',
               'source': {'database': 'KnotInfo',
                          'package_version': importlib.metadata.version('database_knotinfo'),
                          'ranges_sha256': hashlib.sha256(args.table.read_bytes()).hexdigest()},
               'counts': dict(Counter(r['verdict'] for r in out.values())), 'knots': out}
    text = json.dumps(payload, indent=1) + '\n'
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
        print(f'{len(out)} knots: {payload["counts"]}; wrote {args.out}')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
