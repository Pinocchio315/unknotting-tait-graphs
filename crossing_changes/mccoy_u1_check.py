#!/usr/bin/env python
"""Apply McCoy's theorem to reduced alternating reference diagrams (paper Sections 2–4).

A determinant-one crossing neighbor is only a candidate unknot. Every such
neighbor must be decided by Reidemeister simplification or knot Floer genus.
Any failed computation leaves the result undecided, never a lower bound.
Importing this module performs no scan and writes no files.
"""
import argparse
import ast
import csv
import io
import json
import os
from pathlib import Path
import re
import sys
import zipfile
from paths import ROOT, RESULTS, CC
from xtait.graph import from_pd, to_pd, pd_passages
from xtait import moves as mv
from xtait.reduce import greedy_reduce
from flipdet import exact_det


def recognize_unknot(pd):
    """Exact genus-zero detection after diagram simplification; errors propagate."""
    from spherogram import Link
    from knot_floer_homology import pd_to_hfk
    link = Link([tuple(int(x) for x in q) for q in pd])
    link.simplify('global')
    if not link.crossings:
        return True
    return pd_to_hfk(link.PD_code())['seifert_genus'] == 0


def check_diagram(pd, unknot_test=None):
    """Decide the unknotting crossings of a nontrivial alternating knot diagram.

    McCoy supplies the quantifier over all diagrams, but only for alternating
    knots. Check the input hypothesis instead of trusting a database flag.
    """
    if not pd:
        raise ValueError('McCoy test expects a nontrivial alternating knot')
    graph = from_pd([list(q) for q in pd])
    pd_passages(pd)  # reject links; the genus-zero theorem used below is for knots
    if len(set(graph.signs.values())) != 1:
        raise ValueError('McCoy test requires an alternating diagram')
    if exact_det(graph) == 1:
        raise ValueError('input diagram represents the unknot, outside the u=1 test')
    test = unknot_test or recognize_unknot
    unknotting, undecided = [], []
    errors = {}
    for edge in sorted(graph.edges):
        try:
            red, _ = greedy_reduce(mv.crossing_change(graph, edge))
            if red.n_crossings() == 0:
                unknotting.append(edge)
            elif exact_det(red) == 1 and test(to_pd(red)):
                unknotting.append(edge)
        except Exception as exc:
            undecided.append(edge)
            errors[str(edge)] = f'{type(exc).__name__}: {exc}'
    verdict = 'u = 1' if unknotting else ('undecided' if undecided else 'u >= 2 (McCoy)')
    return {'unknotting_crossings': unknotting, 'undecided': undecided,
            'verdict': verdict, 'errors': errors}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all', metavar='DIR', help='use the archived April2026 wheel in DIR, including u=1 controls')
    parser.add_argument('--reference', type=Path, default=Path(RESULTS) / 'paper_v1_1_snapshot.json')
    parser.add_argument('--out', type=Path, help='report path (default: generated/mccoy_alternating_u1.json, or mccoy_alternating_all.json with --all)')
    parser.add_argument('--names', nargs='+', help='optional bounded subset of knot names')
    args = parser.parse_args(argv)
    if args.out is None:
        filename = 'mccoy_alternating_all.json' if args.all else 'mccoy_alternating_u1.json'
        args.out = Path(ROOT) / 'generated' / filename
    import database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    if args.all:
        csv.field_size_limit(sys.maxsize)
        wheel = Path(args.all) / 'database_knotinfo-2026.4.1-py3-none-any.whl'
        ref = {}
        with zipfile.ZipFile(wheel) as archive:
            with archive.open('database_knotinfo/csv_data/knotinfo_data_complete.csv') as raw:
                reader = csv.reader(io.TextIOWrapper(raw, encoding='utf-8'), delimiter='|')
                header = [x.strip() for x in next(reader)]
                next(reader)
                ni, ui = header.index('name'), header.index('unknotting_number')
                for row in reader:
                    bounds = re.findall(r'\d+', row[ui]) if len(row) > ui else []
                    if bounds and row[ni].strip() in rows:
                        ref[row[ni].strip()] = [int(bounds[0]), int(bounds[-1])]
    else:
        with args.reference.open() as handle:
            ref = json.load(handle)
    targets = sorted(name for name, bounds in ref.items() if name in rows and name != '0_1'
                     and rows[name]['alternating'].upper().startswith('Y')
                     and (bounds[0] <= 1 if args.all else bounds[0] == 1 and bounds[1] > 1))
    if args.names:
        targets = [name for name in targets if name in args.names]
    if not targets:
        raise SystemExit('No matching alternating knots; no verification performed.')
    print(f'{len(targets)} alternating knots to check', flush=True)
    out = {}
    for name in targets:
        result = check_diagram(ast.literal_eval(rows[name]['pd_notation']))
        result['reference_range'] = ref[name]
        out[name] = result
        print(name, ref[name], result['verdict'], flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w') as handle:
        json.dump(out, handle, indent=1)
    print(f'Wrote {len(out)} results to {args.out}')
    if any(record['undecided'] for record in out.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
