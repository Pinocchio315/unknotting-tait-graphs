#!/usr/bin/env python3
"""Bounded, single-process v1.3 Greene audit; no full sweep is launched.

Checks all frozen premises against the installed KnotInfo rows, all archived
vector domains/conjugation and Ni--Wu exhaustion counters, then replays a few
representative knots. The three new conflicts also receive an independent
rank-two comparison: a complete rectangular characteristic-vector search,
adjugate labels and explicit cyclic isomorphisms, with no production m_Q or
form_admits calls. Run from any directory with the repository environment.
"""
from __future__ import annotations
import ast
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import importlib.metadata
import itertools
import json
import math
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'lower_bounds' / 'greene'))
import database_knotinfo as dk
import sympy
from greene_sweep import run_one, validate_entry
from greene_ranks import candidate_forms, boundary_invariants
from greene_dinv import Diagram
from linkform import compatible, linking_invariants


def independent_rank2(d, D):
    """Use box minima and a separately enumerated cyclic coordinate map.

    Every characteristic class has a minimum with -Qii <= xi_i < Qii:
    subtracting/adding 2Qe_i reduces the norm when a coordinate lies outside
    this interval; at an upper-end equality it preserves norm and reduces
    1^T Q^-1 xi. Thus the finite box contains a minimum in every class.
    """
    reports = []
    for key, raw in candidate_forms(D, 2):
        Q = sympy.Matrix(raw.tolist())
        adj = Q.adjugate()
        assert Q.det() == D
        assert all(Q[:i, :i].det() > 0 for i in range(1, Q.rows + 1))
        labels = lambda v: tuple(int(x) % D for x in adj * sympy.Matrix(v))
        best = {}
        for xi in itertools.product(*(range(-int(Q[i, i]), int(Q[i, i]), 2)
                                      for i in range(Q.rows))):
            x = sympy.Matrix(xi)
            value = (Fraction(int((x.T * adj * x)[0]), D) - Q.rows) / 4
            label = labels(xi)
            if label not in best or value < best[label]:
                best[label] = value
        assert len(best) == D
        generator = next((g for g in best if math.gcd(D, *g) == 1), None)
        if generator is None:
            reports.append({'form': key, 'reason': 'noncyclic discriminant group'})
            continue
        index = {k: tuple(k * x % D for x in generator) for k in range(D)}
        assert len(set(index.values())) == D
        failures = []
        for a in range(1, D):
            if math.gcd(a, D) != 1:
                continue
            for k in range(D):
                difference = best[index[k]] - d[a * k % D]
                if difference < 0 or difference.denominator != 1 or difference.numerator % 2:
                    failures.append({'unit': a, 'class': k, 'difference': str(difference)})
                    break
            else:
                raise AssertionError(f'Independent comparison admitted {key}, unit {a}')
        reports.append({'form': key, 'isomorphisms_excluded': len(failures),
                        'failure_witnesses': failures})
    return reports


def main():
    started = time.time()
    inputs_path = ROOT / 'results/greene/greene_targets_2026-09-08.json'
    raw_path = ROOT / 'results/greene/sweep_2026-09-08.jsonl.gz'
    inputs = json.loads(inputs_path.read_text())
    records = [json.loads(line) for line in gzip.open(raw_path, 'rt')]
    by_name = {r['name']: r for r in records}
    assert len(by_name) == len(records), 'duplicate raw knot records'
    rows = {r['name']: r for r in dk.link_list()}
    premise_counts = Counter()
    for role in ('targets', 'controls'):
        for name, e in inputs[role].items():
            validate_entry(e)
            row = rows[name]
            assert e['pd'] == ast.literal_eval(row['pd_notation']), name
            assert e['determinant'] == int(row['determinant']), name
            assert e['signature'] == int(row['signature']), name
            assert (ast.literal_eval(e['reduced_mod2_vector_raw'])
                    == ast.literal_eval(row['khovanov_reduced_mod2_vector'])), name
            assert len(dict(ast.literal_eval(row['torsion_numbers']))[2]) == 1, name
            premise_counts[role] += 1
    unresolved = Counter()
    for r in records:
        if r['verdict'] in ('ERROR', 'TIMEOUT'):
            continue
        D = r['det']
        if 'd' in r:
            options = {int(k): {Fraction(v)} for k, v in r['d'].items()}
        else:
            options = {int(k): {Fraction(v)} for k, v in r['pinned'].items()}
            options.update({int(k): {Fraction(v) for v in vs}
                            for k, vs in r['ambiguous'].items()})
        assert set(options) == set(range(D))
        assert all(options[k] == options[-k % D] for k in range(D))
        vectors = math.prod(len(options[k]) for k in range(D) if k <= -k % D)
        assert vectors == r['candidate_vectors'] == len(r['candidate_verdicts'])
        assert r['admitting_vectors'] == sum(v['admits'] for v in r['candidate_verdicts'])
        if r['verdict'] == 'OBSTRUCTED':
            assert r['admitting_vectors'] == 0
        for v in r['candidate_verdicts']:
            detail = v['detail']
            if r['test'] == 'u1':
                ph = sum(math.gcd(a, D) == 1 for a in range(D))
                for c in detail['orientations'].values():
                    assert c['tested'] == D * ph
                    assert c['tested'] == sum(c[k] for k in ('invalid_gap', 'inconsistent_V',
                        'invalid_successive_difference', 'surviving'))
        unresolved[r['verdict']] += 1

    # A bounded representative sample, including all three new conflicts and
    # a known positive control for every implemented surgery rank.
    names = ['13n_111', '13n_142', '13n_196', '11n_120', '11n_180',
             '8_20', '9_43', '10_142', '12n_474']
    replay = {}
    independent = {}
    for name in names:
        entry = inputs['targets'].get(name, inputs['controls'].get(name))
        fresh = run_one(entry)
        stored = by_name[name]
        for field in ('verdict', 'candidate_vectors', 'admitting_vectors', 'unresolved_classes'):
            assert fresh[field] == stored[field], (name, field)
        for field in ('d', 'pinned', 'ambiguous'):
            assert fresh.get(field) == stored.get(field), (name, field)
        replay[name] = fresh
        print(name, fresh['verdict'], fresh['seconds'], flush=True)
        if name in ('13n_111', '13n_142', '13n_196'):
            sign = -1 if entry['signature'] < 0 else 1
            d = {int(k): sign * Fraction(v) for k, v in fresh['d'].items()}
            independent[name] = independent_rank2(d, entry['determinant'])

    # Compare the new boundary pre-filter with the old exact matrix-derived
    # invariants on ten positive definite cyclic Goeritz matrices. The filter
    # accepts both signs, so the actual boundary pairing must be retained.
    filter_checks = []
    filter_names = ['3_1', '4_1', '5_1', '7_1', '7_4', '7_3', '7_5',
                    '8_2', '9_10', '9_13']
    for name in filter_names:
        pd = ast.literal_eval(rows[name]['pd_notation'])
        diagram = Diagram(pd)
        classes = [diagram.analyse(x)['label'] for x in diagram.states()]
        _, pairing = diagram.spinc_index(classes)
        D = diagram.D
        sign = 1 if diagram.G[0, 0] > 0 else -1
        Q = sign * diagram.G
        assert all(Q[:i, :i].det() > 0 for i in range(1, Q.rows + 1))
        matrix_invariants = linking_invariants(Q.tolist())
        wanted = boundary_invariants(D, sign * pairing % 1)
        assert matrix_invariants == wanted, (name, sign, matrix_invariants, wanted)
        assert any(compatible(matrix_invariants, boundary_invariants(D, eps * pairing % 1))
                   for eps in (1, -1))
        filter_checks.append(name)
    report = {'input_sha256': hashlib.sha256(inputs_path.read_bytes()).hexdigest(),
              'raw_sha256': hashlib.sha256(raw_path.read_bytes()).hexdigest(),
              'database_knotinfo_version': importlib.metadata.version('database_knotinfo'),
              'implementation_sha256': {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                  for path in [*(ROOT / 'lower_bounds/greene' / name for name in
                      ('greene_dinv.py', 'pin_dinv.py', 'half_integral.py', 'greene_ranks.py', 'greene_sweep.py')),
                      *(ROOT / 'lower_bounds/owens' / name for name in
                      ('owens_obstruction.py', 'owens_u3.py', 'owens_u4.py', 'linkform.py'))]},
              'premises': dict(premise_counts), 'record_counts': dict(unresolved),
              'replayed': replay, 'independent_rank2': independent,
              'linking_filter_matrix_checks': filter_checks,
              'seconds': round(time.time() - started, 1)}
    (HERE / 'greene_code_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print('All checks completed', report['seconds'], flush=True)


if __name__ == '__main__':
    main()
