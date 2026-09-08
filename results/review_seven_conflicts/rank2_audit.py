#!/usr/bin/env python3
"""Bounded independent replay of the five disputed rank-two exclusions.

Run with the repository's topology Python, from any directory. This writes
rank2_audit.json and rank2_audit_full.json.gz beside this script. The independent calculation uses face-based
checkerboard forms, full characteristic-vector boxes, adjugate class labels and
direct generator images. It does not use the production ellipsoid, Smith-form
coordinates, candidate enumeration or isomorphism search for its verdict.

The full-box bound follows by replacing xi with xi +/- 2Qe_i: a global minimum
has |xi_i| <= Q_ii. Boundary ties can be chosen in the half-open box
-Q_ii <= xi_i < Q_ii. Enumeration of this entire finite box therefore suffices.
The comparison with production functions is an additional consistency check.
"""
from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
from importlib.metadata import version
from itertools import product
import json
from math import gcd, isqrt, lcm, prod
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'lower_bounds' / 'owens'))
import crosscheck as independent
import kinfo
import owens_obstruction as production

NAMES = ['13a_15', '13a_55', '13a_422', '13a_568', '13a_650']


def inertia_signature(matrix):
    """Exact rational congruence elimination, allowing 2x2 hyperbolic pivots."""
    a = [[Fraction(x) for x in row] for row in matrix]
    signature = 0
    while a:
        n = len(a)
        diagonal = next((i for i in range(n) if a[i][i]), None)
        if diagonal is not None:
            order = [diagonal] + [i for i in range(n) if i != diagonal]
            a = [[a[i][j] for j in order] for i in order]
            pivot = a[0][0]
            signature += 1 if pivot > 0 else -1
            a = [[a[i][j] - a[i][0] * a[0][j] / pivot
                  for j in range(1, n)] for i in range(1, n)]
        else:
            pair = next(((i, j) for i in range(n) for j in range(i + 1, n)
                         if a[i][j]), None)
            if pair is None:
                break
            order = list(pair) + [i for i in range(n) if i not in pair]
            a = [[a[i][j] for j in order] for i in order]
            pivot = a[0][1]
            a = [[a[i][j] - (a[i][0] * a[1][j] + a[i][1] * a[0][j]) / pivot
                  for j in range(2, n)] for i in range(2, n)]
    return signature


def independent_candidates(determinant):
    """Search both diagonal entries and test whether the residual is a square.

    Reduction gives 0 <= a < m1 <= m2. Thus D >= 4m1-3 and
    m1 <= (D+3)//4. Also m2 <= (D+1)//2 + (D+3)//4 is a generous
    uniform bound obtained from D=(2m1-1)(2m2-1)-4a^2.
    Both m_i are even for the signature-four, two-change case.
    """
    bound1 = (determinant + 3) // 4
    bound2 = (determinant + 1) // 2 + bound1
    answer = []
    for m1 in range(2, bound1 + 1, 2):
        for m2 in range(m1, bound2 + 1, 2):
            remainder = (2 * m1 - 1) * (2 * m2 - 1) - determinant
            if remainder < 0 or remainder % 4:
                continue
            a = isqrt(remainder // 4)
            if 4 * a * a == remainder and a < m1:
                answer.append((m1, m2, a))
    return answer


def full_box(form):
    """Independent exact minima and direct adjugate labels in (Z/D)^rank."""
    determinant = abs(independent.det_frac(form))
    inverse = independent.inv_frac(form)
    rank = len(form)
    adjugate = [[int(determinant * x) for x in row] for row in inverse]
    best = {}
    correspondence = {}
    # This class map is used only to compare with production, never to label
    # the independent minima or to construct an independent isomorphism.
    class_map = production.ClassMap(form)
    for vector in product(*(range(-form[i][i], form[i][i], 2)
                            for i in range(rank))):
        adj_vector = [sum(adjugate[i][j] * vector[j] for j in range(rank))
                      for i in range(rank)]
        label = tuple(x % determinant for x in adj_vector)
        value = Fraction(sum(x * y for x, y in zip(vector, adj_vector))
                         - rank * determinant, 4 * determinant)
        if label not in best or value < best[label]:
            best[label] = value
        production_label = class_map.coords(vector)
        assert correspondence.setdefault(production_label, label) == label
    assert len(best) == determinant
    exact, _ = production.m_Q(form)
    assert exact == {g: best[label] for g, label in correspondence.items()}
    return best


def order(element, determinant):
    return lcm(*(determinant // gcd(x, determinant) for x in element))


def sum_multiple(a, x, b, y, determinant):
    return tuple((a * u + b * v) % determinant for u, v in zip(x, y))


def direct_coordinates(elements, determinant):
    """Find a direct cyclic product with at most two generators, without SNF."""
    zero = (0,) * len(next(iter(elements)))
    first = max(sorted(elements), key=lambda g: order(g, determinant))
    n1 = order(first, determinant)
    n2 = determinant // n1
    seconds = [zero] if n2 == 1 else [g for g in sorted(elements)
                                     if order(g, determinant) == n2]
    for second in seconds:
        coordinates = {sum_multiple(i, first, j, second, determinant): (i, j)
                       for i in range(n1) for j in range(n2)}
        if len(coordinates) == determinant:
            assert set(coordinates) == set(elements)
            return (first, second), (n1, n2), coordinates
    raise AssertionError('no direct product with two generators found')


def all_matchings(source, target, determinant):
    """Record a failing class for every group isomorphism, with no pairing filter."""
    source_orders = Counter(order(g, determinant) for g in source)
    target_orders = Counter(order(g, determinant) for g in target)
    if source_orders != target_orders:
        return {'verdict': 'OBSTRUCTED', 'reason': 'nonisomorphic groups',
                'source_element_orders': dict(sorted(source_orders.items())),
                'target_element_orders': dict(sorted(target_orders.items())),
                'isomorphisms_tested': 0, 'failure_witnesses': []}
    generators, (n1, n2), coordinates = direct_coordinates(source, determinant)
    possible_first = [g for g in sorted(target) if order(g, determinant) == n1]
    possible_second = [g for g in sorted(target) if order(g, determinant) == n2]
    witnesses = []
    for first, second in product(possible_first, possible_second):
        images = {g: sum_multiple(i, first, j, second, determinant)
                  for g, (i, j) in coordinates.items()}
        if len(set(images.values())) != determinant:
            continue
        # Exact orders ensure the direct product's relations; bijectivity
        # completes the isomorphism test. No possible affine translation is
        # omitted: all group maps fix the unique spin structure at zero.
        for g in sorted(source):
            h = images[g]
            difference = source[g] - target[h]
            if difference < 0 or difference.denominator != 1 or difference.numerator % 2:
                witnesses.append({'generator_images': [first, second],
                                  'source_class': g, 'target_class': h,
                                  'm_surgery': str(source[g]),
                                  'd_cover': str(target[h]),
                                  'difference': str(difference)})
                break
        else:
            return {'verdict': 'PASS', 'generator_images': [first, second]}
    return {'verdict': 'OBSTRUCTED', 'reason': 'all correction-term matches fail',
            'source_generators': generators, 'direct_orders': [n1, n2],
            'isomorphisms_tested': len(witnesses), 'failure_witnesses': witnesses}


def encoded(values):
    return [{'class': list(g), 'value': str(value)} for g, value in sorted(values.items())]


def audit_control(name, expected):
    row = kinfo.row(name)
    determinant, sigma = int(row['determinant']), abs(int(row['signature']))
    assert sigma == 4
    target = None
    for form in independent.checkerboard_laplacians(kinfo.pd_code(name)):
        table = full_box(form)
        if table[(0,) * len(form)] == Fraction(-sigma, 4):
            target = table
    assert target is not None
    verdict = 'OBSTRUCTED'
    for m1, m2, a in independent_candidates(determinant):
        form = [[m1, 1, a, 0], [1, 2, 0, 0], [a, 0, m2, 1], [0, 0, 1, 2]]
        result = all_matchings(full_box(form), target, determinant)
        if result['verdict'] == 'PASS':
            verdict = 'PASS'
            break
    assert verdict == expected
    print('Control', name, verdict, 'as expected', flush=True)
    return {'name': name, 'expected': expected, 'independent_verdict': verdict}


def audit(name):
    started = time.monotonic()
    row = kinfo.row(name)
    pd = kinfo.pd_code(name)
    determinant = int(row['determinant'])
    signature = int(row['signature'])
    assert abs(signature) == 4
    # These named links come from Spherogram's built-in diagram table, not
    # from database_knotinfo. Equality of canonical diagram signatures is
    # exact and permits reflection, which does not change unknotting number.
    # The numerical exterior-isometry comparison is an additional check.
    import snappy  # noqa: F401; installs Link.exterior
    import regina
    from spherogram import Link
    supplied, tabulated = Link(pd), Link(name.replace('_', ''))
    supplied_sig = regina.Link.fromPD(supplied.PD_code(min_strand_index=1)).sig(True, True)
    tabulated_sig = regina.Link.fromPD(tabulated.PD_code(min_strand_index=1)).sig(True, True)
    assert supplied_sig == tabulated_sig
    supplied_exterior, tabulated_exterior = supplied.exterior(), tabulated.exterior()
    numerical_isometry = supplied_exterior.is_isometric_to(tabulated_exterior)
    assert numerical_isometry
    seifert = ast.literal_eval(row['seifert_matrix'])
    symmetric = [[seifert[i][j] + seifert[j][i] for j in range(len(seifert))]
                 for i in range(len(seifert))]
    seifert_signature = inertia_signature(symmetric)
    assert seifert_signature == signature
    assert abs(independent.det_frac(symmetric)) == determinant
    forms = independent.checkerboard_laplacians(pd)
    minima = []
    for form in forms:
        assert independent.det_frac(form) == determinant
        assert all(independent.det_frac([row[:r] for row in form[:r]]) > 0
                   for r in range(1, len(form) + 1))
        minima.append(full_box(form))
    assert sorted(minima[0].values()) == sorted(-v for v in minima[1].values())
    spins = [m[(0,) * len(form)] for form, m in zip(forms, minima)]
    chosen = spins.index(Fraction(-abs(signature), 4))
    target = minima[chosen]
    # Relabelled crossings/arcs and the mirror must give the same selected
    # correction-term multiset, despite a changed face enumeration.
    relabelled = [[(x + 10) % 26 + 1 for x in q] for q in reversed(pd)]
    for variant in (relabelled, [[q[1], q[2], q[3], q[0]] for q in pd]):
        variant_tables = [full_box(f) for f in independent.checkerboard_laplacians(variant)]
        assert sorted(sorted(m.values()) for m in variant_tables) == sorted(
            sorted(m.values()) for m in minima)
    tuples = independent_candidates(determinant)
    assert tuples == production.candidates(determinant, 2)
    candidates = []
    for m1, m2, a in tuples:
        form = [[m1, 1, a, 0], [1, 2, 0, 0], [a, 0, m2, 1], [0, 0, 1, 2]]
        assert independent.det_frac(form) == determinant
        source = full_box(form)
        result = all_matchings(source, target, determinant)
        assert result['verdict'] == 'OBSTRUCTED'
        result.update({'coefficients': [m1, m2, a], 'form': form,
                       'correction_terms': encoded(source)})
        candidates.append(result)
    replay = production.obstruct_u2(pd, signature, determinant)
    assert replay == {'verdict': 'OBSTRUCTED', 'candidates': len(tuples)}
    archived = next(r for r in json.loads((ROOT / 'results' / 'owens_rank2' /
                    'owens_verdicts_sigma4_alternating_2026-09-07.json').read_text())
                    if r['name'] == name)
    assert archived['verdict'] == replay['verdict']
    assert archived['n_candidates'] == len(tuples)
    assert sorted(archived['dvals']) == sorted(str(value) for value in target.values())
    result = {'name': name, 'determinant': determinant, 'signature': signature,
              'seifert_signature': seifert_signature,
              'seifert_determinant_matches': True, 'pd': pd,
              'independent_named_diagram': name.replace('_', ''),
              'canonical_diagram_signature': supplied_sig,
              'exact_named_diagram_signature_match': True,
              'numerical_exterior_isometry': numerical_isometry,
              'goeritz_forms': forms, 'goeritz_spins': [str(v) for v in spins],
              'goeritz_box_sizes': [prod(f[i][i] for i in range(len(f))) for f in forms],
              'goeritz_correction_terms': [encoded(m) for m in minima],
              'selected_goeritz': chosen, 'mirror_and_relabelling_checks': True,
              'every_full_box_matches_production_classwise': True,
              'deposited_verdict_candidate_count_and_d_values_match': True,
              'candidate_count': len(candidates), 'candidates': candidates,
              'production_replay': replay, 'independent_verdict': 'OBSTRUCTED',
              'seconds': round(time.monotonic() - started, 3)}
    print(name, 'OBSTRUCTED', 'forms', len(candidates), 'isomorphisms',
          [c['isomorphisms_tested'] for c in candidates], 'seconds', result['seconds'], flush=True)
    return result


if __name__ == '__main__':
    results = [audit(name) for name in NAMES]
    controls = [audit_control(name, expected) for name, expected in
                [('8_2', 'PASS'), ('9_18', 'PASS'), ('9_10', 'OBSTRUCTED')]]
    files = ['lower_bounds/owens/owens_obstruction.py',
             'lower_bounds/owens/crosscheck.py', 'lower_bounds/owens/kinfo.py',
             'results/owens_rank2/owens_verdicts_sigma4_alternating_2026-09-07.json']
    report = {'scope': 'Five disputed signature-four alternating knots: complete rank-two exclusion',
              'algorithm': 'Independent full boxes and all direct generator-image isomorphisms',
              'inputs': 'Installed database_knotinfo PD and Seifert matrices',
              'package_versions': {p: version(p) for p in
                                   ('database_knotinfo', 'spherogram', 'snappy', 'regina')},
              'limitations': 'This audit proves lower bounds only; upper bounds are separate certificates.',
              'source_sha256': {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in files},
              'independent_controls': controls,
              'knots': results}
    destination = Path(__file__).with_suffix('.json')
    complete = json.dumps(report, separators=(',', ':')).encode()
    full_destination = Path(__file__).with_name('rank2_audit_full.json.gz')
    full_destination.write_bytes(gzip.compress(complete, mtime=0))
    # Keep the human-readable report compact; the compressed companion retains
    # every exact correction term and a failing class for every isomorphism.
    report['full_evidence'] = full_destination.name
    report['full_evidence_sha256'] = hashlib.sha256(full_destination.read_bytes()).hexdigest()
    for knot in report['knots']:
        knot.pop('goeritz_correction_terms')
        for candidate in knot['candidates']:
            candidate.pop('correction_terms')
            failures = candidate.pop('failure_witnesses')
            candidate['first_failure_witness'] = failures[0] if failures else None
    destination.write_text(json.dumps(report, indent=2) + '\n')
    print('Wrote', destination, flush=True)
