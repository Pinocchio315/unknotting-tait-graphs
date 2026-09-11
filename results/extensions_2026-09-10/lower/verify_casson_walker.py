#!/usr/bin/env python3
"""Refine Greene candidates using the Casson--Walker sum of correction terms.

This script certifies u(13n1619)=2 and supplies the sum-constraint helpers
used by the rank-four extensions in manuscript v1.8. Reruns write to generated/,
not to the frozen certificates or sweep. All arithmetic
in the new constraint and the surgery obstruction uses exact fractions.

Dependencies: numpy, sympy, regina; run from any working directory.
The computation modules and frozen inputs are found relative to this file.
See README.md for the hypotheses, normalization, sources, and limitations.
"""
from __future__ import annotations

import ast
import argparse
import hashlib
import itertools
import json
import math
import platform
import sys
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE = HERE.parents[2]
ROOT = CODE
OUTPUT = CODE / 'generated/v1_8_extension_recalculation'
GREENE = CODE / 'lower_bounds' / 'greene'
sys.path[:0] = [str(CODE), str(GREENE)]

import regina
import consolidate_results as consolidate
from greene_sweep import validate_entry
from greene_ranks import oriented, rank_test
from half_integral import surgery_test
from pin_dinv import khovanov_lspace_evidence, pin

ARCHIVE = CODE / 'results/knotinfo_calculation_inputs_2026-09-09.json.gz'
INPUTS = [
    ARCHIVE,
    CODE / 'results/greene/greene_targets_2026-09-08.json',
    CODE / 'results/greene/sweep_2026-09-08.jsonl.gz',
    CODE / 'consolidate_results.py',
    *[GREENE / name for name in (
        'pin_dinv.py', 'greene_dinv.py', 'greene_sweep.py',
        'greene_ranks.py', 'half_integral.py')],
    *sorted((CODE / 'lower_bounds/owens').glob('*.py')),
    Path(__file__).resolve(),
]
SOURCES = {
    'Ozsvath_Szabo_branched_cover_spectral_sequence': 'https://arxiv.org/abs/math/0309170',
    'Greene': 'https://arxiv.org/abs/0805.1381',
    'Rustamov_Theorem_3_3': 'https://arxiv.org/abs/math/0409294',
    'Mullins_formula_in_Greene_Watson_Theorem_13': 'https://arxiv.org/abs/1106.5559',
    'spin_identity_Lin_Ruberman_Saveliev': 'https://arxiv.org/abs/1802.07704',
    'Ni_Wu_Proposition_1_6': 'https://arxiv.org/abs/1009.4720',
}


def require(condition, message):
    # Unlike assert, this guard is not disabled by python -O.
    if not condition:
        raise ValueError(message)


def hashes():
    return {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS}


def dump(name, value):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / name).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def clean(poly):
    return {k: v for k, v in poly.items() if v}


def add(left, right):
    out = dict(left)
    for k, value in right.items():
        out[k] = out.get(k, F(0)) + value
    return clean(out)


def multiply(left, right):
    out = {}
    for a, x in left.items():
        for b, y in right.items():
            out[a + b] = out.get(a + b, F(0)) + x * y
    return clean(out)


def power(poly, exponent):
    if exponent < 0:
        require(len(poly) == 1, 'A Laurent denominator must be a monomial')
        k, value = next(iter(poly.items()))
        return {k * exponent: value ** exponent}
    out = {0: F(1)}
    for _ in range(exponent):
        out = multiply(out, poly)
    return out


def parse_jones(expression):
    """Parse a Laurent polynomial without eval/sympify or executable names.

    The archived expressions use integer coefficients, t, and arithmetic.
    Restricting denominators to monomials gives an exact coefficient map and
    prevents silently treating an arbitrary rational function as a polynomial.
    """
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return clean({0: F(node.value)})
        if isinstance(node, ast.Name) and node.id == 't':
            return {1: F(1)}
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            sign = -1 if isinstance(node.op, ast.USub) else 1
            return {k: sign * v for k, v in visit(node.operand).items()}
        if isinstance(node, ast.BinOp):
            left, right = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(left, right)
            if isinstance(node.op, ast.Sub):
                return add(left, {k: -v for k, v in right.items()})
            if isinstance(node.op, ast.Mult):
                return multiply(left, right)
            if isinstance(node.op, ast.Div):
                return multiply(left, power(right, -1))
            if isinstance(node.op, ast.Pow):
                require(set(right) <= {0}, 'Nonconstant polynomial exponent')
                exponent = right.get(0, F(0))
                require(exponent.denominator == 1 and abs(exponent) <= 1000,
                        'Nonintegral or unreasonable polynomial exponent')
                return power(left, int(exponent))
        raise ValueError(f'Unsupported Jones expression: {ast.dump(node)}')
    return visit(ast.parse(expression.replace('^', '**'), mode='eval').body)


def jones_data(poly, determinant, signature):
    value = sum(v * F(-1) ** k for k, v in poly.items())
    derivative = sum(k * v * F(-1) ** (k - 1) for k, v in poly.items())
    require(abs(value) == determinant, 'Jones determinant does not match the input')
    # Rustamov's Casson normalization: lambda(boundary of negative E8) = -1.
    # Greene--Watson print twice this lambda. Keeping the factor explicit
    # avoids combining their normalization with Rustamov's by accident.
    casson_walker = F(signature, 8) - derivative / (12 * value)
    required_sum = -2 * determinant * casson_walker
    return {
        'V_at_minus_one': str(value), 'V_derivative_at_minus_one': str(derivative),
        'lambda_Casson_normalization': str(casson_walker),
        'required_sum_d': str(required_sum),
        'coefficients_t': {str(k): str(v) for k, v in sorted(poly.items())},
    }


def metadata():
    from knotinfo_inputs import load_metadata
    return load_metadata()


def candidate_sets(record):
    if 'd' in record:
        return {int(k): {F(v)} for k, v in record['d'].items()}
    out = {int(k): {F(v)} for k, v in record['pinned'].items()}
    out.update({int(k): set(map(F, values)) for k, values in record['ambiguous'].items()})
    return out


def vectors(sets, determinant):
    """Enumerate all candidates satisfying conjugation, in frozen sweep order."""
    require(set(sets) == set(range(determinant)), 'Missing Spin^c class')
    require(all(sets[k] and sets[k] == sets[-k % determinant] for k in sets),
            'Candidate sets are empty or violate conjugation')
    pinned = {k: next(iter(values)) for k, values in sets.items() if len(values) == 1}
    free = sorted(k for k, values in sets.items()
                  if len(values) > 1 and k <= -k % determinant)
    for choice in itertools.product(*(sorted(sets[k]) for k in free)):
        vector = dict(pinned)
        for k, value in zip(free, choice):
            vector[k] = vector[-k % determinant] = value
        yield vector


def filter_vectors(sets, determinant, signature, required_sum):
    total, after_spin, survivors = 0, 0, []
    for i, vector in enumerate(vectors(sets, determinant)):
        total += 1
        if vector[0] != F(-signature, 4):
            continue
        after_spin += 1
        if sum(vector.values()) == required_sum:
            survivors.append((i, vector))
    # Empty candidates are an inconsistent computation, never an obstruction.
    require(survivors, 'The invariant constraints eliminated every candidate')
    return total, after_spin, survivors


def regina_jones(pd):
    link = regina.Link.fromPD(pd)
    require(link.countComponents() == 1, 'PD is not a knot')
    polynomial = link.jones()  # Regina uses sqrt(t), not t.
    out = {}
    for exponent in range(polynomial.minExp(), polynomial.maxExp() + 1):
        value = int(str(polynomial[exponent]))
        if value:
            require(exponent % 2 == 0, 'Unexpected half-integral Jones exponent')
            out[exponent // 2] = F(value)
    return out


def independent_upper_bound(pd):
    """Verify two explicit witnesses by actual Reidemeister simplification.

    A Jones polynomial equal to 1 is not an unknot certificate. We require a
    one-component, zero-crossing final diagram. Regina simplify uses type I,
    II and III moves. Because its type III selection can be random, make up
    to 20 independent attempts; failure is inconclusive, never nonexistence.
    """
    witnesses = []
    for indices in ([0, 8], [0, 12]):
        reduced = False
        for attempt in range(1, 21):
            link = regina.Link.fromPD(pd)
            for index in indices:
                link.change(link.crossing(index))
            changed_pd = link.pdData()
            link.simplify()
            if link.size() == 0 and link.countComponents() == 1:
                reduced = True
                break
        require(reduced, f'Regina did not certify the witness {indices}')
        witnesses.append({
            'crossing_indices_zero_based': indices, 'changed_pd': changed_pd,
            'attempt': attempt, 'final_crossings': link.size(),
            'final_components': link.countComponents(),
        })
    return {'upper_bound': 2, 'method': 'Regina simplify: Reidemeister I/II/III',
            'witnesses': witnesses}


def main():
    initial_hashes = hashes()
    meta = metadata()
    frozen, records = consolidate.read_greene_sweep()
    _, table = consolidate.reconstruct(manuscript='v1.3')

    # First check the normalization on every complete vector in the old sweep.
    checked = []
    for name, record in records.items():
        if 'd' not in record:
            continue
        poly = parse_jones(meta[name]['jones_polynomial'])
        data = jones_data(poly, record['det'], record['sigma'])
        actual = sum(F(v) for v in record['d'].values())
        require(actual == F(data['required_sum_d']), f'Sum mismatch on {name}')
        checked.append(name)
    print(f'Existing exact vectors: {len(checked)} sum identities agree', flush=True)

    # Screen only knots still open in the reconstructed paper table. This
    # does not count further information about already determined knots as
    # additional unknotting numbers. Old verdicts are used for prioritization;
    # every claimed new obstruction is separately replayed below.
    screen = []
    for name, record in records.items():
        if ('ambiguous' not in record or record['verdict'] not in ('PASS', 'UNDECIDED')
                or table[name][0] == table[name][1]):
            continue
        data = jones_data(parse_jones(meta[name]['jones_polynomial']),
                          record['det'], record['sigma'])
        total, spin_count, survivors = filter_vectors(
            candidate_sets(record), record['det'], record['sigma'], F(data['required_sum_d']))
        require(total == record['candidate_vectors'], 'Candidate enumeration changed')
        screen.append({
            'name': name, 'range': table[name], 'test': record['test'],
            'old_status': record['verdict'], 'candidate_vectors': total,
            'spin_survivors': spin_count, 'sum_and_spin_survivors': len(survivors),
            'required_sum': data['required_sum_d'],
            'admitting_survivors': sum(record['candidate_verdicts'][i]['admits']
                                       for i, _ in survivors),
        })

    selected = {}
    for name in ('13n_1619', '13n_2731'):
        entry, old = frozen['targets'][name], records[name]
        evidence = validate_entry(entry)
        pd, determinant, signature = entry['pd'], entry['determinant'], entry['signature']
        require(pd == ast.literal_eval(meta[name]['pd_notation']), 'Archived PD changed')
        poly = regina_jones(pd)
        require(poly == parse_jones(meta[name]['jones_polynomial']),
                f'Independent Jones polynomial mismatch on {name}')
        data = jones_data(poly, determinant, signature)
        diagnostics = {}
        fresh, pairing, reference = pin(name, pd, determinant, verbose=False,
                                       diagnostics=diagnostics)
        require(fresh == candidate_sets(old), 'Fresh Greene candidates differ from the sweep')
        total, spin_count, survivors = filter_vectors(
            fresh, determinant, signature, F(data['required_sum_d']))
        require(len(survivors) == 1, 'A unique correction-term vector was expected')
        index, vector = survivors[0]
        if entry['test'] == 'u1':
            # Here there is a shorter proof of uniqueness, requiring neither
            # enumeration nor the separate spin identity: each actual d(k)
            # is at least min(C_k), and these minima already have the required
            # total. Hence every one of the lower bounds is attained.
            require(sum(min(v) for v in fresh.values()) == F(data['required_sum_d']),
                    'The componentwise minima do not attain the required sum')
            require(all(vector[k] == min(fresh[k]) for k in fresh),
                    'The selected vector is not the componentwise minimum')
            outcome = surgery_test(vector, determinant)
            expected = 2 * determinant * sum(math.gcd(a, determinant) == 1
                                             for a in range(1, determinant))
            require(sum(v['tested'] for v in outcome['orientations'].values()) == expected,
                    'Incomplete affine search')
            require(not outcome['fits'], 'New vector admits a surgery fit')
            upper = independent_upper_bound(pd)
            result = {'unknotting_number': 2, 'new_relative_to_paper': True,
                      'lower_bound': 2, 'upper_bound_verification': upper,
                      'uniqueness_proof': 'componentwise minima attain the Casson--Walker sum'}
        else:
            outcome = rank_test(oriented(vector, signature), determinant, 2, pairing)
            require(outcome['verdict'] == 'PASS', '13n2731 changed its surviving verdict')
            result = {'range': table[name], 'new_relative_to_paper': False}
        result.update({
            'name': name, 'input': entry, 'lspace_evidence': evidence,
            'Jones_and_Casson_Walker': data, 'Greene_diagnostics': diagnostics,
            'reference_marking': reference, 'pairing': str(pairing),
            'candidate_sets': {str(k): list(map(str, sorted(v))) for k, v in fresh.items()},
            'old_candidate_vectors': total, 'after_spin_constraint': spin_count,
            'after_spin_and_sum_constraints': 1, 'frozen_vector_index': index,
            'minimum_candidate_sum': str(sum(min(v) for v in fresh.values())),
            'maximum_candidate_sum': str(sum(max(v) for v in fresh.values())),
            'd': {str(k): str(vector[k]) for k in range(determinant)},
            'surgery_test': outcome,
        })
        selected[name] = result
        print(f'{name}: unique vector {index}; {result.get("unknotting_number", result.get("range"))}',
              flush=True)

    # Known u=1 controls must still admit the permissive Ni--Wu test after
    # applying the new constraints. Their Jones polynomials and Greene sets
    # are independently recomputed, rather than copied from the sweep.
    controls = []
    for name in ('3_1', '4_1', '5_2', '8_20', '9_44'):
        row = meta[name]
        determinant, signature = int(row['determinant']), int(row['signature'])
        evidence = khovanov_lspace_evidence(row)
        pd = ast.literal_eval(row['pd_notation'])
        poly = regina_jones(pd)
        require(poly == parse_jones(row['jones_polynomial']), 'Control Jones mismatch')
        data = jones_data(poly, determinant, signature)
        current, _, _ = pin(name, pd, determinant, verbose=False)
        total, _, survivors = filter_vectors(current, determinant, signature,
                                             F(data['required_sum_d']))
        outcomes = [surgery_test(vector, determinant) for _, vector in survivors]
        require(any(outcome['fits'] for outcome in outcomes), f'False obstruction on {name}')
        controls.append({'name': name, 'lspace_evidence': evidence,
                         'candidates': total, 'survivors': len(survivors),
                         'surgery_tests': outcomes})

    require(hashes() == initial_hashes, 'An input or implementation changed during the run')
    new_names = [n for n, r in selected.items() if r['new_relative_to_paper']]
    summary = {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'python': platform.python_version(), 'regina': regina.versionString(),
        'sources': SOURCES, 'sha256': initial_hashes,
        'paper_open_knots_before': sum(a != b for a, b in table.values()),
        'paper_open_range_counts': {str(k): v for k, v in sorted(Counter(
            tuple(interval) for interval in table.values() if interval[0] != interval[1]).items())},
        'new_exact_knots': new_names,
        'paper_open_knots_after_these_results': sum(a != b for a, b in table.values()) - len(new_names),
        'existing_exact_vector_sum_checks': len(checked), 'sum_mismatches': 0,
        'existing_exact_vectors_checked': sorted(checked),
        'open_ambiguous_screen': sorted(screen, key=lambda r: r['name']),
        'known_u1_controls': controls,
        'limitations': [
            'Khovanov ranks are checked from archived mod-2 homology vectors, not recomputed from chains.',
            'Newness is relative to the frozen paper results; this is not a literature-wide priority claim.',
            'PASS is inconclusive and does not give an unknotting sequence.',
            'The original paper, generated counts, and frozen sweep are not modified.',
        ],
    }
    for name, result in selected.items():
        result['sources'] = SOURCES
        result['input_and_code_sha256'] = initial_hashes
        dump(name + '_certificate.json', result)
    dump('summary.json', summary)
    dump('open_candidate_screen.json', sorted(screen, key=lambda r: r['name']))
    print(f'Complete: {new_names}; {len(controls)} known u=1 controls pass', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=OUTPUT)
    OUTPUT = parser.parse_args().outdir
    main()
