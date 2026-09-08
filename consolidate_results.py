#!/usr/bin/env python3
"""Reconstruct the manuscript result tables from explicitly frozen, deposited inputs.

This is bookkeeping, not a rerun of the mathematical obstructions.  A completed
OBSTRUCTED verdict supplies a lower bound; PASS only means that this obstruction
is silent.  In particular, PASS, a missing record, or a timeout never proves an
upper bound. Upper bounds come from separately verifiable diagram certificates
or explicitly cited published constructions.

Run from any directory: python consolidate_results.py --outdir generated
Historical results are read-only. No installed KnotInfo table or live worker log
is consulted, and conflicting bounds stop the run before any output is written.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import gzip
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
PAPER_SCAN = 'open23/priority_u23_paper_v1_3.json'
PRIORITY = ('L', 't', 'a', 'OT', 'O2', 'g', 'c', 'M', 'G', 'K', 'O3', 'O4', 'w', 'b', 'BH')
COMPLETED = {'PASS', 'OBSTRUCTED'}
RANK4_LOGS = ('owens_rank4/results_u4_local_2026-09-07_snapshot.jsonl',
              'owens_rank4/results_u4_local_2026-09-08_snapshot.jsonl')


def paper_scan(consolidated):
    """Keep the reviewed and historical appendices independent of live searches."""
    return ('open23/priority_u23_final.json' if consolidated.get('manuscript') == 'v1.3'
            else PAPER_SCAN)


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def verify_manifest(results=RESULTS, manuscript="v1.3"):
    """Reject missing or altered deposited inputs rather than silently recounting."""
    results = Path(results)
    if manuscript not in {'v1.1', 'v1.2', 'v1.3', 'v1.3-review'}:
        raise ValueError(f'Unsupported manuscript version: {manuscript}')
    version = '1_3_review' if manuscript == 'v1.3-review' else ('1_3' if manuscript == 'v1.3' else '1_1')
    manifest = read_json(results / f'paper_v{version}_manifest.json')
    for name, expected in manifest['sha256'].items():
        if hashlib.sha256((results / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Deposited input differs from the {manifest["manuscript"]} manifest: {name}')
    return manifest


@dataclass(frozen=True)
class Bound:
    name: str
    side: str
    value: int
    tag: str
    source: str


def validate_interval(interval):
    if (not isinstance(interval, (list, tuple)) or len(interval) != 2
            or any(type(x) is not int for x in interval)
            or not 0 <= interval[0] <= interval[1]):
        raise ValueError(f'Invalid unknotting-number interval: {interval!r}')


def completed_rank4(records):
    """Deduplicate resumed jobs; unfinished records do not erase completed work."""
    completed = {}
    for record in records:
        name, verdict = record['name'], record['verdict']
        if verdict not in COMPLETED:
            continue
        if name in completed and completed[name] != verdict:
            raise ValueError(f'Conflicting completed rank-4 verdicts for {name}')
        completed[name] = verdict
    return completed


def greene_obstructed(record):
    """Read a completed surgery test, rejecting empty or partial candidate data.

    A pinned correction-term vector and a finite set of candidate vectors are
    different kinds of output. In the second case *every* conjugation-symmetric
    candidate must have been tested. Zero tested candidates is not an obstruction.
    These checks authenticate the shape of the deposited computation, not its
    mathematical hypotheses or the implementation of Greene's model.
    """
    order = record['det']
    if type(order) is not int or order <= 0 or order % 2 == 0:
        raise ValueError('Greene records require positive odd determinant')
    if 'd' in record:
        choices = {int(k): {Fraction(v)} for k, v in record['d'].items()}
        fits = record['u1_fits']
        if not isinstance(fits, list):
            raise ValueError('Missing completed surgery-test fits')
        obstructed = not fits
    else:
        pinned, ambiguous = record['pinned'], record['ambiguous']
        if pinned.keys() & ambiguous.keys():
            raise ValueError('A Greene class cannot be both pinned and ambiguous')
        choices = {int(k): {Fraction(v)} for k, v in pinned.items()}
        choices.update({int(k): set(map(Fraction, values))
                        for k, values in ambiguous.items()})
        # One choice for each pair t,-t determines a conjugation-symmetric vector.
        tested = record['candidate_vectors']
        admitting = record['admitting_u1']
        expected = math.prod(len(values) for k, values in choices.items()
                             if k <= (-k) % order)
        if (type(tested) is not int or tested <= 0 or tested != expected
                or type(admitting) is not int or not 0 <= admitting <= tested):
            raise ValueError('Incomplete Greene candidate-vector enumeration')
        obstructed = admitting == 0
    if (set(choices) != set(range(order)) or not all(choices.values())
            or any(values != choices[(-k) % order] for k, values in choices.items())):
        raise ValueError('Incomplete or inconsistent Greene Spin^c classes')
    return obstructed


def validate_sweep_record(row, entry, role):
    """Check the deposited premises and exhaustive candidate bookkeeping.

    This does not recompute the Floer complex. It prevents an incomplete page,
    skipped candidate, or failed worker from being reported as a lower bound.
    Inconclusive candidates remain admissible for this purpose.
    """
    name, order = entry['name'], entry['determinant']
    validate_interval(entry['range'])
    length = max(entry['range'][0], 1)
    expected_test = 'u1' if length == 1 else f'rank{length}'
    vector = ast.literal_eval(entry['reduced_mod2_vector_raw'])
    if (type(order) is not int or order < 3 or order % 2 == 0
            or length not in (1, 2, 3, 4) or entry['test'] != expected_test
            or not vector or any(len(e) != 4 or e[0] != 2 or e[1] < 0
                or any(type(x) is not int for x in e) for e in vector)
            or sum(e[1] for e in vector) != order or entry['khovanov_rank'] != order
            or type(entry['signature']) is not int
            or (length == 1 and abs(entry['signature']) > 2)
            or (length > 1 and abs(entry['signature']) != 2 * length)):
        raise ValueError(f'Invalid frozen Greene premise: {name}')
    if role == 'control':
        if (entry['range'] != [length, length]
                or entry.get('known_unknotting_number') != length):
            raise ValueError(f'Invalid Greene control: {name}')
    elif entry['range'][0] == entry['range'][1]:
        raise ValueError(f'A Greene target must have an open range: {name}')
    if row['name'] != name or row.get('role') != role:
        raise ValueError(f'Greene knot or role mismatch: {name}')
    verdict = row['verdict']
    if verdict not in {'OBSTRUCTED', 'PASS', 'UNDECIDED', 'NOT_APPLICABLE', 'ERROR', 'TIMEOUT'}:
        raise ValueError(f'Unknown Greene verdict: {name}: {verdict}')
    if verdict in {'ERROR', 'TIMEOUT', 'NOT_APPLICABLE'}:
        return
    if (row['det'] != order or row['sigma'] != entry['signature']
            or row['test'] != expected_test or row['range'] != entry['range']
            or row['excluded_length'] != length):
        raise ValueError(f'Greene output disagrees with the frozen premise: {name}')
    if 'd' in row:
        choices = {int(k): {Fraction(v)} for k, v in row['d'].items()}
    else:
        pinned, ambiguous = row['pinned'], row['ambiguous']
        if pinned.keys() & ambiguous.keys():
            raise ValueError(f'Overlapping Greene candidate classes: {name}')
        choices = {int(k): {Fraction(v)} for k, v in pinned.items()}
        choices.update({int(k): set(map(Fraction, v)) for k, v in ambiguous.items()})
    if (set(choices) != set(range(order)) or not all(choices.values())
            or any(v != choices[(-k) % order] for k, v in choices.items())):
        raise ValueError(f'Incomplete Greene candidate classes: {name}')
    expected = math.prod(len(v) for k, v in choices.items() if k <= (-k) % order)
    tests = row['candidate_verdicts']
    if (type(row['candidate_vectors']) is not int or row['candidate_vectors'] != expected
            or len(tests) != expected or type(row['admitting_vectors']) is not int):
        raise ValueError(f'Incomplete Greene candidate-vector enumeration: {name}')
    for test in tests:
        if test['test'] != expected_test or type(test['admits']) is not bool:
            raise ValueError(f'Mismatched Greene candidate test: {name}')
        detail = test['detail']
        if length == 1:
            if not isinstance(detail['fits'], list) or test['admits'] != bool(detail['fits']):
                raise ValueError(f'Inconsistent Greene surgery-test fits: {name}')
        elif test['admits'] != (detail['verdict'] != 'OBSTRUCTED'):
            raise ValueError(f'Inconsistent Greene definite-form verdict: {name}')
    admitting = sum(test['admits'] for test in tests)
    if (row['admitting_vectors'] != admitting
            or (verdict == 'OBSTRUCTED') != (admitting == 0)):
        raise ValueError(f'Inconsistent completed Greene verdict: {name}')
    if verdict == 'OBSTRUCTED':
        if role == 'control':
            raise ValueError(f'Greene sweep obstructs the control {name}')
        if row.get('lower_bound') != length + 1:
            raise ValueError(f'Wrong Greene lower bound: {name}')


def read_greene_sweep(results=RESULTS):
    """Read each frozen job once; unreported jobs remain explicitly uncompleted."""
    results = Path(results)
    frozen = read_json(results / 'greene/greene_targets_2026-09-08.json')
    if set(frozen['controls']) & set(frozen['targets']):
        raise ValueError('The Greene target and control sets overlap')
    records = {}
    with gzip.open(results / 'greene/sweep_2026-09-08.jsonl.gz', 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f'Malformed Greene sweep JSON at line {number}') from exc
            name = row['name']
            if name in records:
                raise ValueError(f'Duplicate frozen Greene sweep record: {name}')
            role = 'control' if name in frozen['controls'] else 'target'
            entry = frozen[role + 's'].get(name)
            if entry is None:
                raise ValueError(f'Greene sweep knot outside the frozen list: {name}')
            validate_sweep_record(row, entry, role)
            records[name] = row
    return frozen, records


def collect_bounds(results=RESULTS, include_sweep=True):
    """Translate each deposited theorem application to a bound with provenance.

    Tags select the primary explanation in the appendices; they do not count a
    knot twice when several independent methods establish the same inequality.
    """
    results = Path(results)
    bounds = []
    def data(name): return read_json(results / name)
    def add(name, value, tag, source, side='lo'):
        bounds.append(Bound(name, side, value, tag, source))

    source = 'lickorish_cw_obstructed_815_knotinfo2026.8.1.json'
    for name, _ in data(source): add(name, 2, 'L', source)
    source = 'torsion_lower_bounds_2026-08-24.json'
    for row in data(source)['hits']: add(row['name'], 2, 't', source)
    source = 'owens_u3_2026-08-24.json'
    traczyk = {'12a_634', '13a_2008', '13a_2745', '13a_3607'}
    for row in data(source):
        if row['verdict'] != 'OBSTRUCTED':
            raise ValueError(f'Unexpected verdict in {source}: {row}')
        name = row['name']
        add(name, 3, 'a' if name == '13a_2890' else 'OT' if name in traczyk else 'O2', source)
    source = 'generator_bound_2026-08-24.json'
    for row in data(source): add(row['name'], 3, 'g', source)
    source = 'owens_rank2/owens_verdicts_sigma4_alternating_2026-09-07.json'
    for row in data(source):
        if row['verdict'] == 'OBSTRUCTED': add(row['name'], 3, 'O2', source)
    source = 'owens_rank3/results_u3_2026-09-07.json'
    for part in data(source)['records'].values():
        for name, row in part.items():
            if row['verdict'] == 'OBSTRUCTED': add(name, 4, 'O3', source)

    source = 'owens_rank4/results_u4_partial_2026-09-07.json'
    partial = data(source)
    records = [r for key in ('targets', 'controls') for r in partial[key].values()]
    for log in RANK4_LOGS:
        for number, line in enumerate((results / log).read_text().splitlines(), 1):
            if not line.strip(): continue
            try: records.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f'Malformed JSON in {log}:{number}') from exc
    rank4 = completed_rank4(records)
    for name, verdict in rank4.items():
        if verdict == 'OBSTRUCTED': add(name, 5, 'O4', 'owens_rank4/ (frozen completed records)')

    source = 'montesinos/montesinos_u1_2026-09-07.json'
    mont = data(source)
    for key in ('validate', 'apply'):
        for name, row in mont[key].items():
            if row['verdict'] in {'OBSTRUCTED', 'OBSTRUCTED_HF'}: add(name, 2, 'M', source)
    # The L-space premise is an explicit deposited input, not a consequence of
    # the surgery-test output. These are tabulated F_2 ranks, not a fresh
    # Khovanov chain-complex computation performed by this reporting script.
    kh = data('bernhard_jablan/khovanov_inputs_2026-09-08.json')
    if kh['coefficient_field'] != 'F_2':
        raise ValueError('The L-space rank bound requires the deposited F_2 data')
    for name in ('12n_491', '13n_3370'):
        source = f'bernhard_jablan/greene_d_{name}_2026-09-08.json'
        row = data(source)
        if row.get('schema_version') != 2 or row['name'] != name:
            raise ValueError(f'Mismatched Greene knot label: {source}')
        premise = kh['knots'][name]
        if (premise['determinant'] != row['det']
                or premise['reduced_mod2_rank'] != row['det']
                or premise['l_space_rank_equality'] is not True
                or row['pd'] != premise['pd']
                or row['lspace_evidence']['vector'] != premise['reduced_mod2_vector']):
            raise ValueError(f'Missing L-space rank-equality premise: {name}')
        obstructed = greene_obstructed(row)
        tests = row['candidate_tests']
        if (type(row['candidate_vectors']) is not int or row['candidate_vectors'] <= 0
                or ('d' in row and row['candidate_vectors'] != 1)
                or len(tests) != row['candidate_vectors']
                or sum(bool(test['surgery_test']['fits']) for test in tests)
                    != row['admitting_u1']
                or (row['verdict'] == 'OBSTRUCTED') != obstructed
                or (obstructed and row.get('lower_bound') != 2)):
            raise ValueError(f'Inconsistent completed Greene verdict: {name}')
        if obstructed: add(name, 2, 'G', source)

    # The sweep applies the same two theorems to every knot of the frozen list whose reduced Khovanov
    # homology over F_2 has rank equal to its determinant, which is the L-space premise recorded there.
    # A record contributes only when it is a completed obstruction whose excluded length matches the
    # comparison range and the signature, so a timeout or an inconclusive comparison supplies nothing.
    if include_sweep:
        source = 'greene/sweep_2026-09-08.jsonl.gz'
        _, records = read_greene_sweep(results)
        for name, row in records.items():
            if row['verdict'] == 'OBSTRUCTED':
                add(name, row['lower_bound'], 'G', source)

    # The archived comparison has upper bound three for 13n_3370. Its sharper
    # upper bound is already in Brittenham--Hermiller, Theorem 1.3(a), and is
    # supplied separately rather than relabelled as a new diagram construction.
    source = 'bernhard_jablan/brittenham_hermiller_upper_bounds.json'
    published = data(source)
    if published['evidence_kind'] != 'published_construction':
        raise ValueError('Missing published upper-bound provenance')
    for name, row in published['records'].items():
        if row['upper_bound'] != row['child_upper_bound'] + 1:
            raise ValueError(f'Inconsistent published crossing-change bound: {name}')
        add(name, row['upper_bound'], 'BH', source, 'hi')
    source = 'crossing_changes/mccoy_alternating_u1_2026-09-08.json'
    for name, row in data(source).items():
        if row['verdict'] == 'u >= 2 (McCoy)': add(name, 2, 'K', source)
    source = 'cyclic_cover/cyclic_cover_bound_2026-09-07.json'
    for name, row in data(source).items():
        if row['bound'] >= 2: add(name, row['bound'], 'g' if row['n'] == 2 else 'c', source)

    source = 'summary.json'
    diagrams = data(source)
    if len(diagrams) != 8 or len({r['knot'] for r in diagrams}) != 8:
        raise ValueError('The paper requires eight distinct presentation certificates')
    for row in diagrams:
        # These flags are a record of the original verification, not a replacement
        # for verify_presentation_diagrams.py, which checks the actual diagrams.
        best = row['best']
        if (not best['verified'] or not best['unknot_after_change']
                or len(set(best['change_rows_0based'])) != 2
                or row['new_unknotting_number'] != 2):
            raise ValueError(f'Incomplete presentation certificate: {row["knot"]}')
        add(row['knot'], 2, 'w', source, 'hi')
    # The braid certificate is attributed to Brittenham--Hermiller in the paper.
    source = '13n_1587_u_le_2_braid_certificate.json'
    data(source)  # Require the actual certificate to be present and valid JSON.
    add('13n_1587', 2, 'b', source, 'hi')
    return bounds, rank4


def merge_bounds(snapshot, bounds):
    """Intersect valid intervals; a contradiction is an error, never a result."""
    for interval in snapshot.values(): validate_interval(interval)
    ends = {'lo': defaultdict(dict), 'hi': defaultdict(dict)}
    provenance = defaultdict(set)
    for bound in bounds:
        if bound.name not in snapshot:
            raise ValueError(f'Bound for knot absent from comparison snapshot: {bound.name}')
        if (bound.side not in ends or bound.tag not in PRIORITY
                or type(bound.value) is not int or bound.value < 0):
            raise ValueError(f'Invalid bound: {bound}')
        dest = ends[bound.side][bound.name]
        op = max if bound.side == 'lo' else min
        dest[bound.tag] = op(dest.get(bound.tag, bound.value), bound.value)
        provenance[bound.name].add(bound.source)
    table, changed = {}, {}
    for name, (a, b) in snapshot.items():
        low, high = ends['lo'][name], ends['hi'][name]
        lo, hi = max([a, *low.values()]), min([b, *high.values()])
        if lo > hi:
            raise ValueError(f'Contradictory bounds for {name}: [{lo},{hi}], input [{a},{b}]')
        table[name] = [lo, hi]
        if [lo, hi] != [a, b]:
            changed[name] = {'reference': [a, b], 'new': [lo, hi],
                'lower_tags': [t for t in PRIORITY if low.get(t) == lo] if lo > a else [],
                'upper_tags': [t for t in PRIORITY if high.get(t) == hi] if hi < b else [],
                'sources': sorted(provenance[name])}
    return table, changed


def primary_tag(record):
    return (record['lower_tags'] or record['upper_tags'])[0]


def summarize(table, changed):
    exact = [r for r in changed.values() if r['new'][0] == r['new'][1]]
    improved = [r for r in changed.values() if r['new'][0] != r['new'][1]]
    return {'knots': len(table), 'changed': len(changed), 'exact': len(exact),
        'exact_by_u': dict(sorted(Counter(str(r['new'][0]) for r in exact).items())),
        'improved': len(improved),
        'exact_lower': sum(bool(r['lower_tags']) for r in exact),
        'exact_upper': sum(bool(r['upper_tags']) and not r['lower_tags'] for r in exact),
        'improved_lower': sum(bool(r['lower_tags']) for r in improved),
        'improved_upper': sum(bool(r['upper_tags']) and not r['lower_tags'] for r in improved),
        'exact_by_primary_method': dict(sorted(Counter(map(primary_tag, exact)).items())),
        'improved_by_primary_method': dict(sorted(Counter(map(primary_tag, improved)).items()))}


def reconstruct(results=RESULTS, manuscript="v1.3"):
    """Validate all inputs and counts in memory before the caller writes outputs."""
    results = Path(results)
    manifest = verify_manifest(results, manuscript)
    bounds, rank4 = collect_bounds(results, include_sweep=manuscript.startswith("v1.3"))
    table, changed = merge_bounds(read_json(results / 'paper_v1_1_snapshot.json'), bounds)
    counts = summarize(table, changed)
    for key, expected in manifest['expected'].items():
        if key == 'dichotomy':
            scan = read_json(results / paper_scan({'manuscript': manuscript}))
            actual = sum(heuristic == 0 for _, _, _, heuristic in scan['zero_candidate_knots'])
        elif key in counts:
            actual = counts[key]
        else:
            raise ValueError(f'Unknown expected manuscript count: {key}')
        if actual != expected:
            raise ValueError(f'Manuscript count mismatch: {key}: {actual} != {expected}')
    return {'manuscript': manuscript, 'reference': 'results/paper_v1_1_snapshot.json', 'changed': changed,
            'rank4_status': rank4, 'counts': counts}, table



def manuscript_counts(consolidated, table, results=RESULTS):
    """Read the manuscript's numerical premises without generating LaTeX files."""
    results = Path(results)
    changed = consolidated['changed']
    totals = consolidated['counts']
    exact = [r for r in changed.values() if r['new'][0] == r['new'][1]]
    by_value_method = Counter((r['new'][0], primary_tag(r)) for r in exact)
    improved_methods = totals['improved_by_primary_method']
    scan = read_json(results/paper_scan(consolidated))
    moved = set(scan.get('moved_to_tierB_after_resolution', []))
    zero_names = {name for name, _, _, _ in scan['zero_candidate_knots']}
    candidate_names = set(scan['knots_with_candidates'])
    if (moved & zero_names or not moved <= candidate_names
            or len(candidate_names) != len(scan['knots_with_candidates'])
            or len(zero_names) != len(scan['zero_candidate_knots'])):
        raise ValueError('Candidate scan categories overlap or contain duplicate knots')
    cyclic = read_json(results/'cyclic_cover/cyclic_cover_bound_2026-09-07.json')
    cyclic_degrees = Counter(cyclic[name]['n'] for name, r in changed.items() if primary_tag(r) == 'c')
    snapshot = read_json(results/'paper_v1_1_snapshot.json')
    sig4 = read_json(results/'owens_rank2/owens_verdicts_sigma4_alternating_2026-09-07.json')
    with gzip.open(results/'crossing_changes/dataset_v2.json.gz', 'rt') as stream:
        data_knots = {r['knot'] for r in json.load(stream)}
    sweep = Counter()
    frozen = {'targets': {}, 'controls': {}}
    if consolidated.get('manuscript', 'v1.3').startswith('v1.3'):
        frozen, records = read_greene_sweep(results)
        settled = {'OBSTRUCTED', 'PASS', 'UNDECIDED', 'NOT_APPLICABLE'}
        for row in records.values():
            role = row['role']
            sweep[role, 'settled' if row['verdict'] in settled else 'unsettled'] += 1
            sweep[role, 'verdict_' + row['verdict']] += 1
            if row['verdict'] == 'OBSTRUCTED': sweep[role, row['test']] += 1
    greene_exact = sum(count for (value, tag), count in by_value_method.items() if tag == 'G')
    # Sigma_2(A # B) = Sigma_2(A) # Sigma_2(B), so the generator bound of the double branched cover is
    # additive.  Where it is sharp for both summands, the unknotting number of the connected sum is
    # forced; these are the knots for which that happens.
    import ast as _ast, database_knotinfo as _dk
    kinfo = {r['name']: r for r in _dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    def _g2(name):
        raw = str(kinfo[name].get('torsion_numbers', '')).strip()
        if not raw:
            return None
        try:
            factors = dict((int(m), list(f)) for m, f in _ast.literal_eval(raw))
        except (SyntaxError, ValueError, TypeError):
            return None
        return None if 2 not in factors else sum(1 for x in factors[2] if x != 1)
    sharp = [n for n, v in table.items() if v[0] == v[1] >= 1 and n in kinfo and _g2(n) == v[0]]
    sharp_u1 = sum(1 for n in sharp if table[n][0] == 1)
    u1_scan, u1_open, u1_with, u1_open_with = {}, {}, 0, 0
    if consolidated.get('manuscript') == 'v1.3':
        u1_scan = read_json(results/'crossing_changes/u1_minimal_diagram_scan_2026-09-08.json')
        u1_open = read_json(results/'crossing_changes/u1_open_diagram_scan_2026-09-08.json')
        u1_with = sum(bool(r['unknotting_crossings']) for r in u1_scan.values())
        u1_open_with = sum(bool(r['unknotting_crossings']) for r in u1_open.values())
        if any(r['undecided'] for r in u1_scan.values()) or any(r['undecided'] for r in u1_open.values()):
            raise ValueError('An undecided crossing change leaves the minimal-diagram scan incomplete')
    counts = {
        'nChildren': len(scan['children']),
        'nCyclicNew': sum(cyclic_degrees.values()),
        'nCyclicThree': cyclic_degrees[3], 'nCyclicFour': cyclic_degrees[4],
        'nCyclicFive': cyclic_degrees[5], 'nDataKnots': len(data_knots),
        'nDichotomy': sum(heuristic == 0 for _, _, _, heuristic in scan['zero_candidate_knots']),
        # The seven Jones-only composite identifications moved out of the
        # candidate list remain heuristic; they never join the theorem list.
        # Thus the latest cohort is 959 certified + (39 + 7) heuristic + 22.
        'nDichotomyHeuristic': sum(heuristic != 0 for _, _, _, heuristic in scan['zero_candidate_knots']) + len(moved),
        'nExact': totals['exact'], 'nExactLower': totals['exact_lower'],
        'nImproved': totals['improved'], 'nImprovedL': improved_methods.get('L', 0),
        'nImprovedC': improved_methods.get('c', 0), 'nMcCoyKnots': improved_methods.get('K', 0),
        'nOpenTwoThreeAlt': len(scan['zero_candidate_knots']) + len(scan['knots_with_candidates']),
        'nRefKnots': len(table),
        'nSigFourImproved': sum(snapshot[r['name']] == [2, 4] and r['verdict'] == 'OBSTRUCTED' for r in sig4),
        'nSigFourObstructed': sum(snapshot[r['name']] == [2, 3] and r['verdict'] == 'OBSTRUCTED' for r in sig4),
        'nSigFourOpen': sum(snapshot[r['name']] == [2, 3] for r in sig4),
        'nUtwo': totals['exact_by_u']['2'], 'nUtwoL': by_value_method[2, 'L'],
        'nUtwoM': by_value_method[2, 'M'], 'nUtwoC': by_value_method[2, 'c'],
        'nUtwoG': by_value_method[2, 'G'],
        'nUthree': totals['exact_by_u']['3'], 'nUthreeC': by_value_method[3, 'c'],
        'nUthreeOwens': sum(by_value_method[3, tag] for tag in ('O2', 'a', 'OT', 'g')),
        'nUfour': totals['exact_by_u']['4'], 'nUfive': totals['exact_by_u']['5'],
        'nWithCandidates': len(candidate_names - moved),
        'nSweepTargets': len(frozen['targets']), 'nSweepControls': len(frozen['controls']),
        'nSweepTargetsSettled': sweep['target', 'settled'],
        'nSweepControlsSettled': sweep['control', 'settled'],
        'nSweepControlsPass': sweep['control', 'verdict_PASS'],
        'nSweepControlsUndecided': sweep['control', 'verdict_UNDECIDED'],
        'nSweepControlsUnfinished': len(frozen['controls']) - sweep['control', 'settled'],
        'nSweepObstructed': sum(count for (role, key), count in sweep.items()
                                if role == 'target' and key.startswith(('u1', 'rank'))),
        'nSweepUone': sweep['target', 'u1'], 'nSweepRankTwo': sweep['target', 'rank2'],
        'nSweepRankThree': sweep['target', 'rank3'], 'nSweepRankFour': sweep['target', 'rank4'],
        'nSweepExact': greene_exact, 'nSweepImproved': improved_methods.get('G', 0),
        'nSweepNewExact': sum(1 for r in exact if primary_tag(r) == 'G'
                            and 'greene/sweep_2026-09-08.jsonl.gz' in r['sources']),
        'nSharpGenerator': len(sharp), 'nSharpGeneratorUone': sharp_u1,
        'nSharpGeneratorAdded': sum(1 for n in sharp if snapshot.get(n, [0, 9])[0] != snapshot.get(n, [0, 9])[1]),
        'nSharpGeneratorSums': len(sharp) * (len(sharp) + 1) // 2,
        'nSharpGeneratorSumsNontrivial': (len(sharp) * (len(sharp) + 1) // 2
                                          - sharp_u1 * (sharp_u1 + 1) // 2),
        'nUoneVerified': len(u1_scan), 'nUoneWithCrossing': u1_with,
        'nUoneOpenScanned': len(u1_open), 'nUoneOpenWithCrossing': u1_open_with,
        'nSweepUoneTargets': sum(1 for e in frozen['targets'].values() if e['test'] == 'u1'),
        'nSweepRankTargets': sum(1 for e in frozen['targets'].values() if e['test'] != 'u1'),
    }

    if consolidated.get('manuscript') != 'v1.3':
        counts = {k: v for k, v in counts.items()
                  if not k.startswith('nUone') and k != 'nDichotomyHeuristic'}
    if not consolidated.get('manuscript', 'v1.3').startswith('v1.3'):
        counts = {k: v for k, v in counts.items() if not k.startswith('nSweep')}
    return counts


def knot_sort_key(name):
    """Use table order: crossing number, alternating type, then knot index."""
    match = re.fullmatch(r'(\d+)([an]?)_(\d+)', name)
    if not match:
        raise ValueError(f'Unrecognized knot label: {name}')
    return int(match[1]), match[2], int(match[3])


def tex_knot(name):
    knot_sort_key(name)  # Do not let an arbitrary input string become TeX code.
    left, right = name.split('_')
    return f'{left}{right}' if left[-1] in 'an' else f'${left}_{{{right}}}$'


def appendix_tag(record):
    """Preserve the manuscript's method notation, with distinct G and g.

    For lower bound three, if both the elementary generator bound and a
    correction-term obstruction apply, the appendix highlights the generator
    bound. At lower bound two the linking-pairing method keeps its unmarked
    notation even when a generator bound also applies. The lower-bound summary
    uses the fixed primary-method priority to avoid double counting throughout.
    """
    tag = ('g' if record['new'][0] == 3 and 'g' in record['lower_tags']
           else primary_tag(record))
    if record['new'][0] == record['new'][1]:
        symbols = {'L': '', 't': 't', 'a': 'a', 'OT': 't', 'O2': '',
                   'g': 'g', 'c': 'c', 'M': 'm', 'G': 'G', 'O3': '',
                   'O4': '', 'w': 'w'}
    else:
        symbols = {'L': '', 't': 't', 'g': 'g', 'c': 'c', 'M': 'm',
                   'O2': 'o', 'O3': 'h', 'K': 'k', 'b': 'b', 'G': 'G'}
    if tag not in symbols:
        raise ValueError(f'No appendix notation for method {tag}')
    return symbols[tag]


def tex_entry(name, record):
    tag = appendix_tag(record)
    return tex_knot(name) + (rf'$^{{\mathrm{{{tag}}}}}$' if tag else '')


def column_rows(cells, columns):
    """Read down columns and keep the last three rows with the final legend."""
    height = math.ceil(len(cells)/columns)
    rows = []
    for row in range(height):
        entries = [cells[row + column*height] if row + column*height < len(cells)
                   else '' for column in range(columns)]
        rows.append(' & '.join(entries) + (r' \\*' if row >= height-3 else r' \\'))
    return '\n'.join(rows) + '\n'


def comparison_report(consolidated, table):
    """Recompute release comparisons and require their deposited provenance."""
    from audit_knotinfo_releases import audit
    summary, _ = audit(consolidated, table)
    filename = ('baseline_comparison_v1_3_review.json'
                if consolidated.get('manuscript') == 'v1.3-review' else 'baseline_comparison.json')
    deposited = read_json(RESULTS / 'comparison' / filename)
    if consolidated.get('manuscript', 'v1.3').startswith('v1.3') and summary != deposited:
        raise ValueError('Release comparison differs from the deposited audit; '
                         'review the input provenance before updating it')
    return summary


def comparison_counts(summary):
    keys = {
        'nBaselineOverrides': 'modified_entries',
        'nBaselineLowerRestored': 'restored_lower_entries',
        'nAprilDifferences': 'baseline_differences_2026_4_1',
        'nJuneDifferences': 'baseline_differences_2026_6_1',
        'nAugustDifferences': 'baseline_differences_2026_8_1',
        'nSeptemberDifferences': 'baseline_differences_2026_9_1',
        'nSeptemberExactEntries': 'september_exact_entries',
        'nSeptemberExactAgreements': 'september_exact_agreements',
        'nSeptemberExactDisagreements': 'september_exact_disagreements',
        'nSeptemberUnresolvedExact': 'september_unresolved_exact',
        'nReleaseConflicts': 'september_disjoint_intervals',
        'nRetainedLargerUpperBounds': 'retained_larger_upper_bounds',
        'nStrictlyWeakerRanges': 'september_strictly_weaker_intervals',
        'nIncomparableRanges': 'september_incomparable_intervals',
        'nMcCoyVerified': 'mccoy_total',
        'nMcCoyUnknotting': 'mccoy_u1',
        'nMcCoyExcluded': 'mccoy_excluded',
        'nMcCoyDatasetCovered': 'mccoy_dataset_overlap',
        'nMcCoyDatasetMissing': 'mccoy_dataset_missing',
        'nGPExact': 'gebel_prangley_exact_agreements',
        'nGPImproved': 'gebel_prangley_range_agreements',
    }
    return {macro: summary['counts'][key] for macro, key in keys.items()}


def tex_count_macros(counts):
    return '\n'.join(rf'\newcommand{{\{name}}}{{{value:,}}}'.replace(',', '{,}')
                     for name, value in sorted(counts.items())) + '\n'


def manuscript_tex(consolidated, table, results=RESULTS):
    """Return generated numerical inputs; all layout and prose stay in LaTeX.

    These are regenerated from authenticated records, never extracted from an
    existing manuscript or a historical generated table. Thus the Greene cases
    cannot disappear from a stale appendix while remaining in the headline count.
    """
    changed = consolidated['changed']
    prefix = '% Generated by consolidate_results.py; do not edit by hand.\n'
    counts = manuscript_counts(consolidated, table, results)
    outputs = {'counts.tex': prefix + tex_count_macros(counts)}
    for u in (5, 4, 3, 2):
        names = sorted((name for name, record in changed.items()
                        if record['new'] == [u, u]), key=knot_sort_key)
        outputs[f'rows_u{u}.tex'] = prefix + column_rows(
            [tex_entry(name, changed[name]) for name in names], 7)
    names = sorted((name for name, record in changed.items()
                    if record['new'][0] != record['new'][1]), key=knot_sort_key)
    cells = []
    for name in names:
        record = changed[name]
        before, after = (f'$[{lo},{hi}]$' for lo, hi in (record['reference'], record['new']))
        cells.append(f'{tex_entry(name, record)} & {before} & {after}')
    outputs['rows_improvements.tex'] = prefix + column_rows(cells, 2)
    scan = read_json(Path(results)/paper_scan(consolidated))
    names = sorted((name for name, _, _, heuristic in scan['zero_candidate_knots']
                    if heuristic == 0), key=knot_sort_key)
    outputs['rows_dichotomy.tex'] = prefix + column_rows(list(map(tex_knot, names)), 6)
    totals = consolidated['counts']
    methods = (
        ('Linking pairing', ('L',)),
        ('Homological torsion', ('t',)),
        ("Correction terms, Traczyk's criterion, and homology", ('O2', 'OT', 'a', 'g')),
        ('Obstruction to three crossing changes', ('O3',)),
        ('Obstruction to four crossing changes (completed cases)', ('O4',)),
        ('Higher cyclic covers', ('c',)),
        ('Branched covers of Montesinos knots', ('M',)),
        ("Greene's model for branched covers", ('G',)),
        ("McCoy's criterion for alternating diagrams", ('K',)),
    )
    method_rows = []
    for title, tags in methods:
        exact = sum(totals['exact_by_primary_method'].get(tag, 0) for tag in tags)
        improved = sum(totals['improved_by_primary_method'].get(tag, 0) for tag in tags)
        method_rows.append(f'{title} & {exact} & {improved}' + r' \\')
    outputs['rows_lower_summary.tex'] = prefix + '\n'.join(method_rows) + '\n'
    comparison = comparison_report(consolidated, table)
    outputs['comparison_counts.tex'] = prefix + tex_count_macros(comparison_counts(comparison))
    conflicts = []
    for row in sorted(comparison['disputes'], key=lambda row: knot_sort_key(row['knot'])):
        lo, hi = row['paper']
        paper_range = f'${lo}$' if lo == hi else f'$[{lo},{hi}]$'
        lo, hi = row['september']
        current_range = f'${lo}$' if lo == hi else f'$[{lo},{hi}]$'
        conflicts.append(f"{tex_knot(row['knot'])} & {paper_range} & {current_range}" + r' \\')
    outputs['rows_release_conflicts.tex'] = prefix + '\n'.join(conflicts) + '\n'
    return outputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'generated')
    parser.add_argument('--manuscript', choices=('v1.2', 'v1.3', 'v1.3-review'), default='v1.3')
    args = parser.parse_args()
    consolidated, table = reconstruct(manuscript=args.manuscript)
    tex = manuscript_tex(consolidated, table)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, value in [('consolidated.json', consolidated), ('u_table.json', table),
                        ('counts.json', consolidated['counts'])]:
        (args.outdir / name).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    paper_dir = args.outdir / ('paper_' + args.manuscript.replace('.', '_').replace('-', '_'))
    paper_dir.mkdir(parents=True, exist_ok=True)
    for name, source in tex.items():
        (paper_dir / name).write_text(source)
    print(json.dumps(consolidated['counts'], indent=2))


if __name__ == '__main__':
    main()
