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
PRIORITY = ('L', 't', 'a', 'OT', 'O2', 'g', 'c', 'M', 'G', 'K', 'O3', 'O4', 'w', 'b', 'BH')
COMPLETED = {'PASS', 'OBSTRUCTED'}
RANK4_LOGS = ('owens_rank4/results_u4_local_2026-09-07_snapshot.jsonl',
              'owens_rank4/results_u4_local_2026-09-08_snapshot.jsonl')


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def verify_manifest(results=RESULTS):
    """Reject missing or altered deposited inputs rather than silently recounting."""
    results = Path(results)
    manifest = read_json(results / 'paper_v1_1_manifest.json')
    for name, expected in manifest['sha256'].items():
        if hashlib.sha256((results / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Deposited input differs from the v1.1 manifest: {name}')
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


def collect_bounds(results=RESULTS):
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


def reconstruct(results=RESULTS):
    """Validate all inputs and counts in memory before the caller writes outputs."""
    results = Path(results)
    manifest = verify_manifest(results)
    bounds, rank4 = collect_bounds(results)
    table, changed = merge_bounds(read_json(results / 'paper_v1_1_snapshot.json'), bounds)
    counts = summarize(table, changed)
    for key, expected in manifest['expected'].items():
        if key in counts and counts[key] != expected:
            raise ValueError(f'Manuscript count mismatch: {key}: {counts[key]} != {expected}')
    return {'reference': 'results/paper_v1_1_snapshot.json', 'changed': changed,
            'rank4_status': rank4, 'counts': counts}, table



def manuscript_counts(consolidated, table, results=RESULTS):
    """Read the manuscript's numerical premises without generating LaTeX files."""
    results = Path(results)
    changed = consolidated['changed']
    totals = consolidated['counts']
    exact = [r for r in changed.values() if r['new'][0] == r['new'][1]]
    by_value_method = Counter((r['new'][0], primary_tag(r)) for r in exact)
    improved_methods = totals['improved_by_primary_method']
    scan = read_json(results/'open23/priority_u23_final.json')
    cyclic = read_json(results/'cyclic_cover/cyclic_cover_bound_2026-09-07.json')
    cyclic_degrees = Counter(cyclic[name]['n'] for name, r in changed.items() if primary_tag(r) == 'c')
    snapshot = read_json(results/'paper_v1_1_snapshot.json')
    sig4 = read_json(results/'owens_rank2/owens_verdicts_sigma4_alternating_2026-09-07.json')
    with gzip.open(results/'crossing_changes/dataset_v2.json.gz', 'rt') as stream:
        data_knots = {r['knot'] for r in json.load(stream)}
    return {
        'nChildren': len(scan['children']),
        'nCyclicNew': sum(cyclic_degrees.values()),
        'nCyclicThree': cyclic_degrees[3], 'nCyclicFour': cyclic_degrees[4],
        'nCyclicFive': cyclic_degrees[5], 'nDataKnots': len(data_knots),
        'nDichotomy': sum(heuristic == 0 for _, _, _, heuristic in scan['zero_candidate_knots']),
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
        'nWithCandidates': len(scan['knots_with_candidates']) - len(scan.get('moved_to_tierB_after_resolution', [])),
    }


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
                   'O2': 'o', 'O3': 'h', 'K': 'k', 'b': 'b'}
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
    deposited = read_json(RESULTS / 'comparison/baseline_comparison.json')
    if summary != deposited:
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
    scan = read_json(Path(results)/'open23/priority_u23_final.json')
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
    args = parser.parse_args()
    consolidated, table = reconstruct()
    tex = manuscript_tex(consolidated, table)
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, value in [('consolidated.json', consolidated), ('u_table.json', table),
                        ('counts.json', consolidated['counts'])]:
        (args.outdir / name).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    paper_dir = args.outdir / 'paper_v1_2'
    paper_dir.mkdir(parents=True, exist_ok=True)
    for name, source in tex.items():
        (paper_dir / name).write_text(source)
    print(json.dumps(consolidated['counts'], indent=2))


if __name__ == '__main__':
    main()
