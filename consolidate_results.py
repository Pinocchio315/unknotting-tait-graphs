#!/usr/bin/env python3
"""Reconstruct the v1.1 result tables from explicitly frozen, deposited inputs.

This is bookkeeping, not a rerun of the mathematical obstructions.  A completed
OBSTRUCTED verdict supplies a lower bound; PASS only means that this obstruction
is silent.  In particular, PASS, a missing record, or a timeout never proves an
upper bound.  Upper bounds come from separately verifiable diagram certificates.

Run from any directory: python consolidate_results.py --outdir generated
Historical results are read-only. No installed KnotInfo table or live worker log
is consulted, and conflicting bounds stop the run before any output is written.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
PRIORITY = ('L', 't', 'a', 'OT', 'O2', 'g', 'c', 'M', 'K', 'O3', 'O4', 'w', 'b')
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', type=Path, default=ROOT / 'generated')
    args = parser.parse_args()
    consolidated, table = reconstruct()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, value in [('consolidated.json', consolidated), ('u_table.json', table),
                        ('counts.json', consolidated['counts'])]:
        (args.outdir / name).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')
    print(json.dumps(consolidated['counts'], indent=2))


if __name__ == '__main__':
    main()
