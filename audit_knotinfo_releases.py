#!/usr/bin/env python3
"""Compare the historical research baseline with public KnotInfo package snapshots.

The frozen range extracts are *tabulated inputs*, not results proved here.  The
paper's obstruction outputs are reconstructed separately, and incompatible
intervals are reported explicitly rather than merged into a false interval.

This audit performs no knot search and no live database access.  To reproduce the
public input extraction, download the four SHA-256-identified wheels named in
results/comparison/releases/manifest.json and run --freeze-from DIRECTORY.  The
installed database_knotinfo package is never imported, installed, or modified.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import zipfile

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
OUTPUT = RESULTS / 'comparison'
VERSIONS = ('2026.4.1', '2026.6.1', '2026.8.1', '2026.9.1')
CSV_MEMBER = 'database_knotinfo/csv_data/knotinfo_data_complete.csv'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def interval(raw):
    values = [int(x) for x in re.findall(r'\d+', raw)]
    if not values or len(values) > 2:
        raise ValueError(f'Unrecognized tabulated unknotting range: {raw!r}')
    answer = [values[0], values[-1]]
    if answer[0] > answer[1]:
        raise ValueError(f'Reversed range: {answer}')
    return answer


def json_text(value):
    return json.dumps(value, indent=2, sort_keys=True) + '\n'


def freeze_releases(directory):
    """Extract only the named public fields, preserving source archive hashes.

    Each wheel is verified against its PyPI SHA-256 from the adjacent metadata
    file.  Gzip timestamps are fixed so repeated extraction gives identical bytes.
    The expected knot names come from the archived study population, not from an
    assumed release-dependent ordering or from the installed package.
    """
    directory = Path(directory)
    names = set(read_json(RESULTS / 'paper_v1_1_snapshot.json'))
    dest = OUTPUT / 'releases'
    dest.mkdir(parents=True, exist_ok=True)
    release_manifest = {'schema_version': 1, 'evidence_kind': 'tabulated_input',
                        'package': 'database_knotinfo', 'releases': {}}
    csv.field_size_limit(sys.maxsize)
    for version in VERSIONS:
        metadata_path = directory / f'{version}.metadata.json'
        if metadata_path.exists():
            metadata = read_json(metadata_path)
        else:
            saved = read_json(dest / 'manifest.json')['releases'][version]
            metadata = {key: saved[key] for key in
                        ('version', 'filename', 'url', 'upload_time_iso_8601', 'sha256')}
        archive_path = directory / metadata['filename']
        if digest(archive_path) != metadata['sha256']:
            raise ValueError(f'Public wheel hash mismatch: {version}')
        with zipfile.ZipFile(archive_path) as archive:
            raw = archive.read(CSV_MEMBER)
        reader = csv.DictReader(io.StringIO(raw.decode('utf-8')), delimiter='|')
        rows = {}
        for row in reader:
            name = row['name'].strip()
            if name not in names:
                continue
            rows[name] = {'range': interval(row['unknotting_number']),
                          'raw_unknotting_number': row['unknotting_number'],
                          'alternating': row['alternating'].strip() == 'Y',
                          'crossing_number': int(row['crossing_number'])}
        if set(rows) != names:
            raise ValueError(f'Release population mismatch: {version}')
        source = {**metadata, 'csv_member': CSV_MEMBER,
                  'csv_sha256': hashlib.sha256(raw).hexdigest(),
                  'metadata_url': f'https://pypi.org/pypi/database-knotinfo/{version}/json'}
        document = {'schema_version': 1, 'evidence_kind': 'tabulated_input',
                    'source': source, 'rows': rows}
        filename = f'knotinfo_{version}.json.gz'
        (dest / filename).write_bytes(gzip.compress(json_text(document).encode(), mtime=0))
        release_manifest['releases'][version] = {
            **source, 'extracted_file': filename, 'extracted_sha256': digest(dest / filename)}
    (dest / 'manifest.json').write_text(json_text(release_manifest))


def load_releases():
    manifest = read_json(OUTPUT / 'releases/manifest.json')
    releases = {}
    if set(manifest['releases']) != set(VERSIONS):
        raise ValueError('Expected all four comparison releases')
    for version, source in manifest['releases'].items():
        path = OUTPUT / 'releases' / source['extracted_file']
        if digest(path) != source['extracted_sha256']:
            raise ValueError(f'Frozen release extract changed: {version}')
        with gzip.open(path, 'rt') as stream:
            document = json.load(stream)
        if document['source']['sha256'] != source['sha256']:
            raise ValueError(f'Archive provenance mismatch: {version}')
        releases[version] = document['rows']
    return releases


def relation(paper, release):
    """Compare intervals without calling incompatible values 'improvements'."""
    if paper == release:
        return 'equal'
    if paper[0] > release[1] or paper[1] < release[0]:
        return 'disjoint'
    if paper[0] >= release[0] and paper[1] <= release[1]:
        return 'paper_stronger'
    if paper[0] <= release[0] and paper[1] >= release[1]:
        return 'paper_weaker'
    return 'overlapping_incomparable'


def audit(consolidated=None, table=None):
    """Recompute the comparisons from authenticated frozen inputs, in memory.

    The optional arguments let the reporting pipeline pass its already verified
    reconstruction.  Without them this function obtains exactly that same result
    from consolidate_results.reconstruct().  It never treats a release value as
    new mathematical evidence for a deposited obstruction.
    """
    if consolidated is None or table is None:
        from consolidate_results import reconstruct
        consolidated, table = reconstruct()
    baseline = read_json(RESULTS / 'paper_v1_1_snapshot.json')
    releases = load_releases()
    if any(set(rows) != set(baseline) for rows in releases.values()):
        raise ValueError('Comparison populations differ')
    april, june, august, september = [releases[v] for v in VERSIONS]
    rollback_path = RESULTS / 'dkt/knotinfo_rollback_list_2026-09-07.csv'
    with rollback_path.open(newline='') as stream:
        rollback = list(csv.DictReader(stream))
    prior_column = 'KnotInfo value before adoption (DKT workbook column)'
    adopted_column = 'KnotInfo 2026.8.1 value'
    overrides = {k for k in baseline if baseline[k] != august[k]['range']}
    if overrides != {r['knot'] for r in rollback}:
        raise ValueError('Baseline overrides do not equal the deposited 104-entry list')
    for row in rollback:
        name = row['knot']
        if not (baseline[name] == interval(row[prior_column]) == april[name]['range']):
            raise ValueError(f'Restored range differs from April release: {name}')
        if august[name]['range'] != interval(row[adopted_column]):
            raise ValueError(f'Adopted range differs from August release: {name}')
    changed = consolidated['changed']
    exact = {k: r for k, r in changed.items() if r['new'][0] == r['new'][1]}
    exact_partition = Counter()
    for name, row in exact.items():
        current = september[name]['range']
        category = ('september_unresolved' if current[0] < current[1] else
                    'same_exact_value' if row['new'] == current else 'different_exact_value')
        exact_partition[category] += 1
    # Literature overlap is established by matching the actual theorem lists,
    # not by attributing every database update to a nearby publication date.
    gp_path = OUTPUT / 'gebel_prangley_2026.json'
    gp_check_path = OUTPUT / 'gebel_prangley_source_verification.json'
    gp_check = read_json(gp_check_path)
    if digest(gp_path) != gp_check['extracted_file_sha256']:
        raise ValueError('Gebel--Prangley list changed since primary-source verification')
    gp = read_json(gp_path)
    for key, count in gp_check['theorem_sizes'].items():
        if len(gp[key]) != count or len(set(gp[key])) != count:
            raise ValueError(f'Invalid published theorem list: {key}')
    gp_exact = set(gp['theorem1_u3']) | set(gp['theorem2_u3'])
    exact_agreements = {k for k in exact if exact[k]['new'] == september[k]['range']}
    if gp_exact != exact_agreements or any(table[k] != [3, 3] for k in gp_exact):
        raise ValueError('Published exact list differs from the September agreement set')
    if any(table[k] != [3, 4] for k in gp['theorem3_34']):
        raise ValueError('A Gebel--Prangley Theorem 3 range disagrees')
    disputes, comparison_rows = [], []
    for name in sorted(baseline):
        current = september[name]['range']
        kind = relation(table[name], current)
        if baseline[name] != current or table[name] != current:
            comparison_rows.append({'knot': name, 'baseline': baseline[name],
                                    'paper': table[name], 'september': current, 'relation': kind})
        if kind == 'disjoint':
            row = changed[name]
            disputes.append({'knot': name, 'baseline': baseline[name],
                             'paper': table[name], 'september': current,
                             'computed_lower_bound': table[name][0],
                             'lower_tags': row['lower_tags'], 'sources': row['sources']})
    retained_upper = {k for k in table if table[k][1] > september[k]['range'][1]}
    mccoy_path = RESULTS / 'crossing_changes/mccoy_alternating_all_2026-09-08.json'
    mccoy = read_json(mccoy_path)
    cohort = {k for k, r in april.items() if r['alternating'] and r['range'][0] == 1}
    if set(mccoy) != cohort:
        raise ValueError('McCoy data are not exactly the April lower-bound-one cohort')
    positive, negative = set(), set()
    for name, row in mccoy.items():
        if row['reference_range'] != april[name]['range'] or row['undecided']:
            raise ValueError(f'Incomplete or mismatched McCoy record: {name}')
        if row['unknotting_crossings'] and row['verdict'] == 'u = 1':
            positive.add(name)
        elif not row['unknotting_crossings'] and row['verdict'] == 'u >= 2 (McCoy)':
            negative.add(name)
        else:
            raise ValueError(f'Inconsistent McCoy verdict: {name}')
    if {k for k in cohort if june[k]['range'] == [1, 1]} != positive:
        raise ValueError('Positive controls differ from the June exact-u=1 set')
    if {k for k in cohort if june[k]['range'][0] == 2} != negative:
        raise ValueError('Excluded set differs from June lower-bound increases')
    historical_mccoy = read_json(RESULTS / 'crossing_changes/mccoy_alternating_u1_2026-09-08.json')
    restored_lower = {k for k in overrides if baseline[k][0] != august[k]['range'][0]}
    if set(historical_mccoy) != restored_lower or not restored_lower <= negative:
        raise ValueError('Historical McCoy cases differ from the restored lower-bound entries')
    dataset_path = RESULTS / 'crossing_changes/dataset_v2.json.gz'
    with gzip.open(dataset_path, 'rt') as stream:
        dataset = json.load(stream)
    by_parent = defaultdict(list)
    for row in dataset:
        by_parent[row['knot']].append(row)
    overlap = cohort & set(by_parent)
    for name in overlap:
        rows = by_parent[name]
        if len(rows) != april[name]['crossing_number']:
            raise ValueError(f'Incomplete crossing dataset for {name}')
        unknot_crossings = sorted(r['crossing'] for r in rows if r['result'] == '0_1')
        if unknot_crossings != sorted(mccoy[name]['unknotting_crossings']):
            raise ValueError(f'Dataset and McCoy unknotting crossings differ for {name}')
    counts = {
        'reference_knots': len(baseline), 'modified_entries': len(overrides),
        'restored_lower_entries': sum(baseline[k][0] != august[k]['range'][0] for k in overrides),
        'paper_exact': len(exact), 'september_exact_entries': sum(exact_partition[k] for k in
            ('same_exact_value', 'different_exact_value')),
        'september_exact_agreements': exact_partition['same_exact_value'],
        'september_exact_disagreements': exact_partition['different_exact_value'],
        'september_unresolved_exact': exact_partition['september_unresolved'],
        'september_disjoint_intervals': len(disputes),
        'all_larger_upper_bounds': len(retained_upper),
        'retained_larger_upper_bounds': sum(relation(table[k], september[k]['range']) !=
                                            'disjoint' for k in retained_upper),
        'september_strictly_weaker_intervals': sum(relation(table[k], september[k]['range']) ==
                                                  'paper_weaker' for k in table),
        'september_incomparable_intervals': sum(relation(table[k], september[k]['range']) ==
                                               'overlapping_incomparable' for k in table),
        'gebel_prangley_exact_agreements': len(gp_exact),
        'gebel_prangley_range_agreements': len(gp['theorem3_34']),
        'mccoy_total': len(cohort), 'mccoy_u1': len(positive), 'mccoy_excluded': len(negative),
        'mccoy_historical_subset': len(historical_mccoy),
        'mccoy_dataset_overlap': len(overlap),
        'mccoy_dataset_missing': len(cohort - overlap)}
    for version, rows in releases.items():
        counts['baseline_differences_' + version.replace('.', '_')] = sum(
            baseline[k] != rows[k]['range'] for k in baseline)
    paths = [RESULTS / 'paper_v1_1_snapshot.json', rollback_path, mccoy_path, dataset_path,
             gp_path, gp_check_path,
             RESULTS / 'crossing_changes/mccoy_alternating_u1_2026-09-08.json',
             OUTPUT / 'releases/manifest.json', RESULTS /
             ('paper_v1_3_review_manifest.json' if consolidated.get('manuscript') == 'v1.3-review'
              else 'paper_v1_3_manifest.json' if consolidated.get('manuscript') == 'v1.3'
              else 'paper_v1_1_manifest.json')]
    summary = {
        'schema_version': 1,
        'description': 'Historical research baseline and public package comparison; '
                       'tabulated ranges are not new mathematical evidence.',
        'baseline_construction': 'Use database_knotinfo 2026.8.1 ranges, replacing the 104 '
            'entries in the deposited rollback list by their 2026.4.1 ranges. This also '
            'restores 27 lower bounds from 2 to 1. The result is not an unmodified release.',
        'comparison_release': '2026.9.1', 'counts': counts, 'disputes': disputes,
        'gebel_prangley_overlap': {
            'source_url': gp_check['source_url'],
            'exact_agreement_set_equals_theorems_1_and_2': True,
            'all_theorem_3_ranges_agree': True,
            'exact_knots': sorted(gp_exact),
            'range_knots': sorted(gp['theorem3_34'])},
        'exact_partition': dict(sorted(exact_partition.items())),
        'incomparable_intervals': [r for r in comparison_rows if r['relation'] ==
                                  'overlapping_incomparable'],
        'interpretation': [
            f'The {len(exact)} count refers to exact determinations relative to the historical baseline; '
            'it is not a count of previously unpublished values.',
            f'The September package has exact entries for {counts["september_exact_entries"]} of these knots: '
            f'{counts["september_exact_agreements"]} agree and {counts["september_exact_disagreements"]} disagree. '
            f'The other {counts["september_unresolved_exact"]} entries are unresolved in that package.',
            'The 194 agreeing exact entries coincide exactly with the knots in Theorems 1 and 2 '
            'of Gebel--Prangley; their 30 Theorem 3 ranges also agree with the reconstruction.',
            f'{len(disputes)} disjoint intervals require explicit reconciliation; no contradictory '
            'release range is silently intersected with a proved lower bound.',
            f'Excluding the {len(disputes)} conflicts, {counts["retained_larger_upper_bounds"]} retained upper bounds '
            f'exceed the September bounds: {counts["september_strictly_weaker_intervals"]} intervals are strictly '
            f'weaker and {counts["september_incomparable_intervals"]} have a stronger lower bound but a weaker upper bound.',
            'McCoy verification is a consistency audit of deposited computations and tabulated '
            'release inputs, not a new scan of all diagrams. The crossing dataset covers 3105 '
            'of the 3627 knots and agrees on every stored unknotting-crossing list.'],
        'input_sha256': {str(p.relative_to(ROOT)): digest(p) for p in paths}}
    mccoy_rows = [{'knot': k, 'april': april[k]['range'], 'june': june[k]['range'],
                   'verdict': mccoy[k]['verdict'],
                   'unknotting_crossings': mccoy[k]['unknotting_crossings'],
                   'in_crossing_dataset': k in overlap} for k in sorted(cohort)]
    override_rows = [{'knot': k, 'april': april[k]['range'], 'august': august[k]['range'],
                      'baseline': baseline[k], 'lower_restored': baseline[k][0] !=
                      august[k]['range'][0]} for k in sorted(overrides)]
    return summary, {'baseline_overrides.csv': override_rows,
                     'september_comparison.csv': comparison_rows,
                     'mccoy_release_comparison.csv': mccoy_rows}


def write_csv(rows):
    output = io.StringIO(newline='')
    writer = csv.DictWriter(output, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows({k: json.dumps(v, separators=(',', ':')) if isinstance(v, (list, dict))
                     else v for k, v in row.items()} for row in rows)
    return output.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze-from', type=Path, help='directory of official wheels and metadata')
    parser.add_argument('--verify', action='store_true', help='check deposited comparisons without writes')
    args = parser.parse_args()
    if args.freeze_from:
        freeze_releases(args.freeze_from)
    summary, details = audit()
    outputs = {'baseline_comparison.json': json_text(summary),
               **{name: write_csv(rows) for name, rows in details.items()}}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name, content in outputs.items():
        path = OUTPUT / name
        if args.verify:
            if path.read_text() != content:
                raise ValueError(f'Deposited comparison is stale: {name}')
        else:
            path.write_text(content)
    print(json_text(summary['counts']), end='')


if __name__ == '__main__':
    main()
