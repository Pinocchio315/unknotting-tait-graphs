#!/usr/bin/env python3
"""Compare computed intervals with the frozen official KnotInfo ranges.

This numerical comparison identifies which exact determinations were still open
on 9 September 2026. It does not infer mathematical bounds from database changes.
No manuscript, TeX generator, network request, or live worker output is used.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from consolidate_results import RESULTS, read_json, reconstruct, validate_interval


def load_snapshot():
    path = RESULTS / 'knotinfo_2026-09-09.json'
    manifest = read_json(RESULTS / 'knotinfo_2026-09-09_manifest.json')
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['sha256']:
        raise ValueError('The frozen KnotInfo ranges have changed')
    rows = read_json(path)
    if len(rows) != manifest['rows']:
        raise ValueError('The frozen KnotInfo population has changed')
    for interval in rows.values():
        validate_interval(interval)
    return rows


def compare(consolidated, table, snapshot):
    if set(table) != set(snapshot):
        raise ValueError('Computed and reference knot populations differ')
    changed = consolidated['changed']
    if not set(changed) <= set(table):
        raise ValueError('Changed knot outside the computation population')
    groups = {key: [] for key in ('unresolved_exact', 'exact_agreements',
                                 'range_improvements', 'range_agreements',
                                 'other_changed_ranges', 'disjoint_intervals')}
    for name in sorted(table):
        result, reference = table[name], snapshot[name]
        validate_interval(result)
        validate_interval(reference)
        if name in changed and changed[name]['new'] != result:
            raise ValueError(f'Changed/result interval mismatch: {name}')
        if result[0] > reference[1] or reference[0] > result[1]:
            groups['disjoint_intervals'].append(name)
            continue
        if name not in changed:
            continue
        if result[0] == result[1]:
            key = 'exact_agreements' if reference[0] == reference[1] else 'unresolved_exact'
        elif result == reference:
            key = 'range_agreements'
        elif reference[0] <= result[0] <= result[1] <= reference[1]:
            key = 'range_improvements'
        else:
            key = 'other_changed_ranges'
        groups[key].append(name)
    return {'snapshot_date': '2026-09-09',
            'counts': {'study_knots': len(table), **{k: len(v) for k, v in groups.items()}},
            'knots': groups}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, help='optionally write the per-knot JSON comparison')
    args = parser.parse_args()
    consolidated, table = reconstruct()
    result = compare(consolidated, table, load_snapshot())
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result['counts'], indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
