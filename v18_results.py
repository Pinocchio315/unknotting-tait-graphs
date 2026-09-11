"""Reconstruct v1.8 while keeping the historical v1.3 inputs reproducible.

The lower-bound stage precedes the new upper-bound constructions. In ten
cases both stages are needed: a G tag must not make a new upper bound vanish
from the accounting. This module does bookkeeping, not knot recognition.
Mathematical certificates and their independent replays are deposited with
the inputs authenticated by paper_v1_8_manifest.json.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path

import consolidate_results as legacy

RESULTS = legacy.RESULTS


def inputs(results=RESULTS):
    results = Path(results)
    manifest = legacy.read_json(results / 'paper_v1_8_manifest.json')
    for name, expected in manifest['sha256'].items():
        if hashlib.sha256((results / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Changed v1.8 input: {name}')
    rows = legacy.read_json(results / 'extensions_2026-09-10/bounds.json')
    if len({(r['name'], r['side']) for r in rows}) != len(rows):
        raise ValueError('Duplicate extension bound')
    bounds = [legacy.Bound(r['name'], r['side'], r['value'], r['tag'], r['source']) for r in rows]
    return manifest, rows, bounds


def reconstruct(results=RESULTS):
    results = Path(results)
    baseline, _ = legacy.reconstruct(results, 'v1.3')
    manifest, rows, extra = inputs(results)
    old_bounds, rank4 = legacy.collect_bounds(results)
    snapshot = legacy.read_json(results / 'paper_v1_1_snapshot.json')
    all_bounds = old_bounds + extra
    table, changed = legacy.merge_bounds(snapshot, all_bounds)
    # 'w' denotes the explicit constructions supplied by this paper. Published
    # upper bounds, including Brittenham--Hermiller's, remain in the first stage.
    lower_table, lower_changed = legacy.merge_bounds(snapshot, [b for b in all_bounds if b.tag != 'w'])
    lower = legacy.summarize(lower_table, lower_changed)
    totals = legacy.summarize(table, changed)
    totals['exact_with_lower_tag'] = totals['exact_lower']
    totals['exact_lower'] = lower['exact']
    totals['exact_upper'] = totals['exact'] - lower['exact']
    totals['lower_improved'] = lower['improved']
    totals['joint_exact'] = sum(r['new'][0] == r['new'][1] and bool(r['lower_tags'])
                                and 'w' in r['upper_tags'] for r in changed.values())
    # The same frozen 9 September snapshot is used for both manuscript versions.
    from compare_knotinfo import load_snapshot
    official = load_snapshot(results)
    totals['latest_unresolved_exact'] = sum(r['new'][0] == r['new'][1]
        and official[n][0] != official[n][1] for n, r in changed.items())
    totals['latest_range_improvements'] = sum(r['new'][0] != r['new'][1]
        and r['new'] != official[n]
        and official[n][0] <= r['new'][0] <= r['new'][1] <= official[n][1]
        for n, r in changed.items())
    if any(not official[n][0] <= lo <= hi <= official[n][1]
           for n, (lo, hi) in table.items()):
        raise ValueError('A v1.8 result is incompatible with the frozen official snapshot')
    for key, value in manifest['expected'].items():
        if totals[key] != value:
            raise ValueError(f'v1.8 count mismatch: {key}: {totals[key]} != {value}')
    return {'manuscript': 'v1.8', 'reference': baseline['reference'], 'counts': totals,
            'changed': changed, 'lower_stage_changed': lower_changed,
            'lower_stage_counts': lower, 'extension_rows': rows,
            'rank4_status': rank4}, table


def counts(consolidated, table, results=RESULTS):
    old, old_table = legacy.reconstruct(results, 'v1.3')
    out = legacy.computation_counts(old, old_table, results)
    total, lower = consolidated['counts'], consolidated['lower_stage_counts']
    extra_lower = [r for r in consolidated['extension_rows'] if r['side'] == 'lo']
    tests = Counter(r['test'] for r in extra_lower)
    new_exact = sum(table[r['name']][0] == table[r['name']][1]
                    and old_table[r['name']][0] != old_table[r['name']][1] for r in extra_lower)
    for macro, key in [('nExact', 'exact'), ('nExactLower', 'exact_lower'),
                       ('nExactUpper', 'exact_upper'), ('nImproved', 'improved'),
                       ('nLowerImproved', 'lower_improved'), ('nJointExact', 'joint_exact')]:
        out[macro] = total[key]
    for u, macro in [(2, 'nUtwo'), (3, 'nUthree'), (4, 'nUfour'), (5, 'nUfive')]:
        out[macro] = total['exact_by_u'][str(u)]
    out['nUtwoG'] = sum(r['new'] == [2, 2] and legacy.primary_tag(r) == 'G'
                       for r in consolidated['changed'].values())
    out['nSweepExact'] = lower['exact_by_primary_method']['G']
    out['nSweepImproved'] = lower['improved_by_primary_method']['G']
    out['nSweepNewExact'] += new_exact
    out['nSweepUone'] += tests['u1']
    out['nSweepRankFour'] += tests['rank4']
    out['nSweepObstructed'] += len(extra_lower)
    frozen, records = legacy.read_greene_sweep(results)
    settled = {'OBSTRUCTED', 'PASS', 'UNDECIDED', 'NOT_APPLICABLE'}
    out['nSweepTargetsSettled'] += sum(records.get(r['name'], {}).get('verdict') not in settled
                                      for r in extra_lower)
    # Literature recoveries and existing agreements are not changed by these 15 knots.
    out['nLatestUnresolvedExact'] = total['latest_unresolved_exact']
    out['nLatestRangeImprovements'] = total['latest_range_improvements']
    return out
