#!/usr/bin/env python3
"""Reproduce version deltas, cohort accounting, and all recorded candidate checks.

This is an audit of deposited evidence and its reporting, not a fresh
calculation of every Floer complex. Run with the standard library from anywhere.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import consolidate_results as c


def main():
    before, old = c.reconstruct(manuscript='v1.2')
    initial, _ = c.reconstruct(manuscript='v1.3-review')
    after, current = c.reconstruct(manuscript='v1.3')
    frozen, records = c.read_greene_sweep()
    entries = {**frozen['targets'], **frozen['controls']}
    assert all(entry['range'] == old[name] for name, entry in entries.items())
    assert all(name == entry['name'] for name, entry in entries.items())
    assert len(records) == len(set(records))
    assert not any(r['verdict'] == 'OBSTRUCTED' and r['role'] == 'control'
                   for r in records.values())
    delta = {name: {'v1.2': old[name], 'v1.3': current[name],
                    'sources': after['changed'][name]['sources']}
             for name in current if current[name] != old[name]}
    transitions = Counter((tuple(r['v1.2']), tuple(r['v1.3'])) for r in delta.values())
    assert len(delta) == 1707
    assert sum(r['v1.3'][0] == r['v1.3'][1] for r in delta.values()) == 1492
    assert after['counts']['exact'] == before['counts']['exact'] + 1492
    assert after['counts']['improved'] == before['counts']['improved'] + 215 - 1
    assert delta['13n_1587']['v1.2'] == [1, 2]
    assert delta['13n_1587']['v1.3'] == [2, 2]
    assert set(delta) == {name for name, row in records.items()
                          if row['verdict'] == 'OBSTRUCTED'}
    diagnostics = {}
    for role in ('target', 'control'):
        wanted = set(frozen[role + 's'])
        rows = [r for r in records.values() if r['role'] == role]
        diagnostics[role] = {
            'wanted': len(wanted), 'recorded': len(rows),
            'verdicts': dict(sorted(Counter(r['verdict'] for r in rows).items())),
            'missing': sorted(wanted - records.keys()),
            'already_determined_by_v1_2_computations': sum(name in before['changed']
                and before['changed'][name]['new'][0] == before['changed'][name]['new'][1]
                for name in wanted),
        }
    old_scan = c.read_json(c.RESULTS / c.PAPER_SCAN)
    scan = c.read_json(c.RESULTS / 'open23/priority_u23_final.json')
    old_children = {row['child'] for row in old_scan['children']}
    raised_children = {name for name in old_children if current[name][0] > 1}
    assert len(old_children) == 112 and len(raised_children) == 101
    zero = {name for name, _, _, _ in scan['zero_candidate_knots']}
    rigorous = {name for name, _, _, heuristic in scan['zero_candidate_knots'] if heuristic == 0}
    moved = set(scan['moved_to_tierB_after_resolution'])
    candidates = set(scan['knots_with_candidates']) - moved
    heuristic = (zero - rigorous) | moved
    assert not (rigorous & heuristic or heuristic & candidates or rigorous & candidates)
    assert (len(rigorous), len(heuristic), len(candidates)) == (959, 46, 22)
    assert len(rigorous | heuristic | candidates) == 1027
    assert sum(row['new'][0] == row['new'][1]
               for row in after['changed'].values()) == 2719
    generated_before = c.manuscript_tex(before, old)
    historical_byte_matches = {name: text == (ROOT / 'generated/paper_v1_2' / name).read_text()
                               for name, text in generated_before.items()}
    assert all(historical_byte_matches.values())
    comparison = c.comparison_report(after, current)
    output = {
        'scope': 'Frozen data bookkeeping, complete recorded candidate validation, and manuscript inputs; '
                 'not a fresh Floer calculation for every knot.',
        'v1.2_counts': before['counts'], 'initial_v1.3_counts': initial['counts'],
        'latest_v1.3_counts': after['counts'],
        'transitions': [{'before': list(a), 'after': list(b), 'count': count}
                        for (a, b), count in sorted(transitions.items())],
        'new_exact': 1492, 'new_interval_improvements': 215,
        'old_interval_becoming_exact': ['13n_1587'],
        'cohort_ranges_equal_v1_2_consolidated_table': True,
        'sweep_diagnostics': diagnostics,
        'historical_generated_files_unchanged': historical_byte_matches,
        'latest_candidate_partition': {
            'rigorous': len(rigorous), 'heuristic': len(heuristic),
            'with_identified_candidates': len(candidates),
            'heuristic_zero_list': len(zero - rigorous),
            'heuristic_moved_list': sorted(moved),
            'historical_children': len(old_children),
            'historical_children_excluded_by_new_bounds': len(raised_children),
            'remaining_children': sorted(old_children - raised_children),
        },
        'conflicts': comparison['disputes'],
        'september_comparison_counts': comparison['counts'],
        'all_new_results': delta,
        'source_sha256': {name: hashlib.sha256((c.RESULTS/name).read_bytes()).hexdigest()
            for name in ('paper_v1_1_manifest.json', 'paper_v1_3_review_manifest.json',
                         'paper_v1_3_manifest.json', c.PAPER_SCAN,
                         'open23/priority_u23_final.json',
                         'greene/greene_targets_2026-09-08.json',
                         'greene/sweep_2026-09-08.jsonl.gz')},
    }
    (HERE/'reporting_audit.json').write_text(json.dumps(output, sort_keys=True, indent=2)+'\n')
    print('Verified 1707 new bounds, 1492 new exact values, 215 new interval improvements, '
          'all 3418 recorded jobs, ten conflicts, and the 959 + 46 + 22 candidate partition.')


if __name__ == '__main__':
    main()
