#!/usr/bin/env python3
"""Audit the enlarged dichotomy without machine-learning or search scores.

Every crossing is checked for coverage and a necessary-condition exclusion.
The review compares against the committed pre-update records, and checks that
each changed prime-neighbour interval is supplied by the Greene summary.
This replays classification and data provenance; it does not rerun all
geometric identifications or Greene state calculations.
"""
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

def old_file(path):
    return subprocess.check_output(['git', 'show', '6468494:' + path], cwd=ROOT)

def read_rows(data):
    return {(r['knot'], r['crossing']): r
            for r in map(json.loads, gzip.decompress(data).decode().splitlines())}

def exclusion(row):
    """Independent proof-sensitive decision; probability fields are unused."""
    sigma = row['sigma_c']
    if sigma is not None and abs(sigma) > 2:
        return 'signature'
    how, name, interval = row['result_how'], row['result'], row['u_range']
    if (name is not None and how in {'code', 'trivial', 'isometry',
            'isometry-jones', 'prime-isometry', 'sum-blocks'}
            and interval is not None and interval[0] >= 2):
        return 'identified lower bound'
    return None

def main():
    new_bytes = (ROOT / 'results/open23/rows833.jsonl.gz').read_bytes()
    new = read_rows(new_bytes)
    old = read_rows(old_file('results/open23/rows833.jsonl.gz'))
    assert set(new) == set(old) and len(new) == 13105
    final = json.loads((ROOT / 'results/open23/priority_u23_final.json').read_text())
    old_final = json.loads(old_file('results/open23/priority_u23_final.json'))
    sweep = json.loads((ROOT / 'results/greene/sweep_2026-09-08.json').read_text())['bounds']
    table = json.loads((ROOT / 'generated/u_table.json').read_text())
    by_knot = defaultdict(list)
    for row in new.values():
        by_knot[row['knot']].append(row)
    for knot, rows in by_knot.items():
        n = rows[0]['crossings']
        assert sorted(r['crossing'] for r in rows) == list(range(n)), knot
        assert table[knot] == [2, 3], knot
    rigorous = {k for k, rows in by_knot.items() if all(exclusion(r) for r in rows)}
    claimed = {r[0] for r in final['zero_candidate_knots'] if r[3] == 0}
    old_rigorous = {r[0] for r in old_final['zero_candidate_knots'] if r[3] == 0}
    assert rigorous == claimed and len(rigorous) == 959
    assert old_rigorous <= rigorous and len(rigorous - old_rigorous) == 153
    changes = Counter()
    changed_intervals = []
    for key, before in old.items():
        after = new[key]
        for field in before.keys() | after.keys():
            if before.get(field) != after.get(field):
                changes[field] += 1
                assert field in {'result', 'u_range'}
        if before['result'] != after['result']:
            assert sorted(before['result'].split('#')) == sorted(after['result'].split('#'))
        if before['u_range'] != after['u_range']:
            result = after['result']
            assert result in sweep, result
            assert after['u_range'] == sweep[result]['new_range'], result
            assert after['u_range'][0] == sweep[result]['lower_bound'], result
            changed_intervals.append({'parent': key[0], 'crossing': key[1],
                                      'child': result, 'old': before['u_range'],
                                      'new': after['u_range'],
                                      'identification': after['result_how']})
    dependencies = []
    for parent in sorted(rigorous - old_rigorous):
        decisive = []
        for row in by_knot[parent]:
            previous = old[(parent, row['crossing'])]
            if exclusion(previous) is None:
                assert previous['u_range'][0] == 1
                assert row['result'] in sweep and sweep[row['result']]['lower_bound'] >= 2
                decisive.append({'crossing': row['crossing'], 'child': row['result']})
        assert decisive, parent
        dependencies.append({'parent': parent, 'new_exclusions': decisive})
    child_parents = {p for child in final['children'] for p in child['parents']}
    unresolved = json.loads((ROOT / 'results/open23/unresolved_candidates_resolution.json').read_text())
    moved = set(final['moved_to_tierB_after_resolution'])
    assert len(moved) == 7
    for parent in moved:
        records = [r for r in unresolved if r['knot'] == parent]
        assert records and all(r['verdict'][0] == 'composite(jones)' for r in records)
        assert parent not in rigorous
    flypes = {r['knot']: r for r in map(json.loads,
        (ROOT / 'results/open23/flype_orbit_check.jsonl').read_text().splitlines())}
    assert rigorous <= set(flypes)
    orbit_unresolved = {k: len(flypes[k]['candidates']) for k in sorted(rigorous)
                        if flypes[k]['candidates']}
    report = {
        'rows_sha256': hashlib.sha256(new_bytes).hexdigest(),
        'crossing_rows': len(new), 'parents': len(by_knot),
        'rigorous_parents': len(rigorous), 'previous_rigorous': len(old_rigorous),
        'added_rigorous': len(rigorous - old_rigorous),
        'all_parent_crossings_present_once': True,
        'rigorous_set_independently_reconstructed': True,
        'changed_fields': dict(changes), 'changed_intervals': changed_intervals,
        'new_parent_dependencies': dependencies,
        'new_decisive_identification_methods': dict(Counter(
            new[(p['parent'], r['crossing'])]['result_how']
            for p in dependencies for r in p['new_exclusions'])),
        'named_children': len(final['children']), 'named_child_parents': len(child_parents),
        'candidate_parent_list': len(final['knots_with_candidates']),
        'jones_only_moved_parents': sorted(moved),
        'jones_only_zero_candidate_parents': sum(r[3] != 0 for r in final['zero_candidate_knots']),
        'rigorous_parents_with_unresolved_supplementary_orbit_identifications': orbit_unresolved,
        'limits': ['Existing numerical SnapPy identifications are inherited evidence.',
                   'This audit verifies classification and dependence on the separately audited Greene bounds.',
                   'The seven Jones-only moved parents are not included among the 959 proved dichotomies.',
                   'Five supplementary flype checks retain unidentified neighbours; the reference-diagram exclusions and flyping lemma establish their dichotomies.']}
    (HERE / 'updated_dichotomy_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print({k: v for k, v in report.items() if k not in
           {'changed_intervals', 'new_parent_dependencies', 'limits'}})

if __name__ == '__main__':
    main()
