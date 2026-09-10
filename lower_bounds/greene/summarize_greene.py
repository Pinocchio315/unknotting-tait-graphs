#!/usr/bin/env python
"""Combine the sweep result files, check the controls, and list the bounds obtained.

The sweep writes one JSONL file per task, so a run split between a cluster and a local machine leaves
several files.  This script reads any number of them, keeps one record per knot, reports the coverage
against the frozen target list, and writes the new lower bounds in the form the consolidation reads.

A knot may appear in more than one file, for instance when a shard was recomputed elsewhere.  A
conclusive verdict replaces an error or a timeout; two different conclusive verdicts for the same knot
are a contradiction and are reported rather than silently resolved.

The reader accepts plain JSONL; decompress archived .jsonl.gz inputs first.
For a paper replay, explicitly select the deposited pre-sweep cohort instead
of relying on the local data/greene_targets.json default. From this directory:

    python summarize_greene.py runs/targets/*.jsonl runs/controls/*.jsonl \
        --targets ../../results/greene/greene_targets_2026-09-08.json \
        --out /tmp/greene-summary.json
"""
from __future__ import annotations
import argparse, json, sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
CONCLUSIVE = {'OBSTRUCTED', 'PASS', 'UNDECIDED', 'NOT_APPLICABLE'}


def load(paths):
    """One record per knot: a conclusive verdict wins over an error or a timeout."""
    best, clashes, files = {}, [], 0
    for path in paths:
        if not path.exists():
            continue
        files += 1
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f'{path}:{lineno}: malformed result record') from exc
            name = record.get('name')
            if name is None:
                raise ValueError(f'{path}:{lineno}: result record has no knot name')
            record['source_file'] = path.name
            old = best.get(name)
            if old is None:
                best[name] = record
                continue
            new_ok = record['verdict'] in CONCLUSIVE
            old_ok = old['verdict'] in CONCLUSIVE
            if new_ok and old_ok and record['verdict'] != old['verdict']:
                clashes.append((name, old['verdict'], old['source_file'], record['verdict'], path.name))
            if new_ok and not old_ok:
                best[name] = record
    return best, clashes, files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('files', nargs='+', type=Path)
    ap.add_argument('--targets', type=Path, default=HERE / 'data' / 'greene_targets.json')
    ap.add_argument('--out', type=Path, help='write the combined bounds to this JSON file')
    ap.add_argument('--list-missing', action='store_true', help='print the names still to compute')
    args = ap.parse_args()

    data = json.loads(args.targets.read_text())
    try:
        records, clashes, files = load(args.files)
    except ValueError as exc:
        ap.error(str(exc))
    print(f'{files} result files, {len(records)} knots\n')

    if clashes:
        print('CONTRADICTORY VERDICTS (the same knot settled two ways):')
        for name, a, fa, b, fb in clashes:
            print(f'   {name}: {a} in {fa} vs {b} in {fb}')
        print()

    summary = {}
    for which in ('controls', 'targets'):
        wanted = set(data[which])
        got = {n: r for n, r in records.items() if n in wanted}
        settled = {n: r for n, r in got.items() if r['verdict'] in CONCLUSIVE}
        missing = sorted(wanted - set(settled))
        counts = Counter(r['verdict'] for r in got.values())
        by_test = defaultdict(Counter)
        for r in got.values():
            by_test[r.get('test') or data[which][r['name']]['test']][r['verdict']] += 1
        print(f'{which}: {len(settled)} settled of {len(wanted)}   ' + str(dict(sorted(counts.items()))))
        for test, c in sorted(by_test.items()):
            print(f'   {test:6s} ' + str(dict(sorted(c.items()))))
        if missing:
            print(f'   still to compute: {len(missing)}' + ('  ' + ' '.join(missing) if args.list_missing else ''))
        summary[which] = {'wanted': len(wanted), 'settled': len(settled), 'verdicts': dict(counts),
                          'missing': missing}
        print()

    # The frozen control list, not an optional flag supplied by a worker,
    # determines whether an obstruction contradicts an adopted exact value.
    bad = sorted(n for n, r in records.items() if r.get('control_ok') is False
                 or (n in data['controls'] and r['verdict'] == 'OBSTRUCTED'))
    if bad:
        print('*** CONTROLS OBSTRUCTED: the implementation contradicts a known unknotting number ***')
        print('    ' + ' '.join(bad))
        print()
    else:
        settled_controls = summary['controls']['settled']
        print(f'control check: no obstructed control among the {settled_controls} settled\n')

    # Never publish a bounds file after detecting inconsistent source records.
    # In particular, an earlier OBSTRUCTED record must not survive a later PASS
    # merely because it happened to appear first in the command-line file list.
    if clashes or bad:
        sys.exit('No bounds written: resolve conflicting verdicts or obstructed controls first')

    bounds = {}
    for name, r in sorted(records.items()):
        if name in data['targets'] and r['verdict'] == 'OBSTRUCTED':
            bounds[name] = {'reference_range': data['targets'][name]['range'],
                            'new_range': r['new_range'], 'lower_bound': r['lower_bound'],
                            'test': r['test'], 'determinant': r['det'], 'signature': r['sigma'],
                            'unresolved_classes': r.get('unresolved_classes'),
                            'candidate_vectors': r.get('candidate_vectors')}
    exact = {n: b for n, b in bounds.items() if b['new_range'][0] == b['new_range'][1]}
    print(f'new lower bounds: {len(bounds)}   of which the range becomes exact: {len(exact)}')
    print('   by test:   ' + str(dict(sorted(Counter(b['test'] for b in bounds.values()).items()))))
    print('   new exact values by u: ' + str(dict(sorted(Counter(b['new_range'][0] for b in exact.values()).items()))))
    widened = Counter((tuple(b['reference_range']), tuple(b['new_range'])) for b in bounds.values())
    print('   range changes: ' + str({f'{a} -> {c}': n for (a, c), n in sorted(widened.items())}))

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({'summary': summary, 'contradictions': clashes, 'bounds': bounds},
                                       sort_keys=True, indent=1) + '\n')
        print(f'\nwrote {args.out}')


if __name__ == '__main__':
    main()
