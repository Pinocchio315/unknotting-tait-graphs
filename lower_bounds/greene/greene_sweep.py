#!/usr/bin/env python
"""Sweep the correction-term obstructions over knots whose double branched cover is an L-space.

Every knot in `data/greene_targets.json` has reduced Khovanov homology over F_2 of rank equal to its
determinant, so its double branched cover is an L-space and Greene's spanning-tree model determines the
correction terms from a diagram.  The recorded range then selects the test:

    range [1, b]   the half-integral surgery pattern of Ni and Wu excludes u = 1
    range [n, b]   with |sigma| = 2n and n in {2, 3, 4}: every definite form of rank n is compared with
                   the correction terms, which excludes u = n

A verdict of OBSTRUCTED raises the lower bound to n + 1; PASS, UNDECIDED, TIMEOUT and ERROR supply no
bound.  Controls have a known unknotting number equal to n, so OBSTRUCTED on a control is an
implementation error and is reported as such.

    python greene_sweep.py --shard 0/16                      # one shard of the target list
    python greene_sweep.py --workers 16                      # one node, 16 concurrent knots
    python greene_sweep.py --controls --workers 16           # the validation set
    python greene_sweep.py --one 12n_491                     # a single knot, one JSON line on stdout

Each knot runs in its own subprocess, so a failure or a memory spike is confined to that knot.  Results
are appended to a JSONL file and finished knots are skipped on a rerun.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_TARGETS = HERE / 'data' / 'greene_targets.json'


# --------------------------------------------------------------------------- one knot
def run_one(entry):
    """Compute the verdict for a single knot.  Imports live here so the orchestrator stays light."""
    sys.path.insert(0, str(HERE))
    from pin_dinv import pin
    from greene_dinv import surgery_test
    from greene_ranks import rank_test, oriented
    import itertools

    name, D, sigma = entry['name'], int(entry['determinant']), int(entry['signature'])
    lo = int(entry['range'][0])
    n = max(lo, 1)
    started = time.time()
    diagnostics = {}
    current, pairing, reference = pin(name, entry['pd'], D, verbose=False, diagnostics=diagnostics)
    ambiguous = {k: sorted(v) for k, v in current.items() if len(v) != 1}
    pinned = {k: next(iter(v)) for k, v in current.items() if len(v) == 1}
    free = sorted(k for k in ambiguous if k <= -k % D)

    admitted, total, verdicts = 0, 0, []
    for choice in itertools.product(*(ambiguous[k] for k in free)):
        vector = dict(pinned)
        for k, value in zip(free, choice):
            vector[k] = vector[-k % D] = value
        if set(vector) != set(range(D)):
            raise ValueError('a candidate vector omits Spin^c classes')
        if n == 1:
            outcome = surgery_test(vector, D)
            fits = bool(outcome['fits'])
            verdicts.append({'test': 'u1', 'admits': fits, 'detail': outcome})
        else:
            outcome = rank_test(oriented(vector, sigma), D, n)
            fits = outcome['verdict'] != 'OBSTRUCTED'
            verdicts.append({'test': f'rank{n}', 'admits': fits, 'detail': outcome})
        total += 1
        admitted += fits
    if total == 0:
        raise ValueError('no complete candidate vector was tested')

    inconclusive = [v for v in verdicts if v['admits'] and v['detail'].get('verdict') in ('UNDECIDED', 'ERROR')]
    verdict = ('OBSTRUCTED' if admitted == 0 else
               'UNDECIDED' if inconclusive else 'PASS')
    record = {'name': name, 'det': D, 'sigma': sigma, 'range': entry['range'],
              'crossing_number': entry.get('crossing_number'), 'test': entry['test'],
              'role': 'control' if 'known_unknotting_number' in entry else 'target',
              'verdict': verdict, 'excluded_length': n,
              'pages': diagnostics.get('pages'), 'intersection_passes': diagnostics.get('intersection_passes'),
              'unresolved_classes': len(ambiguous), 'candidate_vectors': total, 'admitting_vectors': admitted,
              'pairing': str(pairing), 'reference_marking': list(reference),
              'candidate_verdicts': verdicts, 'seconds': round(time.time() - started, 1)}
    if verdict == 'OBSTRUCTED':
        record['lower_bound'] = n + 1
        record['new_range'] = [n + 1, entry['range'][1]]
    if 'known_unknotting_number' in entry:
        record['known_unknotting_number'] = entry['known_unknotting_number']
        record['control_ok'] = verdict != 'OBSTRUCTED'
    if ambiguous:
        record['pinned'] = {str(k): str(v) for k, v in pinned.items()}
        record['ambiguous'] = {str(k): [str(x) for x in v] for k, v in ambiguous.items()}
    else:
        record['d'] = {str(k): str(v) for k, v in pinned.items()}
    return record


# --------------------------------------------------------------------------- orchestration
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--targets', type=Path, default=DEFAULT_TARGETS)
    ap.add_argument('--one', help='compute this knot and print one JSON line (used by the workers)')
    ap.add_argument('--controls', action='store_true', help='sweep the control set instead of the targets')
    ap.add_argument('--test', help='restrict to one test kind: u1, rank2, rank3, rank4')
    ap.add_argument('--shard', default='0/1', help='i/n split of the list for this task')
    ap.add_argument('--workers', type=int, default=1, help='concurrent knots in this task')
    ap.add_argument('--out', type=Path, help='JSONL result file (default results_<set><shard>.jsonl)')
    ap.add_argument('--time-per-knot', type=float, default=10800.0,
                    help='seconds before a knot is abandoned; a timeout is not a negative result')
    ap.add_argument('--max-det', type=int, default=0, help='skip knots of larger determinant (0 = no limit)')
    args = ap.parse_args()

    data = json.loads(args.targets.read_text())
    which = 'controls' if args.controls else 'targets'
    entries = data[which]

    if args.one:
        print(json.dumps(run_one(entries[args.one]), sort_keys=True), flush=True)
        return

    order = [n for n in data['control_order' if args.controls else 'target_order'] if n in entries]
    if args.test:
        order = [n for n in order if entries[n]['test'] == args.test]
    if args.max_det:
        order = [n for n in order if entries[n]['determinant'] <= args.max_det]
    i, total_shards = (int(x) for x in args.shard.split('/'))
    order = order[i::total_shards]

    out = args.out or HERE / f'results_{which}_{i}_{total_shards}.jsonl'
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            try:
                done.add(json.loads(line)['name'])
            except Exception:
                pass
    todo = [n for n in order if n not in done]
    print(f'{which}: {len(order)} in shard {i}/{total_shards}, {len(done & set(order))} already recorded, '
          f'{len(todo)} to run, {args.workers} workers -> {out}', flush=True)

    def work(name):
        started = time.time()
        command = [sys.executable, str(Path(__file__).resolve()), '--targets', str(args.targets),
                   '--one', name] + (['--controls'] if args.controls else [])
        try:
            finished = subprocess.run(command, capture_output=True, text=True, timeout=args.time_per_knot)
            lines = [l for l in finished.stdout.splitlines() if l.startswith('{')]
            record = json.loads(lines[-1]) if lines else {
                'name': name, 'verdict': 'ERROR', 'stderr': finished.stderr[-800:],
                'returncode': finished.returncode}
        except subprocess.TimeoutExpired:
            record = {'name': name, 'verdict': 'TIMEOUT', 'seconds': round(time.time() - started, 1)}
        except Exception as exc:                                   # noqa: BLE001 - one knot must not stop the shard
            record = {'name': name, 'verdict': 'ERROR', 'error': repr(exc)[:400]}
        record.setdefault('seconds', round(time.time() - started, 1))
        record.setdefault('role', 'control' if args.controls else 'target')
        with open(out, 'a') as handle:            # one short line per knot; O_APPEND keeps concurrent writes intact
            handle.write(json.dumps(record, sort_keys=True) + '\n')
            handle.flush()
            os.fsync(handle.fileno())
        flag = ''
        if record.get('control_ok') is False:
            flag = '   *** CONTROL OBSTRUCTED: implementation error ***'
        print(f"{name} {record['verdict']} {record.get('seconds')}s{flag}", flush=True)
        return record

    with ThreadPoolExecutor(max(1, args.workers)) as pool:
        results = list(pool.map(work, todo))

    from collections import Counter
    print('done ' + str(dict(sorted(Counter(r['verdict'] for r in results).items()))), flush=True)
    bad = [r['name'] for r in results if r.get('control_ok') is False]
    if bad:
        print('CONTROLS OBSTRUCTED (implementation error): ' + ' '.join(bad), flush=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
