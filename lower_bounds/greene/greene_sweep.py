#!/usr/bin/env python
"""Sweep the correction-term obstructions over knots whose double branched cover is an L-space.

Every knot in `data/greene_targets.json` has reduced Khovanov homology over F_2 of rank equal to its
determinant, so its double branched cover is an L-space and Greene's spanning-tree model gives finite
sets containing the correction terms. Every vector remaining after comparing the marked diagrams is
tested. The recorded input range, frozen before this sweep, then selects the test:

    range [1, b]   the half-integral surgery pattern of Ni and Wu excludes u = 1
    range [n, b]   with |sigma| = 2n and n in {2, 3, 4}: every definite form of rank n is compared with
                   the correction terms, which excludes u = n

A verdict of OBSTRUCTED raises the lower bound to n + 1; PASS, UNDECIDED, TIMEOUT and ERROR supply no
bound. PASS only means that a necessary surgery comparison admits a candidate, not that an unknotting
sequence exists. Controls have a known unknotting number equal to n, so OBSTRUCTED on a control is an
implementation error and is reported as such.

    python greene_sweep.py --shard 0/16                      # one shard of the target list
    python greene_sweep.py --workers 16                      # one node, 16 concurrent knots
    python greene_sweep.py --controls --workers 16           # the validation set
    python greene_sweep.py --one 13n_111                     # a single target, one JSON line on stdout

Each knot runs in its own subprocess, so a failure or a memory spike is confined to that knot.  Results
are appended to a JSONL file and finished knots are skipped on a rerun.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time, threading
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_TARGETS = HERE / 'data' / 'greene_targets.json'


# --------------------------------------------------------------------------- one knot
def validate_entry(entry):
    """Verify the frozen hypotheses before any solitary-state grading is used.

    The input builder performs these checks when it freezes the table. Repeating
    them here prevents an edited or mismatched input file from turning the
    conditional L-space argument into an unconditional obstruction. The PD's
    determinant and cyclic labelling are checked independently by ``pin``.
    """
    from pin_dinv import khovanov_lspace_evidence
    D, sigma = entry['determinant'], entry['signature']
    if type(D) is not int or D < 3 or D % 2 == 0 or type(sigma) is not int:
        raise ValueError('An odd determinant greater than one and an integer signature are required')
    interval = entry['range']
    if (len(interval) != 2 or any(type(x) is not int for x in interval)
            or not 0 <= interval[0] <= interval[1]):
        raise ValueError('The frozen range must be an ordered pair of nonnegative integers')
    n = max(interval[0], 1)
    if n not in (1, 2, 3, 4):
        raise ValueError('Only tests for one through four crossing changes are implemented')
    if ((n == 1 and abs(sigma) > 2)
            or (n >= 2 and abs(sigma) != 2 * n)):
        raise ValueError('The signature does not satisfy the selected surgery test hypotheses')
    if entry['test'] != ('u1' if n == 1 else f'rank{n}'):
        raise ValueError('The recorded test does not agree with the lower end of the range')
    evidence = khovanov_lspace_evidence({
        'determinant': D,
        'khovanov_reduced_mod2_vector': entry['reduced_mod2_vector_raw'],
    })
    if entry.get('khovanov_rank', D) != evidence['rank']:
        raise ValueError('The cached Khovanov rank disagrees with the actual mod-2 vector')
    if 'known_unknotting_number' in entry and (
            type(entry['known_unknotting_number']) is not int
            or interval != [n, n] or entry['known_unknotting_number'] != n):
        raise ValueError('A control must have a recorded exact value equal to the tested length')
    return evidence


def run_one(entry):
    """Compute the verdict for a single knot.  Imports live here so the orchestrator stays light."""
    sys.path.insert(0, str(HERE))
    from pin_dinv import pin
    from greene_dinv import surgery_test
    from greene_ranks import rank_test, oriented
    import itertools

    validate_entry(entry)
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
    # One free entry for each conjugate pair, including the spin class when
    # ambiguous. There is no vector-count cap or early stop: the full Cartesian
    # product is tested. A process time limit supplies no mathematical verdict.
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
            # The signature fixes the cover's orientation; rank_test applies
            # its spin-entry guard only after this choice has been made.
            outcome = rank_test(oriented(vector, sigma), D, n, pairing=pairing)
            fits = outcome['verdict'] != 'OBSTRUCTED'
            verdicts.append({'test': f'rank{n}', 'admits': fits, 'detail': outcome})
        total += 1
        admitted += fits
    if total == 0:
        raise ValueError('no complete candidate vector was tested')

    # In particular, a candidate failing the spin-entry guard remains
    # inconclusive. Excluding only the other candidates cannot establish a bound.
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
    ap.add_argument('--skip-test', help='comma-separated test kinds to leave out, e.g. rank4')
    ap.add_argument('--names', type=Path,
                    help='restrict to the knots in this file (a JSON list, or one name per line); '
                         'useful for finishing a run whose shards were split across machines')
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

    if not args.one:
        sys.path.insert(0, str(HERE))
        try:
            import numpy, sympy                                    # noqa: F401
            import greene_dinv, half_integral, pin_dinv, greene_ranks   # noqa: F401
        except Exception as exc:                                   # noqa: BLE001
            sys.exit(f'the computation modules do not import here: {exc!r}\n'
                     'the interpreter needs numpy and sympy, and this directory needs greene_dinv.py, '
                     'half_integral.py, pin_dinv.py, greene_ranks.py and the owens modules '
                     '(owens_obstruction.py, owens_u3.py, owens_u4.py, linkform.py); '
                     'build the directory with make_server_package.py')

    if args.one:
        print(json.dumps(run_one(entries[args.one]), sort_keys=True), flush=True)
        return

    order = [n for n in data['control_order' if args.controls else 'target_order'] if n in entries]
    if args.test:
        order = [n for n in order if entries[n]['test'] == args.test]
    if args.skip_test:
        left_out = {t.strip() for t in args.skip_test.split(',') if t.strip()}
        order = [n for n in order if entries[n]['test'] not in left_out]
    if args.names:
        text = args.names.read_text().strip()
        wanted = set(json.loads(text)) if text.startswith('[') else {l.strip() for l in text.splitlines() if l.strip()}
        unknown = wanted - set(entries)
        if unknown:
            sys.exit(f'{len(unknown)} names are not in the {which} list, e.g. {sorted(unknown)[:5]}')
        order = [n for n in order if n in wanted]
    if args.max_det:
        order = [n for n in order if entries[n]['determinant'] <= args.max_det]
    i, total_shards = (int(x) for x in args.shard.split('/'))
    order = order[i::total_shards]

    out = args.out or HERE / f'results_{which}_{i}_{total_shards}.jsonl'
    CONCLUSIVE = {'OBSTRUCTED', 'PASS', 'UNDECIDED', 'NOT_APPLICABLE'}
    done, retry, previous_bad_controls = set(), set(), set()
    if out.exists():
        for line in out.read_text().splitlines():
            try:
                record = json.loads(line)
            except Exception:
                continue
            if (record.get('name') in data['controls']
                    and record.get('verdict') == 'OBSTRUCTED'):
                previous_bad_controls.add(record['name'])
            (done if record.get('verdict') in CONCLUSIVE else retry).add(record['name'])
    if previous_bad_controls:
        sys.exit('Previously recorded controls are obstructed; refusing to resume: '
                 + ' '.join(sorted(previous_bad_controls)))
    retry -= done
    todo = [n for n in order if n not in done]
    again = len(retry & set(order))
    print(f'{which}: {len(order)} in shard {i}/{total_shards}, {len(done & set(order))} already settled, '
          f'{again} to retry after an error or timeout, {len(todo)} to run, '
          f'{args.workers} workers -> {out}', flush=True)

    write_lock = threading.Lock()

    def work(name):
        started = time.time()
        command = [sys.executable, str(Path(__file__).resolve()), '--targets', str(args.targets),
                   '--one', name] + (['--controls'] if args.controls else [])
        try:
            finished = subprocess.run(command, capture_output=True, text=True, timeout=args.time_per_knot)
            lines = [l for l in finished.stdout.splitlines() if l.startswith('{')]
            record = json.loads(lines[-1]) if lines and finished.returncode == 0 else {
                'name': name, 'verdict': 'ERROR', 'stderr': finished.stderr[-800:],
                'returncode': finished.returncode}
        except subprocess.TimeoutExpired:
            record = {'name': name, 'verdict': 'TIMEOUT', 'seconds': round(time.time() - started, 1)}
        except Exception as exc:                                   # noqa: BLE001 - one knot must not stop the shard
            record = {'name': name, 'verdict': 'ERROR', 'error': repr(exc)[:400]}
        record.setdefault('seconds', round(time.time() - started, 1))
        record.setdefault('role', 'control' if args.controls else 'target')
        # A record may exceed a stream buffer. Serialize whole records instead
        # of assuming append mode makes several buffered writes atomic.
        with write_lock:
            with open(out, 'a') as handle:
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
    bad = [r['name'] for r in results if r.get('control_ok') is False
           or (r['name'] in data['controls'] and r['verdict'] == 'OBSTRUCTED')]
    if bad:
        print('CONTROLS OBSTRUCTED (implementation error): ' + ' '.join(bad), flush=True)
        sys.exit(1)


if __name__ == '__main__':
    main()
