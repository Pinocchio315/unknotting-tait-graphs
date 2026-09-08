#!/usr/bin/env python
"""Sweep the Owens u=2 obstruction (|sigma| = 4 only) over KnotInfo alternating knots with unknotting
interval [2,3].  Every OBSTRUCTED knot has u = 3 (lower bound 3 from the obstruction + upper bound 3).

    python run_sweep.py --out results_13a.jsonl --prefix 13a            # all 13a [2,3] |sigma|=4
    python run_sweep.py --out results.jsonl --names 12a_824 ...         # explicit list
    python run_sweep.py --out r.jsonl --prefix 13a --shard 0/16         # slurm shard
Resumable: reruns skip names already in --out."""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _compat  # noqa
import kinfo
from owens_obstruction import obstruct_u2


def one(nm):
    r = kinfo.row(nm)
    t0 = time.time()
    try:
        jones = kinfo.parse_jones(nm) if abs(int(r['signature'])) == 2 else None
        res = obstruct_u2(kinfo.pd_code(nm), int(r['signature']), int(r['determinant']), jones=jones)
    except Exception as e:
        res = {'verdict': 'ERROR', 'reason': repr(e)[:200]}
    res.update({'name': nm, 'interval': r['unknotting_number'].replace(' ', ''),
                'sigma': int(r['signature']), 'det': int(r['determinant']),
                'seconds': round(time.time() - t0, 2)})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--prefix', help="knot name prefix filter, e.g. 13a / 12a / '' for all alternating")
    ap.add_argument('--names', nargs='*')
    ap.add_argument('--interval', default='[2,3]')
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--sigma', default='4', choices=['4', '2', '24'], help='which |sigma| rows to scan')
    args = ap.parse_args()
    if args.names:
        names = args.names
    else:
        names = []
        for nm, r in kinfo.rows().items():
            if r['alternating'] != 'Y':
                continue
            if args.prefix and not nm.startswith(args.prefix):
                continue
            if r['unknotting_number'].replace(' ', '') != args.interval:
                continue
            if str(abs(int(r['signature']))) not in args.sigma:
                continue
            names.append(nm)
    i, N = (int(x) for x in args.shard.split('/'))
    names = [nm for k, nm in enumerate(sorted(names)) if k % N == i]
    done = set()
    if os.path.exists(args.out):
        for ln in open(args.out):
            done.add(json.loads(ln)['name'])
    todo = [nm for nm in names if nm not in done]
    print(f'{len(names)} targets in shard {i}/{N}, {len(todo)} to run', flush=True)
    t0 = time.time()
    obst = 0
    with open(args.out, 'a') as h:
        if args.workers > 1:
            from concurrent.futures import ProcessPoolExecutor, as_completed
            with ProcessPoolExecutor(max_workers=args.workers) as ex:
                for k, fu in enumerate(as_completed({ex.submit(one, nm) for nm in todo}), 1):
                    res = fu.result()
                    h.write(json.dumps(res) + '\n'); h.flush()
                    if res['verdict'] == 'OBSTRUCTED':
                        obst += 1
                        print(f"  ** {res['name']}: OBSTRUCTED -> u = 3", flush=True)
                    if k % 50 == 0:
                        print(f'  {k}/{len(todo)} ({time.time()-t0:.0f}s)', flush=True)
        else:
            for k, nm in enumerate(todo, 1):
                res = one(nm)
                h.write(json.dumps(res) + '\n'); h.flush()
                if res['verdict'] == 'OBSTRUCTED':
                    obst += 1
                    print(f"  ** {res['name']}: OBSTRUCTED -> u = 3", flush=True)
                if k % 50 == 0:
                    print(f'  {k}/{len(todo)} ({time.time()-t0:.0f}s)', flush=True)
    rows = [json.loads(ln) for ln in open(args.out)]
    from collections import Counter
    print('done:', dict(Counter(r['verdict'] for r in rows)), f'{time.time()-t0:.0f}s -> {args.out}')


if __name__ == '__main__':
    main()
