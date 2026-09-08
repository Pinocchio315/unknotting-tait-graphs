#!/usr/bin/env python
"""Second-pass settle hunt over ALL stored campaign candidates (stdlib only; no SnapPy).

Reads the raw results_*.jsonl candidate records (which carry the reduced partner PDs),
de-duplicates them by (target, canonical code), and tries to settle each target by two
certificate-grade mechanisms that need no census identification:

  1. CODE MATCH: greedy-pass-reduce the partner diagram and match its canonical code against
     ALL KnotInfo knots (data/all_knot_codes.json; catches non-hyperbolic partners like
     torus knots that the SnapPy census cannot see).  KnotInfo's upper bound `hi` is
     certified, so  u(target) <= j + hi  --- settles when j + hi <= k.
  2. DIRECT UNKNOTTING: search the partner diagram itself for <= k-j crossing changes that
     greedy-reduce to the unknot: target --(j)--> partner --(<= k-j)--> unknot proves
     u(target) <= k with no identification at all (works for connected sums).

Already-settled targets (SUCCESS/CHAIN certificates, or --skip names) are skipped.

    ~/.pyenv/versions/unknot-venv/bin/python resolve_unidentified.py xup_runs
"""
from __future__ import annotations

import argparse
import glob
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd, to_pd  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.reduce import greedy_reduce, search_reduce, ReduceConfig  # noqa: E402


def _monotone(g):
    cur = g
    while True:
        rem = mv.find_r1_removals(cur)
        if rem:
            cur = mv.r1_remove(cur, rem[0][1])
            continue
        rem2 = mv.find_r2_removals(cur)
        if rem2:
            kind, params = rem2[0]
            cur = mv.apply_move(cur, 'r2-' + kind[0], params)
            continue
        return cur


def direct_unknotting(g, max_changes, greedy_gate=9):
    """Smallest crossing-change subset (PD rows of to_pd(g)) whose result reduces to the
    unknot (monotone R1-/R2- first; the pass-greedy runs only when that gets close)."""
    edges = sorted(g.edges)
    for size in range(1, max_changes + 1):
        for rows in itertools.combinations(range(len(edges)), size):
            h = g
            for i in rows:
                h = mv.crossing_change(h, edges[i])
            red = _monotone(h)
            if red.n_crossings() == 0:
                return list(rows)
            if red.n_crossings() <= greedy_gate:
                red2, _ = greedy_reduce(red)
                if red2.n_crossings() == 0:
                    return list(rows)
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--out', default=None)
    ap.add_argument('--skip', nargs='*', default=[], help='targets to skip (already settled)')
    args = ap.parse_args()
    out_path = args.out or os.path.join(args.run_dir, 'unidentified_resolved.jsonl')

    allk = json.load(open(os.path.join(HERE, 'data', 'all_knot_codes.json')))
    code_info = {}
    for x in allk:
        for c in x['codes']:
            code_info[c] = x

    done = {os.path.basename(p)[len(prefix):-5].replace('SUCCESS_', '').replace('CHAIN_', '')
            for prefix in ('SUCCESS_', 'CHAIN_')
            for p in glob.glob(os.path.join(args.run_dir, prefix + '*.json'))}
    done |= set(args.skip)

    jobs = {}
    n_rows = 0
    for path in sorted(glob.glob(os.path.join(args.run_dir, 'results_*.jsonl'))):
        for line in open(path):
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec['name'] in done:
                continue
            for c in rec.get('candidates', []):
                n_rows += 1
                code = repr(canonical_code(from_pd([list(q) for q in c['reduced_pd']])))
                key = (rec['name'], code)
                if key not in jobs or c['j'] < jobs[key][1]['j']:
                    jobs[key] = (rec, c)
    print(f'{n_rows} candidate rows -> {len(jobs)} distinct (target, partner-diagram) pairs',
          flush=True)

    settles, matched, still = [], 0, 0
    t0 = time.time()
    settled_targets = set()
    cache = {}      # canonical code of reduced_pd -> [info, n_greedy, red_pd, hunts{}]
    with open(out_path, 'w') as out:
        for n_done, ((target, code0), (rec, c)) in enumerate(sorted(jobs.items()), 1):
            if target in settled_targets:
                continue
            k, j = rec['k'], c['j']
            hit = cache.get(code0)
            if hit is None:
                g = from_pd([list(q) for q in c['reduced_pd']])
                red, _ = greedy_reduce(g)
                info = code_info.get(repr(canonical_code(red)))
                hit = [info, red.n_crossings(), to_pd(red), {}]
                cache[code0] = hit
            info, n_greedy, red_pd, hunts = hit
            red = from_pd([list(q) for q in red_pd])
            res = {'target': target, 'k': k, 'j': j, 'rows': c['rows'],
                   'n_reduced': c['n_reduced'], 'n_greedy': n_greedy}
            if info is not None:
                matched += 1
                res.update(partner=info['name'], partner_interval=[info['lo'], info['hi']])
                if info['hi'] is not None and j + info['hi'] <= k:
                    res['claim'] = f"u({target}) <= {j + info['hi']} via {info['name']}"
                    res['mechanism'] = 'code-match'
            if 'claim' not in res:
                if (k - j) not in hunts:
                    hunts[k - j] = direct_unknotting(red, k - j)
                changes = hunts[k - j]
                if changes is not None:
                    res.update(mechanism='direct', partner_pd=to_pd(red),
                               partner_changes=changes,
                               claim=f"u({target}) <= {j + len(changes)}: partner unknotted "
                                     f"by {len(changes)} change(s)")
                elif info is None:
                    still += 1
            if 'claim' in res:
                settles.append(res)
                settled_targets.add(target)
            out.write(json.dumps(res) + '\n')
            if n_done % 2000 == 0:
                print(f'  {n_done}/{len(jobs)} (matched {matched}, settles {len(settles)}, '
                      f'unresolved {still}) {time.time()-t0:.0f}s', flush=True)

    print(f'\ncode-matched {matched}; NEW SETTLES: {len(settles)}; unresolved: {still}')
    for r in settles:
        print('  ***', r['claim'])


if __name__ == '__main__':
    main()
