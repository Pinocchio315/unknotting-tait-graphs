#!/usr/bin/env python
"""Summarise enum_search results and project the cost of larger crossing budgets (stdlib only).

    python summarize.py <run_dir> [--budgets 15,16,18,20]

For every knot the results line records, per crossing level, the number of states and the time at
which the level was fully expanded.  The projection for budget B is the time at which all states
with at most B crossings had been expanded (or the knot's total time if the search was cut off
before that), so it is a lower estimate for cut-off knots; the table reports what fraction of the
knots reached each level completely.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import statistics


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('run_dir')
    ap.add_argument('--budgets', default='14,15,16')
    ap.add_argument('--total-targets', type=int, default=224)
    args = ap.parse_args()
    rows = []
    for path in sorted(glob.glob(os.path.join(args.run_dir, 'results_*.jsonl'))):
        for line in open(path):
            if line.strip():
                rows.append(json.loads(line))
    if not rows:
        print('no results')
        return
    print(f'{len(rows)} knots; successes: {sum(1 for r in rows if r["success"])}; exhausted searches: '
          f'{sum(1 for r in rows if r.get("exhausted"))}; median time {statistics.median(r["seconds"] for r in rows):.0f}s; '
          f'total {sum(r["seconds"] for r in rows) / 3600:.1f} CPU-h')
    for r in rows:
        if r['success']:
            print(f"  success {r['name']}: {r['success']} after {r['seconds']}s, states {r['states']}")
    budgets = [int(b) for b in args.budgets.split(',')]
    print(f'\nbudget  knots-complete  median-time  mean-time  projected CPU-h for {args.total_targets} knots')
    for B in budgets:
        times, complete = [], 0
        for r in rows:
            c0 = r['crossings']
            ldt = {int(k): v for k, v in r['level_done_time'].items()}
            done_levels = [lv for lv in ldt if lv <= B]
            if r.get('exhausted') or (r['frontier_min_crossings'] is not None and r['frontier_min_crossings'] > B):
                complete += 1
                times.append(max([ldt[lv] for lv in done_levels], default=r['seconds']) if done_levels else r['seconds'])
            else:
                times.append(r['seconds'])          # cut off: lower estimate
        mean = statistics.mean(times)
        print(f'  {B:2d}      {complete:4d}/{len(rows):<4d}      {statistics.median(times):7.0f}s   {mean:7.0f}s   {mean * args.total_targets / 3600:8.0f}')
    print('\nstates per level (median over knots):')
    levels = sorted({int(k) for r in rows for k in r['states']})
    for lv in levels:
        vals = [r['states'].get(str(lv), r['states'].get(lv, 0)) for r in rows]
        print(f'  {lv:2d} crossings: median {statistics.median(vals):8.0f}  max {max(vals):8d}')


if __name__ == '__main__':
    main()
