#!/usr/bin/env python
"""Upper-bound campaign: expanded-diagram crossing changes + pass-move reduction (stdlib only).

For every target knot K with updated range [k, k+1] (data/targets.json):
  1. expand the KnotInfo diagram to ~N inequivalent diagrams, with the campaign crossing
     cutoff of 50, by R2+/R3/flype random walks (NO Reidemeister I) --- expand.py;
  2. on every diagram, screen ALL j-subsets of crossings, j = 1..k, by the determinant of the
     changed knot, computed in O(1) per subset by FlipDet (matrix determinant lemma mod 2^61-1):
       j = k : det must be 1              (candidate unknot        -> u(K) <= k, settled)
       j < k : det must match a knot with known u <= k-j            (u(K) <= j + u' <= k, settled)
  3. det-survivors are actually changed and reduced with the pass-move engine of the
     crossing-additivity campaign (greedy pass reduction, then a bounded best-first search);
     reaching 0 crossings proves the unknot; reductions to <= 13 crossings are identified
     against the canonical codes of all knots with known u <= 2 (data/known_u.json).
  4. every success records the expansion path from the KnotInfo PD and the changed PD rows;
     verify locally with verify_certificates.py.  When a replayable partner reduction is
     needed, rehunt_certificates.py re-derives it as a full CHAIN certificate.

    python3 run_search.py --shard 0/16 --out-dir runs                  # server shard
    python3 run_search.py --names 13n1587 --diagrams 500 --seed 1      # pilot
Resume-safe: finished knots (a line in results_<shard>.jsonl) are skipped on rerun.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd, to_pd  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.reduce import greedy_reduce, search_reduce, ReduceConfig  # noqa: E402
from flipdet import FlipDet, det_pairs, P  # noqa: E402
from expand import expand_diagrams, jsonable_path  # noqa: E402


def load_data(targets_path=None):
    # An explicit frozen list permits replaying a later cohort without replacing
    # the historical inputs. Diagram generation and subset selection are unchanged.
    targets = json.load(open(targets_path or os.path.join(HERE, 'data', 'targets.json')))
    known = json.load(open(os.path.join(HERE, 'data', 'known_u.json')))
    return targets, known


_CTX = {}


def _ctx():
    """Worker-local lazy load of the known-u tables (avoids pickling them per job)."""
    if not _CTX:
        known = json.load(open(os.path.join(HERE, 'data', 'known_u.json')))
        code_u = {}
        for x in known:
            for c in x['codes']:
                code_u[c] = (x['name'], x['u'])
        dets_le = {1: det_pairs(x['det'] for x in known if x['u'] <= 1),
                   2: det_pairs(x['det'] for x in known if x['u'] <= 2)}
        _CTX.update(code_u=code_u, dets_le=dets_le)
    return _CTX


def monotone_reduce(g):
    """R1-/R2- removals only --- the cheap triage pass (no pass-move search)."""
    cur, path = g, []
    while True:
        rem = mv.find_r1_removals(cur)
        if rem:
            kind, e = rem[0]
            cur = mv.r1_remove(cur, e)
            path.append(('r1-', (e,)))
            continue
        rem2 = mv.find_r2_removals(cur)
        if rem2:
            kind, params = rem2[0]
            cur = mv.apply_move(cur, 'r2-' + kind[0], params)
            path.append(('r2-' + kind[0], params))
            continue
        return cur, path


def escalate_reduce(red, target_crossings, nodes, seconds, seed):
    """Bounded best-first pass search on an already-greedy-reduced diagram."""
    cfg = ReduceConfig(target=target_crossings, max_extra=2, max_nodes=nodes,
                       time_limit=seconds, seed=seed, heap_max=50000, log_every=10 ** 9)
    best, more, _stats = search_reduce(red, cfg)
    return best, more


def scan_knot(job):
    t, opt = job
    t0 = time.time()
    ctx = _ctx()
    k = t['k']
    g0 = from_pd([list(q) for q in t['pd']], t['name'])
    # determinant lookup per change count j (mod-P residue pairs) and code tables
    allowed = {k: det_pairs([1])}
    code_u = ctx['code_u']
    for j in range(1, k):
        allowed[j] = ctx['dets_le'][k - j]
    stats = {'name': t['name'], 'k': k, 'range': t['range'], 'diagrams': 0, 'subsets': 0,
             'det_hits': 0, 'reduced': 0, 'escalated': 0, 'near_miss': 0,
             'candidates': [], 'unknot_near': [], 'success': None}
    esc_left = [opt['escalations']]
    near_seen = set()

    def try_full(g, F, rows):
        """j = k route: changed diagram must be the unknot.  Returns success record or None."""
        h = g
        for i in rows:
            h = mv.crossing_change(h, F.edge_ids[i])
        red, _ = monotone_reduce(h)
        stats['reduced'] += 1
        if red.n_crossings() > opt['greedy_gate']:
            stats['near_miss'] += 1
            return None
        red, _ = greedy_reduce(red)
        if red.n_crossings() > 0 and red.n_crossings() <= opt['escalate_at'] and esc_left[0] > 0:
            esc_left[0] -= 1
            stats['escalated'] += 1
            red, _ = escalate_reduce(red, 0, opt['nodes'], opt['seconds'], seed=len(rows))
        if red.n_crossings() == 0:
            return {'kind': 'unknot', 'u': k, 'j': k}
        stats['near_miss'] += 1
        # det = 1 but not reduced to 0: store small survivors for the offline
        # HFK genus-0 backstop (a true unknot here is proof-grade even without a move path)
        if red.n_crossings() <= opt['near_store'] and len(stats['unknot_near']) < 30:
            code = repr(canonical_code(red))
            if code not in near_seen:
                near_seen.add(code)
                stats['unknot_near'].append({'rows': list(rows),
                                             'n_reduced': red.n_crossings(),
                                             'reduced_pd': to_pd(red)})
        return None

    def try_partial(g, F, rows, j):
        """j < k route: changed diagram must reduce to a knot with known u <= k-j."""
        h = g
        for i in rows:
            h = mv.crossing_change(h, F.edge_ids[i])
        red, _ = monotone_reduce(h)
        stats['reduced'] += 1
        if red.n_crossings() > 16:
            stats['near_miss'] += 1
            return None
        if red.n_crossings() > 13:
            red2, _ = greedy_reduce(red)
            red = red2 if red2.n_crossings() < red.n_crossings() else red
        n = red.n_crossings()
        if not 0 < n <= 13:
            stats['near_miss'] += 1
            return None
        code = repr(canonical_code(red))
        hit = code_u.get(code)
        if hit and hit[1] <= k - j:
            return {'kind': 'known', 'u': j + hit[1], 'j': j,
                    'partner': hit[0], 'partner_u': hit[1], 'reduced_pd': to_pd(red)}
        stats['candidates'].append({'j': j, 'rows': list(rows), 'n_reduced': n,
                                    'reduced_pd': to_pd(red)})
        return None

    success = None
    budget_end = t0 + opt['time_per_knot']
    # source diagrams: the KnotInfo diagram itself, then the expanded ones
    def sources():
        yield g0, []
        yield from expand_diagrams(g0, n_target=opt['n_target'], count=opt['diagrams'],
                                   seed=opt['seed'] + sum(ord(c) for c in t['name']))

    for g, path in sources():
        if success or time.time() > budget_end:
            break
        stats['diagrams'] += 1
        try:
            F = FlipDet(g)
        except ZeroDivisionError:
            continue
        m = len(F.edge_ids)
        rnd = random.Random(stats['diagrams'])
        for j in ((k,) if opt.get('no_partial') else range(1, k + 1)):
            if success:
                break
            pairs = allowed[j]
            # det screen first (cheap), then a sampled, budgeted reduction pass
            hits = []
            for rows in itertools.combinations(range(m), j):
                stats['subsets'] += 1
                if F.det_after(rows) in pairs:
                    hits.append(rows)
            stats['det_hits'] += len(hits)
            cap = opt['partial_per_diagram'] if j < k else opt['full_per_diagram']
            if cap and len(hits) > cap:
                hits = rnd.sample(hits, cap)
            for rows in hits:
                rec = try_partial(g, F, rows, j) if j < k else try_full(g, F, rows)
                if rec is not None:
                    rec.update({'name': t['name'], 'range': t['range'],
                                'start_pd': t['pd'], 'expansion_path': jsonable_path(path),
                                'witness_pd': to_pd(g), 'change_rows': list(rows),
                                'n_witness': m,
                                'det_start': t['det']})
                    success = rec
                    break
                if time.time() > budget_end:
                    break
    stats['seconds'] = round(time.time() - t0, 1)
    stats['success'] = success
    stats['candidates'] = stats['candidates'][:20]
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--targets', help='frozen target JSON; defaults to data/targets.json')
    ap.add_argument('--shard', default='0/1', help='i/n split of the target list')
    ap.add_argument('--control', action='store_true',
                    help='scan the five SETTLED minimal-witness knots instead of the targets; '
                         'each must yield a SUCCESS (u <= 2) or the pipeline is broken')
    ap.add_argument('--k-class', default='all', choices=['1', '2', '3', 'all'])
    ap.add_argument('--names', nargs='*', help='explicit target names (KnotInfo style, e.g. 13n1587)')
    ap.add_argument('--diagrams', type=int, default=10000)
    ap.add_argument('--n-target', type=int, default=50,
                    help='campaign crossing cutoff (paper: 50); R1 is never used, so the '
                         'effective size is the largest value <= this with the start parity')
    ap.add_argument('--time-per-knot', type=float, default=1800.0)
    ap.add_argument('--nodes', type=int, default=1500, help='best-first nodes per escalation')
    ap.add_argument('--seconds', type=float, default=4.0, help='best-first seconds per escalation')
    ap.add_argument('--escalate-at', type=int, default=8,
                    help='run the best-first search when greedy ends at <= this many crossings')
    ap.add_argument('--escalations', type=int, default=150, help='best-first budget per knot')
    ap.add_argument('--partial-per-diagram', type=int, default=6,
                    help='max j<k det survivors tried per diagram (sampled; 0 = no cap)')
    ap.add_argument('--full-per-diagram', type=int, default=15,
                    help='max j=k det-1 survivors tried per diagram (sampled; 0 = no cap)')
    ap.add_argument('--no-partial', action='store_true',
                    help='j = k unknot route only (skip the partial j<k routes)')
    ap.add_argument('--near-store', type=int, default=16,
                    help='store j=k det-1 near-misses reduced to <= this many crossings '
                         '(offline HFK backstop; 30 per knot, deduplicated)')
    ap.add_argument('--greedy-gate', type=int, default=24,
                    help='run the greedy pass reducer only when the monotone pass lands at <= this')
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out-dir', default=os.environ.get('XUP_OUT', 'xup_runs'))
    args = ap.parse_args()

    targets, known = load_data(args.targets)
    if args.control:
        targets = json.load(open(os.path.join(HERE, 'data', 'controls.json')))
    if args.k_class != 'all':
        targets = [t for t in targets if t['k'] == int(args.k_class)]
    if args.names:
        wanted = {n.replace('_', '') for n in args.names}
        targets = [t for t in targets if t['name'].replace('_', '') in wanted]
    i, n = map(int, args.shard.split('/'))
    targets = targets[i::n]

    opt = {'diagrams': args.diagrams,
           'n_target': args.n_target, 'time_per_knot': args.time_per_knot,
           'nodes': args.nodes, 'seconds': args.seconds, 'seed': args.seed,
           'escalate_at': args.escalate_at, 'escalations': args.escalations,
           'partial_per_diagram': args.partial_per_diagram,
           'no_partial': args.no_partial, 'near_store': args.near_store,
           'full_per_diagram': args.full_per_diagram, 'greedy_gate': args.greedy_gate}

    os.makedirs(args.out_dir, exist_ok=True)
    res_path = os.path.join(args.out_dir, f'results_{i}_{n}.jsonl')
    done = set()
    if os.path.exists(res_path):
        with open(res_path) as h:
            for line in h:
                try:
                    done.add(json.loads(line)['name'])
                except Exception:
                    pass
    todo = [t for t in targets if t['name'] not in done]
    print(f'shard {i}/{n}: {len(todo)} targets to scan ({len(done)} already done)', flush=True)

    def record(stats):
        with open(res_path, 'a') as h:
            h.write(json.dumps(stats) + '\n')
        if stats['success']:
            sp = os.path.join(args.out_dir, f"SUCCESS_{stats['name']}.json")
            with open(sp, 'w') as h:
                json.dump(stats['success'], h)
            print(f"*** SUCCESS {stats['name']}  u <= {stats['success']['u']} "
                  f"({stats['success']['kind']}, j={stats['success']['j']}) -> {sp}", flush=True)
        else:
            print(f"    {stats['name']}: diagrams={stats['diagrams']} subsets={stats['subsets']} "
                  f"det_hits={stats['det_hits']} esc={stats['escalated']} "
                  f"near_miss={stats['near_miss']} "
                  f"cand={len(stats['candidates'])} {stats['seconds']}s", flush=True)

    jobs = [(t, opt) for t in todo]
    if args.workers <= 1:
        for job in jobs:
            record(scan_knot(job))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(scan_knot, job) for job in jobs]
            for f in as_completed(futs):
                record(f.result())
    print('shard finished.', flush=True)


if __name__ == '__main__':
    main()
