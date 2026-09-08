#!/usr/bin/env python
"""Level-by-level enumeration search for upper bounds on unknotting numbers (stdlib only).

For each target knot K with updated range [k, k+1] (normally k = 3) the KnotInfo diagram is expanded
into the set of KINK-FREE diagrams reachable by pass moves (gain >= -PASS_EXTRA), Reidemeister III
moves and flypes (on both checkerboard graphs) without exceeding MAX_CROSSINGS crossings.  States
are expanded in order of increasing crossing number (so small diagrams are exhausted first),
deduplicated by the canonical code of the signed embedded Tait graph, and screened as they appear:

  route 1  change ONE crossing whose determinant matches a knot of known u <= k-1, reduce, and
           identify the result among the minimal diagrams of the known knots  ->  u(K) <= 1 + u
  route 2  k = 2: change two crossings with determinant 1, reduce to the trivial diagram (u <= 2)
           k = 3: change two crossings whose determinant matches a knot of known u = 1 (u <= 3)
  route 3  k = 3 only, optional: change three crossings with determinant 1 (u <= 3)

Every success is written as a replayable certificate (moves from the KnotInfo PD, rows to change,
claimed partner); verify locally with verify_enum_certificates.py.  Per-knot statistics (states
and time per crossing level, whether the search was exhausted) go to results_<i>_<n>.jsonl and
feed summarize.py, which projects the cost of larger crossing budgets.  Resume-safe.

    python enum_search.py --k-class 3 --extra-crossings 3 --workers 16 --wall-time 171900
    python enum_search.py --control                      # the nine settled knots must be re-found
"""
from __future__ import annotations

import argparse
import collections
import glob
import heapq
import itertools
import json
import multiprocessing
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from xtait.graph import from_pd, to_pd, dual  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.flypes import find_flypes_all, flype_block  # noqa: E402
from xtait.reduce import greedy_reduce, search_reduce, ReduceConfig  # noqa: E402
from xtait.passmove import iter_graph_pass_candidates  # noqa: E402
from flipdet import FlipDet, det_pairs, P  # noqa: E402
from moves_io import apply_record  # noqa: E402
from lickorish import LickorishFilter  # noqa: E402
from alexander import fingerprint, UNKNOT_FP  # noqa: E402


def det_value(residue):
    return residue if residue < P // 2 else P - residue


def code(g):
    return repr(canonical_code(g))


def turn_over(g):
    h = g.copy()
    for v in h.vertices:
        h.rotation[v] = list(reversed(h.rotation[v]))
    return h


def mirror(g):
    h = g.copy()
    for e in h.signs:
        h.signs[e] = -h.signs[e]
    return h


def variant_codes(g):
    """codes of the diagram, seen from the other side, and their mirrors (KnotInfo names are
    chiral-blind; the mirror is recorded in the certificate as 'mirror' when it matches)"""
    return [(code(g), False), (code(turn_over(g)), False), (code(mirror(g)), True), (code(mirror(turn_over(g))), True)]


# ----------------------------------------------------------------------------- neighbours
def flype_neighbours(g):
    """(record, graph) for all flypes visible in either checkerboard graph."""
    out = []
    for is_dual, base in ((False, g), (True, dual(g))):
        for e, Pc in find_flypes_all(base):
            try:
                h = flype_block(base, e, Pc)
            except ValueError:
                continue
            out.append((['flypeB', {'dual': is_dual, 'e': e, 'edges': sorted(Pc['edges'])}], h))
    return out


def r2_neighbours(g, max_crossings):
    """(record, graph) for every Reidemeister II addition (parallel and series pairs in the
    Tait graph, both strand orders) that stays within the crossing budget."""
    if g.n_crossings() + 2 > max_crossings:
        return
    for fi, i, j in mv.find_r2_parallel_additions(g):
        for sign in (1, -1):
            try:
                h = mv.r2_add_parallel(g, fi, i, j, sign)
            except ValueError:
                continue
            if not mv.find_r1_removals(h):
                yield ['r2P', [fi, i, j, sign]], h
    for v, i, j in mv.find_r2_series_additions(g):
        for sign in (1, -1):
            try:
                h = mv.r2_add_series(g, v, i, j, sign)
            except ValueError:
                continue
            if not mv.find_r1_removals(h):
                yield ['r2S', [v, i, j, sign]], h


def neighbours(g, max_crossings, pass_extra, max_routes, include_flypes=True,
               should_stop=None, include_r2=False):
    """(record, graph) for every kink-free neighbour within the crossing budget."""
    for fi in mv.find_r3_triangles(g):
        yield ['r3D', fi], mv.r3_delta_to_y(g, fi)
    for t in mv.find_r3_stars(g):
        yield ['r3Y', t], mv.r3_y_to_delta(g, t)
    if include_flypes:
        yield from flype_neighbours(g)
    if include_r2:
        yield from r2_neighbours(g, max_crossings)
    for prm in iter_graph_pass_candidates(g, max_extra=pass_extra,
                                          max_routes=max_routes,
                                          should_stop=should_stop):
        try:
            h = mv.apply_move(g, 'pass', (prm,))
        except Exception:
            continue
        if h.n_crossings() <= max_crossings and not mv.find_r1_removals(h):
            yield ['pass', prm], h


# ----------------------------------------------------------------------------- identification
class Known:
    """known knots with u <= 2: determinant residue sets per u, and canonical codes of the
    minimal diagrams (flype orbits of the KnotInfo diagrams) -> (name, u)"""

    def __init__(self, path):
        d = json.load(open(path))
        self.u_of = {}
        self.codes = {}
        dets = {1: set(), 2: set()}
        for x in d['knots']:
            self.u_of[x['name']] = x['u']
            dets[x['u']].add(x['det'])
            for c in x['codes']:
                self.codes[c] = (x['name'], x['u'])
        self.det_pairs = {1: det_pairs(dets[1]), 2: det_pairs(dets[1] | dets[2])}   # u <= 1, u <= 2
        self.by_fp = {}                                    # (det, Alexander fingerprint) -> [(name, u)]
        for x in d['knots']:
            if 'alex' in x:
                self.by_fp.setdefault((x['det'], tuple(x['alex'])), []).append((x['name'], x['u']))
        # how many known knots share a determinant residue (rarer residues are screened first)
        self.det_count = {1: collections.Counter(), 2: collections.Counter()}
        for x in d['knots']:
            for r in det_pairs([x['det']]):
                self.det_count[1][r] += (x['u'] == 1)
                self.det_count[2][r] += 1

    def identify(self, g):
        for c, mirrored in variant_codes(g):
            if c in self.codes:
                nm, u = self.codes[c]
                return nm, u, mirrored
        return None


# ----------------------------------------------------------------------------- search
def reduce_small(g, deep=False):
    """greedy pass/Reidemeister reduction; a small bounded search only for small leftovers"""
    red, path = greedy_reduce(g)
    n = red.n_crossings()
    if deep and 0 < n <= 10:
        try:
            r2, p2, _ = search_reduce(red, ReduceConfig(target=0, max_extra=1, pass_extra=1, inflate_samples=0,
                                                       use_r3=True, use_flype=False, max_nodes=80,
                                                       time_limit=0.06, seed=0))
            if r2.n_crossings() < n:
                red, path = r2, path + p2
        except Exception:
            pass
    return red, path


def search_knot(t, known, args, log):
    t0 = time.monotonic()
    k = t['k']
    g0 = from_pd([list(q) for q in t['pd']])
    c0 = g0.n_crossings()
    budget = args.max_crossings if args.extra_crossings is None else min(args.max_crossings, c0 + args.extra_crossings)
    routes = [int(r) for r in args.routes.split(',')]
    partner_u_max = k - 1
    dp_unknot = det_pairs([1])
    dp1 = known.det_pairs[partner_u_max] if k >= 2 else dp_unknot   # route 1: partner with u <= k-1; k = 1: the unknot
    dp_u1 = known.det_pairs[1]
    stats = {'name': t['name'], 'k': k, 'crossings': c0, 'max_crossings': budget,
             'routes': routes, 'states': collections.Counter(), 'expanded': collections.Counter(),
             'level_done_time': {}, 'screened': collections.Counter(), 'det_hits': collections.Counter(), 'lick_pass': collections.Counter(),
             'candidates': [], 'success': None, 'hits_skipped': 0,
             'skipped_flype_class': 0, 'flype_closure_states': 0}
    codes = {code(g0): 0}
    graphs = [g0]
    parent = [(-1, None)]
    uf = [0]                        # union-find over states: flype moves join classes
    uf_screened = [False]           # whether a member of the class has been screened
    flype_expanded = [False]        # eager flype closure for alternating targets
    alternating = str(t.get('alternating', '')).upper() == 'Y'

    def find(i):
        while uf[i] != i:
            uf[i] = uf[uf[i]]
            i = uf[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            root, child = min(ri, rj), max(ri, rj)
            uf[child] = root
            uf_screened[root] = uf_screened[ri] or uf_screened[rj]
    heap = [(c0, 0, 0)]
    seq = itertools.count(1)
    stats['states'][c0] += 1
    cur_level = c0

    def time_limit_reason():
        now = time.monotonic()
        job_deadline = getattr(args, 'job_deadline', 0.0)
        if job_deadline and now >= job_deadline:
            return 'job_time'
        if args.time_per_knot > 0 and now - t0 >= args.time_per_knot:
            return 'knot_time'
        return None

    class SearchInterrupted(Exception):
        pass

    def check_time():
        if time_limit_reason() is not None:
            raise SearchInterrupted

    def add_state(source, rec, h):
        """Insert one state and preserve a replayable path; return (index, is_new)."""
        c = code(h)
        if c in codes:
            j = codes[c]
            if rec[0] == 'flypeB':
                union(source, j)
            return j, False
        if args.max_states > 0 and len(graphs) >= args.max_states:
            return None, False
        j = len(graphs)
        codes[c] = j
        graphs.append(h)
        parent.append((source, rec))
        uf.append(j)
        uf_screened.append(False)
        flype_expanded.append(False)
        if rec[0] == 'flypeB':
            union(source, j)
        stats['states'][h.n_crossings()] += 1
        heapq.heappush(heap, (h.n_crossings(), next(seq), j))
        return j, True

    def close_flype_class(seed):
        """Eagerly close one flype class, used for alternating target knots.

        Every member remains in the search because it may admit different pass
        moves, but the completed union-find class is screened only once.
        """
        todo = collections.deque([seed])
        queued = {seed}
        while todo:
            reason = time_limit_reason()
            if reason is not None:
                return reason
            i = todo.popleft()
            if flype_expanded[i]:
                continue
            flype_expanded[i] = True
            for rec, h in flype_neighbours(graphs[i]):
                reason = time_limit_reason()
                if reason is not None:
                    return reason
                j, is_new = add_state(i, rec, h)
                if j is None:
                    return 'state_limit'
                if is_new:
                    stats['flype_closure_states'] += 1
                if j not in queued and not flype_expanded[j]:
                    queued.add(j)
                    todo.append(j)
        return None

    def moves_to(idx):
        recs = []
        while idx > 0:
            p, rec = parent[idx]
            recs.append(rec)
            idx = p
        recs.reverse()
        return recs

    def certificate(idx, kind, rows, j, partner=None, partner_u=None, mirrored=False, reduced=None):
        wpd = to_pd(graphs[idx])
        return {'name': t['name'], 'kind': kind, 'j': j, 'u': (j if kind == 'unknot' else j + partner_u),
                'partner': partner, 'partner_u': partner_u, 'partner_mirror': mirrored,
                'start_pd': t['pd'], 'moves': moves_to(idx), 'witness_pd': wpd, 'n_witness': len(wpd),
                'change_rows': list(rows), 'reduced_pd': reduced, 'identification': 'canonical code' if kind == 'known' else 'reduced to the trivial diagram'}

    cand_seen = set()

    def add_candidate(idx, rows, j, red, D, kind, partners):
        """a (det, Alexander) match that could not be identified exactly: keep for local resolution"""
        key = code(red)
        if key in cand_seen:
            return
        cand_seen.add(key)
        rec = {'rows': list(rows), 'j': j, 'n_state': graphs[idx].n_crossings(), 'n_reduced': red.n_crossings(),
               'det': D, 'kind': kind, 'partners': partners, 'reduced_pd': to_pd(red), 'moves': moves_to(idx),
               'witness_pd': to_pd(graphs[idx])}
        # keep the most promising: few partner candidates, small reduced diagram
        stats['candidates'].append(rec)
        stats['candidates'].sort(key=lambda r: (len(r['partners']), r['n_reduced']))
        del stats['candidates'][args.cand_store:]

    def screen(idx, g):
        check_time()
        ts = time.time()
        n = g.n_crossings()
        fd = FlipDet(g)
        edges = sorted(g.edges)
        lick = None
        if (1 in routes and partner_u_max == 1) or (2 in routes and k == 3):
            lick = LickorishFilter(g)
        # route 1: one change to a known partner (rarest determinants first;
        # optionally capped for pilot runs)
        if 1 in routes and k == 1:
            for i in range(n):
                check_time()
                stats['screened'][1] += 1
                if fd.det_after([i]) not in dp_unknot:
                    continue
                stats['det_hits'][1] += 1
                x = mv.crossing_change(g, edges[i])
                red, _ = reduce_small(x, deep=True)
                if red.n_crossings() == 0:
                    return certificate(idx, 'unknot', [i], 1), None
                if red.n_crossings() <= args.cand_max_crossings and fingerprint(to_pd(red)) == UNKNOT_FP:
                    add_candidate(idx, [i], 1, red, 1, 'unknot', ['unknot'])
        if 1 in routes and k >= 2:
            hits = []
            for i in range(n):
                check_time()
                stats['screened'][1] += 1
                d = fd.det_after([i])
                if d in dp1:
                    stats['det_hits'][1] += 1
                    if partner_u_max == 1 and not lick.passes([i], det_value(d)):
                        continue
                    hits.append((known.det_count[partner_u_max][d], i))
            stats['lick_pass'][1] += len(hits)
            hits.sort()
            if args.hits_per_state > 0 and len(hits) > args.hits_per_state:
                stats['hits_skipped'] += len(hits) - args.hits_per_state
                hits = hits[:args.hits_per_state]
            for _, i in hits:
                check_time()
                x = mv.crossing_change(g, edges[i])
                red, _ = reduce_small(x)
                D = det_value(fd.det_after([i]))
                if red.n_crossings() <= 13:
                    hit = known.identify(red)
                    if hit is None:
                        red2, _ = reduce_small(x, deep=True)
                        hit = known.identify(red2)
                        if hit:
                            red = red2
                    if hit and hit[1] <= partner_u_max:
                        return certificate(idx, 'known', [i], 1, hit[0], hit[1], hit[2], to_pd(red)), None
                if red.n_crossings() <= args.cand_max_crossings:
                    partners = [nm for nm, u in known.by_fp.get((D, fingerprint(to_pd(red))), []) if u <= partner_u_max]
                    if partners:
                        add_candidate(idx, [i], 1, red, D, 'known', partners)
        # route 2
        if 2 in routes and k >= 2:
            pair_hits = []
            for a, b in itertools.combinations(range(n), 2):
                check_time()
                stats['screened'][2] += 1
                d = fd.det_after([a, b])
                if k == 2:
                    if d not in dp_unknot:
                        continue
                    stats['det_hits'][2] += 1
                    x = mv.crossing_change(mv.crossing_change(g, edges[a]), edges[b])
                    red, _ = reduce_small(x, deep=True)
                    if red.n_crossings() == 0:
                        return certificate(idx, 'unknot', [a, b], 2), None
                    if red.n_crossings() <= args.cand_max_crossings and fingerprint(to_pd(red)) == UNKNOT_FP:
                        add_candidate(idx, [a, b], 2, red, 1, 'unknot', ['unknot'])
                else:
                    if d not in dp_u1:
                        continue
                    stats['det_hits'][2] += 1
                    if not lick.passes([a, b], det_value(d)):
                        continue
                    stats['lick_pass'][2] += 1
                    pair_hits.append((known.det_count[1][d], a, b))
            if k == 3:
                pair_hits.sort()
                if args.hits_per_state > 0 and len(pair_hits) > args.hits_per_state:
                    stats['hits_skipped'] += len(pair_hits) - args.hits_per_state
                    pair_hits = pair_hits[:args.hits_per_state]
                for _, a, b in pair_hits:
                    check_time()
                    x = mv.crossing_change(mv.crossing_change(g, edges[a]), edges[b])
                    red, _ = reduce_small(x)
                    D = det_value(fd.det_after([a, b]))
                    if red.n_crossings() <= 13:
                        hit = known.identify(red)
                        if hit is None:
                            red2, _ = reduce_small(x, deep=True)
                            hit = known.identify(red2)
                            if hit:
                                red = red2
                        if hit and hit[1] == 1:
                            return certificate(idx, 'known', [a, b], 2, hit[0], 1, hit[2], to_pd(red)), None
                    if red.n_crossings() <= args.cand_max_crossings:
                        partners = [nm for nm, u in known.by_fp.get((D, fingerprint(to_pd(red))), []) if u == 1]
                        if partners:
                            add_candidate(idx, [a, b], 2, red, D, 'known', partners)
        # route 3 (k = 3 only): three changes to the unknot
        if 3 in routes and k == 3:
            for sub in itertools.combinations(range(n), 3):
                check_time()
                stats['screened'][3] += 1
                if fd.det_after(list(sub)) not in dp_unknot:
                    continue
                stats['det_hits'][3] += 1
                x = g
                for i in sub:
                    x = mv.crossing_change(x, edges[i])
                red, _ = reduce_small(x, deep=True)
                if red.n_crossings() == 0:
                    return certificate(idx, 'unknot', list(sub), 3), None
                if red.n_crossings() <= args.cand_max_crossings and fingerprint(to_pd(red)) == UNKNOT_FP:
                    add_candidate(idx, list(sub), 3, red, 1, 'unknot', ['unknot'])
        return None, None

    # screen the start diagram, then expand in order of crossing number
    t_screen = [0.0]
    def screen_timed(idx, g):
        ts = time.monotonic()
        try:
            return screen(idx, g)
        finally:
            t_screen[0] += time.monotonic() - ts
    cert = None
    exhausted = False
    stop_reason = None
    incomplete_level = None
    while heap and cert is None:
        stop_reason = time_limit_reason()
        if stop_reason is not None:
            break
        n, _, idx = heapq.heappop(heap)
        if n > cur_level:
            stats['level_done_time'][cur_level] = round(time.monotonic() - t0, 1)
            log(f"  {t['name']}: level {cur_level} complete after {stats['level_done_time'][cur_level]}s; "
                f"states {dict(sorted(stats['states'].items()))}; candidates {len(stats['candidates'])}")
            cur_level = n
        g = graphs[idx]
        stats['expanded'][n] += 1
        n_expanded_total = sum(stats['expanded'].values())
        if n_expanded_total % 1000 == 0:
            log(f"  {t['name']}: {n_expanded_total} states expanded, {len(graphs)} generated, "
                f"level {n}, {round(time.monotonic() - t0)}s (screening {round(t_screen[0])}s)")
        # For alternating targets, discover the complete flype class before
        # screening.  This makes the flype/crossing-change commutation an actual
        # computational saving rather than an order-dependent heuristic.
        if alternating:
            stop_reason = close_flype_class(idx)
            if stop_reason is not None:
                incomplete_level = n
                break
        r = find(idx)
        if uf_screened[r]:
            stats['skipped_flype_class'] += 1
        else:
            try:
                cert, _ = screen_timed(idx, g)
            except SearchInterrupted:
                stop_reason = time_limit_reason() or 'knot_time'
                incomplete_level = n
                break
            uf_screened[find(idx)] = True
            if cert is not None:
                break
        for rec, h in neighbours(g, budget, args.pass_extra, args.max_routes,
                                 include_flypes=not alternating,
                                 should_stop=lambda: time_limit_reason() is not None,
                                 include_r2=getattr(args, 'r2', False)):
            stop_reason = time_limit_reason()
            if stop_reason is not None:
                incomplete_level = n
                break
            j, _ = add_state(idx, rec, h)
            if j is None:
                stop_reason = 'state_limit'
                incomplete_level = n
                break
        if stop_reason is not None:
            break
    graph_exhausted = not heap and cert is None and stop_reason is None
    exhaustive_settings = (args.max_routes <= 0 and args.hits_per_state <= 0 and
                           args.max_states <= 0)
    if graph_exhausted and exhaustive_settings:
        exhausted = True
        stats['level_done_time'][cur_level] = round(time.monotonic() - t0, 1)
    stats['exhausted'] = exhausted
    stats['stop_reason'] = ('success' if cert is not None else
                            'exhausted' if exhausted else
                            'configured_caps' if graph_exhausted else
                            stop_reason or 'incomplete')
    stats['job_interrupted'] = stats['stop_reason'] == 'job_time'
    stats['alternating_flype_quotient'] = alternating
    frontier = heap[0][0] if heap else None
    if incomplete_level is not None:
        frontier = incomplete_level if frontier is None else min(frontier, incomplete_level)
    stats['frontier_min_crossings'] = frontier
    stats['n_states'] = len(graphs)
    stats['seconds'] = round(time.monotonic() - t0, 1)
    stats['seconds_screening'] = round(t_screen[0], 1)
    stats['states'] = dict(sorted(stats['states'].items()))
    stats['expanded'] = dict(sorted(stats['expanded'].items()))
    stats['screened'] = dict(stats['screened'])
    stats['det_hits'] = dict(stats['det_hits'])
    stats['lick_pass'] = dict(stats['lick_pass'])
    if cert is not None:
        stats['success'] = {'kind': cert['kind'], 'j': cert['j'], 'u': cert['u'], 'partner': cert['partner'],
                            'n_moves': len(cert['moves'])}
    stats['n_candidates'] = len(stats['candidates'])
    log(f"{t['name']}: k={k} states={len(graphs)} {stats['states']} stop={stats['stop_reason']} candidates={len(stats['candidates'])} "
        f"{stats['seconds']}s (screening {stats['seconds_screening']}s, {stats['skipped_flype_class']} flype-duplicates skipped; det hits {stats['det_hits']}, Lickorish pass {stats['lick_pass']}) success={stats['success']}")
    return stats, cert


_WORKER_KNOWN = None
_WORKER_ARGS = None


def _worker_init(known_path, args_dict):
    global _WORKER_KNOWN, _WORKER_ARGS
    _WORKER_KNOWN = Known(known_path)
    _WORKER_ARGS = argparse.Namespace(**args_dict)


def _worker_search(t):
    """Top-level multiprocessing entry point (therefore spawn-safe)."""
    if _WORKER_ARGS.job_deadline and time.monotonic() >= _WORKER_ARGS.job_deadline:
        return None
    # progress lines are printed by the worker as they happen (flushed), so the
    # scheduler log shows levels completing long before the knot finishes
    stats, cert = search_knot(t, _WORKER_KNOWN, _WORKER_ARGS, lambda s: print(s, flush=True))
    return stats, cert, []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--k-class', default='3', choices=['1', '2', '3', 'all'])
    ap.add_argument('--names', nargs='*')
    ap.add_argument('--control', action='store_true', help='run the nine settled knots (data/controls.json)')
    ap.add_argument('--targets', default=None,
                    help='alternative target list (JSON, same fields as data/targets_23_34.json)')
    ap.add_argument('--known', default=None,
                    help='alternative known-u list (JSON, same layout as data/known_u_codes.json)')
    ap.add_argument('--r2', action='store_true',
                    help='also generate Reidemeister II additions (fuller move closure, many more states)')
    ap.add_argument('--max-crossings', type=int, default=16)
    ap.add_argument('--extra-crossings', type=int, default=3,
                    help='budget relative to the knot: crossing number + this (capped by --max-crossings)')
    ap.add_argument('--time-per-knot', type=float, default=10800.0,
                    help='seconds allowed per knot; 0 disables the per-knot cap')
    ap.add_argument('--wall-time', type=float, default=0.0,
                    help='total seconds for this invocation; 0 disables the job deadline')
    ap.add_argument('--workers', type=int, default=1,
                    help='parallel worker processes sharing a dynamic target queue')
    ap.add_argument('--max-states', type=int, default=0,
                    help='states per knot; 0 is unlimited (needed for exhaustive runs)')
    ap.add_argument('--pass-extra', type=int, default=3, help='pass moves may increase the crossing number by this much')
    ap.add_argument('--max-routes', type=int, default=0,
                    help='pass routes per run; 0 enumerates every route')
    ap.add_argument('--routes', default='1,2', help='screening routes, e.g. 1,2 or 1,2,3')
    ap.add_argument('--hits-per-state', type=int, default=0,
                    help='determinant hits reduced per state/route; 0 tests every hit')
    ap.add_argument('--cand-store', type=int, default=30, help='candidates kept per knot for local resolution')
    ap.add_argument('--cand-max-crossings', type=int, default=16)
    ap.add_argument('--out-dir', default=os.environ.get('ENUM_OUT', 'enum_runs'))
    ap.add_argument('--result-tag', default=None, help='stable suffix for the JSONL result file')
    args = ap.parse_args()
    if args.workers < 1:
        ap.error('--workers must be at least 1')
    if args.extra_crossings is not None and args.extra_crossings < 0:
        ap.error('--extra-crossings must be nonnegative')
    args.job_deadline = time.monotonic() + args.wall_time if args.wall_time > 0 else 0.0

    known_path = args.known or os.path.join(HERE, 'data', 'known_u_codes.json')
    if args.control:
        targets = json.load(open(os.path.join(HERE, 'data', 'controls.json')))
    elif args.targets:
        targets = json.load(open(args.targets))
        if args.k_class != 'all':
            targets = [t for t in targets if t['k'] == int(args.k_class)]
    else:
        targets = json.load(open(os.path.join(HERE, 'data', 'targets_23_34.json')))
        if args.k_class != 'all':
            targets = [t for t in targets if t['k'] == int(args.k_class)]
    if args.names:
        want = {n if '_' in n else n.replace('a', 'a_', 1).replace('n', 'n_', 1) for n in args.names}
        targets = [t for t in targets if t['name'] in want]
    i, n = (int(x) for x in args.shard.split('/'))
    targets = targets[i::n]
    os.makedirs(args.out_dir, exist_ok=True)
    tag = args.result_tag or ('control' if args.control else f'{i}_{n}')
    res_path = os.path.join(args.out_dir, f'results_{tag}.jsonl')
    done = set()
    # Scan every result file, not only this worker configuration.  A resumed
    # campaign may therefore change its process count without redoing knots.
    for old_path in glob.glob(os.path.join(args.out_dir, 'results_*.jsonl')):
        for line in open(old_path):
            if not line.strip():
                continue
            try:
                done.add(json.loads(line)['name'])
            except (ValueError, KeyError):
                # A scheduler kill can leave one final, incomplete JSON line.
                continue
    targets = [t for t in targets if t['name'] not in done]
    log = lambda s: print(s, flush=True)  # noqa: E731

    log(f'run {tag}: {len(targets)} pending targets, {len(done)} already recorded; '
        f'class={args.k_class}, extra_crossings={args.extra_crossings}, max_crossings={args.max_crossings}, '
        f'workers={args.workers}, time/knot={args.time_per_knot}s, wall_time={args.wall_time}s, '
        f'routes={args.routes}, max_routes={args.max_routes}, hits/state={args.hits_per_state}, r2={args.r2}')

    def record(result, handle):
        if result is None:
            return False
        stats, cert, messages = result
        for message in messages:
            log(message)
        if stats.get('job_interrupted'):
            # Do not put a deadline-interrupted knot in the resume ledger.  Its
            # partial statistics are retained separately and overwritten next run.
            pp = os.path.join(args.out_dir, f"PARTIAL_{stats['name']}.json")
            with open(pp, 'w') as ph:
                json.dump(stats, ph, indent=1)
            return False
        if cert is not None:
            cp = os.path.join(args.out_dir, f"SUCCESS_{stats['name']}.json")
            with open(cp, 'w') as ch:
                json.dump(cert, ch, indent=1)
            log(f"  *** SUCCESS {stats['name']}: u <= {cert['u']} ({cert['kind']}, "
                f"j={cert['j']}, partner={cert['partner']}) -> {cp}")
        handle.write(json.dumps(stats) + '\n')
        handle.flush()
        return True

    completed = 0
    with open(res_path, 'a') as h:
        if args.workers == 1:
            known = Known(known_path)
            for t in targets:
                if args.job_deadline and time.monotonic() >= args.job_deadline:
                    break
                result = (*search_knot(t, known, args, log), [])
                completed += int(record(result, h))
        else:
            args_dict = vars(args).copy()
            with multiprocessing.Pool(args.workers, initializer=_worker_init,
                                      initargs=(known_path, args_dict)) as pool:
                try:
                    for result in pool.imap_unordered(_worker_search, targets, chunksize=1):
                        completed += int(record(result, h))
                except KeyboardInterrupt:
                    pool.terminate()
                    raise
    log(f'run {tag}: recorded {completed} new knot results; output={res_path}')


if __name__ == '__main__':
    main()
