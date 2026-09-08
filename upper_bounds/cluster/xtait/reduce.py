"""Pass-first crossing minimisation — greedy reducer and best-first search (stdlib only).

greedy_reduce   repeatedly apply R1-/R2- removals and the best crossing-reducing pass move
search_reduce   best-first over diagrams: passes (including sideways/inflating routes), R3, flype,
                R1-/R2- and sampled R2+ inflation, de-duplicated by canonical hash, with crossing budget,
                node/time budgets, per-crossing-level accounting and a bounded heap (OOM-safe).

Both record the full move path (elementary moves + 'pass' moves), replayable with verify.replay.
"""
from __future__ import annotations

import heapq
import itertools
import random
import time
from dataclasses import dataclass, field

from . import moves as mv
from .canonical import canonical_hash
from .graph import EmbeddedTaitGraph, ribbon_faces
from .passmove import graph_pass_candidates


# ----------------------------------------------------------------------------- greedy
def greedy_reduce(g: EmbeddedTaitGraph, max_rounds: int = 10000):
    """Monotone reduction: R1-/R2- whenever available, else the best strictly-reducing pass move.
    Returns (reduced graph, path)."""
    cur = g
    path = []
    for _ in range(max_rounds):
        progressed = False
        while True:
            rem = mv.find_r1_removals(cur)
            if rem:
                kind, e = rem[0]
                cur = mv.r1_remove(cur, e)
                path.append(('r1-', (e,)))
                progressed = True
                continue
            rem2 = mv.find_r2_removals(cur)
            if rem2:
                kind, params = rem2[0]
                cur = mv.apply_move(cur, 'r2-' + kind[0], params)
                path.append(('r2-' + kind[0], params))
                progressed = True
                continue
            break
        best = None
        for p in graph_pass_candidates(cur, max_extra=-1, max_routes=8):
            if p['gain'] > 0 and (best is None or p['gain'] > best['gain']):
                best = p
        if best is not None:
            cur = mv.apply_move(cur, 'pass', (best,))
            path.append(('pass', (best,)))
            progressed = True
        if not progressed:
            break
    return cur, path


# ----------------------------------------------------------------------------- best-first search
@dataclass
class ReduceConfig:
    target: int | None = None          # stop as soon as n_crossings <= target
    max_extra: int = 4                 # crossing budget above the *start* diagram
    pass_extra: int = 1                # allow pass routes up to |run| + pass_extra (sideways/inflating)
    pass_routes: int = 8               # routes per run
    inflate_samples: int = 6           # sampled R2+ moves per node
    max_nodes: int = 20000
    time_limit: float = 600.0
    depth_weight: float = 0.05
    jitter: float = 0.4
    heap_max: int = 500000             # prune the heap to half when exceeded
    per_level_cap: int = 200000
    seed: int = 0
    use_flype: bool = True
    use_r3: bool = True
    log_every: int = 2000


def _expansions(g, budget, cfg, rnd):
    out = []
    for kind, e in mv.find_r1_removals(g):
        out.append(('r1-', (e,)))
    for kind, params in mv.find_r2_removals(g):
        out.append(('r2-' + kind[0], params))
    for p in graph_pass_candidates(g, max_extra=cfg.pass_extra, max_routes=cfg.pass_routes):
        if g.n_crossings() - p['gain'] <= budget:
            out.append(('pass', (p,)))
    if cfg.use_r3:
        for fi in mv.find_r3_triangles(g):
            out.append(('r3D', (fi,)))
        for t in mv.find_r3_stars(g):
            out.append(('r3Y', (t,)))
    if cfg.use_flype:
        for e, pi in mv.find_flypes(g):
            out.append(('flype', (e, pi)))
    if g.n_crossings() + 2 <= budget and cfg.inflate_samples > 0:
        faces = ribbon_faces(g)
        par = [(fi, i, j) for fi, f in enumerate(faces)
               for i, j in itertools.combinations(range(len(f)), 2) if f[i][0] != f[j][0]]
        for _ in range(cfg.inflate_samples):
            if par:
                fi, i, j = par[rnd.randrange(len(par))]
                out.append(('r2+p', (fi, i, j, rnd.choice((1, -1)))))
    return out


def search_reduce(g0: EmbeddedTaitGraph, cfg: ReduceConfig, log=None):
    """Best-first crossing minimisation.  Returns (best_graph, best_path, stats); if cfg.target is set and
    reached, best_* is the first diagram at or below the target (goal hit; stats['goal'] = True)."""
    rnd = random.Random(cfg.seed)
    t0 = time.time()
    n0 = g0.n_crossings()
    budget = n0 + cfg.max_extra
    visited = {canonical_hash(g0)}
    counter = itertools.count()
    heap = []
    stats = {'expanded': 0, 'generated': 0, 'dedup': 0, 'by_crossings': {}, 'heap_prunes': 0,
             'min_crossings': n0, 'goal': False}
    best_g, best_path = g0, []

    def push(g, depth, path):
        pr = g.n_crossings() + cfg.depth_weight * depth + cfg.jitter * rnd.random()
        heapq.heappush(heap, (pr, next(counter), g, depth, path))

    push(g0, 0, [])
    while heap:
        if stats['expanded'] >= cfg.max_nodes or time.time() - t0 > cfg.time_limit:
            break
        if len(heap) > cfg.heap_max:
            keep = heapq.nsmallest(cfg.heap_max // 2, heap)
            heap.clear()
            heap.extend(keep)
            heapq.heapify(heap)
            stats['heap_prunes'] += 1
        _, _, g, depth, path = heapq.heappop(heap)
        n = g.n_crossings()
        lvl = stats['by_crossings'].get(n, 0)
        if lvl >= cfg.per_level_cap:
            continue
        stats['by_crossings'][n] = lvl + 1
        stats['expanded'] += 1
        if n < stats['min_crossings']:
            stats['min_crossings'] = n
            best_g, best_path = g, path
        if cfg.target is not None and n <= cfg.target:
            stats['goal'] = True
            stats['seconds'] = round(time.time() - t0, 1)
            return g, path, stats
        for kind, params in _expansions(g, budget, cfg, rnd):
            try:
                h = mv.apply_move(g, kind, params)
            except ValueError:
                continue
            stats['generated'] += 1
            hh = canonical_hash(h)
            if hh in visited:
                stats['dedup'] += 1
                continue
            visited.add(hh)
            push(h, depth + 1, path + [(kind, params)])
        if log and stats['expanded'] % cfg.log_every == 0:
            log(f'  expanded {stats["expanded"]} heap {len(heap)} min {stats["min_crossings"]} '
                f'levels {dict(sorted(stats["by_crossings"].items()))} {time.time() - t0:.0f}s')
    stats['visited'] = len(visited)
    stats['seconds'] = round(time.time() - t0, 1)
    return best_g, best_path, stats
