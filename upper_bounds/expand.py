"""R1-free diagram expansion: many inequivalent diagrams of the same knot at a target size.

A continuous random walk from the start diagram grows with R2 additions (parallel and series
--- each adds a cancelling opposite-sign pair, so the knot is unchanged and no kinks are ever
created; Reidemeister I is NOT used) and keeps moving at the target size with R3 moves, flypes,
and R2 remove/add churn (an R2- removal followed later by a fresh R2+ addition re-routes
strands, which is what makes the emitted diagrams genuinely different).  Snapshots are taken
whenever the walk sits at the target crossing number, de-duplicated by a cheap structural key
first and by the exact canonical code only on key collisions.

All moves are the exact, validated moves of xtait.moves; the move path of every emitted diagram
is recorded so certificates can be replayed from the KnotInfo start diagram.
"""
from __future__ import annotations

import random

from xtait import moves as mv
from xtait.canonical import canonical_code


def _grow(g, rnd):
    for _ in range(6):
        if rnd.random() < 0.7:
            c = mv.find_r2_parallel_additions(g)
            if c:
                fi, i, j = c[rnd.randrange(len(c))]
                return ('r2+p', (fi, i, j, rnd.choice((1, -1))))
        c = mv.find_r2_series_additions(g)
        if c:
            v, i, j = c[rnd.randrange(len(c))]
            return ('r2+s', (v, i, j, rnd.choice((1, -1))))
    return None


def _shrink(g, rnd):
    c = mv.find_r2_removals(g)
    if c:
        kind, params = c[rnd.randrange(len(c))]
        return ('r2-' + kind[0], params)
    return None


def _shuffle(g, rnd):
    r = rnd.random()
    if r < 0.55:                                        # R3, either direction
        c = [('r3D', (fi,)) for fi in mv.find_r3_triangles(g)]
        c += [('r3Y', (t,)) for t in mv.find_r3_stars(g)]
        if c:
            return c[rnd.randrange(len(c))]
    c = mv.find_flypes(g)                               # flype
    if c:
        return ('flype', c[rnd.randrange(len(c))])
    return None


def _pick_move(g, rnd, n_target):
    """One random knot-preserving move; grows below the target size, churns at it.
    The crossing number never exceeds n_target."""
    n = g.n_crossings()
    r = rnd.random()
    if n < n_target:
        return (_grow(g, rnd) if r < 0.75 else _shuffle(g, rnd)) or _grow(g, rnd)
    # at the target size: churn down (the walk regrows next steps) or shuffle in place
    if r < 0.25:
        m = _shrink(g, rnd)
        if m is not None:
            return m
    return _shuffle(g, rnd) or _shrink(g, rnd)


def expand_diagrams(g0, n_target=50, count=10000, seed=0, snap_every=3,
                    max_stale=60, restart_every=400, walk_cap=10 ** 7,
                    max_dry_restarts=25):
    """Yield (graph, path) pairs at n_target crossings; path replays from g0.

    Emits at most `count` canonically distinct diagrams.  R1 is never used; every move
    preserves the knot.  R2/R3/flype preserve the parity of the crossing number, so the
    effective target is the largest size <= n_target with the parity of the start diagram
    (e.g. 49 for a 13-crossing knot with n_target = 50).
    """
    n_target -= (n_target - g0.n_crossings()) % 2
    rnd = random.Random(seed)
    seen = set()
    emitted = stale = steps = since_restart = 0
    dry_restarts = 0          # consecutive restarts without a single new emission
    cur, path = g0, []
    since_snap = 0
    while emitted < count and steps < walk_cap:
        steps += 1
        if stale > max_stale or since_restart > restart_every:
            if since_restart == 0:
                dry_restarts += 1
                if dry_restarts > max_dry_restarts:
                    return    # the reachable set at this size is exhausted
            else:
                dry_restarts = 0
            cur, path = g0, []
            stale = since_restart = 0
        m = _pick_move(cur, rnd, n_target)
        if m is None:
            cur, path = g0, []
            continue
        try:
            nxt = mv.apply_move(cur, *m)
        except Exception:
            continue
        cur = nxt
        path.append(m)
        if cur.n_crossings() != n_target:
            continue
        since_snap += 1
        if since_snap < snap_every:
            continue
        since_snap = 0
        code = canonical_code(cur)
        if code in seen:
            stale += 1
            continue
        seen.add(code)
        stale = 0
        emitted += 1
        since_restart += 1
        yield cur, list(path)


def jsonable_path(path):
    """Expansion paths contain tuples; make them JSON-round-trippable."""
    def enc(x):
        if isinstance(x, tuple):
            return ['#t'] + [enc(y) for y in x]
        if isinstance(x, list):
            return [enc(y) for y in x]
        return x
    return [[kind, enc(list(params))] for kind, params in path]


def decode_path(jp):
    def dec(x):
        if isinstance(x, list):
            if x and x[0] == '#t':
                return tuple(dec(y) for y in x[1:])
            return [dec(y) for y in x]
        return x
    return [(kind, tuple(dec(params))) for kind, params in jp]
