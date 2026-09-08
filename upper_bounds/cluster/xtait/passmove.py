"""P0: the pass move (strand pick-up) as a recorded, replayable move — stdlib only.

A *pass strand* is a maximal run of consecutive passages of the knot that are all OVER (or all UNDER).
Classical pick-up lemma: erasing such a strand and redrawing it along any other embedded route, again
passing entirely over (under) the rest of the diagram, is an isotopy of the knot.  Surgery on PD codes:

    delete    every crossing of the run disappears (the strands crossed under it are re-joined);
    reroute   a new simple arc from the run's entry stub to its exit stub crossing a chosen sequence of
              arcs of the remaining diagram, passing over (for an over-run) each of them.

Scope (P0): simple routes — pairwise-distinct faces and pairwise-distinct arcs, and the route does not
cross the two stub arcs themselves.  This is what BFS shortest routes produce and matches the power of
spherogram's pickup simplifier in practice.

Determinism: a pass applies to a graph through the deterministic PD `to_pd(g)`; its parameters reference
that PD (run index in `find_runs`, route as (arc label, side) steps), so replaying a recorded move on the
same graph reproduces the same result exactly.  Every application re-validates planarity (`pd_faces`) and
the single-component condition (`pd_passages`). Bounded Kauffman-bracket regressions
are in tests/test_upper_bounds.py.
"""
from __future__ import annotations

from collections import deque

from .graph import (EmbeddedTaitGraph, from_pd, to_pd, pd_faces, pd_passages,
                    pd_normalize_labels)


# =============================================================================
# runs (pass strands)
# =============================================================================
def find_runs(pd) -> list:
    """Maximal same-type passage runs.  [{'over','crossings','entry','entry_arc','exit_arc'}, ...]
    (a fully descending diagram — a single run — returns [])."""
    if not pd:
        return []
    walk = pd_passages(pd)
    n2 = len(walk)
    kinds = [i % 2 for (_c, i) in walk]
    if len(set(kinds)) == 1:
        return []
    runs = []
    start = 0
    while kinds[start] == kinds[(start - 1) % n2]:
        start += 1
    i, seen = start, 0
    while seen < n2:
        cur = []
        j = i
        while kinds[j % n2] == kinds[i % n2] and len(cur) < n2:
            cur.append(walk[j % n2])
            j += 1
        c0, p0 = cur[0]
        cl, pl = cur[-1]
        runs.append({'over': bool(p0 % 2),
                     'crossings': [c for (c, _p) in cur],
                     'entry': (c0, p0),
                     'entry_arc': int(pd[c0][p0]),
                     'exit_arc': int(pd[cl][(pl + 2) % 4])})
        seen += len(cur)
        i = j
    return runs


# =============================================================================
# faces of an open-strand PD (loose arc labels get a U-turn at their tip)
# =============================================================================
def _open_faces(pd, loose: set):
    occ = {}
    for c, q in enumerate(pd):
        for i, l in enumerate(q):
            occ.setdefault(int(l), []).append((c, i))
    partner = {}
    for l, ds in occ.items():
        if len(ds) == 2:
            partner[ds[0]] = ds[1]
            partner[ds[1]] = ds[0]
        elif len(ds) == 1 and l in loose:
            partner[ds[0]] = ds[0]
        else:
            raise ValueError(f'label {l} occurs {len(ds)} times (loose={l in loose})')
    face_of, faces = {}, []
    for startd in partner:
        if startd in face_of:
            continue
        fi = len(faces)
        cur = startd
        while cur not in face_of:
            face_of[cur] = fi
            c2, i2 = partner[cur]
            cur = (c2, (i2 + 1) % 4)
        if cur != startd:
            raise ValueError('open face walk did not close')
        faces.append(fi)
    return face_of, occ


# =============================================================================
# the surgery (over-runs; under-runs via the mirror wrapper)
# =============================================================================
def _delete_run(pd, run):
    run_set = set(run['crossings'])
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for c in run['crossings']:
        union(int(pd[c][0]), int(pd[c][2]))
    static = [[find(int(x)) for x in q] for c, q in enumerate(pd) if c not in run_set]
    A, B = find(int(run['entry_arc'])), find(int(run['exit_arc']))
    if A == B:
        raise ValueError('pass: entry and exit stubs merged')
    return static, A, B


def iter_routes(pd, run, max_len: int, max_routes: int = 64, should_stop=None):
    """Yield simple reroutes for an OVER-run, in increasing route length.

    ``max_routes <= 0`` removes the route-count cap.  ``should_stop`` is an
    optional cheap callback used by long server runs to honour their wall-time
    deadline without waiting for a large route queue to drain.
    """
    static, A, B = _delete_run(pd, run)
    face_of, occ = _open_faces(static, {A, B})
    fa = face_of[occ[A][0]]
    fb = face_of[occ[B][0]]
    sides = {}
    for l, ds in occ.items():
        if len(ds) == 2:                      # stubs (loose arcs) may not be crossed in P0
            sides[l] = ((0, face_of[ds[0]], face_of[ds[1]]),
                        (1, face_of[ds[1]], face_of[ds[0]]))
    yielded = 0
    queue = deque()
    queue.append((fa, (), frozenset((fa,))))
    while queue and (max_routes <= 0 or yielded < max_routes):
        if should_stop is not None and should_stop():
            return
        f, path, used_faces = queue.popleft()
        if f == fb:
            yielded += 1
            yield [list(s) for s in path]
            continue
        if len(path) >= max_len:
            continue
        for l, opts in sides.items():
            for (di, f_from, f_to) in opts:
                if f_from != f:
                    continue
                if f_to in used_faces:
                    continue
                queue.append((f_to, path + ((l, di),), used_faces | {f_to}))


def enumerate_routes(pd, run, max_len: int, max_routes: int = 64) -> list:
    """Materialized compatibility wrapper around :func:`iter_routes`."""
    return list(iter_routes(pd, run, max_len, max_routes=max_routes))


_SIDE_FLIP = False        # global orientation convention of route sides; fixed by the selftest suite


def apply_pass(pd, run, route) -> list:
    """Erase the OVER-run and redraw along `route`.  Returns the new normalized PD (validated)."""
    if not run['over']:
        raise ValueError('apply_pass: over-runs only (see pass_move)')
    static, A, B = _delete_run(pd, run)
    face_of, occ = _open_faces(static, {A, B})
    current_face = face_of[occ[A][0]]
    visited_faces = {current_face}
    for (l, di) in route:
        if int(l) in (A, B) or len(occ.get(int(l), ())) != 2:
            raise ValueError('pass: route must cross ordinary arcs (not the stubs)')
        if di not in (0, 1):
            raise ValueError('pass: route side must be zero or one')
        ds = occ[int(l)]
        if face_of[ds[di]] != current_face:
            raise ValueError('pass: route is not a continuous walk from the entry stub')
        current_face = face_of[ds[1 - di]]
        if current_face in visited_faces:
            raise ValueError('pass: route must be simple in the dual graph')
        visited_faces.add(current_face)
    if current_face != face_of[occ[B][0]]:
        raise ValueError('pass: route does not reach the exit stub')
    # A simple dual path redraws one uniform over-run in the plane. Keeping
    # it above the static diagram is an isotopy; the under-run case is mirrored.
    if len(set(int(l) for (l, _d) in route)) != len(route):
        raise ValueError('pass: route arcs must be distinct')
    out = [list(q) for q in static]
    all_labels = {int(x) for q in out for x in q} | {A, B}
    fresh = max(all_labels) + 1
    prev_s = A
    new_rows = []
    for j, (l, di) in enumerate(route):
        l, di = int(l), int(di)
        ds = occ[l]
        d_other = ds[1 - di]
        m = fresh
        fresh += 1
        out[d_other[0]][d_other[1]] = m           # split arc l: label l stays on the ds[di] side
        s_next = B if j == len(route) - 1 else fresh
        if s_next != B:
            fresh += 1
        # The face the strand crosses FROM walks arc l from the ds[di] crossing towards the other one,
        # i.e. it enters the new crossing along piece l and leaves along piece m; the incoming sigma'
        # piece lies on that side.  With the global face-walk orientation this forces one ccw order,
        # independent of di (the di dependence only decides which piece keeps the label l).
        if _SIDE_FLIP:
            new_rows.append([l, s_next, m, prev_s])
        else:
            new_rows.append([l, prev_s, m, s_next])
        prev_s = s_next
    if not route:
        out = [[A if int(x) == B else int(x) for x in q] for q in out]
    new_pd = pd_normalize_labels(out + new_rows)
    pd_faces(new_pd)                              # planarity (raises if the side convention is violated)
    pd_passages(new_pd)                           # still a one-component knot diagram
    return new_pd


def _mirror_pd(pd) -> list:
    return [[q[1], q[2], q[3], q[0]] for q in pd]


def _mirror_run(mpd, run):
    """The same geometric run inside the mirrored PD (the traced orientation may be reversed there)."""
    want = run['crossings']
    cand = [r2 for r2 in find_runs(mpd) if r2['over'] and
            (r2['crossings'] == want or r2['crossings'] == want[::-1])]
    if not cand:
        cand = [r2 for r2 in find_runs(mpd) if r2['over'] and len(r2['crossings']) == len(want)
                and set(r2['crossings']) == set(want)]
    if len(cand) != 1:
        raise ValueError('pass: could not locate the mirrored run unambiguously')
    return cand[0]


def pass_move(pd, run, route) -> list:
    """Pass move for an over- OR under-run (under-runs via the mirror trick)."""
    if run['over']:
        return apply_pass(pd, run, route)
    mpd = _mirror_pd(pd)
    return _mirror_pd(apply_pass(mpd, _mirror_run(mpd, run), route))


# =============================================================================
# graph-level wrappers (deterministic; usable through moves.apply_move('pass', params))
# =============================================================================
def iter_graph_pass_candidates(g: EmbeddedTaitGraph, max_extra: int = 0, max_routes: int = 16,
                               min_run: int = 1, should_stop=None):
    """Pass-move parameter dicts for g: {'run_index', 'route', 'gain'} with
    gain = |run| - |route| >= -max_extra, routes enumerated shortest-first."""
    if not g.edges:
        return
    pd = to_pd(g)
    for ri, run in enumerate(find_runs(pd)):
        if should_stop is not None and should_stop():
            return
        if len(run['crossings']) < min_run:
            continue
        base, brun = (pd, run) if run['over'] else (_mirror_pd(pd), None)
        try:
            if brun is None:
                brun = _mirror_run(base, run)
            routes = iter_routes(base, brun, max_len=len(run['crossings']) + max_extra,
                                 max_routes=max_routes, should_stop=should_stop)
            for route in routes:
                yield {'run_index': ri, 'route': [list(s) for s in route],
                       'gain': len(run['crossings']) - len(route)}
        except ValueError:
            continue


def graph_pass_candidates(g: EmbeddedTaitGraph, max_extra: int = 0, max_routes: int = 16,
                          min_run: int = 1) -> list:
    """Materialized compatibility wrapper around :func:`iter_graph_pass_candidates`."""
    return list(iter_graph_pass_candidates(g, max_extra=max_extra, max_routes=max_routes,
                                           min_run=min_run))


def apply_pass_params(g: EmbeddedTaitGraph, params) -> EmbeddedTaitGraph:
    if isinstance(params, tuple):
        params = params[0]
    pd = to_pd(g)
    run = find_runs(pd)[int(params['run_index'])]
    route = [(int(a), int(b)) for a, b in params['route']]
    h = from_pd(pass_move(pd, run, route), g.knot_name)
    h.meta = dict(g.meta)
    return h
