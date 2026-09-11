"""The marked-move engine of reduce_witnesses.children, extended by SAFE pass moves.

A pass move erases a run of consecutive passages that are all over (or all under) and redraws it
along another route, again over (under) the rest of the diagram.  Only runs that contain no marked
crossing are used.  For such a run, changing the marked crossings alters no passage of the run, so
the run is still a pass strand of the changed diagram, and the new route crosses the same static
arcs in both diagrams.  Hence the move commutes with the crossing changes at the marks.  A run that
contains a marked crossing is never used: the pass deletes every crossing of its run.

Mark tracking: xtait.passmove.apply_pass keeps the rows outside the run in their original order and
appends the new rows, and tait.from_pd numbers edges by PD row, so a marked row r becomes
r - #{run rows < r}.  Each pass child is also guarded: the determinant of the child with the tracked
marks changed must equal that of the input witness, otherwise the child is discarded and counted.
"""
from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

from reduce_witnesses import children as base_children  # noqa: E402
from tait import graph as tg, moves as tmv, invariants as inv  # noqa: E402
from xtait.passmove import find_runs, pass_move, iter_routes, _mirror_pd, _mirror_run  # noqa: E402

STATS = {'pass_children': 0, 'pass_guard_rejected': 0, 'pass_errors': 0, 'unsafe_runs_skipped': 0}


def changed_det(h, marks):
    for e in marks:
        h = tmv.crossing_change(h, e)
    return inv.determinant(h)


def pass_children(g, marks, det_target, pass_extra=0, pass_routes=4):
    pd = tg.to_pd(g)
    ids = sorted(g.edges)
    markrows = sorted(ids.index(m) for m in marks)
    for ri, run in enumerate(find_runs(pd)):
        if set(run['crossings']) & set(markrows):
            STATS['unsafe_runs_skipped'] += 1
            continue
        try:
            base = pd if run['over'] else _mirror_pd(pd)
            brun = run if run['over'] else _mirror_run(base, run)
            routes = list(iter_routes(base, brun, max_len=len(run['crossings']) + pass_extra,
                                      max_routes=pass_routes))
        except ValueError:
            continue
        runrows = sorted(run['crossings'])
        newrows = [r - sum(1 for x in runrows if x < r) for r in markrows]
        for route in routes:
            try:
                h = tg.from_pd(pass_move(pd, run, route))
            except Exception:
                STATS['pass_errors'] += 1
                continue
            hid = sorted(h.edges)
            S2 = {hid[r] for r in newrows}
            if changed_det(h, S2) != det_target:
                STATS['pass_guard_rejected'] += 1
                continue
            STATS['pass_children'] += 1
            yield h, S2, ('pass', [ri, [list(s) for s in route]])


def children_with_pass(g, marks, det_target, pass_extra=0, pass_routes=4):
    yield from base_children(g, marks)
    yield from pass_children(g, marks, det_target, pass_extra, pass_routes)
