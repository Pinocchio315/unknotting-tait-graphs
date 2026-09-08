#!/usr/bin/env python3
"""Exact diagram identification of both 12a684 neighbours of 14a2539.

Reduce by the recorded pass moves, then search flypes on both checkerboards.
Success requires equality of full canonical signed planar-map codes, allowing
reflection because unknotting number is mirror invariant. No numerical or
polynomial knot recognition is used. This is a bounded search: failure to
find a path before 55 seconds establishes nothing about knot equivalence.
"""
from collections import deque
import json
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'upper_bounds'))
from spherogram import Link
from xtait.graph import from_pd, to_pd, dual
from xtait.canonical import canonical_code
from xtait.reduce import greedy_reduce
from xtait.flypes import find_flypes_all, flype_block, turn_over

def code(g):
    return repr(canonical_code(g))

def serial(value):
    if isinstance(value, set):
        return sorted(value)
    raise TypeError(type(value).__name__)

def main():
    def timed_out(*_):
        raise TimeoutError('55 second bound reached; no negative conclusion')
    signal.signal(signal.SIGALRM, timed_out)
    signal.alarm(55)
    start = time.monotonic()
    raw = json.loads((HERE / 'bjset14a2539_replay.json').read_text())
    named = Link('12a684')
    target = from_pd(named.PD_code())
    variants = {'reference': target, 'turned_over': turn_over(target),
                'mirror': target.mirror(),
                'mirror_turned_over': turn_over(target.mirror())}
    target_codes = {code(g): label for label, g in variants.items()}
    results = []
    for row in (raw['rows'][0], raw['rows'][5]):
        pd = row['factors'][0]['pd']
        reduced, reductions = greedy_reduce(from_pd(pd))
        assert reduced.n_crossings() == 12
        todo = deque([(reduced, [])])
        seen = {code(reduced)}
        while todo:
            graph, path = todo.popleft()
            canonical = code(graph)
            if canonical in target_codes:
                results.append({'crossing_zero_based': row['crossing_zero_based'],
                                'input_pd': pd, 'reductions': reductions,
                                'flypes': path, 'endpoint_pd': to_pd(graph),
                                'target_variant': target_codes[canonical],
                                'canonical_signed_map_code': canonical,
                                'visited_diagrams': len(seen)})
                break
            for colouring, g in (('original', graph), ('dual', dual(graph))):
                for edge, piece in find_flypes_all(g):
                    h = flype_block(g, edge, piece)
                    canonical = code(h)
                    if canonical in seen:
                        continue
                    seen.add(canonical)
                    todo.append((h, path + [{'colouring': colouring,
                                             'edge': edge, 'piece': piece}]))
        else:
            raise RuntimeError('No canonical match in completed flype orbit')
    report = {'reference_pd': named.PD_code(), 'results': results,
              'elapsed_seconds': round(time.monotonic() - start, 3),
              'identification': 'exact canonical signed planar-map match after pass moves and flypes'}
    (HERE / '12a684_exact_identification.json').write_text(
        json.dumps(report, default=serial, indent=2) + '\n')
    signal.alarm(0)
    print([(r['crossing_zero_based'], len(r['reductions']), len(r['flypes']),
            r['target_variant'], r['visited_diagrams']) for r in results])
    print('Both identifications are exact; elapsed', report['elapsed_seconds'])

if __name__ == '__main__':
    main()
