#!/usr/bin/env python
"""13a647 as a counterexample to the Bernhard-Jablan conjecture (run LOCALLY, unknot-venv).

The BJ conjecture asserts every knot K has a MINIMAL diagram with a crossing whose change
produces K' with u(K') = u(K) - 1.  For the alternating knot 13a647 (u = 3, this project)
the check is FINITE: by the Tait flyping theorem (Menasco-Thistlethwaite) together with the
Kauffman-Murasugi-Thistlethwaite theorem that minimal diagrams of a prime alternating knot
are exactly its reduced alternating diagrams, the minimal diagrams form one flype orbit.

Two facts make the computation short and robust:
  * flypes commute with crossing changes (a flype is an isotopy that moves a crossing across a
    tangle and is available whatever the sign of that crossing), so the knot obtained by
    changing a given crossing is the same in every diagram of the orbit --- changing the 13
    crossings of ONE reduced alternating diagram already decides the question;
  * the orbit is nevertheless enumerated completely (xtait/flypes.py: flypes on both
    checkerboard graphs, tangles = unions of consecutive bridges) and every single crossing
    change on every minimal diagram is identified, as a consistency check.

This script enumerates the orbit, applies EVERY single crossing change on EVERY minimal
diagram, identifies each changed knot rigorously (canonical-code match; isometry against the
KnotInfo diagram; for the composite, canonical-code membership in the flype orbit of an
explicit connected-sum diagram), and reports u of each.  Outcome (verified 2026-09-04):

    orbit = 15 diagrams on the sphere (8 up to turning the sphere over), 195 changes,
    the same 5 knot types with the same multiplicities in every diagram:
      11a62 (x3 per diagram, u=3, KnotInfo), 11a63 (x4, u=3, THIS PAPER's Owens sweep; KnotInfo
      has [2,3]), mirror 9_3 (x2, u=3), mirror 11a235 (x2, u=3),
      4_1 # mirror(7_3) (x2, u=3: Scharlemann's theorem --- a u=2 knot is prime or a sum of two
      u=1 knots --- with u(7_3)=2; independently, Owens' obstruction OBSTRUCTS u=2 for it).

    Every changed knot has u = 3, never 2  =>  no minimal-diagram change lowers u:
    13a647 is a CONSTRUCTIVE, ALTERNATING counterexample to Bernhard-Jablan
    (the BJ recursion value of 13a647 is 4, but u = 3).
"""
from __future__ import annotations

import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))
sys.path.insert(0, os.path.join(HERE, '..', 'lower_bounds', 'owens'))

from xtait.graph import from_pd, to_pd, pd_connected_sum  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.reduce import greedy_reduce  # noqa: E402
from xtait.jones import jones  # noqa: E402
from xtait.flypes import flype_orbit, orbit_up_to_turning_over  # noqa: E402


def mirror_pd(pd):
    g = from_pd([list(q) for q in pd])
    for e in g.signs:
        g.signs[e] = -g.signs[e]
    return to_pd(g)


def isometric(pd_a, pd_b, tries=10):
    from spherogram import Link
    Ea = Link([list(q) for q in pd_a]).exterior()
    Eb = Link([list(q) for q in pd_b]).exterior()
    for _ in range(tries):
        try:
            if Ea.is_isometric_to(Eb):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False


def main() -> None:
    import tait  # noqa: F401
    from tait import knotinfo

    checks = []

    def check(label, ok):
        checks.append(ok)
        print(('  ok  ' if ok else '  FAIL') + ' ' + label)

    # ---- 1. all minimal diagrams (= complete flype orbit of the reduced alternating diagram)
    g0 = from_pd([list(q) for q in knotinfo.pd_code('13a_647')])
    orbit = flype_orbit(g0)
    classes = orbit_up_to_turning_over(orbit)
    print(f'minimal diagrams of 13a647 (flype orbit): {len(orbit)} on the sphere, '
          f'{len(classes)} up to turning the sphere over')
    check('every diagram of the orbit is a 13-crossing alternating diagram',
          all(len(h.edges) == 13 and len(set(h.signs.values())) == 1 for h in orbit.values()))

    # ---- 2. every single change on every diagram, identified rigorously -----------------
    allk = json.load(open(os.path.join(HERE, 'data', 'all_knot_codes.json')))
    code_info = {}
    for x in allk:
        for c in x['codes']:
            code_info[c] = x
    per_diagram = []
    unmatched = {}
    n_changes = 0
    for g in orbit.values():
        row = collections.Counter()
        for e in sorted(g.edges):
            n_changes += 1
            red, _ = greedy_reduce(mv.crossing_change(g, e))
            code = repr(canonical_code(red))
            info = code_info.get(code)
            if info is not None:
                row[info['name']] += 1
            else:
                unmatched.setdefault(code, [to_pd(red), 0])
                unmatched[code][1] += 1
                row[('unmatched', code)] += 1
        per_diagram.append(row)
    hits = collections.Counter()
    for row in per_diagram:
        for k, v in row.items():
            if not isinstance(k, tuple):
                hits[k] += v
    print(f'{n_changes} crossing changes on {len(orbit)} diagrams; code-matched partners: '
          + ', '.join(f'{nm} (x{c})' for nm, c in sorted(hits.items())))

    # isometry identification of unmatched prime partners
    known_prime = {}
    composite = []
    name_of_code = {}
    for code, (pd, cnt) in unmatched.items():
        found = None
        for nm in ('11a_63', '11a_235', '11a_62', '9_3'):
            if isometric(pd, knotinfo.pd_code(nm)):
                found = nm
                break
        if found:
            known_prime[found] = known_prime.get(found, 0) + cnt
            name_of_code[code] = found
        else:
            composite.append((pd, cnt))
            name_of_code[code] = '4_1#m7_3'
    print('isometry-identified:', known_prime, '; composite candidates:',
          sum(c for _, c in composite))

    # the composite: canonical-code membership in the flype orbit of 4_1 # mirror(7_3)
    S = pd_connected_sum(knotinfo.pd_code('4_1'), 1, mirror_pd(knotinfo.pd_code('7_3')), 1, 1)
    sum_orbit = set(flype_orbit(from_pd([list(q) for q in S])))
    comp_ok = all(repr(canonical_code(from_pd([list(q) for q in pd]))) in sum_orbit
                  and jones(pd) == jones(S) for pd, _ in composite)
    check('composite changes are 4_1 # mirror(7_3) (flype-orbit code + Jones)', comp_ok)

    # flypes commute with crossing changes: identical multiset of changed knots in every diagram
    def named(row):
        out = collections.Counter()
        for k, v in row.items():
            out[name_of_code[k[1]] if isinstance(k, tuple) else k] += v
        return out
    rows = [named(r) for r in per_diagram]
    check('every minimal diagram yields the same multiset of changed knots: '
          + ', '.join(f'{k} x{v}' for k, v in sorted(rows[0].items())),
          all(r == rows[0] for r in rows))
    check('per diagram: 11a62 x3, 11a63 x4, 9_3 x2, 11a235 x2, 4_1#m7_3 x2 (13 changes)',
          rows[0] == collections.Counter({'11a_62': 3, '11a_63': 4, '9_3': 2, '11a_235': 2, '4_1#m7_3': 2}))

    # ---- 3. u of every changed knot ---------------------------------------------------
    check('11a_62: u = 3 (KnotInfo)', knotinfo.unknotting_interval('11a_62') == (3, 3))
    check('9_3: u = 3 (KnotInfo)', knotinfo.unknotting_interval('9_3') == (3, 3))
    check('11a_235: u = 3 (KnotInfo)', knotinfo.unknotting_interval('11a_235') == (3, 3))
    u3 = {x['name'] for x in json.load(open(os.path.join(
        HERE, '..', 'results', 'owens_u3_2026-08-24.json')))}
    check('11a_63: u = 3 (THIS PAPER, Owens sweep)', '11a_63' in u3)
    check('7_3: u = 2 (KnotInfo) => u(4_1 # mirror(7_3)) = 3 by Scharlemann',
          knotinfo.unknotting_interval('7_3') == (2, 2))
    import _compat  # noqa: F401
    from owens_obstruction import obstruct_u2
    res = obstruct_u2(S, 4, 65)
    check(f"Owens independently obstructs u = 2 for the composite ({res['verdict']})",
          res['verdict'] == 'OBSTRUCTED')

    # ---- verdict ----------------------------------------------------------------------
    every = set(hits) | set(known_prime) | ({'4_1#m7_3'} if composite else set())
    print(f'\ndistinct changed knots: {sorted(every)}')
    if all(checks):
        print(f'ALL CHECKS PASSED: every one of the {n_changes} single crossing changes on the '
              f'{len(orbit)} minimal diagrams of 13a647 yields a knot with u = 3 (never 2).')
        print('=> 13a647 is a constructive, alternating counterexample to Bernhard-Jablan '
              '(recursion value 4, true u = 3).')
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
