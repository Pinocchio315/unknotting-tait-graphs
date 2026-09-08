#!/usr/bin/env python
"""Generator bound sweep: u(K) >= minimal number of generators of H_1(Sigma_2(K)).

For every KnotInfo row with a gapped unknotting interval [lo, hi], H_1(Sigma_2) is computed
from the Goeritz matrix of the Tait graph (Smith normal form), and the minimal number of
generators is max_p dim_{F_p} H_1 tensor F_p (the maximum p-rank over primes p dividing det K,
equivalently the number of nonunit invariant factors in Smith normal form).  Whenever that number
exceeds lo, the lower bound improves; combined with the recorded upper bound it can be exact.
When 3 | det K, dim_3 is cross-checked against the Jones evaluation |V(e^{i pi/3})| = sqrt3^{dim_3}
(Lickorish-Millett).

    python lower_bounds/generator_bound_sweep.py [--out generator_bound.json]

Reproduces results/generator_bound_2026-08-24.json: the bound settles u = 3 for 13a1786,
13a2720, 13a2727 and improves 13a4877 from [2,5] to [3,5]; 179 further hits reprove knots
already obstructed by the linking form (reported with --all).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))
sys.path.insert(0, os.path.join(HERE, 'owens'))


def prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def min_generators(tor: list[int]) -> int:
    """Minimal generators of a direct sum of cyclic groups, allowing Z/0 = Z.

    For arbitrary cyclic summands (e.g. Z/2 + Z/3), count the largest
    p-rank, not the number of summands. For Smith invariant factors these
    counts agree. The free rank contributes at every prime.
    """
    tor = [abs(int(t)) for t in tor]
    free = tor.count(0)
    finite = [t for t in tor if t > 1]
    primes = {p for t in finite for p in prime_factors(t)}
    return free + max((sum(t % p == 0 for t in finite) for p in primes), default=0)


def h1_string(tor: list[int]) -> str:
    return '+'.join('Z' if t == 0 else f'Z/{abs(t)}'
                    for t in sorted(tor) if abs(t) != 1) or '0'


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=None, help='write hits as JSON')
    ap.add_argument('--all', action='store_true',
                    help='also list hits that only reprove an existing lower bound')
    args = ap.parse_args()

    import tait  # noqa: F401
    from spherogram import Link
    from tait import knotinfo, graph as tg, invariants as inv
    from kinfo import parse_jones
    from traczyk import eps_and_dim3

    # knots already obstructed by the linking-form sweep (Section 3, "Linking pairings") are not new here
    lick_path = os.path.join(HERE, '..', 'results',
                             'lickorish_cw_obstructed_815_knotinfo2026.8.1.json')
    already = {nm for nm, _ in json.load(open(lick_path))} if os.path.exists(lick_path) else set()

    hits, reproved = [], 0
    names = knotinfo.gapped_names()
    for i, nm in enumerate(names):
        lo, hi = knotinfo.unknotting_interval(nm)
        g = tg.from_link(Link(knotinfo.pd_code(nm)), nm)
        tor = inv.h1_torsion(g)
        gens = min_generators(tor)
        if gens <= lo and not args.all:
            continue
        det = 1
        for t in tor:
            det *= t
        dim3 = sum(1 for t in tor if t % 3 == 0)
        if det % 3 == 0:  # Jones cross-check of the 3-rank
            _, d3_jones = eps_and_dim3(parse_jones(nm))
            assert d3_jones == dim3, (nm, dim3, d3_jones)
        if gens <= lo or (nm in already and not args.all):
            reproved += 1
            continue
        r = knotinfo.row(nm)
        rec = {'name': nm, 'interval': f'[{lo},{hi}]', 'sigma': int(r['signature']),
               'det': det, 'H1': h1_string(tor), 'generators': gens, 'dim3': dim3,
               'conclusion': (f'u = {gens} (exact)' if gens == hi
                              else f'[{lo},{hi}] -> [{gens},{hi}]')}
        hits.append(rec)
        print(f"[{i+1}/{len(names)}] {nm}: H1 = {rec['H1']}, gens = {gens} > lo = {lo}"
              f" -> {rec['conclusion']}")

    print(f"\n{len(hits)} new improvements; {reproved} further hits are subsumed by the "
          "linking-form sweep (or merely reprove the recorded lower bound)")
    if args.out:
        with open(args.out, 'w') as h:
            json.dump(hits, h, indent=1)
        print('wrote', args.out)


if __name__ == '__main__':
    main()
