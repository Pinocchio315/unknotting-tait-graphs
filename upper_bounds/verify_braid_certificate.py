#!/usr/bin/env python
"""Verify the forward certificate for u(13n1587) <= 2.

The certificate consists of two 7-strand braid words differing in exactly one letter
(sigma_i <-> sigma_i^{-1}).  The first closure is converted to an embedded Tait graph, the
corresponding edge sign is changed, and the result is reduced by Reidemeister and pass moves.
The script checks, from scratch:

  1. the two words differ in exactly one position, by inversion of one generator;
  2. the closure of the first word is 13n1587 and of the second 10_113 (SnapPy isometry
     of the knot exteriors against the KnotInfo PD codes; Gordon-Luecke);
  3. changing the corresponding Tait edge produces the second closure;
  4. one Reidemeister II removal and three pass moves reduce the changed diagram from 22
     to 16 crossings, and the reduced endpoint is still 10_113;
  5. a crossing change in the reference diagram of 10_113 gives the unknot,
     independently verified by knot Floer genus detection.

Together: u(13n1587) <= 1 + u(10_113) = 2. This construction was published by
Brittenham--Hermiller, arXiv:1705.05985v2, Section 3. The unchanged archival
JSON's attribution to Applebaum and its unresolved [1,2] status are historical
descriptive errors, documented in results/README.md. The new, separate Greene
lower bound gives u(13n1587)=2; this verifier establishes only the upper bound.

    python upper_bounds/verify_braid_certificate.py
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))


def isometric(link_a, link_b, tries: int = 12) -> bool:
    """Check exterior isometry numerically (hyperbolic case), retrying with randomisation."""
    Ea, Eb = link_a.exterior(), link_b.exterior()
    for _ in range(tries):
        try:
            if Ea.is_isometric_to(Eb):
                return True
        except RuntimeError:
            pass
        Ea.randomize()
        Eb.randomize()
    return False


def require(condition, message):
    """Certificate validation remains active when Python runs with -O."""
    if not condition:
        raise ValueError(message)


def main() -> None:
    import tait  # noqa: F401
    from spherogram import Link
    from tait import knotinfo, invariants as inv
    from xtait.graph import from_pd, to_pd
    from xtait import moves as mv
    from xtait.reduce import greedy_reduce
    from flipdet import exact_det

    cert = json.load(open(os.path.join(HERE, '..', 'results',
                                       '13n_1587_u_le_2_braid_certificate.json')))
    w1, w2 = cert['braid_13n_1587'], cert['braid_10_113']

    diff = [i for i, (a, b) in enumerate(zip(w1, w2)) if a != b]
    require(all(type(x) is int and x != 0 for x in w1 + w2), 'invalid braid generator')
    require(len(w1) == len(w2) and len(diff) == 1 and w1[diff[0]] == -w2[diff[0]],
            'words do not differ by a single crossing change')
    print(f'1. single crossing change at letter {diff[0]} '
          f'(sigma_{abs(w1[diff[0]])}^{{{ "+" if w1[diff[0]]>0 else "-" }1}} inverted): OK')

    closures = {}
    for word, name in ((w1, '13n_1587'), (w2, '10_113')):
        closure = Link(braid_closure=word)
        reference = Link(knotinfo.pd_code(name))
        require(isometric(closure, reference), f'closure not isometric to KnotInfo {name}')
        closures[name] = closure
        print(f'2. closure of braid_{name} is {name} (SnapPy isometry): OK')

    source = from_pd([list(q) for q in closures['13n_1587'].PD_code()])
    edge = sorted(source.edges)[diff[0]]
    changed = mv.crossing_change(source, edge)
    require(isometric(Link(to_pd(changed)), closures['10_113']),
            'Tait edge-sign change does not produce the recorded 10_113 closure')
    print(f'3. Tait edge-sign change at row {diff[0]} produces the 10_113 closure: OK')

    reduced, reduction_path = greedy_reduce(changed)
    reduction_kinds = [kind for kind, _params in reduction_path]
    require(changed.n_crossings() == 22 and reduced.n_crossings() == 16,
            'unexpected crossing counts in the deterministic reduction')
    require(reduction_kinds == ['r2-s', 'pass', 'pass', 'pass'],
            f'unexpected deterministic reduction path: {reduction_kinds}')
    reference_10_113 = Link(knotinfo.pd_code('10_113'))
    require(isometric(Link(to_pd(reduced)), reference_10_113),
            'reduced endpoint not isometric to KnotInfo 10_113')
    print('4. R2 removal + three pass moves: 22 -> 16 crossings; endpoint 10_113: OK')

    partner = from_pd(knotinfo.pd_code('10_113'))
    unknotting_row = None
    for row, edge in enumerate(sorted(partner.edges)):
        neighbor = mv.crossing_change(partner, edge)
        if exact_det(neighbor) == 1 and inv.is_unknot(Link(to_pd(neighbor))):
            unknotting_row = row
            break
    require(unknotting_row is not None, 'no unknotting crossing of 10_113 verified')
    require(exact_det(partner) > 1, 'nontriviality check of 10_113 failed')
    print(f'5. u(10_113) = 1: reference crossing {unknotting_row} gives the unknot (HFK): OK')

    print('\n=> Published upper bound u(13n1587) <= 1 + u(10_113) = 2 verified. '
          'The separately computed Greene lower bound determines u=2.')


if __name__ == '__main__':
    main()
