#!/usr/bin/env python
"""Owens's obstruction to n crossing changes, driven by correction terms supplied from outside.

The published rank-n tests read the correction terms of the double branched cover off the sharp
positive-definite Goeritz form of an alternating diagram.  When the cover is an L-space, Greene's
spanning-tree model supplies the same correction terms for a knot that need not be alternating, and the
enumeration of definite forms is unchanged.  This module contains only the comparison; the forms come
from `owens_obstruction` (rank two), `owens_u3` and `owens_u4`.

Conventions.  The published tests use the orientation with d(spin) = -|sigma|/4 and place the spin
structure at the origin of the discriminant group.  Greene's model gives d(spin) = -sigma/4 for the
signed signature, so the vector is negated when sigma < 0, which is the mirror image.

    python greene_ranks.py --validate 5_1 7_3 9_10 9_13     # alternating controls against the published test
"""
from __future__ import annotations
import math, os, sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
# The forms come from ../owens in the repository.  A server copy may be flat, or may carry owens/ as a
# subdirectory, so every layout is offered to the import machinery before the modules are loaded.
for _candidate in (HERE.parent / 'owens', HERE / 'owens', HERE):
    if (_candidate / 'owens_obstruction.py').is_file():
        sys.path.insert(0, str(_candidate))
        break
else:
    raise ImportError('owens_obstruction.py not found: copy the owens modules next to this file, '
                      'into an owens/ subdirectory, or run from the repository')
from owens_obstruction import candidates, qtilde, m_Q, det_int, definite_goeritz_candidates  # noqa: E402


def cyclic_index(class_map, D):
    """Index a cyclic discriminant group by Z/D with the identity at 0, or None when it is not cyclic."""
    for generator in class_map.all_elements():
        point = tuple(0 for _ in class_map.d)
        index = {}
        for k in range(D):
            index[k] = point
            point = tuple((x + y) % d for x, y, d in zip(point, generator, class_map.d))
        if point == index[0] and len(set(index.values())) == D:
            return index
    return None


def form_admits(m_form, class_map, d, D):
    """Does some group isomorphism fixing the spin origin satisfy m(g) >= d(phi(g)) and m == d mod 2?

    Returns True, False, or None when the m-function is incomplete.  A form whose discriminant group is
    not cyclic cannot be isomorphic to a cyclic H_1, so it is rejected rather than left undecided.
    """
    index = cyclic_index(class_map, D)
    if index is None:
        return False
    if any(index[k] not in m_form for k in range(D)):
        return None
    for unit in range(1, D):
        if math.gcd(unit, D) != 1:
            continue
        for k in range(D):
            difference = m_form[index[k]] - d[unit * k % D]
            if difference < 0 or difference.denominator != 1 or difference.numerator % 2:
                break
        else:
            return True
    return False


def candidate_forms(D, n):
    """The definite forms allowed by Theorem `montesinos` for an unknotting sequence of n negative changes."""
    if n == 2:
        return [(f'({m1},{m2},{a})', qtilde(m1, m2, a)) for (m1, m2, a) in candidates(D, n_even=2)]
    if n == 3:
        from owens_u3 import candidate_plumbings
        return candidate_plumbings(D, n_even=3)
    if n == 4:
        from owens_u4 import candidate_plumbings
        return candidate_plumbings(D, n_even=4)
    raise ValueError(f'rank {n} is not implemented')


def rank_test(d, D, n):
    """Exclude an unknotting sequence of n crossing changes, all negative, from the correction terms.

    `d` maps Z/D to the correction terms in the orientation with d(0) = -|sigma|/4 = -n/2; a verdict of
    OBSTRUCTED proves u >= n + 1.
    """
    spin = Fraction(-n, 2)
    if d.get(0) != spin:
        return {'verdict': 'ERROR', 'reason': f'd(spin) = {d.get(0)} is not -|sigma|/4 = {spin}'}
    forms = candidate_forms(D, n)
    if not forms:
        return {'verdict': 'OBSTRUCTED', 'candidates': 0, 'reason': 'no form of this rank and determinant'}
    undecided = 0
    for key, Q in forms:
        assert det_int(Q) == D, (key, det_int(Q), D)
        m_form, class_map = m_Q(Q)
        verdict = form_admits(m_form, class_map, d, D)
        if verdict is True:
            return {'verdict': 'PASS', 'candidates': len(forms), 'witness_form': str(key)}
        if verdict is None:
            undecided += 1
    if undecided:
        return {'verdict': 'UNDECIDED', 'candidates': len(forms), 'undecided': undecided}
    return {'verdict': 'OBSTRUCTED', 'candidates': len(forms)}


def oriented(d, sigma):
    """Put a Greene correction-term vector in the orientation used by the published tests."""
    return {k: -v for k, v in d.items()} if int(sigma) < 0 else dict(d)


def goeritz_d(pd, sigma, D):
    """Correction terms of the cover of an alternating knot from its sharp Goeritz form, indexed cyclically."""
    target = Fraction(-abs(int(sigma)), 4)
    for G in definite_goeritz_candidates(pd):
        if det_int(G) != D:
            continue
        m_form, class_map = m_Q(G)
        if m_form.get(class_map.coords([0] * len(G))) != target:
            continue
        index = cyclic_index(class_map, D)
        return None if index is None else {k: m_form[index[k]] for k in range(D)}
    return None


def main():
    import argparse, ast, json
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('knots', nargs='*', default=['5_1', '7_3', '7_5', '8_2', '9_10', '9_13'])
    ap.add_argument('--validate', action='store_true',
                    help='alternating controls: compare with the Goeritz correction terms and the published verdict')
    args = ap.parse_args()
    if not args.validate:
        ap.error('this module is a library; use --validate to run the controls')
    import database_knotinfo as dk
    from owens_obstruction import obstruct_u2
    from owens_u3 import obstruct_u3
    from pin_dinv import pin
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    published = {2: obstruct_u2, 3: obstruct_u3}
    for name in args.knots:
        row = rows[name]
        D, sigma = int(row['determinant']), int(row['signature'])
        n = abs(sigma) // 2
        if n not in published:
            print(f'{name}: |sigma| = {abs(sigma)} has no published comparison, skipped')
            continue
        pd = [list(map(int, c)) for c in ast.literal_eval(row['pd_notation'])]
        reference = goeritz_d(pd, sigma, D)
        current, _, _ = pin(name, pd, D, verbose=False)
        if any(len(values) != 1 for values in current.values()):
            print(f'{name}: the pages leave classes ambiguous, skipped')
            continue
        greene = oriented({k: next(iter(v)) for k, v in current.items()}, sigma)
        same = reference is not None and any(
            all(greene[unit * k % D] == reference[k] for k in range(D))
            for unit in range(1, D) if math.gcd(unit, D) == 1)
        mine = rank_test(greene, D, n)
        theirs = published[n](pd, sigma, D)
        flag = 'OK' if mine['verdict'] == theirs['verdict'] else 'MISMATCH'
        print(f'{name}: det {D} sigma {sigma} rank {n} | Greene d == Goeritz m_G: {same} | '
              f'this module {mine["verdict"]} | published {theirs["verdict"]} | {flag}', flush=True)


if __name__ == '__main__':
    main()
