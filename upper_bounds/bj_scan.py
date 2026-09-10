#!/usr/bin/env python
"""Historical Bernhard-Jablan exploration for prime alternating knots with known u.

This script reads earlier selected u=2/u=3 lists and the archived diagram-data
bounds. It does not select the current v1.6 cohort from the consolidated table.
The frozen Appendix F results are reproduced by consolidate_results.py.

For each given knot K (alternating, u known): enumerate ALL minimal diagrams (the complete
flype orbit on both checkerboard graphs --- Tait/KMT + Menasco-Thistlethwaite), apply every
single crossing change, attempt to identify each changed knot, and decide whether some change
realises u(K) - 1.  (Flypes commute with crossing changes, so every minimal diagram gives the
same multiset of changed knots; the orbit is enumerated anyway as a consistency check.)

Identification ladder: canonical-code match against every KnotInfo diagram (+mirrors);
SnapPy reference knots + numerical isometry comparison; connected-sum hypothesis (factor pairs from
KnotInfo filtered by det multiplicativity, confirmed by Jones product AND canonical-code
membership in the flype orbit of an explicit sum diagram). A nontrivial connected sum has
u >= 2 by Scharlemann; if both factors have u = 1, subadditivity gives equality.
In general, u <= the sum of the factor upper bounds. Numerical isometry results
are computational identifications, not interval-certified homeomorphism proofs.

Verdicts per knot:
  COUNTEREXAMPLE  every changed knot has decided u != u(K)-1
  BJ-HOLDS        some changed knot has decided u = u(K)-1
  UNDECIDED       neither (some changed knot's u is open, or unidentified)

    ~/.pyenv/versions/unknot-venv/bin/python bj_scan.py 11a63 11a64 11a144 11a299 11a320 11a329
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'tait_graphs'))

from xtait.graph import from_pd, to_pd, pd_connected_sum  # noqa: E402
from xtait import moves as mv  # noqa: E402
from xtait.canonical import canonical_code  # noqa: E402
from xtait.flypes import flype_orbit as _flype_orbit  # noqa: E402
from xtait.reduce import greedy_reduce  # noqa: E402
from xtait.jones import jones  # noqa: E402
from flipdet import exact_det  # noqa: E402


def flype_orbit(g0):
    """All minimal diagrams of a prime alternating knot: complete flype orbit on both
    checkerboard graphs (see xtait/flypes.py)."""
    return _flype_orbit(g0)


def mirror_pd(pd):
    g = from_pd([list(q) for q in pd])
    for e in g.signs:
        g.signs[e] = -g.signs[e]
    return to_pd(g)


def isometric(pd_a, pd_b, tries=8):
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


def composite_unknotting_range(u_a, u_b):
    """Bounds for two identified nontrivial summands with known exact values.

    Scharlemann excludes unknotting number one, and concatenating the two
    unknotting sequences gives the upper bound. Neither theorem raises the
    lower bound to three when one summand has unknotting number greater than
    one; in particular, values one and two leave the sum in [2, 3].
    """
    return 2, u_a + u_b


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('names', nargs='+', help='KnotInfo-style names, e.g. 11a63 or 11a_63')
    ap.add_argument('--overrides', default=None,
                    help='JSON {"exact": {name: u}, "open": [names]} applied on top of KnotInfo and the paper lists '
                         '(used for the rank-3 Owens determinations and to disregard DKT-derived KnotInfo values)')
    args = ap.parse_args()
    overrides = json.load(open(args.overrides)) if args.overrides else {'exact': {}, 'open': []}
    ov_exact = {k if '_' in k else k.replace('a', 'a_', 1).replace('n', 'n_', 1): int(v) for k, v in overrides.get('exact', {}).items()}
    ov_open = {k if '_' in k else k.replace('a', 'a_', 1).replace('n', 'n_', 1) for k in overrides.get('open', [])}

    import tait  # noqa: F401
    from tait import knotinfo
    from spherogram import Link

    allk = json.load(open(os.path.join(HERE, 'data', 'all_knot_codes.json')))
    code_info, by_name, buckets = {}, {}, {}
    for x in allk:
        by_name[x['name']] = x
        buckets.setdefault(x['crossings'], []).append(x)
        for c in x['codes']:
            code_info[c] = x

    # exact u after OUR updates (u3 list + u2 list of the paper)
    R = os.path.join(HERE, '..', 'results')
    ours_u3 = {y['name'] for y in json.load(open(os.path.join(R, 'owens_u3_2026-08-24.json')))}
    ours_u3 |= {'13a_1786'}   # 13a_647 removed: rank-3 Owens gives u = 4 (see --overrides)
    lick = json.load(open(os.path.join(R, 'lickorish_cw_obstructed_815_knotinfo2026.8.1.json')))
    ours_u2 = {nm for nm, iv in lick if iv.replace(' ', '') == '[1,2]'}
    ours_u2 |= {'13n_1166', '13n_2504', '13n_30', '13n_45', '13n_80', '13n_2379',
                '13n_2809', '13n_2907', '13n_3033', '13n_3589'}

    def decided_u(nm):
        if nm in ov_exact:
            return ov_exact[nm]
        if nm in ov_open:
            return None
        if nm in ours_u3:
            return 3
        if nm in ours_u2:
            return 2
        lo, hi = by_name[nm]['lo'], by_name[nm]['hi']
        return lo if lo == hi else None

    def identify(red):
        """-> ('prime', name) | ('sum', (A, B)) | (None, data)."""
        code = repr(canonical_code(red))
        if code in code_info:
            return 'prime', code_info[code]['name']
        pd = to_pd(red)
        try:
            ids = [str(m) for m in Link([list(q) for q in pd]).exterior().identify()]
        except Exception:
            ids = []
        import re
        for s0 in ids:
            s0 = re.sub(r'\(.*\)$', '', s0.strip())
            cands = [s0]
            if re.fullmatch(r'K\d+[an]\d+', s0):
                t = s0[1:]
                cands += [t, re.sub(r'([an])', r'\1_', t, count=1)]
            for cand in cands:
                if cand in by_name and isometric(pd, knotinfo.pd_code(cand)):
                    return 'prime', cand
        # connected-sum hypothesis: factor pairs with det multiplicativity + Jones + code
        n, det = red.n_crossings(), exact_det(red)
        VK = jones(pd)
        for A in allk:
            cB = n - A['crossings']
            if A['crossings'] < 3 or cB < A['crossings']:
                continue
            for B in buckets.get(cB, []):
                if A['det'] * B['det'] != det or (cB == A['crossings'] and A['name'] > B['name']):
                    continue
                pa0 = knotinfo.pd_code(A['name'])
                pb0 = knotinfo.pd_code(B['name'])
                for pa in (pa0, mirror_pd(pa0)):
                    for pb in (pb0, mirror_pd(pb0)):
                        for variant in (0, 1):
                            S = pd_connected_sum(pa, 1, pb, 1, variant)
                            if jones(S) != VK:
                                continue
                            if code in {repr(canonical_code(h))
                                        for h in flype_orbit(from_pd([list(q) for q in S])).values()}:
                                return 'sum', (A['name'], B['name'])
        return None, {'n': n, 'det': det}

    for raw in args.names:
        nm = raw if '_' in raw else raw.replace('a', 'a_', 1).replace('n', 'n_', 1)
        uK = decided_u(nm)
        print(f'\n=== {nm}: u = {uK} ===')
        if uK is None:
            print('  u not decided; skipping')
            continue
        g0 = from_pd([list(q) for q in knotinfo.pd_code(nm)], nm)
        orbit = flype_orbit(g0)
        partners = {}
        cache = {}
        unresolved = 0
        for g in orbit.values():
            for e in sorted(g.edges):
                red, _ = greedy_reduce(mv.crossing_change(g, e))
                code = repr(canonical_code(red))
                if code not in cache:
                    cache[code] = identify(red)
                kind, data = cache[code]
                if kind == 'prime':
                    partners.setdefault(('prime', data), 0)
                    partners[('prime', data)] += 1
                elif kind == 'sum':
                    partners.setdefault(('sum', data), 0)
                    partners[('sum', data)] += 1
                else:
                    unresolved += 1
        n_changes = sum(partners.values()) + unresolved
        print(f'  minimal diagrams: {len(orbit)}; changes: {n_changes}; unresolved: {unresolved}')
        found_down, undecided = [], []
        for (kind, data), cnt in sorted(partners.items(), key=str):
            if kind == 'prime':
                u = decided_u(data)
                iv = f'[{by_name[data]["lo"]},{by_name[data]["hi"]}]'
                ki_exact = by_name[data]['lo'] == by_name[data]['hi']
                src = ('KnotInfo' if u is not None and ki_exact
                       else 'ours' if u is not None else 'open')
                print(f'    {data:10s} x{cnt}  u = {u if u is not None else iv} ({src})')
                if u == uK - 1:
                    found_down.append(data)
                elif u is None and by_name[data]['lo'] <= uK - 1 <= by_name[data]['hi']:
                    undecided.append(data)
            else:
                A, B = data
                uA, uB = decided_u(A), decided_u(B)
                t = uK - 1
                print(f'    {A}#{B} x{cnt}  factors u = {uA}, {uB}')
                if uA is None or uB is None:
                    undecided.append(f'{A}#{B}')
                else:
                    lo_s, hi_s = composite_unknotting_range(uA, uB)
                    if lo_s == hi_s == t:
                        found_down.append(f'{A}#{B}')
                    elif lo_s <= t <= hi_s:
                        undecided.append(f'{A}#{B} (u in [{lo_s},{hi_s}])')
        if found_down:
            verdict = f'BJ-HOLDS (change to {found_down[0]} realises u-1)'
        elif unresolved:
            verdict = 'UNDECIDED (unresolved identifications)'
        elif undecided:
            verdict = f'UNDECIDED (open partners: {sorted(set(undecided))})'
        else:
            verdict = 'COUNTEREXAMPLE: no minimal-diagram change realises u-1'
        print(f'  VERDICT: {verdict}')


if __name__ == '__main__':
    main()
