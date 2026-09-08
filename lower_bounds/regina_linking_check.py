#!/usr/bin/env python
"""Third-method re-verification of the Lickorish u=1 obstructions (the "815 list").

Method (independent of the Goeritz/Tait-graph and Seifert-matrix computations):
  KnotInfo PD -> spherogram Link -> SnapPy exterior -> order-two meridional orbifold filling -> cyclic double cover
  = double branched cover Sigma_2(K) -> Regina Triangulation3 -> HomologicalData
  -> torsion rank vector (is H1 cyclic of order det?) and odd-prime Legendre symbol vector
     (Kawauchi–Kojima invariants of the torsion linking pairing).

Lickorish: u(K)=1 => H1(Sigma_2 K)=Z/D and linking pairing ≅ <±2/D>.
For cyclic H1 the Legendre symbols (u_p | p), where the p-part of <a/D> is <a*m_p/p^k>, m_p=D/p^k,
classify the form completely, so we compare Regina's vector with the vectors predicted by a=+2 and a=-2.
The overall orientation of Sigma_2 is a convention (multiplies the form by -1), which only swaps the two
predictions, so the obstruction verdict (neither matches) is convention-free.

Usage:  ~/.pyenv/versions/unknot-venv/bin/python regina_linking_check.py [--workers 6] [--out FILE]
Writes a JSON with per-knot records and prints a comparison with the Goeritz-based verdicts.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def legendre(u: int, p: int) -> int:
    u %= p
    if u == 0:
        return 0
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def factor(n: int) -> dict[int, int]:
    if n < 1:
        raise ValueError('factorisation requires a positive integer')
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def parse_vector(s: str) -> dict[int, list[int]]:
    """'3(1) 5(0 1)' -> {3:[1], 5:[0,1]}"""
    out = {}
    for m in re.finditer(r"(\d+)\(([^)]*)\)", s):
        out[int(m.group(1))] = [int(x) for x in m.group(2).split()]
    return out


def predicted_legendre(a: int, D: int) -> dict[int, int]:
    """Legendre symbol vector of the cyclic form <a/D> (odd D)."""
    if D < 1 or D % 2 == 0 or math.gcd(a, D) != 1:
        raise ValueError('an odd-order nonsingular cyclic pairing is required')
    out = {}
    for p, k in factor(D).items():
        m = D // p ** k
        out[p] = legendre(a * m, p)
    return out


def worker(name: str) -> dict:
    import regina
    from tait_tools import link_from_knotinfo, knotinfo_rows, goeritz_matrix, generator_self_linking, h1_torsion, sigma2_snappy
    import snappy  # noqa: F401

    r = knotinfo_rows()[name]
    D = int(r['determinant'])
    rec = {'knot': name, 'interval': r['unknotting_number'].replace(' ', ''), 'det': D}
    try:
        L = link_from_knotinfo(name)
        # --- Goeritz (Tait graph) side ---
        G, _ = goeritz_matrix(L)
        if abs(int(G.det())) != D:
            raise ValueError('Goeritz and recorded determinants disagree')
        tors = h1_torsion(G)
        rec['goeritz_h1'] = tors
        a = generator_self_linking(G, D) if tors == [D] else None
        rec['goeritz_a'] = a
        if D == 1:
            g_allowed = {1, -1}
        elif tors != [D] or a is None:
            g_allowed = set() if tors != [D] else None
        else:
            sq = {(k * k) % D for k in range(1, D) if math.gcd(k, D) == 1}
            ia = pow(a, -1, D)
            g_allowed = {s for s in (1, -1) if (s * 2 * ia) % D in sq}
        rec['goeritz_allowed'] = sorted(g_allowed) if g_allowed is not None else None
        # --- Regina side ---
        if D == 1:
            rec['regina_h1'] = '0'
            rec['regina_allowed'] = [1, -1]
            rec['regina_cyclic'] = True
            rec['legendre_match_goeritz'] = True
            return rec
        # The meridional orbifold filling is lifted with the cover; the
        # chosen closed manifold is not inferred from its homology order.
        N = sigma2_snappy(L)
        rec['slope'] = [int(x) for x in N.cusp_info(0)['filling']]
        T = regina.Triangulation3(N.filled_triangulation()._to_string())
        T.simplify()
        H = regina.HomologicalData(T)
        h1 = str(T.homology())
        ranks = parse_vector(H.torsionRankVectorString())
        leg = parse_vector(H.torsionLegendreSymbolVectorString())
        rec['regina_h1'] = h1
        rec['regina_ranks'] = {str(p): v for p, v in ranks.items()}
        rec['regina_legendre_full'] = {str(p): v for p, v in leg.items()}
        # regina lists one symbol per level p^1, p^2, ...; for a cyclic p-part Z/p^k the relevant
        # Kawauchi–Kojima symbol is the k-th entry (lower levels are absent and reported as +1)
        rec['regina_legendre'] = {str(p): (v[ranks[p].index(1)] if p in ranks and 1 in ranks[p] and len(v) > ranks[p].index(1) else None)
                                  for p, v in leg.items()}
        fD = factor(D)
        cyclic = (set(ranks) == set(fD)) and all(
            sum(v) == 1 and v[fD[p] - 1] == 1 if len(v) >= fD[p] else False for p, v in ranks.items()
        )
        rec['regina_cyclic'] = cyclic
        # order check
        order = 1
        for p, v in ranks.items():
            for i, rk in enumerate(v):
                order *= (p ** (i + 1)) ** rk
        rec['regina_order_ok'] = (order == D)
        if not cyclic:
            rec['regina_allowed'] = []
        else:
            reg = {int(p): sym for p, sym in rec['regina_legendre'].items()}
            allowed = []
            for s in (1, -1):
                pred = predicted_legendre(s * 2, D)
                if all(reg.get(p) == pred[p] for p in pred):
                    allowed.append(s)
            rec['regina_allowed'] = allowed
            if a is not None:
                mine = predicted_legendre(a, D)
                mine_neg = predicted_legendre(-a, D)
                rec['legendre_match_goeritz'] = (all(reg.get(p) == mine[p] for p in mine)
                                                 or all(reg.get(p) == mine_neg[p] for p in mine_neg))
        return rec
    except Exception as e:  # keep going, record the failure
        rec['error'] = f'{type(e).__name__}: {e}'
        return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 2))
    ap.add_argument('--out', default=os.path.join(HERE, 'verify815_regina.json'))
    ap.add_argument('--only-815', action='store_true', help='only the 815 list, skip the u=1 validation set')
    args = ap.parse_args()

    from tait_tools import knotinfo_rows
    rows = knotinfo_rows()
    obs_path = os.path.join(HERE, '..', 'results', 'lickorish_cw_obstructed_815_knotinfo2026.8.1.json')
    obstructed = [n for n, _ in json.load(open(obs_path))]
    targets = [n for n, r in rows.items() if re.fullmatch(r'\[1,\d+\]', r['unknotting_number'].replace(' ', ''))]
    if args.only_815:
        targets = obstructed
    validation = [] if args.only_815 else [n for n, r in rows.items()
                                             if r['unknotting_number'].strip() == '1' and r['pd_notation']]
    names = list(dict.fromkeys(targets + validation))
    print(f'[1,b] targets: {len(targets)} (815-list: {len(obstructed)}), u=1 validation: {len(validation)}, total {len(names)}', flush=True)

    t0 = time.time()
    results = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(worker, n): n for n in names}
        for i, f in enumerate(as_completed(futs), 1):
            rec = f.result()
            results[rec['knot']] = rec
            if i % 250 == 0:
                print(f'  {i}/{len(names)} done, {time.time() - t0:.0f}s', flush=True)
    json.dump(results, open(args.out, 'w'), indent=1)

    # ---- summary ----
    obs_set = set(obstructed)
    err = [n for n, r in results.items() if 'error' in r]
    print(f'errors: {len(err)} {err[:10]}')
    # validation: u=1 knots must be allowed by regina test
    bad_val = [n for n in validation if n in results and 'error' not in results[n] and results[n].get('regina_allowed') == []]
    print(f'u=1 validation set: {len(validation)}; regina test falsely obstructs: {len(bad_val)} {bad_val[:10]}')
    # targets: compare verdicts
    cmp = {}
    for n in targets:
        r = results.get(n, {})
        if 'error' in r or r.get('regina_allowed') is None:
            key = 'error'
        else:
            g = r.get('goeritz_allowed')
            key = (('G-obs' if g == [] else 'G-allow'), ('R-obs' if r['regina_allowed'] == [] else 'R-allow'),
                   ('in815' if n in obs_set else 'not815'))
        cmp[key] = cmp.get(key, 0) + 1
    for k, v in sorted(cmp.items(), key=str):
        print(' ', k, v)
    h1_mismatch = [n for n in targets if results.get(n, {}).get('regina_order_ok') is False]
    leg_mismatch = [n for n in targets if results.get(n, {}).get('legendre_match_goeritz') is False]
    print(f'H1 order mismatch (regina vs det): {len(h1_mismatch)} {h1_mismatch[:10]}')
    print(f'Legendre vector mismatch (regina vs Goeritz form, up to global sign): {len(leg_mismatch)} {leg_mismatch[:10]}')
    reg_obs_targets = [n for n in targets if results.get(n, {}).get('regina_allowed') == []]
    print(f'Regina-based obstructions among targets: {len(reg_obs_targets)}; 815-list confirmed: '
          f'{sum(1 for n in obstructed if results.get(n, {}).get("regina_allowed") == [])}/{len(obstructed)}')
    print(f'elapsed {time.time() - t0:.0f}s; results -> {args.out}')


if __name__ == '__main__':
    main()
