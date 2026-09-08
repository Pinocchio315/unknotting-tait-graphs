#!/usr/bin/env python
"""Pin the correction terms of Sigma_2(K) for a Khovanov-thin knot K by intersecting Greene's solitary-state pages over
all markings and colourings of the diagram.  For every marking, d(t) is one of the gradings of the solitary states in
the class t (the E_infinity page is a subquotient of E_1, and Sigma_2 is an L-space); the classes of two markings are
matched by the unique group automorphism k -> a k of H^2 = Z/D that preserves the self-linking of c_1 and is
consistent on every class; conjugation symmetry d(t) = d(-t) is imposed as well.
    python pin_dinv.py 12n_491 [8_20 ...]
"""
import ast, math, itertools, json, sys, os, database_knotinfo as dk
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'montesinos'))
from greene_dinv import Diagram, u1_admissible

def pages_of(pd):
    out = {}
    for black in (0, 1):
        for mark in range(1, 2 * len(pd) + 1):
            try: Dg = Diagram(pd, black, mark)
            except AssertionError: continue
            X = [Dg.analyse(x) for x in Dg.states()]; classes = {}
            for x in X: classes.setdefault(x['label'], []).append(x)
            idx, lam = Dg.spinc_index(list(classes))
            if idx is None or len(idx) != Dg.D: continue
            out[(black, mark)] = ({idx[lab]: sorted(set(x['gr'] for x in xs if x['solitary'])) for lab, xs in classes.items()}, lam)
    return out

def pin(name, pd, D, verbose=True):
    pages = pages_of(pd); ref = min(pages, key=lambda k: sum(1 for v in pages[k][0].values() if len(v) > 1))
    candR, lamR = pages[ref]; cur = {k: set(v) for k, v in candR.items()}
    for k in range(D): cur[k] = cur[k] & cur[(-k) % D]
    units = [a for a in range(1, D) if math.gcd(a, D) == 1]; used = 0
    for key, (candM, lamM) in pages.items():
        if key == ref: continue
        good = set()
        for a in units:
            if (a * a * lamM - lamR) % 1 != 0: continue                  # k_M = a k_R preserves the self-linking of c_1
            if all(cur[k] & set(candM[(a * k) % D]) for k in range(D)): good.add(min(a, D - a))
        if len(good) == 1:
            a = good.pop(); used += 1
            for k in range(D): cur[k] &= set(candM[(a * k) % D]) | set(candM[(-a * k) % D])
            for k in range(D): cur[k] = cur[k] & cur[(-k) % D]
        elif verbose: print('   ', key, 'ambiguous alignment', sorted(good))
    amb = {k: sorted(v) for k, v in cur.items() if len(v) != 1}
    if verbose: print(f'{name}: det {D}, {len(pages)} pages, reference {ref} (lam {lamR}), {used} pages used; unresolved classes: {len(amb)}', {k: [str(x) for x in v] for k, v in list(amb.items())[:6]})
    return cur, lamR, ref

if __name__ == '__main__':
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    for NAME in sys.argv[1:]:
        pd = ast.literal_eval(rows[NAME]['pd_notation']); D = int(rows[NAME]['determinant'])
        cur, lam, ref = pin(NAME, pd, D)
        if all(len(v) == 1 for v in cur.values()):
            d = {k: next(iter(v)) for k, v in cur.items()}; vals = sorted(d.values())
            print('   d-invariants:', [str(x) for x in vals], '| sum', sum(vals), '| d(spin) =', d[0])
            fits = u1_admissible(d, lam, D); print('   u = 1 surgery test:', ('PASS ' + str(fits[:3])) if fits else 'OBSTRUCTED')
            json.dump({'name': NAME, 'det': D, 'lam': str(lam), 'reference_marking': ref, 'd': {str(k): str(v) for k, v in d.items()}, 'u1_fits': fits}, open(f'greene_d_{NAME}.json', 'w'), indent=1)
        else:
            amb = {k: sorted(v) for k, v in cur.items() if len(v) != 1}; keys = sorted(k for k in amb if k <= (-k) % D)
            base = {k: next(iter(v)) for k, v in cur.items() if len(v) == 1}; n_pass = n_tot = 0
            for choice in itertools.product(*[amb[k] for k in keys]):
                dd = dict(base)
                for k, g in zip(keys, choice): dd[k] = g; dd[(-k) % D] = g
                n_tot += 1; n_pass += bool(u1_admissible(dd, lam, D))
            print(f'   {len(keys)} free classes: {n_tot} candidate vectors, {n_pass} admit the u = 1 surgery structure')
            json.dump({'name': NAME, 'det': D, 'lam': str(lam), 'reference_marking': ref, 'pinned': {str(k): str(v) for k, v in base.items()},
                       'ambiguous': {str(k): [str(x) for x in v] for k, v in amb.items()}, 'candidate_vectors': n_tot, 'admitting_u1': n_pass}, open(f'greene_d_{NAME}.json', 'w'), indent=1)
