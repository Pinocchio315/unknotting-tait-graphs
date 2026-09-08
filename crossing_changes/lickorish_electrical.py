#!/usr/bin/env python
"""Lickorish's condition for the knot K_c obtained by one crossing change, read off the Tait graph of the
alternating diagram.  With G the positive-definite Goeritz form and x = v_i - v_j the edge vector of the
crossing c, the Goeritz form of the changed diagram is G' = G - 2 x x^T, so H_1(Sigma_2(K_c)) = coker G' and
the linking form is +- G'^{-1}.  By Sherman-Morrison, lambda([x],[x]) = x^T G'^{-1} x = R/(1-2R) = tau_e/(tau - 2 tau_e)
(mod 1), where tau = det K counts the spanning trees and tau_e those containing e.  Lickorish: u(K_c) = 1 forces
coker G' cyclic with a generator of self-linking +-2 k^2 / det K_c.  This script evaluates the condition for every
crossing change of the data set (results/crossing_changes/dataset_v2.json.gz) and records, per row, whether
coker G' is cyclic, the self-linking of [x], whether [x] generates, and the Lickorish verdict.
    python lickorish_electrical.py [--workers 8]
"""
import argparse, ast, gzip, itertools, json, math, os, sys, time
from fractions import Fraction
from multiprocessing import Pool
import numpy as np
from paths import *
sys.path.insert(0, os.path.join(ROOT, 'lower_bounds'))
import database_knotinfo as dk
from xtait.graph import from_pd
from seifert_linking_check import h1_torsion, generator_self_linking

def goeritz(g):
    V = sorted(g.vertices); idx = {v: i for i, v in enumerate(V)}; n = len(V)
    L = np.zeros((n, n), dtype=np.int64)
    for e in g.edges:
        u, w = g.edges[e]; L[idx[u], idx[u]] += 1; L[idx[w], idx[w]] += 1; L[idx[u], idx[w]] -= 1; L[idx[w], idx[u]] -= 1
    return L[1:, 1:], {v: idx[v] - 1 for v in V}   # delete the first vertex; index -1 = deleted

def work(args):
    nm, pd, crossings = args
    from sympy import Matrix
    g = from_pd([list(q) for q in pd]); edges = sorted(g.edges); G, idx = goeritz(g)
    n = len(G); out = []
    for i in crossings:
        u, w = g.edges[edges[i]]; x = np.zeros(n, dtype=np.int64)
        if idx[u] >= 0: x[idx[u]] += 1
        if idx[w] >= 0: x[idx[w]] -= 1
        Gp = G - 2 * np.outer(x, x)
        M = Matrix(Gp.tolist()); Dp = abs(int(M.det()))
        rec = {'knot': nm, 'crossing': i, 'det_changed_check': Dp}
        if Dp == 1:
            rec.update(cyclic=True, self_linking_x='0', x_generates=True, lickorish='pass'); out.append(rec); continue
        tors = h1_torsion(M); rec['h1'] = tors; rec['cyclic'] = (tors == [Dp])
        lam_x = Fraction(int(x @ (np.array(M.adjugate().tolist(), dtype=object) @ x)), int(M.det())) % 1
        rec['self_linking_x'] = str(lam_x); rec['x_generates'] = rec['cyclic'] and lam_x.denominator == Dp
        if not rec['cyclic']: rec['lickorish'] = 'fail:not cyclic'
        else:
            a = lam_x.numerator if rec['x_generates'] else generator_self_linking(M, Dp)
            if a is None: rec['lickorish'] = 'unknown'
            else:
                squares = {(k * k) % Dp for k in range(1, Dp) if math.gcd(k, Dp) == 1}; ia = pow(int(a), -1, Dp)
                rec['lickorish'] = 'pass' if any((s * 2 * ia) % Dp in squares for s in (1, -1)) else 'fail:self-linking'
        out.append(rec)
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=8); args = ap.parse_args()
    t0 = time.time()
    p = os.path.join(CC, 'dataset_v2.json'); rows = json.load(open(p)) if os.path.exists(p) else json.load(gzip.open(p + '.gz', 'rt'))
    ki = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    byk = {}
    for r in rows: byk.setdefault(r['knot'], []).append(r['crossing'])
    tasks = [(nm, ast.literal_eval(ki[nm]['pd_notation']), sorted(cs)) for nm, cs in byk.items()]
    with Pool(args.workers) as pool: res = [x for chunk in pool.imap_unordered(work, tasks, chunksize=8) for x in chunk]
    L = {(r['knot'], r['crossing']): r for r in res}
    bad = sum(1 for r in rows if L[(r['knot'], r['crossing'])]['det_changed_check'] != r['det_changed']); print('det(K_c) mismatches:', bad)
    stat = {}
    for r in rows:
        l = L[(r['knot'], r['crossing'])]; key = (r['label'], l['lickorish'].split(':')[0])
        stat[key] = stat.get(key, 0) + 1
    print('label x Lickorish verdict of K_c from the Tait graph:', {f'{a}/{b}': c for (a, b), c in sorted(stat.items(), key=str)})
    u1 = [r for r in rows if r.get('u_result') == 1]; print('rows with u(K_c) = 1:', len(u1), 'Lickorish pass:', sum(1 for r in u1 if L[(r['knot'], r['crossing'])]['lickorish'] == 'pass'), '(must be all)')
    xg = sum(1 for l in res if l.get('x_generates')); print('[x] generates coker G\' in', xg, 'of', len(res), 'rows; cyclic in', sum(1 for l in res if l.get('cyclic')), 'rows')
    json.dump(res, open(os.path.join(CC, 'lickorish_electrical_2026-09-07.json'), 'w'))
    print(f'wrote results/crossing_changes/lickorish_electrical_2026-09-07.json ({time.time() - t0:.0f}s)')

if __name__ == '__main__':
    main()
