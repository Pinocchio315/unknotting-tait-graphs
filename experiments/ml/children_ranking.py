#!/usr/bin/env python
"""Rank the child knots of the open alternating [2,3] knots by an estimated probability of u = 1.

Training data: every non-alternating knot with at most 13 crossings whose unknotting number is known to be 1 or 2
(reference table + our bounds, results/u_table_2026-09-07.json).  Features: KnotInfo invariants (signature,
determinant, genera, tau, s, nu, epsilon, algebraic unknotting number, Nakanishi index, bridge/braid indices,
volume, polynomial coefficients, Floer and Khovanov ranks, homology of the cyclic covers) and the Lickorish
linking-form test computed from the Seifert matrix.  Two models are fitted: one with every feature and one without
the features that are themselves lower bounds for u ("neutral" invariants only).  The script also counts the u = 2
knots on which every standard obstruction is silent, because only those tell the model anything about a knot whose
range is [1,2].
    python experiments/ml/children_ranking.py [--workers 8]
"""
import argparse, ast, json, math, os, re, sys, time, collections
from multiprocessing import Pool
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE)); RES = os.path.join(ROOT, 'results')
sys.path.insert(0, os.path.join(ROOT, 'lower_bounds'))
import database_knotinfo as dk
from seifert_linking_check import verdict as lickorish_verdict

def num(s, default=np.nan):
    s = str(s).strip()
    m = re.match(r'^-?\d+(\.\d+)?$', s)
    return float(s) if m else default
def rng_lo(s):
    m = re.findall(r'-?\d+', str(s)); return float(m[0]) if m else np.nan
def rng_hi(s):
    m = re.findall(r'-?\d+', str(s)); return float(m[-1]) if m else np.nan
def yn(s):
    s = str(s).strip().upper(); return 1.0 if s.startswith('Y') else (0.0 if s.startswith('N') else np.nan)
def vec(s):
    try: return ast.literal_eval(str(s))
    except Exception: return None

def features(r, cc):
    f = {}
    f['crossings'] = num(r['crossing_number']); f['det'] = num(r['determinant']); f['det_mod4'] = f['det'] % 4; f['det_mod8'] = f['det'] % 8
    f['sigma'] = num(r['signature']); f['abs_sigma'] = abs(f['sigma']); f['arf'] = num(r['arf_invariant'])
    f['genus'] = num(r['three_genus']); f['four_genus'] = rng_lo(r['smooth_four_genus']); f['four_genus_top'] = rng_lo(r['topological_four_genus'])
    f['s'] = num(r['rasmussen_invariant']); f['tau'] = num(r['ozsvath_szabo_tau_invariant']); f['nu_lo'] = rng_lo(r['nu']); f['epsilon'] = num(r['epsilon'])
    f['u_alg'] = rng_lo(r['unknotting_number_algebraic']); f['nakanishi'] = num(r['nakanishi_index'])
    f['bridge'] = num(r['bridge_index']); f['braid'] = num(r['braid_index']); f['braid_len'] = num(r['braid_length']); f['arc'] = num(r['arc_index'])
    f['tunnel'] = num(r['tunnel_number']); f['turaev'] = num(r['turaev_genus']); f['width'] = num(r['width']); f['crosscap'] = num(r['crosscap_number'])
    f['volume'] = num(r['volume']); f['cs'] = num(r['chern_simons_invariant']); f['dsg'] = num(r['double_slice_genus']); f['cg'] = num(r['smooth_concordance_genus'])
    for k in ('fibered', 'positive', 'quasipositive', 'strongly_quasipositive', 'almost_alternating', 'adequate', 'quasi_alternating', 'l_space', 'cosmetic_crossing'):
        f[k] = yn(r[k])
    sym = str(r['symmetry_type']).strip().lower()
    for k in ('chiral', 'reversible', 'fully amphicheiral', 'negative amphicheiral', 'positive amphicheiral'):
        f['sym_' + k.replace(' ', '_')] = 1.0 if sym == k else 0.0
    # cyclic covers
    tn = vec(r['torsion_numbers']) or []
    gens = {n: sum(1 for x in fac if x != 1) for n, fac in tn}
    for n in range(2, 10): f[f'g{n}'] = float(gens.get(n, np.nan))
    f['h1_cyclic'] = 1.0 if gens.get(2, 99) <= 1 else 0.0
    f['cyclic_bound'] = float(cc.get(r['name'], {}).get('bound', np.nan))
    # polynomials
    a = vec(r['alexander_polynomial_vector']) or []
    coeffs = a[2:] if len(a) > 2 else []
    f['alex_deg'] = float(len(coeffs) - 1) if coeffs else np.nan; f['alex_abs_sum'] = float(sum(abs(c) for c in coeffs)) if coeffs else np.nan
    f['alex_max'] = float(max(abs(c) for c in coeffs)) if coeffs else np.nan
    for i in range(6): f[f'alex_c{i}'] = float(coeffs[i]) if i < len(coeffs) else 0.0
    j = vec(r['jones_polynomial_vector']) or []
    jc = j[2:] if len(j) > 2 else []
    f['jones_min'] = float(j[0]) if j else np.nan; f['jones_span'] = float(len(jc) - 1) if jc else np.nan
    f['jones_abs_sum'] = float(sum(abs(c) for c in jc)) if jc else np.nan; f['jones_max'] = float(max(abs(c) for c in jc)) if jc else np.nan
    cw = vec(r['conway_polynomial_vector']) or []
    cwc = cw[2:] if len(cw) > 2 else []
    for i in range(4): f[f'conway_c{i}'] = float(cwc[i]) if i < len(cwc) else 0.0
    # knot Floer: total rank and rank at the top Alexander grading
    h = str(r['hfk_polynomial_vector']).strip('[]')
    try:
        terms = [tuple(int(x) for x in t.split(',')) for t in h.split(';') if t.strip()]
        f['hfk_rank'] = float(sum(t[0] for t in terms)); f['hfk_top'] = float(sum(t[0] for t in terms if t[2] == max(x[2] for x in terms)))
        f['hfk_thin'] = 1.0 if len({t[2] - t[1] for t in terms}) == 1 else 0.0
    except Exception:
        f['hfk_rank'] = f['hfk_top'] = f['hfk_thin'] = np.nan
    kv = vec(r['khovanov_reduced_rational_vector'])
    if kv:
        f['kh_rank'] = float(sum(t[1] for t in kv)); f['kh_width'] = float(len({t[3] - 2 * t[2] for t in kv}))
    else: f['kh_rank'] = f['kh_width'] = np.nan
    kz = vec(r['khovanov_reduced_integral_vector'])
    f['kh_torsion'] = 1.0 if (kz and any(t[0] != 0 for t in kz)) else 0.0
    return f

OBSTRUCTION_FEATURES = {'sigma', 'abs_sigma', 'four_genus', 'four_genus_top', 's', 'tau', 'nu_lo', 'epsilon', 'u_alg', 'nakanishi', 'cyclic_bound',
                        'h1_cyclic', 'lickorish_pass', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'dsg', 'cg'}

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--workers', type=int, default=8); args = ap.parse_args()
    t0 = time.time()
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    table = json.load(open(os.path.join(RES, 'u_table_2026-09-07.json')))
    cc = json.load(open(os.path.join(RES, 'cyclic_cover', 'cyclic_cover_bound_2026-09-07.json')))
    fin = json.load(open(os.path.join(RES, 'open23', 'priority_u23_final.json')))
    children = [d['child'] for d in fin['children'] if d['child'] in rows]
    nonalt = [nm for nm, r in rows.items() if str(r.get('alternating', '')).upper().startswith('N')]
    train = [nm for nm in nonalt if table[nm][0] == table[nm][1] and table[nm][0] in (1, 2)]
    names = sorted(set(train) | set(children))
    print(f'non-alternating knots: {len(nonalt)}; training knots with u in {{1,2}}: {len(train)} '
          f'(u=1: {sum(1 for nm in train if table[nm][0] == 1)}, u=2: {sum(1 for nm in train if table[nm][0] == 2)}); children: {len(children)}', flush=True)
    # Lickorish test from the Seifert matrices (multiprocessing)
    with Pool(args.workers) as pool:
        lick = dict(zip(names, pool.map(lickorish_verdict, [(nm, rows[nm]['seifert_matrix'], int(rows[nm]['determinant'])) for nm in names], chunksize=20)))
    print(f'Lickorish test done in {time.time() - t0:.0f}s', flush=True)
    F = {}
    for nm in names:
        f = features(rows[nm], cc); f['lickorish_pass'] = 0.0 if lick[nm].get('obstructed') else 1.0; F[nm] = f
    keys = sorted(F[names[0]].keys())
    X = lambda L, ks: np.array([[F[nm].get(k, np.nan) for k in ks] for nm in L], dtype=float)
    y = np.array([1 if table[nm][0] == 1 else 0 for nm in train])
    # which u = 2 training knots are silent on every standard obstruction?
    def silent(nm):
        f = F[nm]
        return (f['abs_sigma'] <= 2 and f['lickorish_pass'] == 1 and f['h1_cyclic'] == 1 and (np.isnan(f['nakanishi']) or f['nakanishi'] <= 1)
                and (np.isnan(f['u_alg']) or f['u_alg'] <= 1) and (np.isnan(f['four_genus']) or f['four_genus'] <= 1) and (np.isnan(f['tau']) or abs(f['tau']) <= 1)
                and (np.isnan(f['s']) or abs(f['s']) <= 2) and (np.isnan(f['cyclic_bound']) or f['cyclic_bound'] <= 1))
    u2 = [nm for nm in train if table[nm][0] == 2]; u1 = [nm for nm in train if table[nm][0] == 1]
    silent_u2 = [nm for nm in u2 if silent(nm)]; silent_u1 = [nm for nm in u1 if silent(nm)]
    print(f'u=2 knots on which every standard obstruction is silent: {len(silent_u2)} of {len(u2)}: {silent_u2[:30]}')
    print(f'u=1 knots silent (sanity, should be all): {len(silent_u1)} of {len(u1)}; children silent: {sum(1 for nm in children if silent(nm))} of {len(children)}')
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    from sklearn.inspection import permutation_importance
    out = {'children': {}, 'models': {}}
    for label, ks in (('all features', keys), ('neutral features only', [k for k in keys if k not in OBSTRUCTION_FEATURES])):
        Xt = X(train, ks); pred = np.zeros(len(train))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=0).split(Xt, y):
            m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05).fit(Xt[tr], y[tr]); pred[te] = m.predict_proba(Xt[te])[:, 1]
        auc = roc_auc_score(y, pred)
        idx = {nm: i for i, nm in enumerate(train)}
        sub = [idx[nm] for nm in u1] + [idx[nm] for nm in silent_u2]
        auc_silent = roc_auc_score(y[sub], pred[sub]) if silent_u2 else float('nan')
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05).fit(Xt, y)
        pi = permutation_importance(m, Xt, y, n_repeats=3, random_state=0, scoring='roc_auc')
        order = np.argsort(-pi.importances_mean)[:10]
        pc = m.predict_proba(X(children, ks))[:, 1]
        print(f'\n### {label}: {len(ks)} features; cross-validated AUC {auc:.3f}; AUC restricted to u=1 vs the silent u=2 knots {auc_silent:.3f}')
        print('   importance:', ', '.join(f'{ks[i]}={pi.importances_mean[i]:.3f}' for i in order))
        print(f'   children: mean P(u=1) {pc.mean():.3f}, min {pc.min():.3f}, max {pc.max():.3f}')
        out['models'][label] = {'features': ks, 'cv_auc': auc, 'auc_silent_subset': auc_silent, 'silent_u2': silent_u2,
                                'importance': {ks[i]: float(pi.importances_mean[i]) for i in order}}
        for nm, p in zip(children, pc): out['children'].setdefault(nm, {})[label] = float(p)
    parents = {d['child']: d['parents'] for d in fin['children']}
    rank = sorted(children, key=lambda nm: -out['children'][nm]['neutral features only'])
    print('\nchildren ranked by P(u=1) from the neutral-feature model (all-feature model in brackets), with the number of parents:')
    for nm in rank[:25]:
        print(f"   {nm:9s} {out['children'][nm]['neutral features only']:.3f} [{out['children'][nm]['all features']:.3f}]  parents {len(parents[nm]):2d}  sigma {int(F[nm]['sigma'])} det {int(F[nm]['det'])}")
    out['ranking'] = [{'child': nm, 'p_u1_neutral': out['children'][nm]['neutral features only'], 'p_u1_all': out['children'][nm]['all features'], 'n_parents': len(parents[nm])} for nm in rank]
    json.dump(out, open(os.path.join(RES, 'open23', 'children_ranking_2026-09-07.json'), 'w'), indent=1)
    print(f'\nwritten results/open23/children_ranking_2026-09-07.json ({time.time() - t0:.0f}s)')

if __name__ == '__main__':
    main()
