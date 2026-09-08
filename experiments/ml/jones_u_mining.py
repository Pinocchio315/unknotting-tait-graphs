"""Stage 1: classifiers + permutation importance for u = 1 vs u >= 2.  Stage 2: mining of exact
residue conditions that hold for every knot with u = 1 but fail for some knots with u >= 2."""
import re, json, math, itertools, collections, sys, time
import numpy as np
import database_knotinfo as dk
S = sys.argv[1]
rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
def rng(u):
    m = re.findall(r'\d+', str(u)); return (int(m[0]), int(m[-1])) if m else None
def parse_poly(s, var='t'):
    """KnotInfo polynomial string -> {exp: int coeff}, via sympy (handles t^(-3), spacing, unary minus)"""
    import sympy
    t = sympy.symbols('t')
    P = sympy.Poly(sympy.expand(sympy.sympify(s.replace('^', '**')) * t ** 60), t)
    return {int(e[0]) - 60: int(c) for e, c in P.as_dict().items()}
def cyclotomic(r):
    """coefficients (low->high) of Phi_r(x) as ints"""
    import sympy
    x = sympy.symbols('x'); P = sympy.Poly(sympy.cyclotomic_poly(r, x), x)
    return [int(c) for c in reversed(P.all_coeffs())]
PHI = {r: cyclotomic(r) for r in (5, 8, 12, 16, 20, 24)}
def reduce_mod_phi(coef, r):
    """integer coordinates of V(zeta_r) in the basis 1, zeta, ..., zeta^(deg-1)"""
    phi = PHI[r]; deg = len(phi) - 1
    lo = min(coef); hi = max(coef)
    # powers zeta^e for e in [lo, hi] as vectors of length deg; zeta^-1 = zeta^(r-1)
    vec = np.zeros(deg, dtype=np.int64)
    for e, c in coef.items():
        k = e % r
        # zeta^k reduced mod phi: compute by repeated multiplication (k < r small)
        p = np.zeros(deg, dtype=np.int64); p[0] = 1
        for _ in range(k):
            q = np.zeros(deg + 1, dtype=np.int64); q[1:] = p
            if q[deg]:
                top = q[deg]; q = q[:deg] - top * np.array(phi[:deg], dtype=np.int64)
            else:
                q = q[:deg]
            p = q
        vec += c * p
    return vec
# ---------------- dataset ----------------
data = []
for r in rows:
    x = rng(r['unknotting_number'])
    try:
        J = parse_poly(r['jones_polynomial']); A = parse_poly(r['alexander_polynomial'])
    except Exception:
        continue
    rec = {'name': r['name'], 'u_lo': x[0] if x else None, 'u_hi': x[1] if x else None,
           'exact': bool(x and x[0] == x[1]), 'crossings': int(r['crossing_number']), 'alternating': str(r.get('alternating','')).upper().startswith('Y')}
    def intval(k):
        v = str(r.get(k, '')).strip()
        try: return int(v)
        except Exception: return None
    rec['sigma'] = intval('signature'); rec['det'] = intval('determinant'); rec['arf'] = intval('arf_invariant')
    rec['tau'] = intval('ozsvath_szabo_tau'); rec['s'] = intval('rasmussen_invariant'); rec['g3'] = intval('three_genus')
    rec['braid_index'] = intval('braid_index'); rec['bridge'] = intval('bridge_index')
    # Jones features
    lo, hi = min(J), max(J)
    f = {'span': hi - lo, 'deg_lo': lo, 'deg_hi': hi, 'nterms': len(J), 'l1': sum(abs(c) for c in J.values()),
         'coef_lo': J[lo], 'coef_hi': J[hi], 'max_abs_coef': max(abs(c) for c in J.values())}
    f['V_m1'] = sum(c * (-1) ** (e % 2) for e, c in J.items())
    f['dV_m1'] = sum(c * e * (-1) ** ((e - 1) % 2) for e, c in J.items())
    f['d2V_m1'] = sum(c * e * (e - 1) * (-1) ** ((e - 2) % 2) for e, c in J.items())
    f['d3V_m1'] = sum(c * e * (e - 1) * (e - 2) * (-1) ** ((e - 3) % 2) for e, c in J.items())
    f['dV_1'] = sum(c * e for e, c in J.items()); f['d2V_1'] = sum(c * e * (e - 1) for e, c in J.items())
    f['d3V_1'] = sum(c * e * (e - 1) * (e - 2) for e, c in J.items())
    for rr in (5, 8, 12, 16, 20, 24):
        v = reduce_mod_phi(J, rr)
        for i, c in enumerate(v): f[f'z{rr}_{i}'] = int(c)
    # e^{i pi/3}: eps and d (via t^2 - t + 1)
    a = b = 0
    for e, c in J.items():
        A_, B_ = 1, 0
        if e >= 0:
            for _ in range(e % 6): A_, B_ = -B_, A_ + B_
        else:
            for _ in range((-e) % 6): A_, B_ = A_ + B_, -A_
        a += c * A_; b += c * B_
    f['w6_a'] = a; f['w6_b'] = b
    # Alexander: a2 = Delta''(1)/2 with Delta symmetrised (KnotInfo gives Conway-normalised?), use second derivative at 1 of the symmetric version
    lo2, hi2 = min(A), max(A); mid = (lo2 + hi2) / 2
    Asym = {e - int(mid) if (lo2 + hi2) % 2 == 0 else e: c for e, c in A.items()}
    f['alex_d2_1'] = sum(c * e * (e - 1) for e, c in A.items()); f['alex_m1'] = sum(c * (-1) ** (e % 2) for e, c in A.items())
    f['alex_span'] = hi2 - lo2; f['alex_l1'] = sum(abs(c) for c in A.values())
    rec['f'] = f
    data.append(rec)
json.dump(data, open(S + '/dataset.json', 'w'))
ex = [d for d in data if d['exact']]
print('knots parsed:', len(data), '| exact u:', len(ex), '| u=1:', sum(1 for d in ex if d['u_lo'] == 1), flush=True)
# ---------------- stage 1: classifiers ----------------
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, cross_val_predict, StratifiedKFold
from sklearn.inspection import permutation_importance
fnames = sorted(ex[0]['f'].keys())
def mat(recs, cols):
    return np.array([[float(d['f'][c]) if c in d['f'] else float(d[c] if d[c] is not None else 0) for c in cols] for d in recs])
y = np.array([1 if d['u_lo'] == 1 else 0 for d in ex])
fnames_j = [c for c in fnames if not c.startswith('alex')]
groups = {
  'Jones only': fnames_j,
  'Jones + sigma + det': fnames_j + ['sigma', 'det'],
  'Jones + Alexander + sigma + det + Arf': fnames + ['sigma', 'det', 'arf'],
  'all cheap invariants (+tau, s, g3, braid, bridge)': fnames + ['sigma', 'det', 'arf', 'tau', 's', 'g3', 'braid_index', 'bridge'],
  'non-Jones invariants only (sigma, det, Arf, Alexander, tau, s, g3, braid, bridge)': [c for c in fnames if c.startswith('alex')] + ['sigma', 'det', 'arf', 'tau', 's', 'g3', 'braid_index', 'bridge'],
}
cv = StratifiedKFold(5, shuffle=True, random_state=0)
res = {}
for label, cols in groups.items():
    X = mat(ex, cols)
    clf = GradientBoostingClassifier(n_estimators=400, max_depth=3, learning_rate=0.05, random_state=0)
    auc = cross_val_score(clf, X, y, cv=cv, scoring='roc_auc').mean(); acc = cross_val_score(clf, X, y, cv=cv, scoring='accuracy').mean()
    res[label] = (auc, acc); print(f'[stage 1] {label}: AUC {auc:.3f} acc {acc:.3f}', flush=True)
cols = groups['all cheap invariants (+tau, s, g3, braid, bridge)']; X = mat(ex, cols)
clf = GradientBoostingClassifier(n_estimators=400, max_depth=3, learning_rate=0.05, random_state=0).fit(X, y)
pi = permutation_importance(clf, X, y, n_repeats=5, random_state=0, scoring='roc_auc')
order = np.argsort(-pi.importances_mean)[:20]
print('[stage 1] permutation importance (AUC drop), top 20:', flush=True)
for i in order: print(f'   {cols[i]:14s} {pi.importances_mean[i]:.4f}', flush=True)
# ---------------- stage 2: exact residue conditions ----------------
u1 = [d for d in ex if d['u_lo'] == 1]; u2 = [d for d in ex if d['u_lo'] >= 2]
open_rows = [d for d in data if not d['exact'] and d['u_lo'] == 1]
ours815 = set(nm for nm, iv in json.load(open('/Users/pinocchio/Documents/05_unknotting_number/0001_writing_with_pavel/code/results/lickorish_cw_obstructed_815_knotinfo2026.8.1.json')))
intfeats = [c for c in fnames] + ['sigma', 'det', 'arf']
def val(d, c): return d['f'][c] if c in d['f'] else d[c]
found = []
for c in intfeats:
    vals1 = [val(d, c) for d in u1]; vals2 = [val(d, c) for d in u2]
    if any(v is None for v in vals1 + vals2): continue
    for m in range(2, 13):
        R1 = {v % m for v in vals1}
        if len(R1) == m: continue
        excl = sum(1 for v in vals2 if v % m not in R1)
        if excl == 0: continue
        excl_open = sum(1 for d in open_rows if val(d, c) is not None and val(d, c) % m not in R1)
        new_open = sum(1 for d in open_rows if val(d, c) is not None and val(d, c) % m not in R1 and d['name'] not in ours815)
        found.append({'feature': c, 'mod': m, 'allowed_residues_u1': sorted(R1), 'excludes_known_u>=2': excl, 'of': len(u2), 'settles_open_[1,x]_rows': excl_open, 'beyond_lickorish_cw': new_open})
# pairs of features with small moduli
pairs = [('sigma', 'w6_a'), ('sigma', 'dV_m1'), ('det', 'dV_m1'), ('sigma', 'det'), ('arf', 'dV_m1'), ('sigma', 'd2V_m1'), ('det', 'd2V_m1'), ('sigma', 'alex_d2_1'), ('det', 'alex_d2_1'), ('dV_m1', 'alex_d2_1')]
for c1, c2 in pairs:
    for m1 in (2, 3, 4, 8):
        for m2 in (2, 3, 4, 8, 16):
            R1 = {(val(d, c1) % m1, val(d, c2) % m2) for d in u1 if val(d, c1) is not None and val(d, c2) is not None}
            if len(R1) == m1 * m2: continue
            excl = sum(1 for d in u2 if val(d, c1) is not None and val(d, c2) is not None and (val(d, c1) % m1, val(d, c2) % m2) not in R1)
            if excl == 0: continue
            excl_open = sum(1 for d in open_rows if val(d, c1) is not None and val(d, c2) is not None and (val(d, c1) % m1, val(d, c2) % m2) not in R1)
            new_open = sum(1 for d in open_rows if val(d, c1) is not None and val(d, c2) is not None and (val(d, c1) % m1, val(d, c2) % m2) not in R1 and d['name'] not in ours815)
            found.append({'feature': f'({c1} mod {m1}, {c2} mod {m2})', 'mod': None, 'allowed_residues_u1': sorted(R1), 'excludes_known_u>=2': excl, 'of': len(u2), 'settles_open_[1,x]_rows': excl_open, 'beyond_lickorish_cw': new_open})
found.sort(key=lambda r: (-r['beyond_lickorish_cw'], -r['excludes_known_u>=2']))
json.dump(found, open(S + '/conditions.json', 'w'), indent=1)
print(f'[stage 2] candidate conditions with 0 violations on u=1 knots and >0 exclusions among known u>=2: {len(found)}', flush=True)
for r in found[:25]:
    print(f"   {r['feature']:28s} mod {str(r['mod']):4s} allowed {str(r['allowed_residues_u1'])[:40]:40s} excludes {r['excludes_known_u>=2']:4d}/{r['of']} known u>=2 | settles {r['settles_open_[1,x]_rows']:4d} open rows, {r['beyond_lickorish_cw']:4d} beyond Lickorish/CW", flush=True)
