from paths import *
import json, collections, sys
from fractions import Fraction
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.inspection import permutation_importance
import gzip, os
_p = os.path.join(CC, 'dataset_v2.json')
R = json.load(open(_p)) if os.path.exists(_p) else json.load(gzip.open(_p + '.gz', 'rt'))
for r in R:
    r['Rb'] = Fraction(r['R_black']); r['Rf'] = float(r['Rb'])
    r['dsig'] = abs(r['sigma_c']) - abs(r['sigma']); r['sharp'] = (2 * r['u'] == abs(r['sigma'])); r['ratio'] = r['det_changed'] / r['det']
    r['sig_drop'] = int(r['dsig'] == -2); r['sig_rel'] = r['sign'] * (1 if r['sigma'] > 0 else -1 if r['sigma'] < 0 else 0)
    r['R_gt_half'] = int(r['Rb'] > Fraction(1, 2))
FEATS = ['sign', 'sig_rel', 'dsig', 'sig_drop', 'R_gt_half', 'Rf', 'ratio', 'deg_min', 'deg_max', 'face_min', 'face_max', 'twist_par', 'twist_ser', 'twist_max_diagram', 'reduced_crossings', 'crossings', 'det', 'det_changed', 'u', 'sigma', 'n_black', 'n_white']
def X_of(rows): return np.array([[float(r[f]) for f in FEATS] for r in rows])
def run(rows, name, feats=FEATS):
    rows = [r for r in rows if r['label'] is not None]
    X = np.array([[float(r[f]) for f in feats] for r in rows]); y = np.array([r['label'] == 'good' for r in rows]); g = np.array([r['knot'] for r in rows])
    pred = np.zeros(len(rows))
    for tr, te in GroupKFold(5).split(X, y, g):
        m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05).fit(X[tr], y[tr]); pred[te] = m.predict_proba(X[te])[:, 1]
    auc = roc_auc_score(y, pred)
    idx = collections.defaultdict(list)
    for i, r in enumerate(rows): idx[r['knot']].append(i)
    def hit(score):
        h = n = 0
        for nm, ii in idx.items():
            if not any(rows[i]['label'] == 'good' for i in ii): continue
            n += 1; best = max(ii, key=score); h += rows[best]['label'] == 'good'
        return h, n
    print(f'\n### {name}: {len(rows)} rows, {len(idx)} knots, good={int(y.sum())}, AUC={auc:.3f}')
    for lab, sc in [('GBM', lambda i: pred[i]), ('min det ratio', lambda i: -rows[i]['ratio']), ('min reduced crossings', lambda i: -rows[i]['reduced_crossings']), ('sig-drop then min ratio', lambda i: (rows[i]['dsig'] == -2, -rows[i]['ratio'])), ('first crossing', lambda i: -i)]:
        h, n = hit(sc); print(f'   hit@1 {lab:24s} {h}/{n} = {h / n:.3f}')
    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05).fit(X, y)
    pi = permutation_importance(m, X, y, n_repeats=3, random_state=0, scoring='roc_auc')
    order = np.argsort(-pi.importances_mean)[:8]
    print('   permutation importance:', ', '.join(f'{feats[i]}={pi.importances_mean[i]:.3f}' for i in order))
    return pred, rows
print('labels:', collections.Counter(r['label'] for r in R))
print('composite labels used (sum-blocks rigorous, sum-jones heuristic):', collections.Counter(r['result_how'] for r in R if r['label'] and str(r['result_how']).startswith('sum')))
# theorem cross-tab
c = collections.Counter((r['sign'], r['R_gt_half'], r['sigma_c'] - r['sigma']) for r in R)
print('\n(sign, R>1/2, sigma_c - sigma):', dict(c))
# per-knot candidate reduction for sharp knots
byk = collections.defaultdict(list)
for r in R: byk[r['knot']].append(r)
red = []
for nm, rs in byk.items():
    if rs[0]['sharp']: red.append((len(rs), sum(r['sig_drop'] for r in rs)))
print(f'sharp knots: {len(red)}; mean crossings {np.mean([a for a, b in red]):.2f}, mean sigma-dropping crossings {np.mean([b for a, b in red]):.2f}')
run(R, 'ALL knots, all crossings (good vs not)')
run([r for r in R if r['sharp']], 'sharp knots, all crossings')
run([r for r in R if r['sharp'] and r['dsig'] == -2], 'sharp knots, sigma-dropping crossings only')
run([r for r in R if not r['sharp']], 'non-sharp knots, all crossings')
run([r for r in R if not r['sharp'] and r['u'] == 2 and r['sigma'] == 0], 'u=2, sigma=0 knots')
run([r for r in R if not r['sharp'] and r['u'] == 2 and abs(r['sigma']) == 2], 'u=2, |sigma|=2 knots')
# local-only features (no K_c information): can we predict without changing the crossing?
LOCAL = ['sign', 'sig_rel', 'sig_drop', 'R_gt_half', 'Rf', 'deg_min', 'deg_max', 'face_min', 'face_max', 'twist_par', 'twist_ser', 'twist_max_diagram', 'crossings', 'det', 'u', 'sigma', 'n_black', 'n_white']
run(R, 'ALL knots, LOCAL features only (no K_c data)', LOCAL)
run([r for r in R if r['sharp'] and r['dsig'] == -2], 'sharp sigma-dropping, LOCAL only', LOCAL)
# exact rules on all rows
print('\nexact rules (>= 300 rows, <= 0.5% exceptions), predicting good or not-good:')
lab = [r for r in R if r['label'] is not None]
for f in ['dsig', 'sig_drop', 'R_gt_half', 'ratio', 'Rf', 'reduced_crossings', 'twist_par', 'twist_ser', 'det_changed']:
    vals = sorted({r[f] for r in lab})
    if len(vals) > 60: vals = vals[::max(1, len(vals) // 60)]
    for v in vals:
        for op in ('<=', '>='):
            sel = [r for r in lab if (r[f] <= v if op == '<=' else r[f] >= v)]
            if len(sel) < 300: continue
            ng = sum(r['label'] == 'good' for r in sel)
            if ng >= 0.995 * len(sel): print(f'   {f} {op} {v:.4g}: n={len(sel)} good={ng}')
            if ng <= 0.005 * len(sel): print(f'   {f} {op} {v:.4g}: n={len(sel)} good={ng} (never good)')
