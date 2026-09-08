"""A: u=1 vs u>=2 with extended features (classifier, importances, residue mining).
   B: u=2 vs u=3 for alternating knots; ranking of the Owens-silent [2,3] alternating knots."""
import json, sys, re, glob, collections
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold, cross_val_predict
from sklearn.inspection import permutation_importance
S = sys.argv[1]; R = '/Users/pinocchio/Documents/05_unknotting_number/0001_writing_with_pavel/code'
base = {d['name']: d for d in json.load(open(S + '/dataset.json'))}
ext = json.load(open(S + '/ext_features.json'))
dkt = set(json.load(open(R + '/../server_code/data/dkt_104.json')))
for nm, d in base.items():
    d['f'].update({k: v for k, v in ext.get(nm, {}).items() if k != 'error'})
    for k in ('sigma', 'det', 'arf', 'tau', 's', 'g3', 'braid_index', 'bridge'):
        d['f'][k] = d[k] if d[k] is not None else 0
    if nm in dkt: d['exact'] = False        # DKT-derived KnotInfo values are not trusted
# our determinations
u3 = json.load(open(R + '/results/owens_u3_2026-08-24.json'))
for x in u3: base[x['name'] if isinstance(x, dict) else x].update({'exact': True, 'u_lo': 3, 'u_hi': 3})
u4 = json.load(open(R + '/lower_bounds/owens/results_u3_2026-09-07.json'))
for nm in u4['application']['exact_u4']: base[nm].update({'exact': True, 'u_lo': 4, 'u_hi': 4})
for nm in ['12a_107', '13a_660']: base[nm].update({'exact': True, 'u_lo': 4, 'u_hi': 4})
for nm in ['13a_15', '13a_55', '13a_422', '13a_568']: base[nm].update({'exact': True, 'u_lo': 3, 'u_hi': 3})
allf = sorted({k for d in base.values() for k in d['f']})
def mat(recs, cols): return np.array([[float(d['f'].get(c, 0)) for c in cols] for d in recs])
cv = StratifiedKFold(5, shuffle=True, random_state=0)
def gb(): return GradientBoostingClassifier(n_estimators=200, max_depth=3, learning_rate=0.05, random_state=0)
ex = [d for d in base.values() if d['exact']]
# ---------------- A ----------------
print('=== A: u = 1 vs u >= 2, extended features ===', flush=True)
y = np.array([1 if d['u_lo'] == 1 else 0 for d in ex])
J = [c for c in allf if not (c.startswith('kh') or c.startswith('hf_') or c.startswith('kf_') or c.startswith('tor') or c in ('nu_plus','nu_minus','epsilon','tb_max','tb_min','g4','ribbon','sigma','det','arf','tau','s','g3','braid_index','bridge') or c.startswith('alex'))]
KH = [c for c in allf if c.startswith('kh')]; HF = [c for c in allf if c.startswith('hf_')]; KF = [c for c in allf if c.startswith('kf_')]
TOR = [c for c in allf if c.startswith('tor')]; CL = ['sigma','det','arf','tau','s','g3','braid_index','bridge','nu_plus','nu_minus','epsilon','tb_max','tb_min','g4'] + [c for c in allf if c.startswith('alex')]
groups = {'Jones': J, 'Jones+classical': J + CL, 'Khovanov only': KH, 'HOMFLY only': HF, 'Kauffman only': KF, 'cyclic-cover torsion only': TOR,
          'classical only': CL, 'Jones+Khovanov+classical': J + KH + CL, 'everything': allf}
for label, cols in groups.items():
    X = mat(ex, cols)
    auc = cross_val_score(gb(), X, y, cv=cv, scoring='roc_auc').mean(); acc = cross_val_score(gb(), X, y, cv=cv, scoring='accuracy').mean()
    print(f'  {label:28s} AUC {auc:.3f} acc {acc:.3f}', flush=True)
X = mat(ex, allf); clf = gb().fit(X, y)
pi = permutation_importance(clf, X, y, n_repeats=5, random_state=0, scoring='roc_auc'); order = np.argsort(-pi.importances_mean)[:20]
print('  permutation importance top 20:', flush=True)
for i in order: print(f'     {allf[i]:16s} {pi.importances_mean[i]:.4f}', flush=True)
# residue mining on the new integer features
u1 = [d for d in ex if d['u_lo'] == 1]; u2 = [d for d in ex if d['u_lo'] >= 2]
open_rows = [d for d in base.values() if not d['exact'] and d['u_lo'] == 1]
ours815 = set(nm for nm, iv in json.load(open(R + '/results/lickorish_cw_obstructed_815_knotinfo2026.8.1.json')))
found = []
newint = KH + HF + KF + TOR + ['nu_plus','nu_minus','epsilon','tb_max','tb_min','g4']
for c in newint:
    v1 = [d['f'].get(c) for d in u1]; v2 = [d['f'].get(c) for d in u2]
    if any(x is None for x in v1 + v2): continue
    for m in range(2, 13):
        R1 = {int(x) % m for x in v1}
        if len(R1) == m: continue
        excl = sum(1 for x in v2 if int(x) % m not in R1)
        if excl == 0: continue
        eo = [d['name'] for d in open_rows if c in d['f'] and int(d['f'][c]) % m not in R1]
        found.append({'feature': c, 'mod': m, 'allowed_u1': sorted(R1), 'excl_known': excl, 'settles_open': len(eo), 'beyond_lickorish_cw': sum(1 for n in eo if n not in ours815)})
# plain inequalities: feature <= max over u=1 (bound-type conditions)
for c in newint + ['s','tau','g4','sigma']:
    v1 = [d['f'].get(c) for d in u1]; v2 = [d['f'].get(c) for d in u2]
    if any(x is None for x in v1 + v2): continue
    lo, hi = min(v1), max(v1)
    excl = sum(1 for x in v2 if x < lo or x > hi)
    if excl:
        eo = [d['name'] for d in open_rows if c in d['f'] and (d['f'][c] < lo or d['f'][c] > hi)]
        found.append({'feature': c, 'mod': f'range [{lo},{hi}]', 'allowed_u1': None, 'excl_known': excl, 'settles_open': len(eo), 'beyond_lickorish_cw': sum(1 for n in eo if n not in ours815)})
found.sort(key=lambda r: (-r['beyond_lickorish_cw'], -r['excl_known']))
json.dump(found, open(S + '/conditions_ext.json', 'w'), indent=1)
print(f'  residue/range conditions with 0 violations on u=1 (n={len(u1)}) and >0 exclusions: {len(found)}', flush=True)
for r in found[:20]:
    print(f"     {r['feature']:16s} {str(r['mod']):14s} allowed {str(r['allowed_u1'])[:36]:36s} excl {r['excl_known']:5d}/{len(u2)} | open {r['settles_open']:4d} | beyond L/CW {r['beyond_lickorish_cw']:4d}", flush=True)
# ---------------- B ----------------
print('=== B: u = 2 vs u = 3, alternating knots ===', flush=True)
alt = [d for d in ex if d['alternating'] and d['u_lo'] in (2, 3)]
yb = np.array([1 if d['u_lo'] == 3 else 0 for d in alt])
print('  training knots:', len(alt), '| u=3:', int(yb.sum()), '| by |sigma|:', dict(collections.Counter((abs(d['sigma']), d['u_lo']) for d in alt)), flush=True)
for label, cols in {'Jones+classical': J + CL, 'everything': allf, 'everything minus sigma/det': [c for c in allf if c not in ('sigma','det')]}.items():
    X = mat(alt, cols)
    auc = cross_val_score(gb(), X, yb, cv=cv, scoring='roc_auc').mean(); acc = cross_val_score(gb(), X, yb, cv=cv, scoring='accuracy').mean()
    print(f'  {label:28s} AUC {auc:.3f} acc {acc:.3f} (majority {max(yb.mean(),1-yb.mean()):.3f})', flush=True)
# within |sigma| = 4 only (the regime of the open knots)
s4 = [d for d in alt if abs(d['sigma']) == 4]; y4 = np.array([1 if d['u_lo'] == 3 else 0 for d in s4])
X4 = mat(s4, allf)
print(f'  |sigma|=4 subset: {len(s4)} knots, u=3: {int(y4.sum())}', flush=True)
pred = cross_val_predict(gb(), X4, y4, cv=cv, method='predict_proba')[:, 1]
from sklearn.metrics import roc_auc_score
print(f'  |sigma|=4 subset CV AUC (everything): {roc_auc_score(y4, pred):.3f}', flush=True)
clf4 = gb().fit(X4, y4)
pi4 = permutation_importance(clf4, X4, y4, n_repeats=5, random_state=0, scoring='roc_auc'); order = np.argsort(-pi4.importances_mean)[:12]
print('  |sigma|=4 importance top 12:', ', '.join(f'{allf[i]}({pi4.importances_mean[i]:.3f})' for i in order), flush=True)
# targets: Owens-PASS alternating [2,3] knots
passed = set()
for f in glob.glob(R + '/results/owens_sweep/*.jsonl'):
    for l in open(f):
        rr = json.loads(l)
        if rr.get('verdict') == 'PASS': passed.add(rr['name'])
targets = [d for d in base.values() if d['name'] in passed and not d['exact'] and d['u_lo'] == 2 and d['u_hi'] == 3 and d['alternating']]
Xt = mat(targets, allf); pt = clf4.predict_proba(Xt)[:, 1] if len(targets) else []
rank = sorted(zip([d['name'] for d in targets], pt), key=lambda x: x[1])
print(f'  Owens-PASS alternating [2,3] targets: {len(targets)}; P(u=3) distribution: min {min(pt):.2f} median {np.median(pt):.2f} max {max(pt):.2f}', flush=True)
print('  most likely u = 2 (witness-search priority):', [(n, round(p, 2)) for n, p in rank[:20]], flush=True)
print('  most likely u = 3:', [(n, round(p, 2)) for n, p in rank[-10:]], flush=True)
json.dump([{'name': n, 'P_u3': float(p)} for n, p in rank], open(S + '/u23_ranking.json', 'w'), indent=1)
