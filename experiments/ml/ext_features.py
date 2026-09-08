"""Extended feature table: Jones (from dataset.json), Khovanov (unreduced/reduced/odd integral vectors),
HOMFLY and Kauffman polynomials, torsion numbers of cyclic covers, nu, epsilon, tb, four-genus."""
import json, re, sys, ast, collections
import sympy as sp, numpy as np
import database_knotinfo as dk
S = sys.argv[1]
base = {d['name']: d for d in json.load(open(S + '/dataset.json'))}
rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
v, z, a = sp.symbols('v z a')
def poly2(s, x, y):
    P = sp.Poly(sp.expand(sp.sympify(s.replace('^', '**')) * x ** 40 * y ** 40), x, y)
    return {(e[0] - 40, e[1] - 40): int(c) for e, c in P.as_dict().items()}
def kh_feats(vec, prefix):
    V = ast.literal_eval(vec) if vec else []
    f = {}
    free = [(m, t, q) for tor, m, t, q in V if tor == 0]; tors = [(tor, m, t, q) for tor, m, t, q in V if tor != 0]
    f[prefix + 'rank'] = sum(m for m, t, q in free); f[prefix + 'ntors'] = sum(m for tor, m, t, q in tors)
    f[prefix + 'ntors_odd'] = sum(m for tor, m, t, q in tors if tor % 2 == 1)   # torsion of odd order
    f[prefix + 'maxtors'] = max([tor for tor, m, t, q in tors], default=0)
    ts = [t for m, t, q in free] or [0]; qs = [q for m, t, q in free] or [0]
    f[prefix + 'tmin'] = min(ts); f[prefix + 'tmax'] = max(ts); f[prefix + 'qmin'] = min(qs); f[prefix + 'qmax'] = max(qs)
    diags = collections.Counter(q - 2 * t for m, t, q in free for _ in range(m))
    f[prefix + 'width'] = (max(diags) - min(diags)) // 2 + 1 if diags else 0    # homological width (thin = 1 ... for unreduced thin = 2 diagonals)
    f[prefix + 'ndiag'] = len(diags); f[prefix + 'maxdiag'] = max(diags.values()) if diags else 0
    f[prefix + 'nhom'] = len({t for m, t, q in free})
    return f
out = {}
t0 = 0
for i, r in enumerate(rows):
    nm = r['name']
    if nm not in base: continue
    f = {}
    try:
        f.update(kh_feats(r['khovanov_unreduced_integral_vector'], 'khU_'))
        f.update(kh_feats(r['khovanov_reduced_integral_vector'], 'khR_'))
        f.update(kh_feats(r['khovanov_odd_integral_vector'], 'khO_'))
        H = poly2(r['homfly_polynomial'], v, z)
        f['hf_nterms'] = len(H); f['hf_l1'] = sum(abs(c) for c in H.values()); f['hf_vmin'] = min(e[0] for e in H); f['hf_vmax'] = max(e[0] for e in H)
        f['hf_zmax'] = max(e[1] for e in H); f['hf_maxc'] = max(abs(c) for c in H.values())
        for (vv, zz) in ((2, 1), (1, 2), (2, 2), (-1, 2), (3, 1)):
            f[f'hf_eval_{vv}_{zz}'] = int(sum(c * sp.Integer(vv) ** e[0] * sp.Integer(zz) ** e[1] for e, c in H.items()))
        f['hf_z0_sum'] = sum(c for e, c in H.items() if e[1] == 0); f['hf_top_z'] = sum(c for e, c in H.items() if e[1] == f['hf_zmax'])
        K = poly2(r['kauffman_polynomial'], a, z)
        f['kf_nterms'] = len(K); f['kf_l1'] = sum(abs(c) for c in K.values()); f['kf_amin'] = min(e[0] for e in K); f['kf_amax'] = max(e[0] for e in K)
        f['kf_zmax'] = max(e[1] for e in K); f['kf_maxc'] = max(abs(c) for c in K.values())
        for (aa, zz) in ((2, 1), (1, 2), (-1, 2), (2, 2)):
            f[f'kf_eval_{aa}_{zz}'] = int(sum(c * sp.Integer(aa) ** e[0] * sp.Integer(zz) ** e[1] for e, c in K.items()))
        T = json.loads(r['torsion_numbers']) if r['torsion_numbers'] else []
        for n, orders in T:
            f[f'tor{n}_ncyc'] = len(orders); f[f'tor{n}_order'] = int(np.prod(orders)) if all(o > 0 for o in orders) else 0; f[f'tor{n}_b1'] = sum(1 for o in orders if o == 0)
        nu = json.loads(r['nu']) if str(r['nu']).strip() not in ('', 'None') else [0, 0]
        if not isinstance(nu, list): nu = [nu, -nu]
        f['nu_plus'] = nu[0]; f['nu_minus'] = nu[1]
        f['epsilon'] = int(r['epsilon']) if str(r['epsilon']).strip() not in ('', 'None') else 0
        tb = re.findall(r'-?\d+', r['thurston_bennequin_number']); f['tb_max'] = max(int(x) for x in tb) if tb else 0; f['tb_min'] = min(int(x) for x in tb) if tb else 0
        f['g4'] = int(r['smooth_four_genus']) if str(r['smooth_four_genus']).strip().isdigit() else -1
        f['ribbon'] = int(r['ribbon_number']) if str(r['ribbon_number']).strip().lstrip('-').isdigit() else -1
    except Exception as e:
        f['error'] = repr(e)[:80]
    out[nm] = f
    if i % 2000 == 0: print(i, flush=True)
json.dump(out, open(S + '/ext_features.json', 'w'))
print('done', len(out), 'errors', sum(1 for f in out.values() if 'error' in f), flush=True)
