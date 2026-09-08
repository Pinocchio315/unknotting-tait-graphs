#!/usr/bin/env python
"""Obstruction to unknotting number one for Montesinos knots from the correction terms of the double branched
cover.

If u(K) = 1, the Montesinos trick gives Sigma_2(K) = S^3_{+-D/2}(C) for a knot C in S^3, D = det K.  By the
surgery formula of Ni and Wu, the correction terms of S^3_{D/2}(C) are
    d(S^3_{D/2}(C), i) = d(L(D,2), i) - 2 max(V_{floor(i/2)}, V_{floor((D+1-i)/2)}),   0 <= i < D,
where V_0 >= V_1 >= ... >= 0 are integers with V_j - V_{j+1} in {0, 1}; negative surgery is covered by the
mirror image (both orientations of Sigma_2 are tested).  For a Montesinos knot, Sigma_2(K) is a Seifert fibred
space, the boundary of a negative-definite star-shaped plumbing with at most one bad vertex, whose correction
terms are computed exactly from the plumbing form (Ozsvath-Szabo).  The test enumerates every affine
identification of Spin^c(Sigma_2) = Z/D with the labels of the surgery formula; if no identification and no
orientation admits a valid sequence V, then u(K) >= 2.  This is the obstruction in Section 3, "Montesinos knots", of
How to Compute the Unknotting Number: Theory and Computations.
The correction-term calculation uses ../owens. Passing a test only means that
these necessary conditions do not exclude unknotting number one.

    python montesinos_u1.py --validate --out /private/tmp/montesinos_validation.json  # known u = 1 (must pass) and u = 2 Montesinos knots, lens-space checks
    python montesinos_u1.py --apply --out /private/tmp/montesinos_apply.json  # every Montesinos knot with range [1, b]
    python montesinos_u1.py 11n_54 12n_309      # named knots
"""
import argparse, ast, itertools, json, math, os, re, sys, time
from fractions import Fraction
import numpy as np
import sympy
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, 'lower_bounds', 'owens'))
sys.path.insert(0, os.path.dirname(HERE))
from linking_pairing import cyclic_generator
from owens_obstruction import m_Q, ClassMap, det_int, definite_goeritz_candidates

def parse(notation):
    m = re.fullmatch(r'K\((.*)\)', notation.strip())
    if m is None or not m.group(1).strip():
        raise ValueError('expected Montesinos notation K(p/q; ...)')
    return [Fraction(t.strip()) for t in m.group(1).split(';')]

def hirzebruch_jung(alpha, beta):
    """alpha/beta = a_1 - 1/(a_2 - 1/(...)), all a_j >= 2, for 0 < beta < alpha."""
    if not isinstance(alpha, int) or not isinstance(beta, int) or not 0 < beta < alpha:
        raise ValueError('continued fractions require integers 0 < beta < alpha')
    out = []
    while beta:
        a = -(-alpha // beta); out.append(a); alpha, beta = beta, a * beta - alpha
    return out

def star_plumbing(fracs):
    """Negative-definite star plumbing of the Seifert space with unnormalised invariants (0; p_i/q_i).
    Returns (Q_neg, flipped): flipped = True if the orientation had to be reversed (e > 0)."""
    e = sum(fracs)
    if e == 0: raise ValueError('not a rational homology sphere')
    flipped = e > 0
    if flipped: fracs = [-f for f in fracs]
    e0 = sum(math.floor(f) for f in fracs)
    arms = []
    for f in fracs:
        b = f - math.floor(f)
        if b: arms.append(hirzebruch_jung(b.denominator, b.numerator))
    n = 1 + sum(len(a) for a in arms); Q = np.zeros((n, n), dtype=np.int64); Q[0, 0] = e0; k = 1
    for arm in arms:
        prev = 0
        for a in arm:
            Q[k, k] = -a; Q[k, prev] = Q[prev, k] = 1; prev = k; k += 1
    return Q, flipped

def lens_d(p, q, i):
    """d(L(p,q), i) with L(p,q) = p/q surgery on the unknot (Ni-Wu convention), recursive formula."""
    if p == 1: return Fraction(0)
    i %= p
    return Fraction(-1, 4) + Fraction((2 * i + 1 - p - q) ** 2, 4 * p * q) - lens_d(q, p % q, i % q)

def correction_terms(Q_neg):
    """d(Y, t) for Y = boundary of the negative-definite plumbing Q_neg, keyed by ClassMap coordinates."""
    validate_plumbing(Q_neg)
    Qpos = -np.array(Q_neg, dtype=object)
    mq, cm = m_Q(Qpos)
    return {g: -v for g, v in mq.items()}, cm

def validate_plumbing(Q_neg):
    """Check the graph hypotheses needed for the sharp plumbing formula.

    Ozsvath--Szabo, On the Floer homology of plumbed three-manifolds,
    Theorem 1.2, requires a negative-definite tree with at most one bad
    vertex. Negative definiteness is checked by m_Q (or exact LDL for HF).
    """
    matrix = sympy.Matrix(Q_neg)
    n = matrix.rows
    if not n or matrix.cols != n or matrix != matrix.T or any(x.q != 1 for x in matrix):
        raise ValueError('an integral symmetric plumbing matrix is required')
    edges = [(i, j) for i in range(n) for j in range(i)
             if matrix[i, j] == 1]
    if any(matrix[i, j] not in (0, 1) for i in range(n) for j in range(i)):
        raise ValueError('plumbing edges must have intersection +1')
    adjacent = [set() for _ in range(n)]
    for i, j in edges:
        adjacent[i].add(j); adjacent[j].add(i)
    reached = {0}; pending = [0]
    while pending:
        for j in adjacent[pending.pop()] - reached:
            reached.add(j); pending.append(j)
    if len(edges) != n - 1 or len(reached) != n:
        raise ValueError('the plumbing graph must be a tree')
    if any(matrix[i, i] >= 0 for i in range(n)):
        raise ValueError('the plumbing weights must be negative')
    if sum(int(matrix[i, i]) > -len(adjacent[i]) for i in range(n)) > 1:
        raise ValueError('the formula requires at most one bad vertex')


def cyclic_labelling(cm, D, Q_neg):
    """Enumerate a cyclic discriminant group and its +Q_neg**(-1) pairing.

    ClassMap stores primary factors, which may repeat the same prime. A
    unit in every factor need not generate their product: Z/3 + Z/3 is
    not cyclic. The exact inverse-column/CRT construction checks the full
    element order and returns None for every noncyclic group.
    """
    if cm.order != D:
        raise ValueError('the determinant does not match the discriminant group')
    matrix = sympy.Matrix(Q_neg)
    vector = cyclic_generator(matrix, D)
    if vector is None:
        return None, None
    g = cm.coords([int(t) for t in vector])
    labels = {j: tuple((j * int(x)) % int(d) for x, d in zip(g, cm.d))
              for j in range(D)}
    if len(set(labels.values())) != D:
        raise ArithmeticError('the proposed generator does not enumerate every class')
    pairing = Fraction((vector.T * matrix.inv() * vector)[0]) % 1
    return labels, pairing


def mapping_cone_ranks_ok(ranks_by_class, lab, D, a, b, V):
    """Second test (Lidman's argument): rank HF_red(Y, phi(i)) must equal sum_n rank A^red_{floor((i+nD)/2)} + T_i,
    where T_i is the torsion part determined by the sequence V (Gainullin, Cor. 14 and Prop. 15), and the A^red part has
    the form c_i + c_{i-1} with c_j = sum of rank A^red_s over s = j/2 (mod D), c_j >= 0.  True if consistent."""
    if D < 1 or D % 2 == 0:
        raise ValueError('the knot determinant must be positive and odd')
    if any(not isinstance(v, (int, np.integer)) or v < 0 for v in V) or any(V[j] - V[j + 1] not in (0, 1) for j in range(len(V) - 1)):
        raise ValueError('V must be a nonnegative integer sequence with steps zero or one')
    if not V or V[-1] != 0: return True                      # unknown tail of V: cannot refute
    def Vat(s):                                             # V_s for s >= 0 (0 beyond the observed range); V_{-s} = V_s + s
        if s < 0: return Vat(-s) - s
        return V[s] if s < len(V) else 0
    # T_i is the finite-dimensional tower-kernel contribution after the
    # infinite U-tower has been removed, not the tower itself. Every reduced
    # summand A_s contributes at residues 2s and 2s+1; hence r_i=c_i+c_(i-1).
    # For odd D this cyclic linear system has a unique rational solution.
    r = []
    for i in range(D):
        Ti = 0
        caseA = (i // 2) <= -((i - D) // 2)
        nV, nH = (1, 1) if caseA else (0, 2)
        n = nV
        while True:
            sv = (i + n * D) // 2
            if sv >= len(V): break
            Ti += Vat(sv); n += 1
        n = nH
        while True:
            sh = -((i - n * D) // 2)                        # H_{floor((i-nD)/2)} = V_{sh}
            if sh >= len(V): break
            Ti += Vat(sh); n += 1
        r.append(ranks_by_class[lab[(a * i + b) % D]] - Ti)
    if any(x < 0 for x in r): return False
    for i in range(D):                                      # c_i = (1/2) sum_j (-1)^j r_{i-j}  (D odd)
        c2 = sum((-1) ** j * r[(i - j) % D] for j in range(D))
        if c2 < 0 or c2 % 2: return False
    return True

def u1_admissible(d, cm, D, Q_neg, ranks_by_class=None):
    """True if some orientation and some affine identification of Z/D with Spin^c(Y) fits the surgery formula.
    Y is the boundary of the negative-definite plumbing.  Either Y = S^3_{D/2}(C) (eps = +1) or -Y = S^3_{D/2}(C)
    (eps = -1).  The linear part of the identification is the isomorphism H_1(S^3_{D/2}(C)) -> H_1(+-Y) induced by
    the surgery description, which preserves the linking pairing; the meridian class of S^3_{D/2}(C) has self-linking
    2/D, so the image a*g must have self-linking eps * 2/D on Y (the two signs are Lickorish's two allowed signs).
    The pairing of the sign with the orientation was calibrated on the lens spaces L(p,2) = S^3_{p/2}(unknot)."""
    lab, lam = cyclic_labelling(cm, D, Q_neg)
    if lab is None: return None            # H_1 not cyclic: u = 1 is excluded by Lickorish already
    dL = [lens_d(D, 2, i) for i in range(D)]
    idx = [min(i // 2, (D + 1 - i) // 2) for i in range(D)]
    for eps in (1, -1):
        units = [a for a in range(D) if math.gcd(a, D) == 1 and (a * a * lam) % 1 == Fraction(eps * 2, D) % 1]
        if not units: continue
        dY = [eps * d[lab[k]] for k in range(D)]
        for a in units:
            for b in range(D):
                V = {}; ok = True
                for i in range(D):
                    w = dL[i] - dY[(a * i + b) % D]
                    if w < 0 or w.denominator != 1 or w.numerator % 2: ok = False; break
                    j = idx[i]; v = w.numerator // 2
                    if V.setdefault(j, v) != v: ok = False; break
                if not ok: continue
                seq = [V[j] for j in sorted(V)]
                if all(0 <= seq[j] - seq[j + 1] <= 1 for j in range(len(seq) - 1)):
                    if ranks_by_class is None or mapping_cone_ranks_ok(ranks_by_class, lab, D, a, b, seq): return True
    return False

def plumbing_for(notation, det_k):
    """KnotInfo writes a two-bridge knot as K(beta/alpha) with det = alpha, whereas for three or more tangles
    K(p_1/q_1; ...) has det = |prod q_i * sum p_i/q_i|; the two conventions differ by inverting the fraction.
    Both are tried and the one whose plumbing has the right determinant is kept."""
    fr = parse(notation)
    for cand in (fr, [1 / f for f in fr] if len(fr) == 1 else None):
        if cand is None: continue
        try: Q, flipped = star_plumbing(cand)
        except ValueError: continue
        if abs(det_int(-Q)) == det_k: return Q, flipped
    return None, None

def test_knot(notation, det_k, with_hf=True):
    Q, flipped = plumbing_for(notation, det_k)
    if Q is None: return {'verdict': 'ERROR', 'note': f'no plumbing with determinant {det_k}'}
    d, cm = correction_terms(Q)
    adm_d = u1_admissible(d, cm, det_k, Q)
    if adm_d is None: return {'verdict': 'OBSTRUCTED', 'note': 'H_1 not cyclic', 'rank': len(Q)}
    rec = {'verdict_d': 'PASS' if adm_d else 'OBSTRUCTED', 'rank': len(Q), 'flipped': flipped, 'n_classes': cm.order}
    if with_hf and adm_d:
        try:
            from .hf_red import hf_red_all
        except ImportError:
            from hf_red import hf_red_all
        ranks, _ = hf_red_all(Q)
        if ranks is None:
            rec['verdict'] = rec['verdict_d']; rec['note'] = 'Floer ranks not computed (enumeration too large)'
        else:
            rec['hf_red_ranks'] = sorted(int(v) for v in ranks.values())
            rec['verdict'] = 'PASS' if u1_admissible(d, cm, det_k, Q, ranks_by_class=ranks) else 'OBSTRUCTED_HF'
    else:
        rec['verdict'] = rec['verdict_d']
    return rec

def goeritz_check(pd, Q_neg):
    """For alternating Montesinos knots: the correction terms from the star plumbing must equal those of the
    positive-definite Goeritz form (up to the overall orientation), as multisets."""
    d, _ = correction_terms(Q_neg); star = sorted(d.values())
    for G in definite_goeritz_candidates(pd):
        mq, _ = m_Q(np.array(G, dtype=np.int64)); goer = sorted(mq.values())
        if goer == star or sorted(-x for x in goer) == star: return True
    return False

def _job(job):
    nm, nt, det_k, lo = job[0], job[1], job[2], job[3]
    t0 = time.time(); res = test_knot(nt, det_k); res['seconds'] = round(time.time() - t0, 2)
    if len(job) == 5: res['range'] = [job[3], job[4]]
    else: res['u'] = lo
    return nm, res

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('names', nargs='*'); ap.add_argument('--validate', action='store_true'); ap.add_argument('--apply', action='store_true'); ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--out', help='JSON destination; required for --validate or --apply'); args = ap.parse_args()
    if (args.validate or args.apply) and not args.out:
        ap.error('--out is required for batch runs; deposited results are not output defaults')
    import database_knotinfo as dk
    rows = {r['name']: r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()}
    table = json.load(open(os.path.join(ROOT, 'results', 'reference_table_2026-09-07.json')))   # fixed reference intervals used in the paper
    mont = {nm: r['montesinos_notation'].strip() for nm, r in rows.items() if r['montesinos_notation'].strip() and not r['montesinos_notation'].strip().lower().startswith('not')}
    out = {'lens_check': None, 'validate': {}, 'apply': {}, 'named': {}}
    if args.validate:
        # 1. lens spaces: plumbing correction terms of L(p,q) = p/q surgery on the unknot vs the recursive formula
        bad = []
        for p in range(3, 40, 2):
            for q in range(1, p):
                if math.gcd(p, q) != 1: continue
                Q, fl = star_plumbing([Fraction(p, q)])
                if len(Q) > 8: continue                       # long chains are slow to enumerate and add nothing to the check
                d, cm = correction_terms(Q); a = sorted(d.values()); b = sorted(lens_d(p, q, i) for i in range(p))
                if a != b and sorted(-x for x in a) != b: bad.append((p, q))
        out['lens_check'] = {'pairs_tested': 'p<40 odd, plumbing rank <= 8', 'mismatches': bad}; print('lens-space check, mismatches:', bad, flush=True)
        twist = {}
        for p in range(3, 60, 2):
            Q, fl = star_plumbing([Fraction(p, 2)]); d, cm = correction_terms(Q); twist[p] = u1_admissible(d, cm, p, Q)
        out['twist_knot_check'] = {p: v for p, v in twist.items()}; print('L(p,2) = double covers of twist knots (u = 1), all admissible:', all(twist.values()), [p for p, v in twist.items() if not v], flush=True)
        # 2. alternating Montesinos knots: star plumbing vs Goeritz form
        alt_bad = []; n_alt = 0
        for nm, nt in mont.items():
            if not rows[nm]['alternating'].upper().startswith('Y') or len(parse(nt)) < 3: continue
            n_alt += 1
            Q, _ = plumbing_for(nt, int(rows[nm]['determinant']))
            if Q is None or not goeritz_check(ast.literal_eval(rows[nm]['pd_notation']), Q): alt_bad.append(nm)
        out['goeritz_check'] = {'alternating_montesinos_knots': n_alt, 'mismatches': alt_bad}; print(f'Goeritz check on {n_alt} alternating Montesinos knots, mismatches: {alt_bad}', flush=True)
        # 3. known u = 1 (must PASS) and known u = 2, in parallel
        from multiprocessing import Pool
        jobs = [(nm, nt, int(rows[nm]['determinant']), table[nm][0]) for nm, nt in sorted(mont.items()) if table[nm][0] == table[nm][1] and table[nm][0] in (1, 2)]
        with Pool(args.workers) as pool:
            for nm, res in pool.imap_unordered(_job, jobs, chunksize=4):
                out['validate'][nm] = res; print(f"  {nm} u={res['u']} {res['verdict']} {res['seconds']}s", flush=True)
        v = out['validate']
        print(f"known u=1: {sum(1 for r in v.values() if r['u']==1)} knots, PASS {sum(1 for r in v.values() if r['u']==1 and r['verdict']=='PASS')}, "
              f"OBSTRUCTED (false!) {[nm for nm, r in v.items() if r['u']==1 and r['verdict']!='PASS']}")
        print(f"known u=2: {sum(1 for r in v.values() if r['u']==2)} knots, OBSTRUCTED {sum(1 for r in v.values() if r['u']==2 and r['verdict'].startswith('OBSTRUCTED'))} (by HF_red ranks only: {[nm for nm, r in v.items() if r['u']==2 and r['verdict']=='OBSTRUCTED_HF']}), PASS {[nm for nm, r in v.items() if r['u']==2 and r['verdict']=='PASS']}", flush=True)
    if args.apply:
        from multiprocessing import Pool
        jobs = [(nm, nt, int(rows[nm]['determinant']), table[nm][0], table[nm][1]) for nm, nt in sorted(mont.items()) if table[nm][0] == 1 and table[nm][1] > 1]
        with Pool(args.workers) as pool:
            for nm, res in pool.imap_unordered(_job, jobs, chunksize=2):
                out['apply'][nm] = res; print(nm, res, flush=True)
        a = out['apply']; print(f"applied to {len(a)} Montesinos knots with range [1,b]: OBSTRUCTED {sum(1 for r in a.values() if r['verdict'].startswith('OBSTRUCTED'))} (by HF_red ranks only: {[nm for nm, r in a.items() if r['verdict']=='OBSTRUCTED_HF']}), PASS {[nm for nm, r in a.items() if r['verdict']=='PASS']}, ERROR {[nm for nm, r in a.items() if r['verdict']=='ERROR']}")
    for nm in args.names:
        res = test_knot(mont[nm], int(rows[nm]['determinant'])); out['named'][nm] = res; print(nm, mont[nm], res)
    if args.out and (args.validate or args.apply or args.names):
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        json.dump(out, open(args.out, 'w'), indent=1); print('wrote', args.out)

if __name__ == '__main__':
    main()
