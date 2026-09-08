"""Owens verdict + closed-form number-theoretic features for every alternating knot with |sigma| = 4."""
import sys, json, re, ast, time
import database_knotinfo as dk
from fractions import Fraction
from owens_obstruction import obstruct_u2, candidates, qtilde, definite_goeritz_candidates, det_int, m_Q, matching_exists
from linkform import linking_invariants, compatible
out, shard = sys.argv[1], sys.argv[2]; si, sn = (int(x) for x in shard.split('/'))
rows = [r for r in dk.link_list() if str(r.get('crossing_number', '')).strip().isdigit()]
alt = lambda r: str(r.get('alternating', '')).upper().startswith('Y')
sel = [r for r in rows if alt(r) and abs(int(r['signature'])) == 4][si::sn]
with open(out, 'a') as h:
    for r in sel:
        nm = r['name']; pd = ast.literal_eval(r['pd_notation']); sig = int(r['signature']); D = int(r['determinant']); t0 = time.time()
        rec = {'name': nm, 'u': str(r['unknotting_number']).strip(), 'sigma': sig, 'det': D, 'crossings': int(r['crossing_number'])}
        try:
            res = obstruct_u2(pd, sig, D); rec['verdict'] = res['verdict']; rec['witness_form'] = res.get('witness_form')
            cands = candidates(D, 2); rec['n_candidates'] = len(cands)
            # sharp Goeritz side and its linking-form invariants
            G = None
            for g in definite_goeritz_candidates(pd):
                if det_int(g) != D: continue
                mG, cmG = m_Q(g)
                if mG.get(cmG.coords([0] * len(g))) == Fraction(-4, 4): G = g; mGd = mG; cmGd = cmG; break
            invG = linking_invariants(G); rec['linkform_G'] = {str(k): v for k, v in invG.items()}
            compat = [c for c in cands if compatible(linking_invariants(qtilde(*c)), invG)]
            rec['n_linkform_compatible'] = len(compat)
            # d-invariant multiset of the sharp side (for later analysis)
            rec['dvals'] = sorted([str(v) for v in mGd.values()])
        except Exception as e:
            rec['error'] = repr(e)[:100]
        rec['seconds'] = round(time.time() - t0, 1)
        h.write(json.dumps(rec) + '\n'); h.flush()
