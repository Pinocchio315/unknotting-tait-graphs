#!/usr/bin/env python3
"""Compare the half-integral correction-term obstructions used in Section 3.9.

P/E/S/M/B are positivity, evenness, symmetry, monotonicity, and boundedness
in Owens--Strle, arXiv:1311.6702, Theorem 3. P/E/S is the comparison with
the older symmetric matching obstruction; the full test adds M and B.
This is not a claim that OS's paper contains no stronger result: see its
Theorem 8.4. A comparison producing no improvement over P/E/S also produces
no improvement over a strengthening of P/E/S on the same inputs.

Two independently indexed routines are used. The first follows Theorem 3
with spin at zero. The second uses the lens-space recursion and every affine
map in the Ni--Wu surgery coordinates, as in the deposited implementation.
All arithmetic is integral after multiplying correction terms by 4D.
No new spanning-tree enumeration or HF_red mapping-cone test is performed.
An ambiguous knot passes a profile if ANY whole candidate vector passes;
discarding an individual vector is not an obstruction for that knot.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import cache
import gzip
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
REPO = ROOT
OUTPUT = ROOT / 'generated/obstruction_comparison'
OUTPUT.mkdir(parents=True, exist_ok=True)
PROFILES = ('PE', 'PES', 'PEM', 'PESM', 'PESMB')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')


@cache
def units(D):
    # The trivial group (D=1) has one automorphism. It occurs among the
    # Montesinos inputs; its symmetry and successive-step conditions are vacuous.
    return (0,) if D == 1 else tuple(a for a in range(1, D) if gcd(a, D) == 1)


@cache
def lens(p, q, i):
    """Positive surgery convention, computed by the lens-space recursion."""
    if p == 1:
        return F(0)
    return -F(1, 4) + F((2*i+1-p-q)**2, 4*p*q) - lens(q, p % q, i % q)


def scaled(d, D):
    assert D >= 1 and D % 2 == 1 and set(d) == set(range(D))
    dd = [4*D*F(d[i]) for i in range(D)]
    assert all(x.denominator == 1 for x in dd)
    assert all(dd[i] == dd[-i % D] for i in range(D)), 'Spin must be at zero'
    return list(map(int, dd))


def spin_profiles(d, D):
    """Literal Theorem 3, all units and both orientations, spin at zero.

    The returned vectors are normalized gaps divided by two. Thus symmetry
    is equality, monotonicity is a nonnegative successive difference, and
    boundedness says the difference is at most one. Both orientations and
    no linking-form filter make a negative result conservative.
    """
    dd = scaled(d, D)
    half = (D-1)//2
    normal = [2*i*i - (0 if i % 2 == half % 2 else 2*D)
              for i in range(half+1)]
    counts = dict.fromkeys(PROFILES, 0)
    first = {}
    extra = []
    sympairs = ([(2*j, 2*j+1) for j in range((D-5)//4+1)] if D % 4 == 1
                else [(2*j-1, 2*j) for j in range(1, (D-3)//4+1)])
    for sign in (1, -1):
        for a in units(D):
            gap = []
            for i, r in enumerate(normal):
                v = r-sign*dd[a*i % D]
                if v < 0 or v % (8*D):
                    break
                gap.append(v//(8*D))
            else:
                s = all(gap[i] == gap[j] for i, j in sympairs)
                steps = [v-u for u, v in zip(gap, gap[1:])]
                m = all(t >= 0 for t in steps)
                b = all(t <= 1 for t in steps)
                passed = dict(PE=True, PES=s, PEM=m, PESM=s and m,
                              PESMB=s and m and b)
                fit = dict(orientation=sign, unit=a, normalized_gaps_over_two=gap)
                for profile, ok in passed.items():
                    if ok:
                        counts[profile] += 1
                        first.setdefault(profile, fit)
                if s and not (m and b):
                    extra.append(dict(**fit, monotone=m, bounded=b,
                                      failed_steps=[(i,t) for i,t in enumerate(steps)
                                                    if t not in (0,1)]))
    return dict(counts=counts, first_fits=first, symmetric_but_not_full=extra)


@cache
def lens_data(D):
    scaled_lens = [4*D*lens(D, 2, i) for i in range(D)]
    assert all(v.denominator == 1 for v in scaled_lens)
    idx = [min(i//2, (D+1-i)//2) for i in range(D)]
    return tuple(map(int, scaled_lens)), tuple(idx)


def affine_profiles(d, D):
    """Independent replay in surgery coordinates; includes every translation.

    P/E/S here means the same-V-index equality, in addition to nonnegative
    even gaps. M/B then constrains the V sequence. No spin restriction is
    imposed; spin_profiles is independently compared to this enlarged test.
    """
    dd = scaled(d, D)
    ll, idx = lens_data(D)
    counts = dict(PE=0, PES=0, PESMB=0)
    first = {}
    for sign in (1, -1):
        for a in units(D):
            for b in range(D):
                gap = []
                for i, r in enumerate(ll):
                    v = r-sign*dd[(a*i+b) % D]
                    if v < 0 or v % (8*D):
                        break
                    gap.append(v//(8*D))
                else:
                    counts['PE'] += 1
                    first.setdefault('PE', dict(orientation=sign, unit=a, translation=b))
                    vv = {}
                    for j, v in zip(idx, gap):
                        if j in vv and vv[j] != v:
                            break
                        vv[j] = v
                    else:
                        seq = [vv[j] for j in range(max(idx)+1)]
                        fit = dict(orientation=sign, unit=a, translation=b, V=seq)
                        counts['PES'] += 1
                        first.setdefault('PES', fit)
                        if all(u-v in (0,1) for u,v in zip(seq,seq[1:])):
                            counts['PESMB'] += 1
                            first.setdefault('PESMB', fit)
    return dict(counts=counts, first_fits=first)


def vectors(row):
    if row.get('d'):
        yield {int(k): F(v) for k,v in row['d'].items()}
        return
    D = row['det']
    fixed = {int(k): F(v) for k,v in row['pinned'].items()}
    amb = {int(k): sorted(map(F, v)) for k,v in row['ambiguous'].items()}
    assert all(amb[k] == amb[-k % D] for k in amb)
    free = sorted(k for k in amb if k <= -k % D)
    for choices in product(*(amb[k] for k in free)):
        d = dict(fixed)
        for k,v in zip(free,choices):
            d[k] = d[-k % D] = v
        yield d


def compare_record(row, replay_affine=True):
    start = time.monotonic()
    result = {k:row[k] for k in ('name','det','role','verdict','sigma','range') if k in row}
    admitted = dict.fromkeys(PROFILES, 0)
    affine_admitted = dict.fromkeys(('PE','PES','PESMB'),0)
    extras = []
    details = []
    for i,d in enumerate(vectors(row)):
        spin = spin_profiles(d, row['det'])
        aa = affine_profiles(d,row['det']) if replay_affine else None
        for p,n in spin['counts'].items():
            admitted[p] += bool(n)
        if aa:
            for p,n in aa['counts'].items():
                affine_admitted[p] += bool(n)
                assert not spin['counts'][p] or n, (row['name'],i,p,'spin fits absent from affine')
            # Independent profile outcome agrees with the saved production
            # test on every historical vector where its verdict was stored.
            old = row.get('candidate_verdicts', row.get('candidate_tests', []))
            if i < len(old):
                saved = old[i].get('admits')
                if saved is None and 'surgery_test' in old[i]:
                    saved = bool(old[i]['surgery_test']['fits'])
                if saved is not None:
                    assert bool(aa['counts']['PESMB']) == bool(saved), (row['name'], i, 'production disagreement')
        if spin['symmetric_but_not_full']:
            extras.append(dict(candidate_index=i, d={str(k):str(v) for k,v in d.items()},
                               fits=spin['symmetric_but_not_full']))
        details.append(dict(candidate_index=i, spin_counts=spin['counts'],
                            affine_counts=aa['counts'] if aa else None))
    result.update(candidate_vectors=len(details), admitted_vectors=admitted,
                  affine_admitted_vectors=affine_admitted if replay_affine else None,
                  strict_improvement=bool(admitted['PES'] and not admitted['PESMB']),
                  candidate_level_extra_rejections=extras, candidates=details,
                  seconds=round(time.monotonic()-start,4))
    if row.get('role') == 'control' and row.get('known_unknotting_number') == 1:
        assert admitted['PESMB'], (row['name'],'known u=1 control rejected')
    return result


def summarize(rows):
    cyclic = [r for r in rows if 'admitted_vectors' in r]
    return dict(records=len(rows), tested_cyclic=len(cyclic),
                skipped_noncyclic=[r['name'] for r in rows if r.get('noncyclic')],
                by_role=dict(Counter(r['role'] for r in cyclic)),
                candidates=sum(r['candidate_vectors'] for r in cyclic),
                obstructed={p:sum(not r['admitted_vectors'][p] for r in cyclic) for p in PROFILES},
                strict_improvements=[r['name'] for r in cyclic if r['strict_improvement']],
                symmetry_vs_monotonicity_differences=[r['name'] for r in cyclic
                  if bool(r['admitted_vectors']['PES']) != bool(r['admitted_vectors']['PEM'])],
                candidate_extra_rejection_knots=[r['name'] for r in cyclic if r['candidate_level_extra_rejections']],
                candidate_extra_rejections=sum(len(r['candidate_level_extra_rejections']) for r in cyclic),
                spin_affine_outcome_differences=[(r['name'],p) for r in cyclic for p in ('PE','PES','PESMB')
                  if r['affine_admitted_vectors'] is not None
                  and bool(r['admitted_vectors'][p]) != bool(r['affine_admitted_vectors'][p])])


def greene():
    source = REPO/'results/greene/sweep_2026-09-08.jsonl.gz'
    with gzip.open(source,'rt') as f:
        rows = {r['name']:r for line in f if (r:=json.loads(line)).get('test') == 'u1'}
    hashes = {str(source.relative_to(ROOT)):sha(source)}
    # Replace the historical 6561-vector 13n1619 input by the proved exact
    # vector added to v1.8. It is not counted as a second knot.
    update = REPO/'results/extensions_2026-09-10/lower/13n_1619_certificate.json'
    cert = json.loads(update.read_text())
    old = rows['13n_1619']
    rows['13n_1619'] = {k:v for k,v in old.items() if k not in
                        ('d','pinned','ambiguous','candidate_verdicts','candidate_vectors')}
    assert not cert['surgery_test']['fits']
    rows['13n_1619'].update(d=cert['d'],verdict='OBSTRUCTED',candidate_vectors=1,
                           candidate_verdicts=[dict(admits=False)])
    hashes[str(update.relative_to(ROOT))] = sha(update)
    for name in ('12n_491','13n_3370'):
        src = REPO/f'results/bernhard_jablan/greene_d_{name}_2026-09-08.json'
        record = json.loads(src.read_text())
        assert name not in rows
        record['role'] = 'target'
        rows[name] = record
        hashes[str(src.relative_to(ROOT))] = sha(src)
    out = []
    start = time.monotonic()
    with (OUTPUT/'greene_comparison.jsonl').open('w') as f:
        for i,row in enumerate(rows.values(),1):
            result = compare_record(row)
            f.write(json.dumps(result,sort_keys=True)+'\n'); f.flush()
            out.append(result)
            if i % 200 == 0:
                print(f'Greene {i}/{len(rows)}; {time.monotonic()-start:.1f}s',flush=True)
    result = summarize(out)
    result.update(input_sha256=hashes,script_sha256=sha(Path(__file__)),seconds=time.monotonic()-start)
    write(OUTPUT/'greene_summary.json',result)
    print(json.dumps(result,indent=2),flush=True)


def montesinos():
    # Heavy plumbing routines run sequentially; BLAS thread counts should
    # also be set to one by the caller. Cache every completed exact vector.
    sys.path[:0] = [str(REPO), str(REPO/'lower_bounds/montesinos')]
    from knotinfo_inputs import load_metadata
    import montesinos_u1 as mont
    source = REPO/'results/montesinos/montesinos_u1_2026-09-07.json'
    old = json.loads(source.read_text())
    archive = REPO/'results/knotinfo_calculation_inputs_2026-09-09.json.gz'
    meta = load_metadata()
    rows=[]
    for group,role in (('apply','target'),('validate','control')):
        for name,record in old[group].items():
            rows.append(dict(name=name,role=role,**record))
    cachepath=OUTPUT/'montesinos_correction_terms.jsonl'
    saved={r['name']:r for line in cachepath.read_text().splitlines()
           if (r:=json.loads(line))} if cachepath.exists() else {}
    output=[]
    start=time.monotonic()
    with cachepath.open('a') as cf, (OUTPUT/'montesinos_comparison.jsonl').open('w') as f:
        for i,row in enumerate(rows,1):
            name=row['name']; md=meta[name]; D=int(F(str(md['determinant'])))
            row.update(det=D,sigma=int(F(str(md['signature']))))
            if name not in saved:
                began=time.monotonic()
                q,flipped=mont.plumbing_for(md['montesinos_notation'],D)
                assert q is not None, name
                # Check cyclicity before potentially expensive lattice enumeration.
                from owens_obstruction import ClassMap
                cm=ClassMap(-q)
                labels,pairing=mont.cyclic_labelling(cm,D,q)
                rec=dict(name=name,det=D,notation=md['montesinos_notation'],
                         plumbing=q.tolist(),flipped=flipped,noncyclic=labels is None)
                if labels is not None:
                    dd,cm=mont.correction_terms(q)
                    labels,pairing=mont.cyclic_labelling(cm,D,q)
                    rec.update(d={str(j):str(dd[labels[j]]) for j in range(D)},pairing=str(pairing))
                rec['seconds']=time.monotonic()-began
                saved[name]=rec
                cf.write(json.dumps(rec,sort_keys=True)+'\n');cf.flush()
            rec=saved[name]
            assert rec['det']==D and rec['notation']==md['montesinos_notation']
            if rec['noncyclic']:
                result=dict(name=name,det=D,role=row['role'],noncyclic=True)
            else:
                row['d']=rec['d']
                row['known_unknotting_number']=row.get('u')
                result=compare_record(row)
                result['known_unknotting_number']=row.get('u')
                # The archived test also imposes a linking filter. Keep any
                # outcome discrepancy explicit instead of masking it.
                result['archived_verdict_d']=row.get('verdict_d')
                result['archived_disagreement']=bool(result['affine_admitted_vectors']['PESMB']) != (row.get('verdict_d')=='PASS')
            output.append(result)
            f.write(json.dumps(result,sort_keys=True)+'\n');f.flush()
            if i % 50 == 0:
                print(f'Montesinos {i}/{len(rows)}; {time.monotonic()-start:.1f}s',flush=True)
    result=summarize(output)
    result.update(input_sha256={str(p.relative_to(ROOT)):sha(p) for p in (source,archive)},
                  script_sha256=sha(Path(__file__)),seconds=time.monotonic()-start,
                  archived_disagreements=[r['name'] for r in output if r.get('archived_disagreement')])
    write(OUTPUT/'montesinos_summary.json',result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cohort',choices=('greene','montesinos'))
    args=ap.parse_args()
    globals()[args.cohort]()
