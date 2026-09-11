#!/usr/bin/env python3
"""Independent arithmetic and regression checks for compare.py.

This does not recompute the underlying Heegaard Floer invariants. It checks
the separation of the obstruction conditions, indexing conventions, exact
integer scaling, and agreement with the Fraction-based production routine.
The constructed non-geometric vectors are unit tests, not knot examples.
"""
from fractions import Fraction as F
import gzip
import json
import sys
import time
import compare as c


def main():
    start=time.monotonic()
    sys.path.insert(0,str(c.REPO/'lower_bounds/greene'))
    from half_integral import surgery_test,lens_half_integral_d
    counts=dict(lens_entries=0,lens_vectors=0,production_vectors=0)
    # Determinant-one Montesinos inputs have a trivial discriminant group.
    # There is still one identification; later conditions are vacuous.
    for v in (0,2,-2):
        spin=c.spin_profiles({0:F(v)},1)
        affine=c.affine_profiles({0:F(v)},1)
        assert all(spin['counts'].values())
        assert spin['counts']['PESMB']==affine['counts']['PESMB']
    for D in range(3,100,2):
        spin=(D+1)//2
        for i in range(D):
            assert c.lens(D,2,i)==lens_half_integral_d(D,i)
            counts['lens_entries']+=1
        d={(i-spin)%D:c.lens(D,2,i) for i in range(D)}
        half=(D-1)//2
        for i in range(half+1):
            independent=F(i*i,2*D)-(0 if i%2==half%2 else F(1,2))
            assert independent==d[i]
        result=c.spin_profiles(d,D)
        assert all(result['counts'].values()),D
        counts['lens_vectors']+=1

    # Distinguish the conditions using controlled normalized-gap sequences.
    # D=17 has symmetry pairs (0,1),(2,3),(4,5),(6,7).
    D=17; half=8
    sequences={
        'all_conditions':[0,0,0,0,1,1,1,1,2],
        'only_boundedness_fails':[0,0,2,2,2,2,2,2,2],
        'monotonicity_fails':[2,2,1,1,1,1,1,1,1],
        'symmetry_fails':[0,1,1,1,1,1,1,1,1],
    }
    synthetics={}
    for name,gap in sequences.items():
        d={}
        for i,v in enumerate(gap):
            d[i]=d[-i%D]=F(i*i,2*D)-(0 if i%2==half%2 else F(1,2))-2*v
        result=c.spin_profiles(d,D)
        synthetics[name]=result
        if name=='all_conditions':
            assert result['counts']['PESMB']
        elif name=='only_boundedness_fails':
            assert result['counts']['PESM'] and not result['counts']['PESMB']
        elif name=='monotonicity_fails':
            assert result['counts']['PES'] and not result['counts']['PESM']
        else:
            assert result['counts']['PEM'] and not result['counts']['PES']

    with gzip.open(c.REPO/'results/greene/sweep_2026-09-08.jsonl.gz','rt') as f:
        rows=[r for line in f if (r:=json.loads(line)).get('test')=='u1']
    rows.sort(key=lambda r:r['det'])
    samples=[(rows[i],0) for i in range(0,len(rows),max(1,len(rows)//18))]
    # Recheck all five vectors that pass symmetry and fail a later condition.
    extra={r['name']:r for r in rows}
    samples.extend((extra[name],i) for name,i in
                   [('12n_275',0),('12n_340',8),('13n_1465',2),
                    ('13n_1598',2),('13n_2822',0)])
    records=[]
    for row,i in samples:
        d=next(v for j,v in enumerate(c.vectors(row)) if i==j)
        affine=c.affine_profiles(d,row['det'])
        rational=surgery_test(d,row['det'])
        assert affine['counts']['PESMB']==len(rational['fits']),(row['name'],i)
        expected=sum(x['surviving']+x['invalid_successive_difference']
                     for x in rational['orientations'].values())
        assert affine['counts']['PES']==expected,(row['name'],i)
        records.append(dict(name=row['name'],det=row['det'],candidate_index=i,
                            symmetric_fits=expected,full_fits=len(rational['fits'])))
        counts['production_vectors']+=1
    result=dict(checks=counts,production_samples=records,
                synthetic_profiles={k:v['counts'] for k,v in synthetics.items()},
                seconds=time.monotonic()-start,
                compare_sha256=c.sha(c.HERE/'compare.py'),
                production_sha256=c.sha(c.REPO/'lower_bounds/greene/half_integral.py'))
    c.write(c.OUTPUT/'verification.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
