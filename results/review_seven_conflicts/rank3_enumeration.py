"""Independent ternary surgery-form enumeration for the seven-knot audit.

Uses exact scalar arithmetic, solves for f from the determinant equation,
and uses binary unimodular lifts rather than production {-1,0,1} lifts.
No correction terms are computed and no deposited results are modified.
"""
from __future__ import annotations
from itertools import product, permutations
from pathlib import Path
import json, sys, time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'lower_bounds/owens'))


def det3(A):
    return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
            - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
            + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))


def transpose(A):
    return list(map(list,zip(*A)))


def mul(A,B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def congruence(F,P):
    return mul(mul(transpose(P),F),P)


def inverse_unimodular(A):
    determinant = det3(A)
    assert abs(determinant)==1
    C=[]
    for i in range(3):
        row=[]
        for j in range(3):
            M=[[A[r][c] for c in range(3) if c!=j] for r in range(3) if r!=i]
            row.append((-1)**(i+j)*(M[0][0]*M[1][1]-M[0][1]*M[1][0]))
        C.append(row)
    return [[x//determinant for x in row] for row in transpose(C)]


def flat(A):
    return tuple(x for row in A for x in row)


def binary_lifts():
    out=[]
    for vals in product((0,1),repeat=9):
        A=[list(vals[3*i:3*i+3]) for i in range(3)]
        if det3(A)%2:
            assert abs(det3(A))==1
            out.append(A)
    assert len(out)==168
    return out


def reduced_forms(D):
    # D=(ad-b^2)f + 2bce - ae^2 - dc^2.
    out=[]
    a=1
    while a*a*a<=2*D:
        d=a
        while a*d*d<=2*D:
            for b,c in product(range(-(a//2),a//2+1),repeat=2):
                denominator=a*d-b*b
                assert denominator>0
                for e in range(-(d//2),d//2+1):
                    numerator=D-2*b*c*e+a*e*e+d*c*c
                    if numerator%denominator:
                        continue
                    f=numerator//denominator
                    if f<d or a*d*f>2*D:
                        continue
                    F=[[a,b,c],[b,d,e],[c,e,f]]
                    assert det3(F)==D
                    out.append(F)
            d+=1
        a+=1
    assert len(out)==len(set(map(flat,out)))
    return out


PERMS=list(permutations(range(3)))
SIGNS=list(product((-1,1),repeat=3))


def canonical_key(Q):
    # Encode the rank-three Q itself, not the six-dimensional trace.
    return min((tuple(Q[p[i]][p[i]] for i in range(3)),
                tuple(s[i]*s[j]*Q[p[i]][p[j]] for i in range(3) for j in range(i+1,3)))
               for p in PERMS for s in SIGNS)


def trace_matrix(Q):
    T=[[0]*6 for _ in range(6)]
    for i in range(3):
        T[2*i][2*i]=(Q[i][i]+1)//2
        T[2*i+1][2*i+1]=2
        T[2*i][2*i+1]=T[2*i+1][2*i]=1
        for j in range(3):
            if i!=j:
                T[2*i][2*j]=Q[i][j]//2
    return T


def admissible(Q):
    return (all(Q[i][i]%4==3 for i in range(3))
            and all(Q[i][j]%2==0 for i in range(3) for j in range(i+1,3)))


def audit(D):
    started=time.monotonic()
    import owens_u3 as production
    F_list=reduced_forms(D)
    prod_F={tuple(map(int,F.flat)) for F in production.reduced_ternary_forms(D)}
    assert set(map(flat,F_list))==prod_F
    lifts=binary_lifts()
    prod_lifts={tuple(int(x)%2 for x in P.flat):P.tolist() for P in production.gl3f2_lifts()}
    assert set(map(flat,lifts))==set(prod_lifts)
    transitions={}
    for P in lifts:
        R=mul(inverse_unimodular(P),prod_lifts[flat(P)])
        assert abs(det3(R))==1
        assert all(R[i][j]%2==int(i==j) for i in range(3) for j in range(3))
        transitions[flat(P)]=R
    candidates={}
    raw_admissible=0
    parity_equal=0
    candidate_frame_witnesses=[]
    for F in F_list:
        for P in lifts:
            Q=congruence(F,P)
            Qprod=congruence(F,prod_lifts[flat(P)])
            assert congruence(Q,transitions[flat(P)])==Qprod
            assert admissible(Q)==admissible(Qprod)
            parity_equal+=1
            if not admissible(Q):
                continue
            raw_admissible+=1
            key=canonical_key(Q)
            if key not in candidates:
                candidates[key]={'Q':Q,'Qt':trace_matrix(Q),
                                 'reduced_form':F,'binary_lift':P,
                                 'production_Q_in_same_Gamma2_class':Qprod,
                                 'Gamma2_transition':transitions[flat(P)]}
    prod_candidates=production.candidate_plumbings(D,n_even=3)
    output={
        'determinant':D,
        'arithmetic':'Exact Python integers; determinant solved for the last diagonal entry.',
        'reduced_form_count':len(F_list),
        'reduced_form_sets_equal_to_production':True,
        'binary_unimodular_lifts':len(lifts),
        'framewise_Gamma2_transitions_verified':parity_equal,
        'admissible_frames':raw_admissible,
        'admissibility_equal_to_production_for_every_frame':True,
        'alternative_candidate_count_up_to_signed_permutations':len(candidates),
        'production_candidate_count_up_to_signed_permutations':len(prod_candidates),
        'candidate_count_note':'Different integral lifts can produce different numbers of representatives after only signed-permutation deduplication. Framewise Gamma(2) transitions certify identical required class coverage.',
        'elapsed_seconds':round(time.monotonic()-started,3),
        'candidates':list(candidates.values()),
    }
    path=Path(__file__).with_name(f'rank3_enumeration_D{D}.json')
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({k:v for k,v in output.items() if k!='candidates'}),flush=True)
    print(path,flush=True)


if __name__=='__main__':
    for D in (239,155):
        audit(D)
