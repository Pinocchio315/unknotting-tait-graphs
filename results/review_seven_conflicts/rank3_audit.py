"""Independent arithmetic witnesses excluding every cyclic affine-free matching.

For each surgery form and group isomorphism, an explicit characteristic vector
has value (xi^T Q^-1 xi-r)/4 failing the required comparison with d_G. This value
is an upper bound for m_Q in its class. A negative difference excludes a matching;
a difference not in 2Z also excludes it because all characteristic vectors in one
class have congruent values modulo two. No production m_Q or matching code is used
for this exclusion. Exact Goeritz correction terms use full-box enumeration.
"""
from pathlib import Path
import ast,itertools,json,math,sys,time,hashlib
from fractions import Fraction
import gzip
from rank2_audit import inertia_signature
import regina
from spherogram import Link
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'lower_bounds/owens'))
import _compat,kinfo,crosscheck as cc

def matrix_data(Q):
    D=cc.det_frac(Q); n=len(Q)
    inv=cc.inv_frac(Q)
    adj=[[int(x*D) for x in row] for row in inv]
    assert D>0 and all(Fraction(adj[i][j],D)==inv[i][j] for i in range(n) for j in range(n))
    row=None
    for weights in itertools.product((0,1),repeat=n):
        v=[sum(weights[i]*adj[i][j] for i in range(n)) for j in range(n)]
        if math.gcd(D,*v)==1:
            row=v;break
    assert row is not None
    assert all(sum(row[i]*Q[i][j] for i in range(n))%D==0 for j in range(n))
    return D,adj,row

def val_label(xi,D,adj,row):
    label=sum(x*y for x,y in zip(row,xi))%D
    norm_num=sum(xi[i]*sum(adj[i][j]*xi[j] for j in range(len(xi))) for i in range(len(xi)))
    value=Fraction(norm_num-len(xi)*D,4*D)
    return label,value

def full_box(Q):
    D,adj,row=matrix_data(Q)
    best={}
    for xi in itertools.product(*(range(-Q[i][i],Q[i][i],2) for i in range(len(Q)))):
        label,value=val_label(xi,D,adj,row)
        if label not in best or value<best[label]:best[label]=value
    assert len(best)==D
    return best,row

def equivalent(A,B,D,negative=False):
    return [a for a in range(1,D) if math.gcd(a,D)==1 and all(A[j]==(-1 if negative else 1)*B[a*j%D] for j in range(D))]

def witnesses(Q,mg):
    D,adj,row=matrix_data(Q)
    assert D == len(mg)
    assert all(cc.det_frac([r[:i] for r in Q[:i]])>0 for i in range(1,len(Q)+1))
    units=[a for a in range(1,D) if math.gcd(a,D)==1]
    live=set(units); failures={}; tested=0
    for radius in (1,2,3):
        coords=[0]+[v for k in range(1,radius+1) for v in (2*k,-2*k)]
        for xi in itertools.product(coords,repeat=len(Q)):
            if radius>1 and all(abs(x)<2*radius for x in xi):continue
            tested+=1
            label,value=val_label(xi,D,adj,row)
            for a in tuple(live):
                delta=value-mg[a*label%D]
                if delta<0 or delta.denominator!=1 or delta.numerator%2:
                    failures[a]={'xi':list(xi),'class':label,'target_class':a*label%D,
                                 'covector_value':str(value),'d_target':str(mg[a*label%D]),
                                 'difference':str(delta)}
                    live.remove(a)
            if not live:
                return {'Qt':Q,'label_row':row,'isomorphisms':len(units),'vectors_examined':tested,
                        'max_coordinate_radius':2*radius,'failure_witnesses':failures,'verdict':'OBSTRUCTED'}
    return {'Qt':Q,'label_row':row,'isomorphisms':len(units),'remaining_units':sorted(live),'verdict':'INCOMPLETE'}

def main():
    report={'method':__doc__,'source_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'knots':{}}
    for name,D in [('12a_107',239),('13a_660',155)]:
        start=time.monotonic(); r=kinfo.row(name);pd=kinfo.pd_code(name)
        source_link=Link(pd); named_link=Link(name.replace('_',''))
        assert source_link.is_alternating()
        sig_source=regina.Link.fromPD(source_link.PD_code(min_strand_index=1)).sig(True,True)
        sig_named=regina.Link.fromPD(named_link.PD_code(min_strand_index=1)).sig(True,True)
        assert sig_source == sig_named
        V=ast.literal_eval(r['seifert_matrix'])
        S=[[V[i][j]+V[j][i] for j in range(len(V))] for i in range(len(V))]
        computed_signature=inertia_signature(S)
        assert computed_signature == int(r['signature'])
        assert abs(cc.det_frac(S)) == D
        Gs=cc.checkerboard_laplacians(pd)
        tables=[full_box(G) for G in Gs]
        assert all(cc.det_frac(G)==D for G in Gs)
        assert equivalent(tables[0][0],tables[1][0],D,negative=True)
        selected=[i for i,(mg,row) in enumerate(tables) if mg[0]==Fraction(-3,2)]
        assert len(selected)==1
        mg,labelrow=tables[selected[0]]
        source=Path(__file__).with_name(f'rank3_enumeration_D{D}.json')
        candidates=json.loads(source.read_text())['candidates']
        records=[]
        print(name,'independent Goeritz complete;',len(candidates),'candidate forms',flush=True)
        for i,entry in enumerate(candidates):
            Q=entry['Qt']
            assert all(Q[i][i]%2==0 for i in range(6))
            result=witnesses(Q,mg);records.append(result)
            assert result['verdict']=='OBSTRUCTED',result
            print(name,i+1,'/',len(candidates),'excluded;',result['vectors_examined'],'vectors',flush=True)
        report['knots'][name]={'determinant':D,'signature_tabulated':r['signature'],'signature_recomputed_from_Seifert':computed_signature,'pd':pd,'Seifert_matrix':V,'canonical_diagram_signature':sig_source,'matches_named_diagram_exactly':True,
            'Goeritz_forms':Gs,'Goeritz_cyclic_label_rows':[x[1] for x in tables],
            'Goeritz_d_vectors':[[str(x[0][j]) for j in range(D)] for x in tables],
            'selected_Goeritz_index':selected[0],'d_spin':str(mg[0]),
            'opposite_colourings_negate_d_up_to_isomorphism':True,
            'candidate_enumeration_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'candidate_count':len(records),'isomorphisms_excluded':sum(x['isomorphisms'] for x in records),
            'witnesses':records,'seconds':round(time.monotonic()-start,3)}
        out=ROOT/'results/review_seven_conflicts/rank3_independent.json'
        out.with_suffix('.json.gz').write_bytes(gzip.compress((json.dumps(report,sort_keys=True)+'\n').encode(),mtime=0))
        summary={**report,'knots':{k:{a:b for a,b in v.items() if a not in ('witnesses','Goeritz_d_vectors')} for k,v in report['knots'].items()}}
        summary['full_evidence_file']=out.with_suffix('.json.gz').name
        summary['full_evidence_sha256']=hashlib.sha256(out.with_suffix('.json.gz').read_bytes()).hexdigest()
        out.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print('All candidate forms excluded by explicit rational witnesses.',flush=True)

if __name__=='__main__':
    main()
