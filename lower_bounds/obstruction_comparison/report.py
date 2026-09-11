#!/usr/bin/env python3
"""Deduplicate the two comparison cohorts and export the per-knot comparison table."""
from collections import Counter
import csv
import json
import compare as c


def main():
    summaries={group:json.loads((c.OUTPUT/f'{group}_summary.json').read_text())
               for group in ('greene','montesinos')}
    combined={}
    skipped=[]
    for group,summary in summaries.items():
        assert summary['script_sha256']==c.sha(c.HERE/'compare.py')
        for source,digest in summary['input_sha256'].items():
            assert c.sha(c.ROOT/source)==digest,source
        for line in (c.OUTPUT/f'{group}_comparison.jsonl').read_text().splitlines():
            row=json.loads(line)
            if row.get('noncyclic'):
                skipped.append(row['name'])
                continue
            name=row['name']
            profile={p:bool(row['admitted_vectors'][p]) for p in c.PROFILES}
            if name in combined:
                assert combined[name]['passes']==profile,(name,'cohort disagreement')
                assert combined[name]['role']==row['role']
                combined[name]['cohorts'].append(group)
                combined[name]['candidate_vectors'][group]=row['candidate_vectors']
            else:
                combined[name]=dict(name=name,det=row['det'],role=row['role'],
                                    cohorts=[group],passes=profile,
                                    candidate_vectors={group:row['candidate_vectors']})
    values=list(combined.values())
    result=dict(unique_cyclic_knots=len(values),
                overlap=sum(len(r['cohorts'])>1 for r in values),
                excluded_noncyclic=sorted(skipped),
                roles=dict(Counter(r['role'] for r in values)),
                obstructed={p:sum(not r['passes'][p] for r in values) for p in c.PROFILES},
                strict_improvements=[r['name'] for r in values if r['passes']['PES'] and not r['passes']['PESMB']],
                symmetry_vs_monotonicity_differences=[r['name'] for r in values if r['passes']['PES']!=r['passes']['PEM']],
                compare_sha256=c.sha(c.HERE/'compare.py'),
                report_sha256=c.sha(c.HERE/'report.py'),
                source_summaries_sha256={f'{g}_summary.json':c.sha(c.OUTPUT/f'{g}_summary.json') for g in summaries})
    c.write(c.OUTPUT/'summary.json',result)
    with (c.OUTPUT/'knots.csv').open('w',newline='') as f:
        writer=csv.writer(f, lineterminator="\n")
        writer.writerow(['name','determinant','role','cohorts','Greene_candidate_vectors',
                         'Montesinos_candidate_vectors']+[p+'_passes' for p in c.PROFILES])
        for row in sorted(values,key=lambda r:r['name']):
            writer.writerow([row['name'],row['det'],row['role'],';'.join(row['cohorts']),
                             row['candidate_vectors'].get('greene',''),
                             row['candidate_vectors'].get('montesinos','')]+
                            [int(row['passes'][p]) for p in c.PROFILES])
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
