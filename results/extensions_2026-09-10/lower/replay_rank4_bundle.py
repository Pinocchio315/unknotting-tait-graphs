#!/usr/bin/env python3
"""Replay the three rank-four certificates and their five-change witnesses.

The replay regenerates the complete surgery-form lists, checks every basis
change, and evaluates the saved characteristic covectors. It never calls
m_Q or a shortest-vector search. The positive controls use the same direct
covector test and must remain inconclusive.
"""
import gzip
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import regina
from sympy import Matrix
import rank4_extension as r
import certify_rank4 as cert

HERE = Path(__file__).resolve().parent
NAMES = ('13n_4876', '13n_4973', '13n_5102')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    frozen, _ = r.audit.consolidate.read_greene_sweep()
    uppers = json.loads((HERE / 'rank4_upper_bound_witnesses.json').read_text())
    live = {v['name']: v for v in json.loads((HERE / 'rank4_live_retrieval.json').read_text())}
    summary = {'results': [], 'source_and_input_sha256': r.audit.hashes(), 'files_sha256': {}}
    for name in NAMES:
        if name == '13n_4876':
            source = HERE / (name + '_rank4.json')
            report = json.loads(source.read_text())
            proof_file = HERE / (name + '_rank4_covectors.json.gz')
            with gzip.open(proof_file, 'rt') as f:
                certificate = json.load(f)
            r.audit.require(certificate['input_report_sha256'] == digest(source), 'Changed report')
            indices, proof_pool = certificate['form_to_proof'], certificate['proofs']
        else:
            source = HERE / (name + '_rank4_direct.json.gz')
            with gzip.open(source, 'rt') as f:
                report = json.load(f)
            proof_file = source
            indices = [item['proof_index'] for item in report['comparisons']]
            proof_pool = report['proofs']
        summary['files_sha256'][source.name] = digest(source)
        summary['files_sha256'][proof_file.name] = digest(proof_file)
        r.audit.require(report['status'] == 'OBSTRUCTED', 'No completed obstruction')
        r.audit.require(report['input'] == frozen['targets'][name], 'Different frozen input')
        r.audit.validate_entry(report['input'])
        D, sigma = report['input']['determinant'], report['input']['signature']
        d = {int(k): F(v) for k, v in report['oriented_d'].items()}
        choices = {int(k): set(map(F, v)) for k, v in report['candidate_sets'].items()}
        _, _, survivors = r.audit.filter_vectors(choices, D, sigma,
            F(report['Jones_and_Casson_Walker']['required_sum_d']))
        r.audit.require(len(survivors) == 1 and r.oriented(survivors[0][1], sigma) == d,
                        'The correction terms do not follow from the candidate sets')

        forms = r.forms_for(D)
        matrices = [Q.tolist() for _, Q in forms]
        r.audit.require(matrices == [item['Q'] for item in report['comparisons']],
                        'The saved comparisons do not cover the regenerated form list')
        expected_hash = hashlib.sha256(json.dumps(matrices, separators=(',', ':')).encode()).hexdigest()
        r.audit.require(expected_hash == report['forms_sha256'], 'Changed surgery forms')
        r.audit.require(len(indices) == len(forms) == report['candidate_forms'], 'Missing forms')
        for proof in proof_pool:
            reduced = proof.get('Q', proof.get('reduced_Q'))
            cert.verify_certificate(np.array(reduced, dtype=np.int64), d, D, proof['proof'])
        noncyclic = 0
        for item, index in zip(report['comparisons'], indices):
            Q, P, R = (Matrix(item[k]) for k in ('Q', 'P', 'reduced_Q'))
            r.audit.require(abs(int(P.det())) == 1 and R == P.T * Q * P, 'Invalid congruence')
            proof = proof_pool[index]
            r.audit.require(R.tolist() == proof.get('Q', proof.get('reduced_Q')), 'Wrong proof matrix')
            noncyclic += 'noncyclic_minor_gcd' in proof['proof']

        witness = uppers[name]
        r.audit.require(witness['input_pd'] == report['input']['pd'], 'Wrong upper-bound diagram')
        r.audit.require(len(set(witness['crossing_indices_zero_based'])) == 5, 'Need five distinct crossings')
        success = False
        for _ in range(20):
            link = regina.Link.fromPD(witness['input_pd'])
            for i in witness['crossing_indices_zero_based']:
                link.change(link.crossing(i))
            link.simplify()
            if link.size() == 0 and link.countComponents() == 1:
                success = True
                break
        r.audit.require(success, 'Upper witness did not simplify to the unknot')
        r.audit.require(live[name]['unknotting_number'] == '[4,5]', 'Unexpected live comparison')
        summary['results'].append({
            'name': name, 'previous_range': [4, 5], 'unknotting_number': 5,
            'determinant': D, 'candidate_forms': len(forms), 'noncyclic_forms': noncyclic,
            'cyclic_form_isomorphisms_excluded': (len(forms) - noncyclic) * sum(
                math.gcd(a, D) == 1 for a in range(1, D)),
            'distinct_covector_proofs': len(proof_pool),
            'upper_crossings_zero_based': witness['crossing_indices_zero_based'],
        })
        print(name, 'replay complete:', len(forms), 'forms, u=5', flush=True)

    controls = json.loads((HERE / 'rank4_positive_controls.json').read_text())
    for row in controls:
        entry = frozen['controls'][row['name']]
        sets, _, _ = r.audit.pin(row['name'], entry['pd'], entry['determinant'], verbose=False)
        jones = r.audit.jones_data(r.audit.regina_jones(entry['pd']), entry['determinant'], entry['signature'])
        _, _, survivors = r.audit.filter_vectors(sets, entry['determinant'], entry['signature'],
                                                F(jones['required_sum_d']))
        vectors = dict(survivors)
        for match in row['compatible_forms']:
            vector = r.oriented(vectors[match['candidate_index']], entry['signature'])
            proof = cert.small_covector_certificate(np.array(match['R'], dtype=np.int64), vector,
                                                    entry['determinant'], radius=1)
            r.audit.require(proof is None, 'False direct obstruction on a known u=4 control')
    summary['known_u4_controls_pass'] = [row['name'] for row in controls]
    summary['additional_exact_values_this_run'] = 3
    summary['additional_exact_values_including_13n1619'] = 4
    summary['proposed_paper_exact_count_vs_2026_09_09_snapshot'] = 2529
    summary['proposed_open_knots_after_all_four_results'] = 2462
    summary['manuscript_and_frozen_results_modified'] = False
    for filename in ('rank4_extension.py', 'certify_rank4.py', 'direct_rank4_extension.py',
                     'replay_rank4_bundle.py', 'rank4_upper_bound_witnesses.json',
                     'rank4_positive_controls.json', 'rank4_live_retrieval.json'):
        summary['files_sha256'][filename] = digest(HERE / filename)
    (HERE / 'rank4_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('All three lower and upper bounds, form lists, and controls verified', flush=True)


if __name__ == '__main__':
    main()
