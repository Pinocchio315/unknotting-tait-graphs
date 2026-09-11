#!/usr/bin/env python3
"""Exclude surgery candidates directly with short characteristic covectors.

This can establish an obstruction without computing any m-function. A failure
to find short witnesses is explicitly inconclusive; it is not an obstruction
and does not imply the existence of a compatible surgery form.
"""
import argparse
import gzip
import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import rank4_extension as r
import certify_rank4 as cert

HERE = Path(__file__).resolve().parent


def run(name, frozen, metadata):
    began = time.monotonic()
    entry = frozen['targets'][name]
    evidence = r.audit.validate_entry(entry)
    D, sigma = entry['determinant'], entry['signature']
    r.audit.require(entry['test'] == 'rank4' and abs(sigma) == 8, 'Wrong rank-four hypotheses')
    diagnostics = {}
    sets, pairing, reference = r.audit.pin(name, entry['pd'], D, verbose=False,
                                          diagnostics=diagnostics)
    poly = r.audit.regina_jones(entry['pd'])
    r.audit.require(poly == r.audit.parse_jones(metadata[name]['jones_polynomial']), 'Jones mismatch')
    jones = r.audit.jones_data(poly, D, sigma)
    total, spin, surviving = r.audit.filter_vectors(sets, D, sigma, F(jones['required_sum_d']))
    r.audit.require(len(surviving) == 1, 'An unresolved vector requires a larger comparison')
    d = r.oriented(surviving[0][1], sigma)
    r.audit.require(d[0] == -2, 'Wrong oriented spin term')
    report = {
        'name': name, 'input': entry, 'lspace_evidence': evidence,
        'Greene_diagnostics': diagnostics, 'pairing': str(pairing), 'reference': reference,
        'candidate_sets': {str(k): list(map(str, sorted(v))) for k, v in sets.items()},
        'candidate_vectors': total, 'after_spin': spin, 'after_sum': len(surviving),
        'Jones_and_Casson_Walker': jones,
        'oriented_d': {str(k): str(v) for k, v in d.items()},
        'd': {str(k): str(v) for k, v in surviving[0][1].items()},
        'method': 'explicit characteristic covectors, without computing m-functions',
        'status': 'RUNNING', 'comparisons': [], 'proofs': [],
    }
    forms = r.forms_for(D)
    report['candidate_forms'] = len(forms)
    report['forms_sha256'] = hashlib.sha256(json.dumps(
        [Q.tolist() for _, Q in forms], separators=(',', ':')).encode()).hexdigest()
    cache = {}
    missing = []
    def save():
        report['seconds'] = round(time.monotonic() - began, 3)
        # Compact progress file; the full proof is written once at completion.
        progress = {k: report[k] for k in ('name', 'status', 'candidate_forms', 'seconds')}
        progress['compared'] = len(report['comparisons'])
        progress['without_short_witnesses'] = len(missing)
        (HERE / (name + '_direct_progress.json')).write_text(json.dumps(progress, indent=2) + '\n')
    for number, (key, Q) in enumerate(forms, 1):
        reduced, P = r.reduce_gram(Q)
        key_reduced = tuple(map(int, reduced.flatten()))
        if key_reduced not in cache:
            proof = cert.small_covector_certificate(np.array(reduced, dtype=np.int64), d, D, 1)
            if proof is not None:
                # Separate exact arithmetic replay, including every unit.
                cert.verify_certificate(np.array(reduced, dtype=np.int64), d, D, proof)
                cache[key_reduced] = len(report['proofs'])
                report['proofs'].append({'Q': reduced.tolist(), 'proof': proof})
            else:
                cache[key_reduced] = None
        index = cache[key_reduced]
        if index is None:
            missing.append(number)
        report['comparisons'].append({'number': number, 'key': str(key), 'Q': Q.tolist(),
                                      'reduced_Q': reduced.tolist(), 'P': P.tolist(),
                                      'proof_index': index})
        if number % 200 == 0:
            save()
            print(name, 'direct witnesses', number, '/', len(forms), 'missing', len(missing),
                  'seconds', report['seconds'], flush=True)
    report['status'] = 'OBSTRUCTED' if not missing else 'INCONCLUSIVE'
    report['without_short_witnesses'] = missing
    if not missing:
        report['lower_bound'] = 5
    report['distinct_reduced_forms'] = len(cache)
    report['sha256'] = r.audit.hashes()
    for filename in ('rank4_extension.py', 'certify_rank4.py', 'direct_rank4_extension.py'):
        report['sha256'][str((HERE / filename).relative_to(r.audit.ROOT))] = hashlib.sha256(
            (HERE / filename).read_bytes()).hexdigest()
    save()
    with gzip.open(HERE / (name + '_rank4_direct.json.gz'), 'wt') as handle:
        json.dump(report, handle, separators=(',', ':'))
    print(name, report['status'], 'all', len(forms), 'forms, seconds', report['seconds'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('names', nargs='+')
    args = parser.parse_args()
    frozen, _ = r.audit.consolidate.read_greene_sweep()
    metadata = r.audit.metadata()
    for name in args.names:
        run(name, frozen, metadata)
