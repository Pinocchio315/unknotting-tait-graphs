#!/usr/bin/env python3
"""Complete the three previously unfinished rank-four comparisons.

Before enumerating characteristic covectors, change the Gram matrix by an
exact unimodular congruence. This changes only the lattice basis, not the
set of surgery candidates or their correction-term bounds. The stored
transformation allows the change of basis to be checked independently.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from fractions import Fraction as F
from pathlib import Path

import numpy as np
import verify_casson_walker as audit
from greene_ranks import cyclic_index, oriented
from owens_obstruction import det_int, m_Q, positive_definite
from owens_u4 import candidate_plumbings

HERE = Path(__file__).resolve().parent


def gram_schmidt(gram):
    """Exact Gram--Schmidt coefficients, from the Gram matrix alone."""
    n = len(gram)
    mu = [[F(0) for _ in range(n)] for _ in range(n)]
    norms = []
    for i in range(n):
        for j in range(i):
            mu[i][j] = (F(gram[i, j]) - sum(
                mu[i][l] * mu[j][l] * norms[l] for l in range(j))) / norms[j]
        norm = F(gram[i, i]) - sum(mu[i][j] ** 2 * norms[j] for j in range(i))
        audit.require(norm > 0, 'Nonpositive Gram--Schmidt norm')
        norms.append(norm)
    return mu, norms


def reduce_gram(matrix):
    """Exact LLL basis changes; return R, P with R = P^T Q P and det(P)=+-1.

    No floating-point decision is used. More importantly, the downstream
    obstruction only needs the verified integral congruence, not the claim
    that the output satisfies any particular lattice reduction condition.
    """
    original = np.array(matrix, dtype=object)
    gram = original.copy()
    n = len(gram)
    transform = np.eye(n, dtype=object)
    k, steps = 1, 0
    while k < n:
        steps += 1
        audit.require(steps < 10000, 'Basis reduction did not terminate')
        for j in range(k - 1, -1, -1):
            mu, _ = gram_schmidt(gram)
            value = mu[k][j]
            quotient = (2 * value.numerator + value.denominator) // (2 * value.denominator)
            if quotient:
                elementary = np.eye(n, dtype=object)
                elementary[j, k] = -quotient
                gram = elementary.T @ gram @ elementary
                transform = transform @ elementary
        mu, norms = gram_schmidt(gram)
        if norms[k] >= (F(3, 4) - mu[k][k - 1] ** 2) * norms[k - 1]:
            k += 1
        else:
            gram[[k - 1, k], :] = gram[[k, k - 1], :]
            gram[:, [k - 1, k]] = gram[:, [k, k - 1]]
            transform[:, [k - 1, k]] = transform[:, [k, k - 1]]
            k = max(k - 1, 1)
    audit.require(abs(det_int(transform)) == 1, 'Non-unimodular transformation')
    audit.require(np.array_equal(gram, transform.T @ original @ transform),
                  'Incorrect Gram congruence')
    audit.require(positive_definite(gram), 'Reduction lost positive definiteness')
    return gram, transform


def matching_witness(m, class_map, d, determinant):
    """Return a full compatible group isomorphism, or a rejection witness.

    The spin origin is zero on both sides. Enumerating every unit includes
    every group isomorphism. No linking-form prefilter is applied here.
    """
    index = cyclic_index(class_map, determinant)
    if index is None:
        return {'verdict': 'OBSTRUCTED', 'reason': 'noncyclic discriminant group'}
    audit.require(len(m) == determinant and all(index[k] in m for k in range(determinant)),
                  'Incomplete m-function cannot be compared')
    failed = []
    for unit in range(1, determinant):
        if math.gcd(unit, determinant) != 1:
            continue
        for k in range(determinant):
            difference = m[index[k]] - d[unit * k % determinant]
            if difference < 0 or difference.denominator != 1 or difference.numerator % 2:
                failed.append({'unit': unit, 'class': k,
                               'difference': str(difference)})
                break
        else:
            return {'verdict': 'PASS', 'unit': unit,
                    'cyclic_labels': {str(k): list(v) for k, v in index.items()},
                    'm_in_cyclic_labels': {str(k): str(m[v]) for k, v in index.items()}}
    return {'verdict': 'OBSTRUCTED', 'failed_isomorphisms': failed}


def forms_for(determinant):
    # Re-enumerate rather than trusting an unversioned temporary cache.
    raw = candidate_plumbings(determinant, n_even=4)
    # Simple bases first: if a candidate passes, that suffices to show this
    # obstruction is inconclusive, without checking the remaining candidates.
    raw.sort(key=lambda item: (max(int(item[1][i, i]) for i in range(8)),
                               sum(int(item[1][i, i]) for i in range(8)), str(item[0])))
    return raw


def run(name, frozen, meta, max_forms=None):
    started = time.monotonic()
    entry = frozen['targets'][name]
    evidence = audit.validate_entry(entry)
    D, signature = entry['determinant'], entry['signature']
    audit.require(signature == -8 and entry['test'] == 'rank4', 'Unexpected target')
    diagnostics = {}
    sets, pairing, reference = audit.pin(name, entry['pd'], D, verbose=False,
                                        diagnostics=diagnostics)
    poly = audit.regina_jones(entry['pd'])
    audit.require(poly == audit.parse_jones(meta[name]['jones_polynomial']), 'Jones mismatch')
    jones = audit.jones_data(poly, D, signature)
    total, spin_count, survivors = audit.filter_vectors(
        sets, D, signature, F(jones['required_sum_d']))
    audit.require(len(survivors) == 1, 'Expected one exact correction-term vector')
    d = survivors[0][1]
    wanted = oriented(d, signature)
    audit.require(wanted[0] == -2, 'Incorrect oriented spin correction term')
    report = {
        'name': name, 'input': entry, 'lspace_evidence': evidence,
        'Jones_and_Casson_Walker': jones, 'pairing': str(pairing),
        'reference': reference, 'Greene_diagnostics': diagnostics,
        'candidate_sets': {str(k): list(map(str, sorted(v))) for k, v in sets.items()},
        'candidate_vectors': total, 'after_spin': spin_count, 'after_sum': len(survivors),
        'd': {str(k): str(v) for k, v in d.items()},
        'oriented_d': {str(k): str(v) for k, v in wanted.items()},
        'status': 'RUNNING', 'comparisons': [],
    }
    dest = HERE / (name + '_rank4.json')
    def save():
        report['seconds'] = round(time.monotonic() - started, 3)
        dest.write_text(json.dumps(report, indent=2) + '\n')
    forms = forms_for(D)
    report['candidate_forms'] = len(forms)
    report['forms_sha256'] = hashlib.sha256(json.dumps(
        [Q.tolist() for _, Q in forms], separators=(',', ':')).encode()).hexdigest()
    print(name, 'D', D, 'forms', len(forms), flush=True)
    save()
    cache = {}
    for number, (key, Q) in enumerate(forms, 1):
        if max_forms and number > max_forms:
            report['status'] = 'INCOMPLETE'
            break
        begin = time.monotonic()
        audit.require(det_int(Q) == D and positive_definite(Q), 'Invalid candidate form')
        reduced, transform = reduce_gram(Q)
        cache_key = tuple(map(int, reduced.flatten()))
        cache_hit = cache_key in cache
        if cache_hit:
            m, cm = cache[cache_key]
        else:
            m, cm = m_Q(reduced)
            cache[cache_key] = (m, cm)
        outcome = matching_witness(m, cm, wanted, D)
        item = {'number': number, 'key': str(key), 'Q': Q.tolist(),
                'reduced_Q': reduced.tolist(), 'P': transform.tolist(),
                'outcome': outcome, 'cache_hit': cache_hit,
                'seconds': round(time.monotonic() - begin, 3)}
        report['comparisons'].append(item)
        if outcome['verdict'] == 'PASS':
            report['status'] = 'PASS'
            report['range'] = entry['range']
            save()
            print(name, 'PASS', 'form', number, 'of', len(forms),
                  'seconds', report['seconds'], flush=True)
            break
        if number % 20 == 0:
            save()
            print(name, 'compared', number, 'of', len(forms),
                  'seconds', report['seconds'], flush=True)
    else:
        report['status'] = 'OBSTRUCTED'
        report['lower_bound'] = 5
        report['range'] = [5, entry['range'][1]]
    report['distinct_reduced_m_functions'] = len(cache)
    save()
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--knots', nargs='+', default=['13n_4876', '13n_4973', '13n_5102'])
    parser.add_argument('--max-forms', type=int,
                        help='diagnostic cutoff; yields INCOMPLETE, never an obstruction')
    args = parser.parse_args()
    frozen, _ = audit.consolidate.read_greene_sweep()
    meta = audit.metadata()
    for name in args.knots:
        run(name, frozen, meta, args.max_forms)


if __name__ == '__main__':
    main()
