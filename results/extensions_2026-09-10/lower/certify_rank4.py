#!/usr/bin/env python3
"""Certify rank-four rejections with explicit characteristic covectors.

This verification does not use m_Q, shortest-vector enumeration, ClassMap,
cyclic_index, or the matching code. It only evaluates explicit covectors.
One violating covector for each group isomorphism is sufficient: m_Q is a
minimum, so a representative gives an upper bound on it. Its congruence
modulo two is independent of the representative in the characteristic class.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import itertools
import json
import math
import time
from fractions import Fraction as F
from functools import reduce
from pathlib import Path

import numpy as np
from sympy import Matrix

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prime_divisors(number):
    result, p = [], 2
    while p * p <= number:
        if number % p == 0:
            result.append(p)
            while number % p == 0:
                number //= p
        p += 1
    if number > 1:
        result.append(number)
    return result


def cyclic_functional(adjugate, D):
    """An explicit surjection coker(Q) -> Z/D, or a noncyclic certificate.

    The gcd of the (r-1)-minors detects cyclicity. If it is one, for each
    prime p|D choose a row of adj(Q) nonzero mod p. CRT combines these rows
    into a functional primitive modulo every prime dividing D.
    """
    divisor = math.gcd(D, *map(int, adjugate.flatten()))
    if divisor != 1:
        return None, divisor
    primes = prime_divisors(D)
    radical = math.prod(primes)
    row = np.zeros(len(adjugate), dtype=object)
    for p in primes:
        chosen = next(r for r in adjugate if any(int(v) % p for v in r))
        cofactor = radical // p
        coefficient = cofactor * pow(cofactor, -1, p)
        row += coefficient * chosen
    row %= D
    require(math.gcd(D, *map(int, row)) == 1, 'Nonprimitive cyclic functional')
    return row, 1


def small_covector_certificate(Q, d, D, radius):
    """Search a finite cube for sufficient witnesses, never infer from failure."""
    rank = len(Q)
    # Sympy's exact inverse is independent of the production Bareiss helpers.
    q = Matrix(Q.tolist())
    require(int(q.det()) == D, 'Wrong determinant')
    adj = np.array((q.inv() * D).tolist(), dtype=object)
    require(all(v == int(v) for v in adj.flatten()), 'Nonintegral adjugate')
    adj = np.array(adj, dtype=np.int64)
    functional, divisor = cyclic_functional(adj, D)
    if functional is None:
        return {'noncyclic_minor_gcd': divisor}
    require(all(v % D == 0 for v in functional @ Q), 'Functional does not vanish on QZ^r')
    parity = np.array([int(Q[i, i]) % 2 for i in range(rank)], dtype=np.int64)
    lattice = np.array(list(itertools.product(range(-radius, radius + 1), repeat=rank)),
                       dtype=np.int64)
    vectors = 2 * lattice + parity
    euclidean_norms = np.einsum('ni,ni->n', vectors, vectors)
    require(max(map(abs, adj.flatten())) * rank ** 2 * (2 * radius + 1) ** 2 < 2 ** 62,
            'Arithmetic might overflow int64')
    norms = np.einsum('ni,ij,nj->n', vectors, adj, vectors)
    labels = (vectors @ np.array(functional, dtype=np.int64)) % D
    scaled_d = [4 * D * d[k] for k in range(D)]
    require(all(v.denominator == 1 for v in scaled_d), 'Unexpected d denominator')
    scaled_d = np.array(list(map(int, scaled_d)), dtype=np.int64)
    selected, selections, pool = {}, {}, []
    for unit in range(1, D):
        if math.gcd(unit, D) != 1:
            continue
        differences = norms - rank * D - scaled_d[(unit * labels) % D]
        invalid = np.flatnonzero((differences < 0) | (differences % (8 * D) != 0))
        if not len(invalid):
            return None
        # Smallest Euclidean witness is easy to read and fast to replay.
        k = int(invalid[np.argmin(euclidean_norms[invalid])])
        vector = tuple(map(int, vectors[k]))
        if vector not in selected:
            selected[vector] = len(pool)
            pool.append(list(vector))
        selections[str(unit)] = selected[vector]
    return {'cyclic_functional': list(map(int, functional)), 'covectors': pool,
            'unit_witness_indices': selections, 'radius': radius}


def verify_certificate(Q, d, D, proof):
    """Verify each explicit exclusion using only exact matrix arithmetic."""
    q = Matrix(Q.tolist())
    require(int(q.det()) == D, 'Wrong certificate determinant')
    inverse = q.inv()
    adj = inverse * D
    if 'noncyclic_minor_gcd' in proof:
        divisor = math.gcd(D, *map(int, adj))
        require(divisor == proof['noncyclic_minor_gcd'] and divisor > 1,
                'Invalid noncyclic certificate')
        return
    row = Matrix([proof['cyclic_functional']])
    require(math.gcd(D, *map(int, row)) == 1 and all(int(v) % D == 0 for v in row * q),
            'Invalid cyclic coordinates')
    units = {str(a) for a in range(1, D) if math.gcd(a, D) == 1}
    require(set(proof['unit_witness_indices']) == units, 'Missing group isomorphism')
    # Many units use the same covector. Evaluate each norm and characteristic
    # parity only once; the remaining checks are scalar rational arithmetic.
    values = []
    for covector in proof['covectors']:
        require(all((int(covector[i]) - int(q[i, i])) % 2 == 0 for i in range(len(Q))),
                'Non-characteristic covector')
        norm = sum(F(inverse[i, j]) * covector[i] * covector[j]
                   for i in range(len(Q)) for j in range(len(Q)))
        label = sum(int(row[i]) * covector[i] for i in range(len(Q))) % D
        values.append(((norm - len(Q)) / 4, label))
    for a, number in proof['unit_witness_indices'].items():
        value, label = values[number]
        difference = value - d[int(a) * label % D]
        require(difference < 0 or difference.denominator != 1 or difference.numerator % 2,
                'Covector does not exclude the assigned isomorphism')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('names', nargs='+')
    args = parser.parse_args()
    for name in args.names:
        source = HERE / (name + '_rank4.json')
        report = json.loads(source.read_text())
        require(report['status'] == 'OBSTRUCTED', 'Comparison has not finished with an obstruction')
        require(len(report['comparisons']) == report['candidate_forms'], 'Missing surgery forms')
        D = report['input']['determinant']
        d = {int(k): F(v) for k, v in report['oriented_d'].items()}
        cache, proofs = {}, []
        start = time.monotonic()
        for number, item in enumerate(report['comparisons'], 1):
            Q, R, P = (Matrix(item[k]) for k in ('Q', 'reduced_Q', 'P'))
            require(abs(int(P.det())) == 1 and R == P.T * Q * P,
                    'Invalid reduction congruence')
            key = tuple(R)
            if key not in cache:
                array = np.array(R.tolist(), dtype=np.int64)
                proof = small_covector_certificate(array, d, D, radius=1)
                require(proof is not None, 'The cube does not yet certify this form')
                verify_certificate(array, d, D, proof)
                cache[key] = len(proofs)
                proofs.append({'reduced_Q': R.tolist(), 'proof': proof})
            item['independent_proof_index'] = cache[key]
            if number % 200 == 0:
                print(name, 'certified', number, 'of', report['candidate_forms'], flush=True)
        payload = {'name': name, 'input_report_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                   'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'candidate_forms': report['candidate_forms'], 'distinct_proofs': len(proofs),
                   'form_to_proof': [r['independent_proof_index'] for r in report['comparisons']],
                   'proofs': proofs, 'seconds': round(time.monotonic() - start, 3),
                   'method': 'explicit characteristic covectors; no shortest-vector computation'}
        # Convert Sympy integers to Python integers for this exact integer proof.
        with gzip.open(HERE / (name + '_rank4_covectors.json.gz'), 'wt') as f:
            json.dump(payload, f, default=int, separators=(',', ':'))
        print(name, 'independent certificate complete', payload['seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
