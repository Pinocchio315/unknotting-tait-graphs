#!/usr/bin/env python3
"""Review-only audit of an excluded connected-sum claim.

This experiment is excluded from the revised v1.3 manuscript and reporting
pipeline. It records a valid sufficient condition while documenting why the
original generator-additivity argument and its pair counts fail.

This is a standalone review experiment, excluded from manuscript reporting.

For a prime p, r_p(K) = dim H_1(Sigma_2(K); F_p) is additive under
connected sum and is at most u(K). Its maximum over primes need not be
additive. We certify a pair only when the SAME prime attains the known
unknotting number of both summands. Repeated summands are allowed; pairs
of tabulated names are unordered and are not counted as distinct knot types.

The reporting function reads only the deposited homology extract. The
optional --freeze-csv command authenticates the complete source CSV against
the archived public package metadata before making that extract.
"""
from __future__ import annotations
import argparse
import ast
from collections import Counter
import csv
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
INPUT = 'connected_sums/homology_inputs_2026-09-08.json.gz'


def prime_ranks(factors):
    """Compute ranks over each prime field from any finite cyclic decomposition."""
    ranks = Counter()
    for factor in factors:
        if type(factor) is not int or factor < 1:
            raise ValueError('Finite positive integral cyclic factors are required')
        n, p = factor, 2
        while p * p <= n:
            if n % p == 0:
                ranks[p] += 1
                while n % p == 0:
                    n //= p
            p += 1
        if n > 1:
            ranks[n] += 1
    return dict(ranks)


def sharp_primes(factors, unknotting_number):
    """Primes at which the homology bound reaches a supplied exact value."""
    if type(unknotting_number) is not int or unknotting_number < 1:
        raise ValueError('A positive exact unknotting number is required')
    ranks = prime_ranks(factors)
    if max(ranks.values(), default=0) > unknotting_number:
        raise ValueError('Homology lower bound contradicts the supplied exact value')
    return {p for p, rank in ranks.items() if rank == unknotting_number}


def freeze_csv(path, version='2026.8.1'):
    """Extract homology only after matching the archived official CSV checksum."""
    manifest = json.loads((RESULTS/'comparison/releases/manifest.json').read_text())
    source = manifest['releases'][version]
    raw = Path(path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != source['csv_sha256']:
        raise ValueError('Source CSV differs from the authenticated public snapshot')
    names = set(json.loads((RESULTS/'paper_v1_1_snapshot.json').read_text()))
    rows = {}
    csv.field_size_limit(sys.maxsize)
    for row in csv.DictReader(io.StringIO(raw.decode()), delimiter='|'):
        name = row['name'].strip()
        if name not in names:
            continue
        covers = {int(n): list(factors) for n, factors in ast.literal_eval(row['torsion_numbers'])}
        factors = covers[2]
        determinant = int(row['determinant'])
        product = 1
        for factor in factors:
            product *= factor
        if product != determinant:
            raise ValueError(f'Double-cover homology order differs from determinant: {name}')
        prime_ranks(factors)
        rows[name] = {'determinant': determinant, 'cyclic_factors': factors}
    if set(rows) != names:
        raise ValueError('The homology extract does not cover the full study population')
    document = {'schema_version': 1, 'evidence_kind': 'tabulated_homology_input',
                'source': {key: source[key] for key in
                    ('version', 'url', 'sha256', 'csv_member', 'csv_sha256')},
                'rows': rows}
    dest = RESULTS/INPUT
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(gzip.compress((json.dumps(document, sort_keys=True, indent=1)+'\n').encode(), mtime=0))


def certify_pairs(table, results=RESULTS):
    """Count formal pairs and retain a common prime witnessing each certificate."""
    results = Path(results)
    with gzip.open(results/INPUT, 'rt') as stream:
        frozen = json.load(stream)
    release = json.loads((results/'comparison/releases/manifest.json').read_text())['releases'][frozen['source']['version']]
    if any(frozen['source'][key] != release[key] for key in frozen['source']):
        raise ValueError('Frozen homology provenance differs from the public package record')
    if set(frozen['rows']) != set(table):
        raise ValueError('Homology input and unknotting table populations differ')
    sharp = {}
    for name, interval in table.items():
        if interval[0] != interval[1] or interval[0] < 1:
            continue
        u = interval[0]
        primes = sharp_primes(frozen['rows'][name]['cyclic_factors'], u)
        if primes:
            sharp[name] = {'u': u, 'primes': sorted(primes)}
    by_prime, by_sum = Counter(), Counter()
    count = nontrivial = 0
    # Assign each pair to its least common sharp prime, so multiple possible
    # certificates cannot inflate the total. The certificate applies to either
    # mirroring of either factor; those choices are not counted again.
    names = sorted(sharp)
    prime_sets = {name: set(sharp[name]['primes']) for name in names}
    for left, right in itertools.combinations_with_replacement(names, 2):
        common = prime_sets[left] & prime_sets[right]
        if not common:
            continue
        total_u = sharp[left]['u'] + sharp[right]['u']
        count += 1
        nontrivial += total_u > 2
        by_prime[min(common)] += 1
        by_sum[total_u] += 1
    return {'pair_convention': 'Unordered pairs of tabulated knot names with repetition; '
                'mirror choices are not counted separately; not a count of isotopy classes.',
            'criterion': 'A common prime p satisfies r_p(A)=u(A) and r_p(B)=u(B).',
            'counts': {'sharp_summands': len(sharp), 'common_prime_pairs': count,
                       'common_prime_pairs_with_summand_u_ge_2': nontrivial},
            'by_least_witness_prime': {str(p): n for p, n in sorted(by_prime.items())},
            'by_sum_unknotting_number': {str(u): n for u, n in sorted(by_sum.items())},
            'sharp_summands': sharp,
            'homology_input_sha256': hashlib.sha256((results/INPUT).read_bytes()).hexdigest(),
            'unknotting_table_sha256': hashlib.sha256(json.dumps(table, sort_keys=True,
                    separators=(',', ':')).encode()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze-csv', type=Path)
    parser.add_argument('--out', type=Path, default=RESULTS/'connected_sums/common_prime_pairs_2026-09-08.json')
    args = parser.parse_args()
    if args.freeze_csv:
        freeze_csv(args.freeze_csv)
    from consolidate_results import reconstruct
    _, table = reconstruct()
    result = certify_pairs(table)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps(result['counts'], indent=2))


if __name__ == '__main__':
    main()
