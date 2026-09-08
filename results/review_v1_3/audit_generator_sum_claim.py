#!/usr/bin/env python3
"""Audit an incoming, invalid additivity argument; not paper result generation.

For a finite abelian group G, its generator number is max_p dim(G/pG).
The p-ranks add under direct sum; the maximum generally does not. All
corrected pair counts below are labelled experimental and are not inserted
into manuscript totals. Input homology is the installed pinned KnotInfo row,
which is saved in the resulting review artifact for provenance.
"""
import ast
from collections import Counter
from importlib.metadata import version
from itertools import combinations_with_replacement
import json
from pathlib import Path
import database_knotinfo as dk

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def prime_divisors(value):
    value = abs(value)
    factors = set()
    p = 2
    while p * p <= value:
        if value % p == 0:
            factors.add(p)
            while value % p == 0:
                value //= p
        p += 1
    if value > 1:
        factors.add(value)
    return factors

def ranks(factors):
    if any(x == 0 for x in factors):
        raise ValueError('Double-cover homology of a knot must be finite')
    out = Counter()
    for factor in factors:
        out.update(prime_divisors(factor))
    return dict(out)

def main():
    table = json.loads((ROOT / 'generated/u_table.json').read_text())
    rows = {r['name']: r for r in dk.link_list()}
    sharp = []
    disagreements = []
    for name, bounds in table.items():
        if name not in rows or bounds[0] != bounds[1] or bounds[0] < 1:
            continue
        raw = str(rows[name].get('torsion_numbers', '')).strip()
        if not raw:
            continue
        homology = dict(ast.literal_eval(raw))
        if 2 not in homology:
            continue
        factors = homology[2]
        by_prime = ranks(factors)
        generator = max(by_prime.values(), default=0)
        naive = sum(x != 1 for x in factors)
        if generator != naive:
            disagreements.append(name)
        if bounds[0] == generator:
            sharp.append({'knot': name, 'u': bounds[0], 'h1_factors': factors,
                          'p_ranks': by_prime,
                          'sharp_primes': sorted(p for p, rank in by_prime.items()
                                                 if rank == bounds[0])})
    count = Counter()
    unsupported_examples = []
    for a, b in combinations_with_replacement(sharp, 2):
        count['all_unordered_pairs_with_repetition'] += 1
        common = set(a['sharp_primes']) & set(b['sharp_primes'])
        both_one = a['u'] == b['u'] == 1
        count['both_u1_pairs'] += both_one
        count['common_prime_certified_pairs'] += bool(common)
        count['common_prime_with_some_u_at_least_two'] += bool(common) and not both_one
        if not common and not both_one:
            count['not_certified_by_this_generator_argument_or_u1_prime_theorem'] += 1
            if len(unsupported_examples) < 20:
                primes = set(a['p_ranks']) | set(b['p_ranks'])
                joint = max(a['p_ranks'].get(p, 0) + b['p_ranks'].get(p, 0) for p in primes)
                unsupported_examples.append({'A': a['knot'], 'B': b['knot'],
                    'u_A_plus_u_B': a['u'] + b['u'], 'g_of_direct_sum': joint})
    lookup = {r['knot']: r for r in sharp}
    examples = []
    for an, bn in [('3_1', '4_1'), ('8_18', '9_40')]:
        a, b = lookup[an], lookup[bn]
        joint = max(a['p_ranks'].get(p, 0) + b['p_ranks'].get(p, 0)
                    for p in set(a['p_ranks']) | set(b['p_ranks']))
        examples.append({'A': a, 'B': b, 'g_direct_sum': joint,
                         'g_A_plus_g_B': a['u'] + b['u']})
    report = {'scope': 'Review of a rejected incoming additivity argument; not released paper counts',
              'database_package_version': version('database_knotinfo'),
              'sharp_knots': len(sharp), 'sharp_u1': sum(r['u'] == 1 for r in sharp),
              'factor_count_vs_actual_generator_disagreements': disagreements,
              'experimental_pair_counts': dict(count), 'examples': examples,
              'unsupported_examples': unsupported_examples, 'inputs': sharp}
    (HERE / 'generator_sum_claim_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print({k: report[k] for k in ['sharp_knots', 'sharp_u1', 'experimental_pair_counts']})
    print('Examples', [(r['A']['knot'], r['B']['knot'], r['g_direct_sum'], r['g_A_plus_g_B']) for r in examples])

if __name__ == '__main__':
    main()
