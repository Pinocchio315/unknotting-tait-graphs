#!/usr/bin/env python3
"""Independent arithmetic replay of the 13n1619 half-integral obstruction.

This script uses only the standard library and the saved correction terms.
It does not import the production surgery test. The lens-space values are
obtained from their recursion, and the necessary conditions are checked in
two whole-vector stages instead of the production test's early-exit loop.
Monotonicity of V_j is unnecessary here: every map already fails the even
nonnegative-gap condition or equality of entries with the same V index.
"""
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import gcd
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


@cache
def lens_d(p, q, i):
    if p == 1:
        return F(0)
    return (-F(1, 4) + F((2 * i + 1 - p - q) ** 2, 4 * p * q)
            - lens_d(q, p % q, i % q))


def main():
    source = HERE / '13n_1619_certificate.json'
    record = json.loads(source.read_text())
    D = record['input']['determinant']
    d = {int(k): F(v) for k, v in record['d'].items()}
    if D != 33 or set(d) != set(range(D)):
        raise ValueError('This replay expects the complete 13n1619 vector')
    indices = [min(i // 2, (D + 1 - i) // 2) for i in range(D)]
    result = {'tested': 0, 'invalid_gaps': 0, 'gap_admissible_maps': [], 'fits': []}
    for sign in (1, -1):
        for a in range(1, D):
            if gcd(a, D) != 1:
                continue
            for b in range(D):
                result['tested'] += 1
                gaps = [(lens_d(D, 2, i) - sign * d[(a * i + b) % D]) / 2
                        for i in range(D)]
                if any(v < 0 or v.denominator != 1 for v in gaps):
                    result['invalid_gaps'] += 1
                    continue
                contradictions = [(i, j) for i, j in combinations(range(D), 2)
                                  if indices[i] == indices[j] and gaps[i] != gaps[j]]
                item = {'orientation': sign, 'a': a, 'b': b}
                if contradictions:
                    i, j = contradictions[0]
                    item['contradiction'] = {
                        'indices': [i, j], 'required_common_V_index': indices[i],
                        'incompatible_values': [str(gaps[i]), str(gaps[j])],
                    }
                else:
                    result['fits'].append(item)
                result['gap_admissible_maps'].append(item)
    if result['tested'] != 1320 or result['fits']:
        raise ValueError('The independent replay did not prove the claimed obstruction')
    result['certificate_sha256'] = hashlib.sha256(source.read_bytes()).hexdigest()
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE / 'independent_half_integral_replay.json').write_text(
        json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
