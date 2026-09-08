#!/usr/bin/env python
"""Section 3, ``Homology of cyclic branched covers'': Nakanishi's bound.

A sequence of u crossing changes gives (n-1)u generators for H1(Sigma_n),
so u >= ceil(g_n/(n-1)). KnotInfo's torsion_numbers includes zero factors
for free summands, which must be counted. A successful bound is necessary
for unknotting; equality requires a separately established upper bound.

    python lower_bounds/cyclic_cover_bound.py --out /private/tmp/cyclic_cover.json

Importing this module never reads the database or writes deposited results.
"""
import argparse
import ast
import json
import os
import re

try:
    from .generator_bound_sweep import min_generators
except ImportError:
    from generator_bound_sweep import min_generators


def rng(value):
    digits = re.findall(r'\d+', str(value))
    return [int(digits[0]), int(digits[-1])] if digits else None


def cover_bound(covers):
    """Return (bound, attaining cover, generator counts), breaking ties by n.

    Integer ceiling arithmetic is exact even for large generator counts.
    Using p-ranks also handles cyclic decompositions not already in Smith
    invariant-factor order. Empty cover data contributes no lower bound.
    """
    generators = {}
    for n, factors in covers:
        if not isinstance(n, int) or n < 2:
            raise ValueError('a branched-cover degree must be an integer at least two')
        generators[n] = min_generators(factors)
    best = max((((g + n - 2) // (n - 1), -n)
                for n, g in generators.items()), default=(0, -2))
    return best[0], -best[1], generators


def scan_rows(rows):
    output = {}
    for row in rows:
        if not str(row.get('crossing_number', '')).strip().isdigit():
            continue
        notation = row.get('torsion_numbers', '').strip()
        if not notation:
            continue
        bound, degree, generators = cover_bound(ast.literal_eval(notation))
        output[row['name']] = {'bound': bound, 'n': degree, 'generators': generators,
                               'knotinfo_range': rng(row['unknotting_number'])}
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, help='destination for this run; deposited files are not defaults')
    args = parser.parse_args()
    import database_knotinfo as dk
    output = scan_rows(dk.link_list())
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, 'w') as handle:
        json.dump(output, handle)
    improvements = {name: row for name, row in output.items()
                    if row['knotinfo_range'] and row['bound'] > row['knotinfo_range'][0]}
    print(f'{len(output)} knots; {len(improvements)} bounds above the installed KnotInfo lower end')
    print('wrote', args.out)


if __name__ == '__main__':
    main()
