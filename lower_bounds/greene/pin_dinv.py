#!/usr/bin/env python
"""Bound correction terms using Greene's solitary states (arXiv:0805.1381).

For an L-space, its correction term in each Spin^c structure occurs among
that structure's solitary-state gradings on every marked diagram. Intersect
these necessary candidate sets over both colourings and all marked strands.
No identification of Spin^c structures is chosen arbitrarily: all compatible
linking-form isometries are retained, including conjugation.
The candidate sets can be larger than the actual correction-term data. A surgery
obstruction requires failure for every remaining vector; a fitting candidate
does not prove that the knot has an unknotting sequence of the tested length.

The command line verifies the L-space hypothesis from the archived reduced
Khovanov homology over F_2 before reporting a surgery obstruction. It prints
results by default; --out explicitly requests an output file for one knot.
"""
import argparse
from collections import Counter
import ast
import itertools
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

from greene_dinv import Diagram, surgery_test


def pages_of(pd):
    """Return every marked page, without hiding failed diagram assertions.

    Labels are first Chern classes in the cokernel of the Goeritz matrix.
    An odd cyclic cokernel identifies Spin^c structures with Z/D and sends
    the unique spin structure to zero. The current labelling implementation
    requires an explicit cyclic generator found by ``spinc_index``; if that
    search fails, no page or obstruction is returned.
    """
    labels = sorted({strand for crossing in pd for strand in crossing})
    counts = Counter(strand for crossing in pd for strand in crossing)
    if len(labels) != 2 * len(pd) or any(n != 2 for n in counts.values()):
        raise ValueError('Each strand label must occur twice in the PD code')
    out = {}
    for black in (0, 1):
        for mark in labels:
            diagram = Diagram(pd, black, mark)
            if diagram.D < 3 or diagram.D % 2 == 0:
                raise ValueError('This labelling requires odd determinant greater than one')
            classes = {}
            for state in diagram.states():
                row = diagram.analyse(state)
                if any((row['v'][i] - diagram.G[i, i]) % 2
                       for i in range(diagram.m)):
                    raise ValueError('A degree vector is not characteristic')
                classes.setdefault(row['label'], []).append(row)
            indices, pairing = diagram.spinc_index(list(classes))
            if indices is None or len(indices) != diagram.D:
                raise ValueError('The cyclic labelling does not cover every class')
            candidates = {
                indices[label]: sorted({row['gr'] for row in states
                                        if row['solitary']})
                for label, states in classes.items()
            }
            if any(not values for values in candidates.values()):
                raise ValueError('A Spin^c class has no solitary state')
            out[(black, mark)] = candidates, pairing
    return out


def pin(name, pd, D, verbose=True, diagnostics=None):
    """Intersect all pages, retaining every possible compatible alignment.

    A page with several allowable alignments contributes their union, not
    one arbitrarily selected alignment. Intersecting that union with the
    current candidates can only remove impossible values. Repeat until no
    set changes. Global correlations between alignments are deliberately
    relaxed; the surviving vectors therefore form a safe superset.
    """
    pages = pages_of(pd)
    if not pages or any(set(candidates) != set(range(D))
                        for candidates, _ in pages.values()):
        raise ValueError('The supplied determinant disagrees with the pages')
    reference = min(pages, key=lambda key: sum(
        len(values) > 1 for values in pages[key][0].values()))
    initial, pairing = pages[reference]
    current = {k: set(values) for k, values in initial.items()}
    units = [a for a in range(1, D) if math.gcd(a, D) == 1]
    alignments = {}
    passes = 0
    while True:
        before = {k: set(values) for k, values in current.items()}
        # The spin origin is canonical for odd determinant, so conjugation
        # is k -> -k and class identifications have no translation term.
        current = {k: values & current[-k % D]
                   for k, values in current.items()}
        for key, (candidates, other_pairing) in pages.items():
            possible = [a for a in units
                        if (a * a * other_pairing - pairing) % 1 == 0
                        and all(current[k] & set(candidates[a * k % D])
                                for k in range(D))]
            if not possible:
                raise ValueError(f'No compatible Spin^c alignment for page {key}')
            alignments[key] = possible
            for k in range(D):
                allowed = set().union(*(set(candidates[a * k % D])
                                        for a in possible))
                current[k] &= allowed
        current = {k: values & current[-k % D]
                   for k, values in current.items()}
        if any(not values for values in current.values()):
            raise ValueError('Page intersections are inconsistent')
        passes += 1
        if current == before:
            break
    ambiguous = {k: sorted(values) for k, values in current.items()
                 if len(values) != 1}
    if diagnostics is not None:
        diagnostics.update({
            'pages': len(pages), 'intersection_passes': passes,
            'alignment_units': {f'{colour},{mark}': values
                                for (colour, mark), values in alignments.items()},
            'initial_ambiguous_classes': sum(len(v) != 1 for v in initial.values()),
            'unresolved_classes': len(ambiguous),
        })
    if verbose:
        print(f'{name}: det {D}, {len(pages)} pages, reference {reference} '
              f'(pairing {pairing}), {passes} intersection passes; '
              f'unresolved classes: {len(ambiguous)}', flush=True)
    return current, pairing, reference


def khovanov_lspace_evidence(row):
    """Check the F_2 rank equality used by the branched-cover spectral sequence.

    Rational thinness alone is insufficient for this rank argument: the
    actual reduced mod-2 vector must be present. Its entries are
    [field, rank, homological grading, quantum grading] in KnotInfo.
    """
    raw = row.get('khovanov_reduced_mod2_vector', '')
    if not str(raw).strip():
        raise ValueError('Reduced Khovanov homology over F_2 is unavailable')
    entries = ast.literal_eval(raw)
    if not entries or any(len(entry) != 4 or entry[0] != 2
                          or any(type(value) is not int for value in entry)
                          or entry[1] < 0
                          for entry in entries):
        raise ValueError('Malformed reduced mod-2 Khovanov vector')
    rank = sum(entry[1] for entry in entries)
    determinant = int(row['determinant'])
    if rank != determinant:
        raise ValueError(f'The mod-2 rank {rank} is not determinant {determinant}; '
                         'the L-space hypothesis has not been certified')
    return {
        'method': 'reduced Khovanov homology over F_2 and branched-cover spectral sequence',
        'source_field': 'khovanov_reduced_mod2_vector',
        'rank': rank, 'determinant': determinant,
        'vector': entries,
        'diagonals_q_minus_2h': sorted({entry[3] - 2 * entry[2] for entry in entries}),
    }


def result_record(name, row, pd=None, input_provenance=None):
    """Produce a reproducible necessary surgery test, conditional only on cited inputs."""
    evidence = khovanov_lspace_evidence(row)
    if pd is None:
        pd = ast.literal_eval(row['pd_notation'])
    D = int(row['determinant'])
    if input_provenance is not None:
        evidence['provenance'] = input_provenance
    diagnostics = {}
    candidates, pairing, reference = pin(name, pd, D, diagnostics=diagnostics)
    ambiguous = {k: sorted(values) for k, values in candidates.items()
                 if len(values) != 1}
    pinned = {k: next(iter(values)) for k, values in candidates.items()
              if len(values) == 1}
    free = sorted(k for k in ambiguous if k <= -k % D)
    total, admitted = 0, 0
    fits_for_unique = None
    candidate_tests = []
    for choice in itertools.product(*(ambiguous[k] for k in free)):
        vector = dict(pinned)
        for k, value in zip(free, choice):
            vector[k] = vector[-k % D] = value
        if set(vector) != set(range(D)):
            raise ValueError('A candidate vector omits Spin^c classes')
        tested = surgery_test(vector, D)
        fits = tested['fits']
        candidate_tests.append({'choice': {str(k): str(v) for k, v in zip(free, choice)},
                                'surgery_test': tested})
        total += 1
        admitted += bool(fits)
        if not ambiguous:
            fits_for_unique = fits
    if total == 0:
        raise ValueError('No complete candidate vector was tested')
    record = {
        'schema_version': 2,
        'name': name, 'det': D, 'pd': pd, 'lam': str(pairing),
        'reference_marking': reference, 'lspace_evidence': evidence,
        'page_intersection': diagnostics,
        'implementation_sha256': {filename: hashlib.sha256(
            Path(__file__).with_name(filename).read_bytes()).hexdigest()
            for filename in ('greene_dinv.py', 'pin_dinv.py', 'half_integral.py')},
        'candidate_vectors': total, 'admitting_u1': admitted,
        'candidate_tests': candidate_tests,
        'verdict': 'OBSTRUCTED' if admitted == 0 else 'NOT_OBSTRUCTED',
    }
    if admitted == 0:
        record['lower_bound'] = 2
    if ambiguous:
        record['pinned'] = {str(k): str(value) for k, value in pinned.items()}
        record['ambiguous'] = {str(k): list(map(str, values))
                               for k, values in ambiguous.items()}
    else:
        record['d'] = {str(k): str(value) for k, value in pinned.items()}
        record['u1_fits'] = fits_for_unique
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('knots', nargs='+', help='KnotInfo names, e.g. 12n_491')
    parser.add_argument('--inputs', type=Path,
                        default=Path(__file__).resolve().parents[2] / 'results' /
                        'bernhard_jablan' / 'khovanov_inputs_2026-09-08.json',
                        help='Frozen PD and reduced mod-2 Khovanov input JSON')
    parser.add_argument('--out', type=Path,
                        help='Write one knot result to this explicit JSON path')
    args = parser.parse_args()
    if args.out and len(args.knots) != 1:
        parser.error('--out requires exactly one knot')
    import database_knotinfo as dk          # only the command line needs the database; the frozen
    rows = {row['name']: row for row in dk.link_list()}     # input file carries everything else
    frozen = json.loads(args.inputs.read_text())
    provenance = dict(frozen['source'],
                      input_file=args.inputs.name,
                      sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest())
    for name in args.knots:
        if name not in rows:
            parser.error(f'Unknown KnotInfo name: {name}')
        if name in frozen['knots']:
            entry = frozen['knots'][name]
            row = {'pd_notation': repr(entry['pd']),
                   'determinant': entry['determinant'],
                   'khovanov_reduced_mod2_vector': entry['reduced_mod2_vector_raw']}
            record = result_record(name, row, input_provenance=provenance)
        else:
            record = result_record(name, rows[name], input_provenance={
                'database': 'KnotInfo', 'package': 'database_knotinfo',
                'package_version': __import__('importlib.metadata').metadata.version('database_knotinfo'),
                'input_type': 'installed database row'})
        print(json.dumps(record, indent=2))
        if args.out:
            args.out.write_text(json.dumps(record, indent=2) + '\n')


if __name__ == '__main__':
    main()
