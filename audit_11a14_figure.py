#!/usr/bin/env python3
"""Replay the diagram-identification diagnostic for DKT Equation (6D.2).

Default execution reads all PD codes and dated reference values from the deposited
JSON, without importing a live KnotInfo database.  The three exterior comparisons
use SnapPy floating-point canonical isometry calculations, with at most 12 attempts
each; they are not interval-certified proofs.  The Jones-polynomial comparison is
exact and is deliberately restricted to the small 8- and 10-crossing references.
A matching Jones polynomial is never used to identify a knot.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
DEPOSIT = ROOT / 'results/comparison/11a14_figure_verification.json'
FIGURE = ROOT / 'upper_bounds/cluster/data/dkt_11a14_figure_pd.json'
sys.path.insert(0, str(ROOT / 'upper_bounds/cluster'))
sys.path.insert(0, str(ROOT / 'tait_graphs'))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def validate_pd(pd):
    if not isinstance(pd, list) or not pd or any(
            not isinstance(q, list) or len(q) != 4 or
            any(type(x) is not int for x in q) for q in pd):
        raise ValueError('Invalid PD code')
    if any(n != 2 for n in Counter(x for q in pd for x in q).values()):
        raise ValueError('Each PD arc label must occur twice')
    from xtait.graph import pd_passages
    pd_passages(pd)  # Reject disconnected traversals: all inputs must be knots.


def freeze_inputs():
    """One-time extraction; the public replay does not need this local source file."""
    import database_knotinfo
    from tait import knotinfo
    version = importlib.metadata.version('database_knotinfo')
    if version != '2026.8.1':
        raise ValueError('The provenance reconstruction requires the August 2026 package')
    csv_path = Path(database_knotinfo.__file__).parent / 'csv_data/knotinfo_data_complete.csv'
    released = json.loads((ROOT / 'results/comparison/releases/manifest.json').read_text())
    public_source = released['releases'][version]
    if sha256(csv_path) != public_source['csv_sha256']:
        raise ValueError('Installed CSV differs from the authenticated public August archive')
    figure = json.loads(FIGURE.read_text())
    references = {name: {'pd': knotinfo.pd_code(name),
                         'tabulated_unknotting_range': list(knotinfo.unknotting_interval(name))}
                  for name in ('11a_14', '8_8', '10_129')}
    inputs = {'figure': {'left_pd': figure['left_pd'], 'right_pd': figure['right_pd'],
                         'boxed_index_0based': figure['boxed_index']},
              'reference_knots': references}
    return {
        'schema_version': 1,
        'purpose': 'Identify the source, changed source, and claimed target of one published figure.',
        'sources': {
            'figure': {'paper': 'Dranowski--Kabkov--Tubbenhauer, arXiv:2603.07955v3',
                       'url': 'https://arxiv.org/abs/2603.07955v3',
                       'equation': '(6D.2)',
                       'figure_files': ['figs/example1.png', 'figs/example2.png'],
                       'extraction_note': figure['source'],
                       'original_local_record_sha256': sha256(FIGURE)},
            'reference_knots': {'package': 'database_knotinfo', 'version': version,
                                'archive_url': public_source['url'],
                                'archive_sha256': public_source['sha256'],
                                'csv_member': public_source['csv_member'],
                                'csv_sha256': public_source['csv_sha256'],
                                'evidence_kind': 'tabulated_input'}},
        'inputs': inputs, 'inputs_sha256': canonical_hash(inputs)}


def replay(document):
    # Import SnapPy before constructing Spherogram links to register .exterior().
    import snappy  # noqa: F401
    from spherogram import Link
    from upper_bounds.verify_braid_certificate import isometric
    from xtait.jones import jones, jones_mirror
    from xtait.graph import from_pd
    from flipdet import exact_det

    inputs = document['inputs']
    if canonical_hash(inputs) != document['inputs_sha256']:
        raise ValueError('Deposited PD inputs or reference values changed')
    figure, references = inputs['figure'], inputs['reference_knots']
    left, right = figure['left_pd'], figure['right_pd']
    index = figure['boxed_index_0based']
    if type(index) is not int or not 0 <= index < len(left):
        raise ValueError('Invalid boxed crossing index')
    for pd in [left, right, *(r['pd'] for r in references.values())]:
        validate_pd(pd)
    if references['8_8']['tabulated_unknotting_range'] != [2, 2] or \
            references['10_129']['tabulated_unknotting_range'] != [1, 1]:
        raise ValueError('Unexpected dated reference values')
    changed = [q[:] for q in left]
    # A cyclic shift by one position exchanges the over- and under-strands while
    # preserving the four incident arcs and the connectivity of the diagram.
    changed[index] = changed[index][1:] + changed[index][:1]
    validate_pd(changed)
    pairs = [('left_diagram', left, '11a_14'),
             ('left_after_boxed_change', changed, '8_8'),
             ('right_diagram', right, '10_129')]
    identifications = {}
    for label, pd, target in pairs:
        if not isometric(Link(pd), Link(references[target]['pd']), tries=12):
            raise ValueError(f'Exterior comparison inconclusive for {label} and {target}')
        identifications[label] = {'reference_knot': target,
                                  'exterior_isometry_found': True,
                                  'mirror_allowed': True,
                                  'determinant': exact_det(from_pd(pd))}
    # Never evaluate the 24-crossing figure by the exponential state sum.
    if len(references['8_8']['pd']) != 8 or len(references['10_129']['pd']) != 10:
        raise ValueError('Jones verification is restricted to the small reference diagrams')
    j8, j10 = jones(references['8_8']['pd']), jones(references['10_129']['pd'])
    if jones_mirror(j8) != j10:
        raise ValueError('The expected exact mirror Jones-polynomial equality failed')
    return {'identifications': identifications,
            'crossing_change': {'boxed_index_0based': index,
                                'original_pd_row': left[index],
                                'changed_pd_row': changed[index]},
            'jones': {'variable': 't', 'normalization': 'V(unknot)=1',
                      'coefficient_maps_by_integer_exponent': {
                          '8_8': {str(k): v for k, v in sorted(j8.items())},
                          '10_129': {str(k): v for k, v in sorted(j10.items())}},
                      'V_8_8_at_t_inverse_equals_V_10_129_at_t': True,
                      'arithmetic': 'exact integer state sum'},
            'identification_method': 'SnapPy floating-point canonical exterior isometry; '
                'at most 12 attempts per pair, mirrors allowed; not interval-certified.',
            'scope': 'This diagnostic concerns the displayed crossing change only. '
                'It does not disprove every proposed upper bound in the 104-entry list.',
            'conclusion': 'The recorded left-hand crossing change leads to 8_8 up to mirror, '
                'whereas the right-hand diagram represents 10_129. Their tabulated '
                'unknotting numbers are 2 and 1, respectively; their Jones polynomials '
                'agree after mirroring 8_8.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--freeze-inputs', action='store_true',
                        help='extract from the original local record and authenticated August CSV')
    parser.add_argument('--write', action='store_true', help='save replayed verification after success')
    args = parser.parse_args()
    document = freeze_inputs() if args.freeze_inputs else json.loads(DEPOSIT.read_text())
    verification = replay(document)
    if not args.write:
        if verification != document.get('verification'):
            raise ValueError('Deposited verification is absent or differs from the replay')
    else:
        document['verification'] = verification
        document['implementation_sha256'] = {
            str(p.relative_to(ROOT)): sha256(p) for p in [Path(__file__).resolve(),
                ROOT / 'upper_bounds/verify_braid_certificate.py',
                ROOT / 'upper_bounds/cluster/xtait/jones.py',
                ROOT / 'upper_bounds/cluster/xtait/graph.py',
                ROOT / 'upper_bounds/cluster/flipdet.py']}
        document['software'] = {name: importlib.metadata.version(name)
                                for name in ('snappy', 'spherogram')}
        DEPOSIT.write_text(json.dumps(document, indent=2, sort_keys=True) + '\n')
    print('Three exterior comparisons reproduced; exact mirror Jones-polynomial equality verified.')
    print('SnapPy identifications use floating-point canonical isometries, not interval certification.')


if __name__ == '__main__':
    main()
