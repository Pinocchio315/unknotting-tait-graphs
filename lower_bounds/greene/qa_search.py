#!/usr/bin/env python
"""Find and replay a quasi-alternating resolution certificate.

A certificate records planar resolutions, seeded Spherogram diagram
simplifications, and exact determinants. Each node stores the PD code before
and after simplification, both free-circle counts, and the determinant. All
unlinked components removed by Spherogram are retained in these counts.
Leaves are the unknot or connected alternating diagrams of nonzero determinant.
Failure of the search gives no obstruction to being quasi-alternating.
"""
import argparse
import ast
import json
import os
from pathlib import Path
import sys
import time
import random

import spherogram
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'upper_bounds'))
from xtait.graph import from_pd


def relabel(pd):
    labels = sorted({label for crossing in pd for label in crossing})
    mapping = {label: index + 1 for index, label in enumerate(labels)}
    return [[mapping[label] for label in crossing] for crossing in pd]


def projection_connected(pd):
    """Disconnected nonempty projections represent split links."""
    if not pd:
        return False
    incidence = {}
    for i, crossing in enumerate(pd):
        for label in crossing:
            incidence.setdefault(label, []).append(i)
    if any(len(ends) != 2 for ends in incidence.values()):
        raise ValueError('Each strand label must occur twice')
    adjacency = {i: set() for i in range(len(pd))}
    for left, right in incidence.values():
        adjacency[left].add(right)
        adjacency[right].add(left)
    reached, pending = set(), [0]
    while pending:
        vertex = pending.pop()
        if vertex not in reached:
            reached.add(vertex)
            pending.extend(adjacency[vertex] - reached)
    return len(reached) == len(pd)


def tau_of_pd(pd):
    """Signed matrix-tree determinant for a connected nonempty projection."""
    if not pd:
        return 1
    if not projection_connected(pd):
        return 0
    graph = from_pd([list(crossing) for crossing in pd])
    vertices = list(graph.vertices)
    indices = {vertex: i for i, vertex in enumerate(vertices)}
    laplacian = sympy.zeros(len(vertices), len(vertices))
    for edge, (u, v) in graph.edges.items():
        if u == v:
            continue
        sign = graph.signs[edge]
        i, j = indices[u], indices[v]
        laplacian[i, i] += sign
        laplacian[j, j] += sign
        laplacian[i, j] -= sign
        laplacian[j, i] -= sign
    return int(laplacian[1:, 1:].det())


def determinant(pd, free_circles=0):
    """Keep zero-crossing components: a disjoint unknot makes a link split."""
    if not pd:
        return 1 if free_circles == 1 else 0
    if free_circles:
        return 0
    return abs(tau_of_pd(pd))


def smoothings(pd, crossing_index):
    """Return (resolved PD, new free circles) for each planar resolution.

    A union-find joins the two pairs of arc labels at the resolved crossing.
    A joined class absent from every remaining crossing is a free circle;
    retaining it is essential for distinguishing the unknot from an unlink.
    """
    a, b, c, d = pd[crossing_index]
    rest = [list(crossing) for i, crossing in enumerate(pd)
            if i != crossing_index]
    outputs = []
    for pairs in (((a, b), (c, d)), ((a, d), (b, c))):
        parent = {label: label for label in (a, b, c, d)}
        def find(label):
            while parent[label] != label:
                parent[label] = parent[parent[label]]
                label = parent[label]
            return label
        for left, right in pairs:
            parent[find(right)] = find(left)
        mapping = {label: find(label) for label in parent}
        resolved = [[mapping.get(label, label) for label in crossing]
                    for crossing in rest]
        remaining = {label for crossing in resolved for label in crossing}
        free = len(set(mapping.values()) - remaining)
        outputs.append((relabel(resolved), free))
    return outputs


def simplify(pd, free_circles=0):
    """Seeded diagram simplification, with all discarded circles retained.

    Spherogram counts components removed during simplification separately
    in unlinked_unknot_components. Omitting that counter can turn a split
    unlink into an apparent unknot and invalidate a QA certificate.
    """
    pd = relabel(pd)
    if not pd:
        return pd, free_circles
    link = spherogram.Link(pd)
    components_before = len(link.link_components) + link.unlinked_unknot_components
    # Reset the random source for each node so verification can replay the
    # same Reidemeister and strand-pickup moves. Restore the caller's RNG.
    random_state = random.getstate()
    random.seed(0)
    try:
        link.simplify('global')
    finally:
        random.setstate(random_state)
    components_after = len(link.link_components) + link.unlinked_unknot_components
    if components_before != components_after:
        raise ValueError('Simplification changed the component count')
    result = relabel(link.PD_code())
    circles = free_circles + link.unlinked_unknot_components
    if determinant(pd, free_circles) != determinant(result, circles):
        raise ValueError('Simplification changed the determinant')
    return result, circles


def qa(pd, stats=None, depth=0, free_circles=None, _memo=None):
    """Return a complete resolution-and-simplification certificate or None."""
    if stats is None:
        stats = {'nodes': 0}
    if _memo is None:
        _memo = {}
    if free_circles is None:
        free_circles = int(not pd)
    stats['nodes'] += 1
    input_pd, input_free = relabel(pd), free_circles
    pd, free_circles = simplify(input_pd, input_free)
    key = tuple(map(tuple, pd)), free_circles
    incoming = {'input_pd': input_pd, 'input_free_circles': input_free}
    if key in _memo:
        return None if _memo[key] is None else dict(_memo[key], **incoming)
    value = determinant(pd, free_circles)
    if value == 0:
        _memo[key] = None
        return None
    common = dict(incoming, pd=pd, free_circles=free_circles, determinant=value)
    if not pd:
        return dict(common, kind='unknot')
    if spherogram.Link(pd).is_alternating():
        return dict(common, kind='alternating')
    for crossing in range(len(pd)):
        resolutions = smoothings(pd, crossing)
        values = [determinant(child, free_circles + circles)
                  for child, circles in resolutions]
        if min(values) == 0 or sum(values) != value:
            continue
        children = [qa(child, stats, depth + 1, free_circles + circles, _memo)
                    for child, circles in resolutions]
        if all(child is not None for child in children):
            certificate = dict(common, kind='resolution', crossing=crossing,
                               resolution_determinants=values,
                               resolutions=children)
            _memo[key] = certificate
            return certificate
    _memo[key] = None
    return None


def verify_certificate(certificate, pd=None, free_circles=None):
    """Replay each resolution, retaining every component and checking all leaves."""
    if not isinstance(certificate, dict):
        raise ValueError('The certificate must record every PD node')
    if pd is None:
        pd = certificate['input_pd']
    if free_circles is None:
        free_circles = certificate['input_free_circles']
    pd = relabel(pd)
    if pd != certificate['input_pd'] or free_circles != certificate['input_free_circles']:
        raise ValueError('A certificate node differs from the actual resolution')
    pd, free_circles = simplify(pd, free_circles)
    if pd != certificate['pd'] or free_circles != certificate['free_circles']:
        raise ValueError('Stored PD disagrees with seeded Spherogram simplification')
    value = determinant(pd, free_circles)
    if value == 0 or value != certificate['determinant']:
        raise ValueError('Incorrect or zero determinant in certificate')
    kind = certificate['kind']
    if kind == 'unknot':
        if pd or free_circles != 1:
            raise ValueError('The leaf is not a single zero-crossing circle')
    elif kind == 'alternating':
        if free_circles or not pd or not spherogram.Link(pd).is_alternating():
            raise ValueError('The leaf is not a nonsplit alternating link')
    elif kind == 'resolution':
        crossing = certificate['crossing']
        if not isinstance(crossing, int) or not 0 <= crossing < len(pd):
            raise ValueError('Invalid resolution crossing')
        children = certificate['resolutions']
        if len(children) != 2:
            raise ValueError('A resolution needs exactly two children')
        values = []
        for child, (resolved, circles) in zip(children, smoothings(pd, crossing)):
            values.append(verify_certificate(child, resolved, free_circles + circles))
        if min(values) <= 0 or sum(values) != value:
            raise ValueError('The determinant sum condition fails')
        if values != certificate['resolution_determinants']:
            raise ValueError('Stored resolution determinants disagree')
    else:
        raise ValueError('Unknown certificate node kind')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('knots', nargs='*')
    parser.add_argument('--out', type=Path)
    parser.add_argument('--verify', type=Path, help='Replay an existing JSON certificate')
    args = parser.parse_args()
    if args.verify:
        record = json.loads(args.verify.read_text())
        value = verify_certificate(record['certificate'], record['pd'], 0)
        print(f'Certificate verified: determinant {value}')
        return
    if not args.knots or (args.out and len(args.knots) != 1):
        parser.error('Provide knots; --out requires exactly one knot')
    import database_knotinfo as dk
    rows = {row['name']: row for row in dk.link_list()}
    for name in args.knots:
        pd = ast.literal_eval(rows[name]['pd_notation'])
        stats, start = {'nodes': 0}, time.time()
        certificate = qa(pd, stats)
        record = {'schema_version': 2, 'name': name, 'pd': relabel(pd),
                  'determinant': int(rows[name]['determinant']),
                  'certificate': certificate,
                  'method': 'exact PD resolutions and seeded Spherogram global simplification',
                  'simplification_seed_per_node': 0,
                  'verification_status': 'verified' if certificate else 'no_certificate_found',
                  'spherogram_version': __import__('importlib.metadata').metadata.version('spherogram')}
        if certificate:
            verify_certificate(certificate, pd, 0)
        print(f'{name}: {"QUASI-ALTERNATING" if certificate else "no certificate"}; '
              f'{stats["nodes"]} nodes, {time.time() - start:.2f}s')
        if args.out:
            args.out.write_text(json.dumps(record, indent=2) + '\n')
        elif certificate:
            print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
