#!/usr/bin/env python3
"""Compare the separately supplied v1.1 LaTeX source to the frozen computations.

This checks numerical consistency and transcription, not the validity of the
mathematical proofs. The parser expects the manuscript's current table syntax;
a format change should produce a failure requiring review, not a silent pass.
No manuscript file is written.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import json
from pathlib import Path
import re
from consolidate_results import RESULTS, primary_tag, read_json, reconstruct


def manuscript_counts(consolidated, table):
    """Read the manuscript's numerical premises without generating LaTeX files."""
    changed = consolidated['changed']
    totals = consolidated['counts']
    exact = [r for r in changed.values() if r['new'][0] == r['new'][1]]
    by_value_method = Counter((r['new'][0], primary_tag(r)) for r in exact)
    improved_methods = totals['improved_by_primary_method']
    scan = read_json(RESULTS/'open23/priority_u23_final.json')
    cyclic = read_json(RESULTS/'cyclic_cover/cyclic_cover_bound_2026-09-07.json')
    cyclic_degrees = Counter(cyclic[name]['n'] for name, r in changed.items() if primary_tag(r) == 'c')
    snapshot = read_json(RESULTS/'paper_v1_1_snapshot.json')
    sig4 = read_json(RESULTS/'owens_rank2/owens_verdicts_sigma4_alternating_2026-09-07.json')
    with gzip.open(RESULTS/'crossing_changes/dataset_v2.json.gz', 'rt') as stream:
        data_knots = {r['knot'] for r in json.load(stream)}
    return {
        'nChildren': len(scan['children']),
        'nCyclicNew': sum(cyclic_degrees.values()),
        'nCyclicThree': cyclic_degrees[3], 'nCyclicFour': cyclic_degrees[4],
        'nCyclicFive': cyclic_degrees[5], 'nDataKnots': len(data_knots),
        'nDichotomy': sum(heuristic == 0 for _, _, _, heuristic in scan['zero_candidate_knots']),
        'nExact': totals['exact'], 'nExactLower': totals['exact_lower'],
        'nImproved': totals['improved'], 'nImprovedL': improved_methods.get('L', 0),
        'nImprovedC': improved_methods.get('c', 0), 'nMcCoyKnots': improved_methods.get('K', 0),
        'nOpenTwoThreeAlt': len(scan['zero_candidate_knots']) + len(scan['knots_with_candidates']),
        'nRefKnots': len(table),
        'nSigFourImproved': sum(snapshot[r['name']] == [2, 4] and r['verdict'] == 'OBSTRUCTED' for r in sig4),
        'nSigFourObstructed': sum(snapshot[r['name']] == [2, 3] and r['verdict'] == 'OBSTRUCTED' for r in sig4),
        'nSigFourOpen': sum(snapshot[r['name']] == [2, 3] for r in sig4),
        'nUtwo': totals['exact_by_u']['2'], 'nUtwoL': by_value_method[2, 'L'],
        'nUtwoM': by_value_method[2, 'M'], 'nUtwoC': by_value_method[2, 'c'],
        'nUtwoG': by_value_method[2, 'G'],
        'nUthree': totals['exact_by_u']['3'], 'nUthreeC': by_value_method[3, 'c'],
        'nUthreeOwens': sum(by_value_method[3, tag] for tag in ('O2', 'a', 'OT', 'g')),
        'nUfour': totals['exact_by_u']['4'], 'nUfive': totals['exact_by_u']['5'],
        'nWithCandidates': len(scan['knots_with_candidates']) - len(scan.get('moved_to_tierB_after_resolution', [])),
    }


def require(condition, message):
    # Unlike assert, verification remains active under python -O.
    if not condition:
        raise ValueError(message)


def knotinfo_name(label):
    return re.sub(r'(\d+[an])(\d+)', r'\1_\2', label)


def names_in_table(section):
    body = section[section.index(r'\endhead'):section.index(r'\end{longtable}')]
    names = [knotinfo_name(label) for label in re.findall(r'\b\d+[an]\d+\b', body)]
    names.extend(f'{a}_{b}' for a, b in re.findall(r'\$(\d+)_\{(\d+)\}\$', body))
    require(len(set(names)) == len(names), 'Repeated knot in an appendix table')
    return set(names)


def verify(tex):
    consolidated, table = reconstruct()
    counts = manuscript_counts(consolidated, table)
    changed = consolidated['changed']
    macros = {name: int(value.replace('{,}', '')) for name, value in re.findall(
        r'\\newcommand\{\\(n\w+)\}\{(\d+(?:\{,\}\d+)*)\}', tex)}
    require(macros == counts, f'Manuscript count macros disagree: '
            f'{[(k, macros.get(k), counts.get(k)) for k in macros.keys() | counts.keys() if macros.get(k) != counts.get(k)]}')
    sections = re.split(r'\\section\{', tex[tex.index(r'\appendix'):])[1:]
    require(len(sections) == 8, 'Expected Appendices A--H')
    for u, section in zip((5, 4, 3, 2), sections[:4]):
        names = names_in_table(section)
        expected = {name for name, record in changed.items() if record['new'] == [u, u]}
        require(names == expected, f'Appendix u={u}: mismatching knots {names ^ expected}')

    # The improved-range appendix must match both interval ends, not just names.
    entries = {}
    pattern = (r'(\d+[an]\d+)(?:\$\^\{\\mathrm\{[a-z]+\}\}\$)?'
               r'\s*&\s*\$\[(\d+),(\d+)\]\$\s*&\s*\$\[(\d+),(\d+)\]\$')
    for match in re.finditer(pattern, sections[4]):
        name = knotinfo_name(match[1])
        require(name not in entries, f'Repeated improved-range entry: {name}')
        entries[name] = [list(map(int, match.group(2, 3))), list(map(int, match.group(4, 5)))]
    expected = {name: [record['reference'], record['new']] for name, record in changed.items()
                if record['new'][0] != record['new'][1]}
    require(entries == expected, 'Appendix E: missing knots or mismatching intervals')

    # These conditional statements must remain separate from exact values.
    scan = read_json(RESULTS/'open23/priority_u23_final.json')
    rigorous = {name for name, _, _, heuristic in scan['zero_candidate_knots'] if heuristic == 0}
    require(names_in_table(sections[5]) == rigorous, 'Appendix F disagrees with the certified candidate analysis')

    # Verify every four-tuple and every bold row of the actual manuscript.
    certificates = {r['knot']: r['best'] for r in read_json(RESULTS/'summary.json')}
    seen = set()
    for line in sections[6].splitlines():
        label = re.search(r'\\kn\{(\d+[an])\}\{(\d+)\}\s*&\s*(\d+)\s*&', line)
        if not label:
            continue
        name = f'{label[1]}_{label[2]}'
        require(name in certificates and name not in seen, f'Unexpected or repeated PD certificate: {name}')
        tuples = list(re.finditer(r'(\\mathbf\{)?\[(\d+),(\d+),(\d+),(\d+)\]', line))
        pd = [list(map(int, t.group(2, 3, 4, 5))) for t in tuples]
        marked = [i for i, t in enumerate(tuples) if t[1]]
        require(pd == certificates[name]['pd_knot'], f'PD mismatch: {name}')
        require(marked == certificates[name]['change_rows_0based'], f'Marked crossings mismatch: {name}')
        require(len(pd) == int(label[3]), f'Wrong crossing count: {name}')
        seen.add(name)
    require(seen == set(certificates), f'Missing PD certificates: {set(certificates)-seen}')

    # Appendix H prints one class from each conjugate pair. Compare its exact
    # rational values and unresolved candidate sets, not only the row labels.
    floer_tables = re.findall(r'\\begin\{longtable\}.*?\\end\{longtable\}',
                             sections[7], flags=re.S)
    require(len(floer_tables) == 2, 'Expected two correction-term tables in Appendix H')
    representative_classes = 0
    for name, section in zip(('12n_491', '13n_3370'), floer_tables):
        record = read_json(RESULTS/f'bernhard_jablan/greene_d_{name}_2026-09-08.json')
        values = record.get('d', record.get('pinned', {}))
        expected = {int(k): {Fraction(v)} for k, v in values.items()}
        expected.update({int(k): set(map(Fraction, v))
                         for k, v in record.get('ambiguous', {}).items()})
        expected = {k: v for k, v in expected.items() if k <= (record['det'] - 1)//2}
        body = section[section.index(r'\endlastfoot') + len(r'\endlastfoot'):]
        printed = {}
        for label, cell in re.findall(r'(\d+)\s*&\s*\$([^$]+)\$', body):
            k = int(label)
            require(k not in printed, f'Appendix H: repeated class {k} for {name}')
            cell = cell.replace(r'\{', '').replace(r'\}', '').replace(r'\,', '')
            printed[k] = set(map(Fraction, cell.split(',')))
        require(printed == expected, f'Appendix H: incorrect correction terms for {name}')
        representative_classes += len(printed)
    return {'count_macros': len(counts), 'appendix_entries': len(changed)+len(rigorous),
            'pd_certificates': len(seen), 'correction_term_classes': representative_classes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tex', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.tex.read_text())
    print(f"All {result['count_macros']} count macros, {result['appendix_entries']} Appendix A--F entries/ranges, "
          f"{result['pd_certificates']} Appendix G PD certificates, and "
          f"{result['correction_term_classes']} Appendix H correction-term classes agree.")


if __name__ == '__main__':
    main()
