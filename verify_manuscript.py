#!/usr/bin/env python3
"""Compare separately supplied v1.1, v1.2 or v1.3 LaTeX to its frozen computations.

This checks numerical consistency and transcription, not the validity of the
mathematical proofs. The parser expects the manuscript's current table syntax;
a format change should produce a failure requiring review, not a silent pass.
No manuscript file is written.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
from consolidate_results import (RESULTS, paper_scan, appendix_tag, comparison_counts, comparison_report,
                                 manuscript_counts, manuscript_tex, read_json, reconstruct)


def manuscript_version(path):
    """Select the authenticated release from explicit generated-input paths."""
    path = Path(path)
    versions = set(re.findall(r'paper_v1_(2|3(?:_review)?)/', path.read_text()))
    if len(versions) > 1:
        raise ValueError('A manuscript mixes numerical inputs from different versions')
    if versions:
        return 'v1.' + versions.pop().replace('_', '-')
    match = re.search(r'v1\.([123])', path.name)
    if not match:
        raise ValueError('Cannot determine the manuscript version')
    return 'v1.' + match[1]


def load_manuscript(path):
    """Resolve local TeX inputs without executing LaTeX or shell commands.

    The paper's numerical inputs are authenticated by exact regeneration before
    expansion. A missing or stale file is an error, even if its visible numbers
    would happen to match a hand-edited copy elsewhere. Inclusion cycles and
    dynamic paths are rejected rather than silently omitting part of the paper.
    """
    path = Path(path).resolve()
    version = manuscript_version(path)
    consolidated, table = reconstruct(manuscript=version)
    expected = manuscript_tex(consolidated, table)
    inputs = {}

    def expand(source_path, stack):
        if source_path in stack:
            raise ValueError(f'Cyclic manuscript input: {source_path}')
        source = source_path.read_text()
        # An escaped percent is text, not the beginning of a TeX comment.
        source = re.sub(r'(?<!\\)%[^\n]*', '', source)
        def replace(match):
            name = match[1]
            if re.search(r'[\\{}]', name):
                raise ValueError(f'Dynamic manuscript input is unsupported: {name}')
            target = (path.parent / name).resolve()
            if not target.suffix:
                target = target.with_suffix('.tex')
            if target.name in expected and any(re.fullmatch(r'paper_v1_(?:2|3(?:_review)?)', p) for p in target.parts):
                require('paper_' + version.replace('.', '_').replace('-', '_') in target.parts,
                        f'Generated input belongs to a different manuscript version: {target}')
                require(target.read_text() == expected[target.name],
                        f'Stale generated manuscript input: {target}')
                inputs[target.name] = target
            return expand(target, stack + (source_path,))
        return re.sub(r'\\input\s*\{([^}]+)\}', replace, source)

    return expand(path, ()), inputs


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


def verify(tex, manuscript='v1.3'):
    consolidated, table = reconstruct(manuscript=manuscript)
    counts = manuscript_counts(consolidated, table)
    changed = consolidated['changed']
    definitions = re.findall(
        r'\\newcommand\{\\(n\w+)\}\{(\d+(?:\{,\}\d+)*)\}', tex)
    require(len(definitions) == len({name for name, _ in definitions}),
            'Repeated count macro: remove the copied block and retain its generated input')
    macros = {name: int(value.replace('{,}', '')) for name, value in definitions}
    # The historical standalone v1.1 source has no release-comparison macros.
    # If any v1.2 comparison is present, require the complete generated set.
    if 'nBaselineOverrides' in macros:
        counts.update(comparison_counts(comparison_report(consolidated, table)))
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
    pattern = (r'(\d+[an]\d+)(?:\$\^\{\\mathrm\{[A-Za-z]+\}\}\$)?'
               r'\s*&\s*\$\[(\d+),(\d+)\]\$\s*&\s*\$\[(\d+),(\d+)\]\$')
    for match in re.finditer(pattern, sections[4]):
        name = knotinfo_name(match[1])
        require(name not in entries, f'Repeated improved-range entry: {name}')
        entries[name] = [list(map(int, match.group(2, 3))), list(map(int, match.group(4, 5)))]
    expected = {name: [record['reference'], record['new']] for name, record in changed.items()
                if record['new'][0] != record['new'][1]}
    require(entries == expected, 'Appendix E: missing knots or mismatching intervals')

    # Upper-case G denotes Greene's model throughout A--E, including the new
    # higher lower bounds; lower-case g remains the homology generator bound.
    # v1.1 predates this notation. Later versions must match even if all G
    # superscripts were accidentally deleted, rather than skipping that check.
    if manuscript != 'v1.1':
        for index, section in enumerate(sections[:5]):
            marked = {knotinfo_name(name) for name in re.findall(
                r'(\d+[an]\d+)\$\^\{\\mathrm\{G\}\}\$', section)}
            expected_marked = {name for name, record in changed.items()
                if appendix_tag(record) == 'G'
                and (record['new'] == [[5, 5], [4, 4], [3, 3], [2, 2]][index]
                     if index < 4 else record['new'][0] != record['new'][1])}
            require(marked == expected_marked,
                    f'Appendix {chr(65 + index)}: Greene superscripts disagree with the records')

    # These conditional statements must remain separate from exact values.
    scan = read_json(RESULTS/paper_scan(consolidated))
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
    source, inputs = load_manuscript(args.tex)
    result = verify(source, manuscript=manuscript_version(args.tex))
    print(f"All {result['count_macros']} count macros, {result['appendix_entries']} Appendix A--F entries/ranges, "
          f"{result['pd_certificates']} Appendix G PD certificates, and "
          f"{result['correction_term_classes']} Appendix H correction-term classes agree.")
    if inputs:
        print(f'All {len(inputs)} generated numerical inputs agree with exact regeneration.')


if __name__ == '__main__':
    main()
