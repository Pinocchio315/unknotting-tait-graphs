"""Reject manuscript inputs that would conceal stale or incomplete results.

These tests exercise the boundary between authenticated mathematical records
and LaTeX inclusion. They do not rerun an obstruction or certify its proof.
"""
from pathlib import Path
import copy
import gzip
import json
import tempfile
import unittest
from unittest.mock import patch

import consolidate_results as reporting
from consolidate_results import (appendix_tag, manuscript_counts, manuscript_tex, reconstruct,
                                 read_greene_sweep, validate_sweep_record)
from verify_manuscript import load_manuscript, verify


class ManuscriptInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        consolidated, table = reconstruct(manuscript="v1.3-review")
        cls.generated = manuscript_tex(consolidated, table)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.inputs = self.root / 'code' / 'generated' / 'paper_v1_3_review'
        self.inputs.mkdir(parents=True)
        for name, text in self.generated.items():
            (self.inputs / name).write_text(text)
        self.manuscript = self.root / 'main_v1.3.tex'

    def include(self, *names):
        self.manuscript.write_text('\n'.join(
            rf'\input{{code/generated/paper_v1_3_review/{name}}}' for name in names))

    def test_authentic_inputs_include_both_greene_knots(self):
        self.include('counts.tex', 'rows_u2.tex')
        source, included = load_manuscript(self.manuscript)
        self.assertEqual(set(included), {'counts.tex', 'rows_u2.tex'})
        self.assertIn(r'12n491$^{\mathrm{G}}$', source)
        self.assertIn(r'13n3370$^{\mathrm{G}}$', source)

    def test_stale_count_or_missing_greene_case_is_rejected(self):
        self.include('counts.tex', 'rows_u2.tex')
        changes = {
            'counts.tex': lambda text: text.replace('2{,}719', '2{,}717'),
            'rows_u2.tex': lambda text: text.replace(r'13n3370$^{\mathrm{G}}$', ''),
        }
        for name, change in changes.items():
            with self.subTest(name=name):
                path = self.inputs / name
                original = path.read_text()
                changed = change(original)
                self.assertNotEqual(original, changed)
                path.write_text(changed)
                with self.assertRaisesRegex(ValueError, 'Stale generated'):
                    load_manuscript(self.manuscript)
                path.write_text(original)

    def test_missing_input_and_inclusion_cycle_are_rejected(self):
        self.include('absent.tex')
        with self.assertRaises(FileNotFoundError):
            load_manuscript(self.manuscript)
        self.manuscript.write_text(r'\input{main_v1.3.tex}')
        with self.assertRaisesRegex(ValueError, 'Cyclic manuscript input'):
            load_manuscript(self.manuscript)

    def test_commented_input_does_not_hide_missing_file(self):
        self.manuscript.write_text('% \\input{absent.tex}\nText.\n')
        source, included = load_manuscript(self.manuscript)
        self.assertEqual(source.strip(), 'Text.')
        self.assertFalse(included)

    def test_copied_macro_block_is_not_silently_accepted(self):
        source = self.generated['counts.tex'] + r'\newcommand{\nExact}{1{,}227}'
        with self.assertRaisesRegex(ValueError, 'Repeated count macro'):
            verify(source, manuscript="v1.3-review")

    def test_secondary_generator_bound_keeps_the_primary_pairing_notation(self):
        # The linking pairing remains unmarked at lower bound two, both for
        # exact values and ranges. Only the established lower-bound-three
        # convention highlights a secondary generator-bound explanation.
        for interval in ([2, 2], [2, 3]):
            with self.subTest(interval=interval):
                self.assertEqual(appendix_tag({'new': interval,
                    'lower_tags': ['L', 'g'], 'upper_tags': []}), '')
        for interval in ([3, 3], [3, 5]):
            with self.subTest(interval=interval):
                self.assertEqual(appendix_tag({'new': interval,
                    'lower_tags': ['O2', 'g'], 'upper_tags': []}), 'g')

    def test_historical_inputs_reconstruct_without_the_new_sweep(self):
        consolidated, table = reconstruct(manuscript='v1.2')
        self.assertEqual(consolidated['counts']['exact'], 1227)
        self.assertEqual(consolidated['counts']['improved'], 176)
        old = manuscript_tex(consolidated, table)
        self.assertNotIn(r'\nSweep', old['counts.tex'])
        self.assertNotIn(r'13n1587$^{\mathrm{G}}$', old['rows_u2.tex'])

    def test_inputs_from_different_versions_cannot_mix(self):
        self.manuscript.write_text(
            r'\input{code/generated/paper_v1_2/counts.tex}' + '\n' +
            r'\input{code/generated/paper_v1_3_review/rows_u2.tex}')
        with self.assertRaisesRegex(ValueError, 'mixes numerical inputs'):
            load_manuscript(self.manuscript)

    def test_jones_only_moved_parents_remain_outside_the_theorem(self):
        consolidated, table = reconstruct(manuscript='v1.3')
        counts = manuscript_counts(consolidated, table)
        self.assertEqual((counts['nDichotomy'], counts['nDichotomyHeuristic'],
                          counts['nWithCandidates']), (959, 46, 22))
        original_read = reporting.read_json
        def poisoned(path):
            data = original_read(path)
            if Path(path).name == 'priority_u23_final.json':
                data['zero_candidate_knots'].append(
                    [data['moved_to_tierB_after_resolution'][0], 0, 0, 0])
            return data
        with patch.object(reporting, 'read_json', side_effect=poisoned):
            with self.assertRaisesRegex(ValueError, 'categories overlap'):
                manuscript_counts(consolidated, table)


class SweepBookkeepingTests(unittest.TestCase):
    """Toy records test rejection rules, not any mathematical obstruction."""

    def setUp(self):
        self.entry = {'name': 'toy', 'determinant': 3, 'signature': 0,
                      'range': [1, 2], 'test': 'u1', 'khovanov_rank': 3,
                      'reduced_mod2_vector_raw': '[[2,3,0,0]]'}
        self.record = {'name': 'toy', 'role': 'target', 'det': 3, 'sigma': 0,
                       'range': [1, 2], 'test': 'u1', 'excluded_length': 1,
                       'd': {'0': '0', '1': '0', '2': '0'}, 'candidate_vectors': 1,
                       'admitting_vectors': 0, 'lower_bound': 2, 'verdict': 'OBSTRUCTED',
                       'candidate_verdicts': [{'test': 'u1', 'admits': False,
                                               'detail': {'fits': []}}]}

    def test_changed_premise_or_incomplete_class_set_is_rejected(self):
        validate_sweep_record(self.record, self.entry, 'target')
        entry = copy.deepcopy(self.entry)
        entry['reduced_mod2_vector_raw'] = '[[2,5,0,0]]'
        with self.assertRaisesRegex(ValueError, 'premise'):
            validate_sweep_record(self.record, entry, 'target')
        record = copy.deepcopy(self.record)
        del record['d']['2']
        with self.assertRaisesRegex(ValueError, 'candidate classes'):
            validate_sweep_record(record, self.entry, 'target')

    def test_every_candidate_must_have_a_completed_comparison(self):
        for mutation in ('omit_test', 'hidden_fit', 'ambiguous_pair'):
            record = copy.deepcopy(self.record)
            if mutation == 'omit_test':
                record['candidate_verdicts'] = []
            elif mutation == 'hidden_fit':
                record['candidate_verdicts'][0]['detail']['fits'] = [{'a': 1}]
            else:
                del record['d']
                record['pinned'] = {'0': '0'}
                record['ambiguous'] = {'1': ['0', '2'], '2': ['0', '2']}
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                validate_sweep_record(record, self.entry, 'target')

    def test_timeout_is_not_a_completed_obstruction_and_control_cannot_be_obstructed(self):
        validate_sweep_record({'name': 'toy', 'role': 'target', 'verdict': 'TIMEOUT'},
                              self.entry, 'target')
        entry, record = copy.deepcopy(self.entry), copy.deepcopy(self.record)
        entry.update(range=[1, 1], known_unknotting_number=1)
        record.update(role='control', range=[1, 1])
        with self.assertRaisesRegex(ValueError, 'obstructs the control'):
            validate_sweep_record(record, entry, 'control')

    def test_duplicate_jobs_cannot_inflate_the_sweep(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / 'greene'
            base.mkdir()
            (base / 'greene_targets_2026-09-08.json').write_text(json.dumps(
                {'targets': {'toy': self.entry}, 'controls': {}}))
            with gzip.open(base / 'sweep_2026-09-08.jsonl.gz', 'wt') as stream:
                stream.write((json.dumps(self.record) + '\n') * 2)
            with self.assertRaisesRegex(ValueError, 'Duplicate frozen'):
                read_greene_sweep(directory)


if __name__ == '__main__':
    unittest.main()
