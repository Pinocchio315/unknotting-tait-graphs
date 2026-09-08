"""Regression tests for mathematical-premise and result-publication failures.

Run with the repository Python environment:
python -m unittest discover -s results/review_v1_3 -p test_greene_audit.py -v
"""
import contextlib
import copy
import io
import json
from fractions import Fraction
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
GREENE = ROOT / 'lower_bounds' / 'greene'
sys.path.insert(0, str(GREENE))
import greene_sweep
import summarize_greene
from greene_ranks import rank_test


class PremiseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.entry = json.loads((ROOT / 'results/greene/greene_targets_2026-09-08.json')
                              .read_text())['targets']['13n_111']

    def test_actual_frozen_entry_passes(self):
        self.assertEqual(greene_sweep.validate_entry(self.entry)['rank'], 125)

    def assert_bad(self, field, value):
        entry = copy.deepcopy(self.entry)
        entry[field] = value
        with self.assertRaises((ValueError, KeyError)):
            greene_sweep.validate_entry(entry)

    def test_wrong_mod2_rank_rejected(self):
        self.assert_bad('reduced_mod2_vector_raw', '[[2, 124, 0, 0]]')

    def test_rational_khovanov_evidence_rejected(self):
        self.assert_bad('reduced_mod2_vector_raw', '[[0, 125, 0, 0]]')

    def test_missing_khovanov_evidence_rejected(self):
        self.assert_bad('reduced_mod2_vector_raw', '')

    def test_wrong_sign_pattern_rejected(self):
        self.assert_bad('signature', 2)

    def test_mislabeled_test_rejected(self):
        self.assert_bad('test', 'rank3')

    def test_control_requires_exact_tested_value(self):
        self.assert_bad('known_unknotting_number', 2)

    def test_partial_correction_terms_rejected(self):
        with self.assertRaises(ValueError):
            rank_test({0: Fraction(-1)}, 5, 2)

    def test_non_generator_pairing_rejected(self):
        with self.assertRaises(ValueError):
            rank_test({k: Fraction(-1) for k in range(5)}, 5, 2, pairing=0)

    def test_spin_guard_is_inconclusive(self):
        result = rank_test({k: Fraction(0) for k in range(5)}, 5, 2)
        self.assertEqual(result['verdict'], 'ERROR')


class PublicationTests(unittest.TestCase):
    def invoke_summary(self, records, controls=None):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            inputs = tmp / 'inputs.json'
            raw = tmp / 'results.jsonl'
            out = tmp / 'summary.json'
            inputs.write_text(json.dumps({'targets': {'K': {'range': [1, 2], 'test': 'u1'}},
                                         'controls': controls or {}}))
            raw.write_text('\n'.join(json.dumps(r) for r in records) + '\n')
            args = ['summarize_greene.py', str(raw), '--targets', str(inputs), '--out', str(out)]
            with patch.object(sys, 'argv', args), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as error:
                    summarize_greene.main()
            self.assertNotEqual(error.exception.code, 0)
            self.assertFalse(out.exists(), 'inconsistent inputs must not publish any bounds')

    def test_conflicting_completed_verdicts_do_not_publish(self):
        self.invoke_summary([{'name': 'K', 'verdict': 'OBSTRUCTED'},
                             {'name': 'K', 'verdict': 'PASS'}])

    def test_obstructed_control_cannot_hide_by_omitting_flag(self):
        self.invoke_summary([{'name': 'C', 'verdict': 'OBSTRUCTED'}],
                            controls={'C': {'range': [1, 1], 'test': 'u1'}})

    def test_malformed_result_is_not_silently_dropped(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'result.jsonl'
            path.write_text('{"name": "K", "verdict":')
            with self.assertRaisesRegex(ValueError, 'malformed result'):
                summarize_greene.load([path])

    def test_recorded_obstructed_control_blocks_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            inputs, out = tmp / 'inputs.json', tmp / 'results.jsonl'
            inputs.write_text(json.dumps({'targets': {}, 'controls': {'C': {}},
                                         'target_order': [], 'control_order': ['C']}))
            out.write_text(json.dumps({'name': 'C', 'verdict': 'OBSTRUCTED'}) + '\n')
            args = ['greene_sweep.py', '--controls', '--targets', str(inputs), '--out', str(out)]
            with patch.object(sys, 'argv', args):
                with self.assertRaisesRegex(SystemExit, 'refusing to resume'):
                    greene_sweep.main()

    def test_failed_worker_cannot_publish_stdout_obstruction(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            inputs, out = tmp / 'inputs.json', tmp / 'results.jsonl'
            inputs.write_text(json.dumps({'targets': {'K': {}}, 'controls': {},
                                         'target_order': ['K'], 'control_order': []}))
            failed = subprocess.CompletedProcess([], 1,
                stdout=json.dumps({'name': 'K', 'verdict': 'OBSTRUCTED'}) + '\n',
                stderr='calculation failed after partial output')
            args = ['greene_sweep.py', '--targets', str(inputs), '--out', str(out)]
            with patch.object(sys, 'argv', args), patch.object(subprocess, 'run', return_value=failed), \
                    contextlib.redirect_stdout(io.StringIO()):
                greene_sweep.main()
            self.assertEqual(json.loads(out.read_text())['verdict'], 'ERROR')


if __name__ == '__main__':
    unittest.main()
