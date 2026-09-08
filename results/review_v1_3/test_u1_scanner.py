"""Focused regressions for the chosen-diagram scanner and evidence handling."""
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'crossing_changes'))
import u1_minimal_diagram_scan as scanner
import database_knotinfo as dk


class ScannerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = {r['name']: r for r in dk.link_list()}

    def test_positive_witness_has_evidence(self):
        result = scanner.scan_row(self.rows['3_1'])
        self.assertEqual(result['verdict'], 'unknotting crossing')
        self.assertEqual(result['scope'], 'one tabulated minimal diagram')
        self.assertEqual(len(result['crossing_evidence']), 3)
        for edge in result['unknotting_crossings']:
            self.assertEqual(result['crossing_evidence'][edge]['determinant'], 1)

    def test_determinant_one_does_not_mean_unknot(self):
        result = scanner.scan_row(self.rows['11n_116'])
        self.assertEqual(result['verdict'], 'none in this diagram')
        evidence = result['crossing_evidence'][8]
        self.assertEqual(evidence['determinant'], 1)
        self.assertGreater(evidence['seifert_genus'], 0)

    def test_hfk_failure_is_undecided(self):
        with patch('knot_floer_homology.pd_to_hfk', side_effect=RuntimeError('unfinished HFK')):
            result = scanner.scan_row(self.rows['11n_116'])
        self.assertEqual(result['verdict'], 'undecided')
        self.assertIn(8, result['undecided'])

    def test_nonminimal_input_is_rejected(self):
        row = dict(self.rows['3_1'], crossing_number='4')
        with self.assertRaisesRegex(ValueError, 'minimal crossing number'):
            scanner.scan_row(row)

    def test_wrong_determinant_is_rejected(self):
        row = dict(self.rows['3_1'], determinant='5')
        with self.assertRaisesRegex(ValueError, 'determinant disagree'):
            scanner.scan_row(row)

    def test_malformed_hfk_is_rejected(self):
        valid = {'seifert_genus': 0, 'modulus': 2, 'ranks': {(0, 0): 1}, 'total_rank': 1}
        self.assertEqual(scanner.hfk_genus(valid), 0)
        for edit in ({'seifert_genus': None}, {'seifert_genus': False},
                     {'ranks': {(0, 0): 3}, 'total_rank': 3},
                     {'seifert_genus': 1}, {'modulus': 0}, {'total_rank': 0}):
            with self.subTest(edit=edit), self.assertRaises(ValueError):
                scanner.hfk_genus(dict(valid, **edit))

    def test_named_cli_does_not_overwrite_published_cohort(self):
        archived = ROOT / 'results/crossing_changes/u1_minimal_diagram_scan_2026-09-08.json'
        before = hashlib.sha256(archived.read_bytes()).hexdigest()
        with patch.object(sys, 'argv', ['scan', '3_1', '3_1']), contextlib.redirect_stdout(io.StringIO()) as out:
            scanner.main()
        payload = json.loads(out.getvalue())
        self.assertEqual(list(payload['knots']), ['3_1'])
        self.assertEqual(hashlib.sha256(archived.read_bytes()).hexdigest(), before)

    def test_open_cohort_has_its_own_selection_and_explicit_output(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            table, out = tmp / 'table.json', tmp / 'out.json'
            table.write_text(json.dumps({'3_1': [1, 1], '11n_116': [1, 2], '5_1': [2, 2]}))
            with patch.object(sys, 'argv', ['scan', '--cohort', 'open', '--table', str(table), '--out', str(out)]), \
                    contextlib.redirect_stdout(io.StringIO()):
                scanner.main()
            payload = json.loads(out.read_text())
            self.assertEqual(list(payload['knots']), ['11n_116'])
            self.assertEqual(payload['counts'], {'none in this diagram': 1})


if __name__ == '__main__':
    unittest.main()
