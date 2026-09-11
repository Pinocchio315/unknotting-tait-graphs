"""Check the mathematical boundary between recorded evidence and new bounds."""
from pathlib import Path
import copy
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from consolidate_results import Bound, completed_rank4, merge_bounds, reconstruct
from compare_knotinfo import compare, load_snapshot


class EvidenceTests(unittest.TestCase):
    def test_contradictory_intervals_are_rejected(self):
        with self.assertRaises(ValueError):
            merge_bounds({'K': [1, 2]}, [Bound('K', 'lo', 3, 'O2', 'completed')])

    def test_multiple_bounds_are_intersected_without_double_counting(self):
        table, changed = merge_bounds({'K': [1, 3]}, [
            Bound('K', 'lo', 2, 'L', 'pairing'),
            Bound('K', 'lo', 2, 'M', 'plumbing'),
            Bound('K', 'hi', 2, 'w', 'diagram')])
        self.assertEqual(table, {'K': [2, 2]})
        self.assertEqual(len(changed), 1)
        self.assertEqual(changed['K']['sources'], ['diagram', 'pairing', 'plumbing'])

    def test_unfinished_rank4_record_does_not_erase_completed_verdict(self):
        records = [{'name': 'K', 'verdict': 'OBSTRUCTED'},
                   {'name': 'K', 'verdict': 'TIMEOUT'}]
        self.assertEqual(completed_rank4(records), {'K': 'OBSTRUCTED'})

    def test_incompatible_completed_rank4_verdicts_are_rejected(self):
        with self.assertRaises(ValueError):
            completed_rank4([{'name': 'K', 'verdict': 'OBSTRUCTED'},
                             {'name': 'K', 'verdict': 'PASS'}])

    def test_conflicting_snapshot_is_not_counted_as_an_improvement(self):
        result = compare({'changed': {'K': {'new': [3, 3]}}},
                         {'K': [3, 3]}, {'K': [1, 2]})
        self.assertEqual(result['knots']['disjoint_intervals'], ['K'])
        self.assertEqual(result['counts']['unresolved_exact'], 0)

    def test_snapshot_population_must_match(self):
        with self.assertRaises(ValueError):
            compare({'changed': {}}, {'K': [1, 2]}, {'L': [1, 2]})


class FrozenComputationTests(unittest.TestCase):
    def test_completed_sweep_recovers_deposited_values_and_snapshot_partition(self):
        result, table = reconstruct()
        self.assertEqual(result['counts']['exact'], 2734)
        self.assertEqual(result['counts']['improved'], 380)
        for name in ('12n_288', '12n_491', '12n_501', '13n_3370', '13n_1587'):
            self.assertEqual(table[name], [2, 2], name)
        counts = compare(result, table, load_snapshot())['counts']
        self.assertEqual(counts['unresolved_exact'], 2540)
        self.assertEqual(counts['exact_agreements'], 194)
        self.assertEqual(counts['range_improvements'], 323)
        self.assertEqual(counts['range_agreements'], 57)
        self.assertEqual(counts['disjoint_intervals'], 0)

    def test_pre_sweep_profile_does_not_use_the_enlarged_sweep_bounds(self):
        before, _ = reconstruct(manuscript='v1.2')
        after, _ = reconstruct(manuscript='v1.3')
        self.assertEqual(before['counts']['exact'], 1227)
        self.assertEqual(after['counts']['exact'] - before['counts']['exact'], 1492)


if __name__ == '__main__':
    unittest.main()
