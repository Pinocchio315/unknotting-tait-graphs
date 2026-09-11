"""Regression tests for the new bounds, joint attribution, and dated comparison.

These tests read the authenticated mathematical records. They do not substitute
for the separate topological/covector replays documented with those records.
"""
from collections import Counter
import json
from pathlib import Path
import sys
import unittest

CODE = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(CODE), str(CODE/'upper_bounds')]
import consolidate_results as c
import v18_results as v18
from verify_presentation_diagrams import load_paper_diagrams, load_additional_diagrams


class V18Results(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old, cls.old_table = c.reconstruct(manuscript='v1.3')
        cls.new, cls.table = c.reconstruct(manuscript='v1.8')

    def test_exactly_fifteen_changes_from_the_preserved_profile(self):
        changed = {n for n in self.table if self.table[n] != self.old_table[n]}
        lower = {'13n_1619':2, '13n_4876':5, '13n_4973':5, '13n_5102':5}
        upper = {n for n, _, _, _ in load_additional_diagrams()}
        self.assertEqual(changed, set(lower) | upper)
        for name in changed:
            self.assertLess(*self.old_table[name])
            value = lower.get(name, 2)
            self.assertEqual(self.table[name], [value, value])
        self.assertEqual((self.old['counts']['exact'], self.old['counts']['improved']), (2719,390))

    def test_joint_contributions_are_not_counted_as_lower_only(self):
        total = self.new['counts']
        lower_changed = self.new['lower_stage_changed']
        joint = {n for n,r in self.new['changed'].items()
                 if r['new'][0] == r['new'][1] and r['lower_tags'] and 'w' in r['upper_tags']}
        expected = {n for n, _, _, _ in load_additional_diagrams()} - {'13n_3733'}
        self.assertEqual(joint, expected)
        for name in joint:
            self.assertEqual(lower_changed[name]['new'], [2,3])
        self.assertEqual(total['exact'], total['exact_lower'] + total['exact_upper'])
        self.assertEqual((total['exact_lower'],total['exact_upper'],total['joint_exact']), (2715,19,10))
        self.assertEqual(total['lower_improved']-total['joint_exact'],total['improved'])
        counts = c.computation_counts(self.new,self.table)
        self.assertEqual((counts['nSweepExact'],counts['nSweepImproved']), (1498,215))

    def test_september_ninth_comparison_includes_13n3733(self):
        from compare_knotinfo import load_snapshot
        snapshot = load_snapshot()
        self.assertEqual(snapshot['13n_3733'],[2,3])
        for row in self.new['extension_rows']:
            lo,hi=snapshot[row['name']]
            self.assertLess(lo,hi)
        self.assertEqual(self.new['counts']['latest_unresolved_exact'],2540)
        self.assertEqual(self.new['counts']['latest_range_improvements'],323)

    def test_presentation_counts_and_partner_exception(self):
        diagrams = [(n,pd,marks,'unknot') for n,pd,marks in load_paper_diagrams()]
        diagrams += load_additional_diagrams()
        self.assertEqual(Counter(len(pd) for _,pd,_,_ in diagrams),{13:15,14:3,18:1})
        self.assertEqual([(n,len(m),end) for n,_,m,end in diagrams if end!='unknot'], [('13n_447',1,'10_91')])
        self.assertTrue(all(len(m)==2 for _,_,m,end in diagrams if end=='unknot'))
        for name, _, _, _ in diagrams:
            self.assertEqual(self.table[name],[2,2])

    def test_reference_neighbours_use_verified_exact_values(self):
        evidence=json.loads((CODE/'results/extensions_2026-09-10/reverification.json').read_text())
        for name in ('13n_2251','13n_4025'):
            rows=evidence['reference_neighbors'][name]
            self.assertEqual({r['crossing_zero_based'] for r in rows},set(range(13)))
            for row in rows:
                self.assertTrue(row['identification']['extends_to_link'])
                self.assertEqual(self.table[row['target']],[2,2])


if __name__ == '__main__':
    unittest.main()
