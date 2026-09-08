"""Independent small controls for the lower-bound methods in Section 3.

Run from code/: python -m unittest discover -s tests -p test_lower_bounds.py -v
The tests use exact presentations, known knots and explicit rational expected
values. Deposited result files are neither inputs nor overwritten outputs.
"""
import itertools
import math
import unittest
from fractions import Fraction
from unittest.mock import patch

import numpy as np
from sympy import Matrix, diag

from lower_bounds.cyclic_cover_bound import cover_bound, scan_rows
from lower_bounds.generator_bound_sweep import min_generators
from lower_bounds.linking_pairing import cyclic_generator, generator_self_linking
from lower_bounds.scan_torsion_lower_bounds import latest_records, summarise
from lower_bounds.montesinos.hf_red import _ExactEllipsoid, hf_red_all, lattice_points
from lower_bounds.montesinos.montesinos_u1 import (
    ClassMap, correction_terms, cyclic_labelling, hirzebruch_jung, lens_d,
    mapping_cone_ranks_ok, plumbing_for, star_plumbing, test_knot,
    u1_admissible, validate_plumbing,
)


class AbelianGroupTests(unittest.TestCase):
    def test_free_and_coprime_summands(self):
        self.assertEqual(min_generators([0, 0]), 2)
        self.assertEqual(min_generators([0, 2, 3]), 2)
        self.assertEqual(min_generators([2, 3]), 1)
        self.assertEqual(min_generators([3, 3, 9]), 3)
        self.assertEqual(min_generators([1, -1]), 0)

    def test_cover_free_rank_and_exact_ceiling(self):
        self.assertEqual(cover_bound([(3, [0, 0, 0])])[:2], (2, 3))
        self.assertEqual(cover_bound([(4, [0] * 7)])[:2], (3, 4))
        self.assertEqual(cover_bound([(4, [3] * 3), (2, [7])])[:2], (1, 2))
        with self.assertRaises(ValueError):
            cover_bound([(1, [3])])

    def test_cyclic_scan_is_pure(self):
        rows = [{'name': 'test', 'crossing_number': '3', 'unknotting_number': '[1,2]',
                 'torsion_numbers': '[[3, [0,0,0]]]'}]
        with patch('builtins.open', side_effect=AssertionError('no file writes')):
            result = scan_rows(rows)
        self.assertEqual(result['test']['bound'], 2)

    def test_generator_needing_four_primary_components(self):
        # In this presentation every vector supported on <=3 coordinates
        # misses a primary component. The former short-vector search failed.
        matrix = diag(3, 5, 7, 11)
        order = int(matrix.det())
        vector = cyclic_generator(matrix)
        denominators = [int(t.q) for t in matrix.inv() * vector]
        self.assertEqual(math.lcm(*denominators), order)
        self.assertEqual(math.gcd(generator_self_linking(matrix, order), order), 1)

    def test_noncyclic_primary_factors_are_rejected(self):
        matrix = diag(3, 3)
        self.assertIsNone(cyclic_generator(matrix))
        labels, pairing = cyclic_labelling(ClassMap(np.array(matrix).astype(int)), 9, -matrix)
        self.assertIsNone(labels)
        self.assertIsNone(pairing)

    def test_trivial_group_and_invalid_order(self):
        self.assertEqual(generator_self_linking(Matrix([[1]]), 1), 0)
        with self.assertRaises(ValueError):
            cyclic_generator(Matrix([[3]]), 5)


class ExactPlumbingTests(unittest.TestCase):
    def test_rational_ellipsoid_boundary(self):
        self.assertEqual(set(lattice_points([[2]], [0], 2)), {(-1,), (0,), (1,)})

    def test_ill_conditioned_unimodular_shear(self):
        # chi is y1(y1-1)/2 + y2(y2-1)/2 after y=(x1+m*x2,x2).
        # The exact zero sublevel has four points, however large m is.
        m = 10 ** 7
        matrix = [[1, m], [m, m * m + 1]]
        expected = {(a - m * b, b) for a in (0, 1) for b in (0, 1)}
        self.assertEqual(set(lattice_points(matrix, [1, m + 1], 0)), expected)

    def test_exact_enumeration_against_box(self):
        matrix, k = [[3, -1], [-1, 2]], [1, 0]
        expected = {x for x in itertools.product(range(-6, 7), repeat=2)
                    if 3*x[0]**2 - 2*x[0]*x[1] + 2*x[1]**2 - x[0] <= 12}
        self.assertEqual(set(lattice_points(matrix, k, 12)), expected)
        with self.assertRaises(ValueError):
            lattice_points(matrix, [0, 0], 12)

    def test_plumbing_hypotheses(self):
        with self.assertRaises(ValueError):
            validate_plumbing([[-3, 0], [0, -3]])
        with self.assertRaises(ValueError):
            correction_terms(np.array([[-1, 1], [1, -1]]))
        self.assertEqual(hirzebruch_jung(7, 2), [4, 2])
        with self.assertRaises(ValueError):
            hirzebruch_jung(2, 3)

    def test_lens_values_and_twist_knot_admissibility(self):
        self.assertEqual([lens_d(3, 2, i) for i in range(3)],
                         [Fraction(1, 6), Fraction(1, 6), Fraction(-1, 2)])
        for p in (3, 5, 7, 9, 11, 15):
            with self.subTest(p=p):
                qneg, _ = star_plumbing([Fraction(p, 2)])
                values, classes = correction_terms(qneg)
                expected = sorted(lens_d(p, 2, i) for i in range(p))
                self.assertTrue(sorted(values.values()) == expected or
                                sorted(-v for v in values.values()) == expected)
                self.assertTrue(u1_admissible(values, classes, p, qneg))
                ranks, _ = hf_red_all(qneg)
                self.assertEqual(set(ranks.values()), {0})

    def test_mapping_cone_adjacent_support(self):
        degree = 5
        labels = {i: (i,) for i in range(degree)}
        c = [0, 2, 1, 0, 3]
        ranks = {(i,): c[i] + c[(i - 1) % degree] for i in range(degree)}
        self.assertTrue(mapping_cone_ranks_ok(ranks, labels, degree, 1, 0, [0]))
        isolated = {(i,): int(i == 0) for i in range(degree)}
        self.assertFalse(mapping_cone_ranks_ok(isolated, labels, degree, 1, 0, [0]))
        # An unobserved nonzero V-tail supplies no rank obstruction.
        self.assertTrue(mapping_cone_ranks_ok(isolated, labels, degree, 1, 0, [1]))
        with self.assertRaises(ValueError):
            mapping_cone_ranks_ok(ranks, labels, degree, 1, 0, [2, 0])

    def test_lidman_known_obstruction(self):
        # Published independent control: Lidman, arXiv:2606.12431.
        result = test_knot('K(2/3;-1/3;-2/7)', 3)
        self.assertEqual(result['verdict_d'], 'PASS')
        self.assertEqual(result['hf_red_ranks'], [0, 0, 2])
        self.assertEqual(result['verdict'], 'OBSTRUCTED_HF')

    def test_12n457_exact_ranks_and_correction_terms(self):
        notation = 'K(2/3;-1/4;-2/7)'
        qneg, _ = plumbing_for(notation, 11)
        values, _ = correction_terms(qneg)
        expected = [11, 7, -5, 19, -9, -1, -1, -9, 19, -5, 7]
        self.assertEqual(sorted(values.values()), sorted(Fraction(n, 22) for n in expected))
        result = test_knot(notation, 11)
        self.assertEqual(result['verdict_d'], 'PASS')
        self.assertEqual(result['hf_red_ranks'], [0] * 9 + [1, 1])
        self.assertEqual(result['verdict'], 'OBSTRUCTED_HF')

    def test_12n457_all_affine_matches_without_pairing_filter(self):
        # Independent weaker test: enumerate *every* affine map, without the
        # production linking-pairing filter or mapping_cone_ranks_ok routine.
        qneg, _ = plumbing_for('K(2/3;-1/4;-2/7)', 11)
        values, classes = correction_terms(qneg)
        ranks, _ = hf_red_all(qneg)
        labels, _ = cyclic_labelling(classes, 11, qneg)
        survivors = []
        for orientation in (-1, 1):
            for a in range(1, 11):
                for b in range(11):
                    sequence = {}
                    for i in range(11):
                        difference = lens_d(11, 2, i) - orientation * values[labels[(a*i+b) % 11]]
                        if difference < 0 or difference.denominator != 1 or difference.numerator % 2:
                            break
                        j, value = min(i//2, (12-i)//2), int(difference/2)
                        if sequence.setdefault(j, value) != value:
                            break
                    else:
                        vector = [sequence[j] for j in sorted(sequence)]
                        if any(vector[j]-vector[j+1] not in (0, 1) for j in range(len(vector)-1)):
                            continue
                        self.assertEqual(vector, [1, 0, 0, 0])
                        support = [i for i in range(11) if ranks[labels[(a*i+b) % 11]]]
                        self.assertEqual(support, [4, 8])
                        # T_i=0 for this V sequence. Any positive c_i would
                        # force two adjacent positive ranks, which are absent.
                        self.assertNotIn((support[1]-support[0]) % 11, (1, 10))
                        survivors.append((orientation, a, b))
        self.assertEqual(len(survivors), 2)


class KnotControls(unittest.TestCase):
    def test_known_u1_linking_and_branched_cover_homology(self):
        from lower_bounds.tait_tools import (link_from_knotinfo, lickorish_u1_allowed,
                                               sigma2_homology)
        from lower_bounds.regina_linking_check import worker
        for name, determinant in [('3_1', 3), ('4_1', 5), ('5_2', 7), ('6_1', 9)]:
            with self.subTest(name=name):
                link = link_from_knotinfo(name)
                self.assertTrue(lickorish_u1_allowed(link))
                self.assertEqual(sigma2_homology(link).order(), determinant)
                record = worker(name)
                self.assertNotIn('error', record)
                self.assertTrue(record['regina_allowed'])

    def test_known_linking_obstruction_controls(self):
        from lower_bounds.tait_tools import knotinfo_rows, link_from_knotinfo, lickorish_u1_allowed
        from lower_bounds.seifert_linking_check import verdict
        from lower_bounds.regina_linking_check import worker
        # 7_4 has a cyclic pairing that fails both square classes; 8_18
        # instead fails cyclicity. Both have independently known u=2.
        for name in ('7_4', '8_18'):
            with self.subTest(name=name):
                row = knotinfo_rows()[name]
                self.assertEqual(lickorish_u1_allowed(link_from_knotinfo(name)), set())
                seifert = verdict((name, row['seifert_matrix'], int(row['determinant'])))
                self.assertTrue(seifert['obstructed'])
                regina = worker(name)
                self.assertNotIn('error', regina)
                self.assertEqual(regina['regina_allowed'], [])
        wrong = verdict(('wrong label', '[[1, 0], [1, 1]]', 5))
        self.assertIsNone(wrong['obstructed'])

    def test_hfk_independent_known_orders(self):
        from lower_bounds.scan_torsion_lower_bounds import one
        # Trefoil and figure-eight have order1; T(3,4)=8_19 has order2.
        for name, order in [('3_1', 1), ('4_1', 1), ('8_19', 2)]:
            with self.subTest(name=name):
                result = one((name, 2))
                self.assertNotIn('error', result)
                self.assertEqual(result['ord_u'], order)

    def test_failed_torsion_rows_are_retryable_and_not_bounds(self):
        records = [{'name': 'a', 'error': 'limit'}, {'name': 'a', 'ord_u': 2},
                   {'name': 'b', 'error': 'limit'}]
        current = latest_records(records)
        self.assertNotIn('error', current['a'])
        summary = summarise(records, 1)
        self.assertEqual(summary['rows'], 2)
        self.assertEqual(summary['failed_rows'], 1)
        self.assertEqual(summary['ord_distribution'], {2: 1})
        self.assertEqual(summary['new_lower_bounds'], [{'name': 'a', 'ord_u': 2}])


if __name__ == '__main__':
    unittest.main()
