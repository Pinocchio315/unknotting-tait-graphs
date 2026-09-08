"""Independent exact controls for Greene pages, L-space premises and surgery.

Small diagrams and deliberately ambiguous identifications are used; deposited
result files are neither required nor overwritten.
"""
import copy
from collections import Counter
from fractions import Fraction
from math import gcd
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import spherogram

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'lower_bounds' / 'greene'))
import pin_dinv
from greene_dinv import Diagram
from half_integral import lens_half_integral_d, surgery_test
from qa_search import determinant, qa, smoothings, simplify, verify_certificate

class SolitaryStateTests(unittest.TestCase):
    def test_solitary_equals_unique_overstrand_side_word(self):
        # Greene's Lemma 6.1 gives a characterization independent of cycle detection.
        for name in ('3_1', '4_1', '8_20'):
            pd = spherogram.Link(name).PD_code()
            for colour in (0, 1):
                for mark in sorted({label for crossing in pd for label in crossing}):
                    diagram = Diagram(pd, colour, mark)
                    states = diagram.states()
                    words = Counter(tuple(diagram.side(corner) for corner in state) for state in states)
                    for state in states:
                        row = diagram.analyse(state)
                        word = tuple(diagram.side(corner) for corner in state)
                        self.assertEqual(row['solitary'], words[word] == 1)
                        self.assertTrue(all((row['v'][i] - diagram.G[i, i]) % 2 == 0 for i in range(diagram.m)))

    def test_multiple_compatible_linking_isometries_are_retained(self):
        # Z/15 has four square roots of one. Two identifications modulo
        # conjugation survive; selecting either arbitrarily would be unsound.
        D = 15
        reference = {k: [0] for k in range(D)}
        for k in (1, 4, 11, 14): reference[k] = [0, 1]
        other = {k: [0, 1] for k in range(D)}
        other[0] = [0]
        other[1] = other[14] = [0]
        other[4] = other[11] = [1]
        pages = {(0, 1): (reference, Fraction(1, D)), (0, 2): (other, Fraction(1, D))}
        with patch.object(pin_dinv, 'pages_of', return_value=pages):
            values, _, ref = pin_dinv.pin('ambiguous', [], D, verbose=False)
        self.assertEqual(ref, (0, 1))
        self.assertEqual(values[1], {0, 1})
        self.assertEqual(values[4], {0, 1})

    def test_inconsistent_pages_fail_instead_of_fabricating_obstruction(self):
        pages = {(0, 1): ({0: [0], 1: [0], 2: [0]}, Fraction(1, 3)),
                 (0, 2): ({0: [1], 1: [1], 2: [1]}, Fraction(1, 3))}
        with patch.object(pin_dinv, 'pages_of', return_value=pages):
            with self.assertRaises(ValueError): pin_dinv.pin('inconsistent', [], 3, verbose=False)

class LSpacePremiseTests(unittest.TestCase):
    def test_rational_thinness_does_not_substitute_for_mod2_input(self):
        with self.assertRaises(ValueError):
            pin_dinv.khovanov_lspace_evidence({'determinant': 3, 'khovanov_reduced_rational_vector': '[[1,3,0,0]]'})

    def test_nonminimal_mod2_rank_is_rejected(self):
        with self.assertRaises(ValueError):
            pin_dinv.khovanov_lspace_evidence({'determinant': 3, 'khovanov_reduced_mod2_vector': '[[2,5,0,0]]'})
        evidence = pin_dinv.khovanov_lspace_evidence({'determinant': 3, 'khovanov_reduced_mod2_vector': '[[2,3,0,0]]'})
        self.assertEqual(evidence['rank'], 3)

class SurgeryTests(unittest.TestCase):
    def test_lens_spaces_survive_affine_relabelling(self):
        # Surgery on the unknot itself is a positive control regardless of
        # the affine cyclic labels or orientation assigned to the d-vector.
        D = 9
        lens = {i: lens_half_integral_d(D, i) for i in range(D)}
        for a in range(1, D):
            if gcd(a, D) != 1: continue
            for b in (0, 2):
                d = {(a * i + b) % D: -value for i, value in lens.items()}
                result = surgery_test(d, D)
                self.assertTrue(result['fits'])
                self.assertFalse(result['linking_filter'])
                self.assertEqual(sum(v['tested'] for v in result['orientations'].values()), 2 * D * 6)

    def test_invalid_order_never_returns_an_empty_obstruction(self):
        for order in (0, -1, -3, 1, 2, 4, 3.0, '3', None, True):
            with self.subTest(order=order):
                with self.assertRaises(ValueError):
                    surgery_test({}, order)

    def test_incomplete_and_inexact_inputs_are_rejected(self):
        with self.assertRaises(ValueError): surgery_test({0: Fraction(0)}, 3)
        with self.assertRaises(ValueError): surgery_test({0: 0.0, 1: 0, 2: 0}, 3)

class QuasiAlternatingTests(unittest.TestCase):
    def test_resolution_retains_unlinked_circles(self):
        # A curl resolves to a two-component unlink and a one-component unknot.
        children = smoothings([[1, 1, 2, 2]], 0)
        self.assertEqual(sorted(circles for _, circles in children), [1, 2])
        self.assertEqual(sorted(determinant(pd, circles) for pd, circles in children), [0, 1])
        self.assertEqual(determinant(spherogram.Link('3_1').PD_code(), 1), 0)

    def test_simplification_keeps_discarded_unlink_components(self):
        pd, circles = simplify([[1, 1, 2, 2], [3, 3, 4, 4]], 0)
        self.assertEqual(pd, [])
        self.assertEqual(circles, 2)
        self.assertEqual(determinant(pd, circles), 0)

    def test_certificate_replay_rejects_tampered_leaf(self):
        pd = spherogram.Link('3_1').PD_code()
        certificate = qa(pd)
        self.assertEqual(verify_certificate(certificate, pd, 0), 3)
        forged = copy.deepcopy(certificate); forged['determinant'] = 5
        with self.assertRaises(ValueError): verify_certificate(forged, pd, 0)

if __name__ == '__main__': unittest.main()
