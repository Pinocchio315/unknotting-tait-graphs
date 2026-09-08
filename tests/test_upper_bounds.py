"""Deterministic mathematical regressions for the Section 4 code.

Run with: python -m unittest discover -s tests -p test_upper_bounds.py -v
Requires the supplied topology environment; no search campaign or raw-data writes.
"""
from pathlib import Path
import copy
from fractions import Fraction
import importlib.util
import itertools
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT, ROOT / 'tait_graphs', ROOT / 'upper_bounds' / 'cluster',
             ROOT / 'upper_bounds', ROOT / 'crossing_changes'):
    sys.path.insert(0, str(path))
import tait
from tait import knotinfo, invariants as inv
from tait.graph import to_link
from spherogram import Link
from sympy import Matrix
from xtait.graph import from_pd, to_pd, dual, pd_crossing_change, pd_passages
from xtait import moves as mv
from xtait.canonical import canonical_code
from xtait.flypes import flype_orbit, find_flypes_all, flype_block
from xtait.passmove import iter_graph_pass_candidates, apply_pass_params
from flipdet import FlipDet, P, exact_det
from lickorish import LickorishFilter
from sig import symmetric_signature, signature, goeritz_data
from candidates import candidate_status, connected_sum_interval
import mccoy_u1_check as mccoy
import verify_certificates as certificates
from certificate_checks import validate_marks, check_partner_bound
from verify_presentation_diagrams import PAPER_KNOTS, load_paper_diagrams
from verify_enum_certificates import verify_enum_upper_bound


def graph(name):
    return from_pd(knotinfo.pd_code(name))


class ExactFormTests(unittest.TestCase):
    def test_inertia_handles_tiny_eigenvalue_without_tolerance(self):
        n = 10 ** 20
        # Determinant -1: one positive and one negative eigenvalue. Floating
        # eigensolvers lose the small eigenvalue of this integral form.
        self.assertEqual(symmetric_signature([[n + 1, n], [n, n - 1]]), 0)
        self.assertEqual(symmetric_signature([[Fraction(1, 10**50), 0], [0, 2]]), 2)
        self.assertEqual(symmetric_signature([[0, 2, 0], [2, 0, 0], [0, 0, 0]]), 0)
        self.assertEqual(symmetric_signature([]), 0)
        with self.assertRaises(ValueError):
            symmetric_signature([[1, 2], [0, 1]])

    def test_signature_reference_values_and_mirror(self):
        for name in ('3_1', '4_1', '6_2', '7_4', '8_19', '8_20', '9_35', '10_113'):
            with self.subTest(name=name):
                pd = knotinfo.pd_code(name)
                expected = int(knotinfo.row(name)['signature'])
                self.assertEqual(signature(pd), expected)
                self.assertEqual(signature(Link(pd).mirror().PD_code()), -expected)

    def test_subset_determinants_against_independent_goeritz(self):
        for name in ('6_2', '8_20'):
            base = graph(name)
            # Include a loop and a pendant edge: zero incidence vectors and
            # changes of Goeritz rank must both obey the determinant lemma.
            variants = [base, mv.r1_add_loop(base, base.vertices[0], 0, -1),
                        mv.r1_add_pendant(base, base.vertices[0], 0, 1)]
            for g in variants:
                filt = FlipDet(g)
                for size in (1, 2, 3):
                    for rows in itertools.combinations(range(g.n_crossings()), size):
                        changed = g
                        for row in rows:
                            changed = mv.crossing_change(changed, filt.edge_ids[row])
                        determinant = abs(int(inv.goeritz_matrix(changed).det()))
                        self.assertIn(filt.det_after(rows), (determinant % P, -determinant % P))
                        self.assertEqual(exact_det(changed), determinant)

    def test_exact_linking_updates_and_necessary_condition(self):
        for name in ('3_1', '6_2', '8_20'):
            base = graph(name)
            filt = LickorishFilter(base)
            for rows in itertools.chain([()], itertools.combinations(range(base.n_crossings()), 1)):
                changed = base
                for row in rows:
                    changed = mv.crossing_change(changed, filt.edge_ids[row])
                # Build the same deleted-vertex presentation independently.
                vertices = sorted(changed.vertices, key=str)
                idx = {v: i for i, v in enumerate(vertices[:-1])}
                G = Matrix.zeros(len(idx))
                for e, (u, v) in changed.edges.items():
                    x = Matrix([int(w == u) - int(w == v) for w in vertices[:-1]])
                    G += changed.signs[e] * x * x.T
                self.assertEqual(Matrix(filt.inverse_after(rows)), G.inv())
                allowed = inv.lickorish_allowed(changed)
                self.assertEqual(filt.passes(rows, exact_det(changed)), bool(allowed))
        # Determinant one cannot be rejected by the linking pairing.
        self.assertTrue(LickorishFilter(graph('3_1')).passes((0,), 1))
        self.assertFalse(LickorishFilter(graph('7_4')).passes((), 15))

    def test_signature_and_electrical_crossing_formula(self):
        for name in ('3_1', '4_1', '6_2', '7_4', '8_3'):
            g = graph(name)
            if next(iter(g.signs.values())) < 0:
                g = dual(g)
            G = inv.goeritz_matrix(g)
            vertices = list(g.vertices)
            D, s0 = abs(int(G.det())), signature(Link(to_pd(g)).PD_code())
            for edge, (u, v) in g.edges.items():
                x = Matrix([int(w == u) - int(w == v) for w in vertices[:-1]])
                resistance = (x.T * G.inv() * x)[0]
                changed = mv.crossing_change(g, edge)
                self.assertEqual(exact_det(changed), D * abs(1 - 2 * resistance))
                # A single crossing change changes signature by at most two.
                # The preceding reference/mirror test fixes the GL convention.
                s1 = signature(Link(to_pd(changed)).PD_code())
                self.assertIn(s1 - s0, (-2, 0, 2))
                Gp = G - 2 * x * x.T
                self.assertEqual((x.T * Gp.inv() * x)[0], resistance / (1 - 2 * resistance))


class DiagramMoveTests(unittest.TestCase):
    def assert_same_invariants(self, original, changed):
        changed.validate()
        pd_passages(to_pd(changed))
        self.assertEqual(inv.jones_regina(to_pd(changed)), inv.jones_regina(to_pd(original)))
        self.assertEqual(exact_det(changed), exact_det(original))

    def test_empty_graph_converts_to_one_unknot(self):
        empty = from_pd([])
        link = to_link(empty)
        self.assertEqual(len(link.link_components) + link.unlinked_unknot_components, 1)
        self.assertTrue(inv.is_unknot(link))

    def test_reidemeister_additions_and_removals(self):
        g = graph('6_2')
        for add in (mv.r1_add_loop, mv.r1_add_pendant):
            for sign in (-1, 1):
                h = add(g, g.vertices[0], 0, sign)
                self.assert_same_invariants(g, h)
                red = mv.r1_remove(h, max(h.edges))
                self.assertEqual(canonical_code(g), canonical_code(red))
        for candidates, add in ((mv.find_r2_parallel_additions, mv.r2_add_parallel),
                                (mv.find_r2_series_additions, mv.r2_add_series)):
            for parameters in candidates(g)[:6]:
                for sign in (-1, 1):
                    h = add(g, *parameters, sign)
                    self.assert_same_invariants(g, h)
                    inverses = []
                    for kind, parameters2 in mv.find_r2_removals(h):
                        try:
                            inverses.append(mv.apply_move(h, 'r2-' + kind[0], parameters2))
                        except ValueError:
                            pass
                    self.assertTrue(any(canonical_code(r) == canonical_code(g) for r in inverses))

    def test_r3_pass_and_flype_invariants(self):
        base = graph('6_2')
        changed = mv.crossing_change(base, sorted(base.edges)[1])
        count_r3 = count_pass = count_flype = 0
        for g in (base, changed):
            for triangle in mv.find_r3_triangles(g):
                self.assert_same_invariants(g, mv.r3_delta_to_y(g, triangle))
                count_r3 += 1
            for star in mv.find_r3_stars(g):
                self.assert_same_invariants(g, mv.r3_y_to_delta(g, star))
                count_r3 += 1
            for params in itertools.islice(iter_graph_pass_candidates(g, max_extra=1, max_routes=4), 20):
                self.assert_same_invariants(g, apply_pass_params(g, params))
                count_pass += 1
            for h in (g, dual(g)):
                for edge, piece in itertools.islice(find_flypes_all(h), 12):
                    self.assert_same_invariants(g, flype_block(h, edge, piece))
                    count_flype += 1
        self.assertGreater(count_r3, 0)
        self.assertGreater(count_pass, 0)
        self.assertGreater(count_flype, 0)

    def test_flype_orbit_preserves_neighbor_polynomials(self):
        base = graph('8_3')
        def neighbors(g):
            return {tuple(sorted(inv.jones_regina(to_pd(mv.crossing_change(g, e))).items())) for e in g.edges}
        expected = neighbors(base)
        for g in flype_orbit(base).values():
            self.assertEqual(neighbors(g), expected)

    def test_branched_cover_tracks_meridian_not_homology_guess(self):
        for name in ('3_1', '4_1', '6_1', '9_35'):
            cover = inv.sigma2_snappy(Link(knotinfo.pd_code(name)))
            actual = sorted(int(x) for x in cover.homology().elementary_divisors())
            expected = sorted(inv.h1_torsion(graph(name)))
            self.assertEqual(actual, expected)
            self.assertTrue(all(tuple(map(int, slope)) != (0, 0) for slope in cover.cusp_info('filling')))


class CertificateTests(unittest.TestCase):
    def make_certificate(self):
        g = graph('3_1')
        pd = to_pd(g)
        return {'name': '3_1', 'start_pd': knotinfo.pd_code('3_1'),
                'expansion_path': [], 'witness_pd': pd, 'n_witness': len(pd),
                'kind': 'unknot', 'j': 1, 'change_rows': [0], 'range': [1, 1]}

    def test_valid_small_certificate_and_false_source(self):
        cert = self.make_certificate()
        self.assertEqual(certificates.verify_upper_bound(cert, {})[0], 1)
        cert['name'] = '4_1'
        with self.assertRaises(ValueError):
            certificates.verify_upper_bound(cert, {})

    def test_uncounted_crossing_change_is_not_isotopy(self):
        with self.assertRaises(ValueError):
            certificates.replay(knotinfo.pd_code('3_1'), [['cc', [0]]])

    def test_marked_count_and_bound_are_authenticated(self):
        for rows, n, count in (([0, 0], 3, 2), ([-1], 3, 1), ([3], 3, 1), ([True], 3, 1), ([0], 3, 0)):
            with self.assertRaises(ValueError):
                validate_marks(rows, n, count)
        with self.assertRaises(ValueError):
            check_partner_bound(1, [2, 3])
        self.assertEqual(check_partner_bound(4, [2, 3]), 3)
        cert = self.make_certificate()
        cert['j'] = 0
        with self.assertRaises(ValueError):
            certificates.verify_upper_bound(cert, {})

    def test_paper_manifest_requires_all_eight(self):
        diagrams = load_paper_diagrams()
        self.assertEqual({name for name, _, _ in diagrams}, set(PAPER_KNOTS))
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(FileNotFoundError):
                load_paper_diagrams(directory)

    def test_enumeration_does_not_sort_away_crossing_information(self):
        cert = self.make_certificate()
        cert['moves'] = []
        self.assertEqual(verify_enum_upper_bound(cert, {}), 1)
        cert['witness_pd'] = pd_crossing_change(cert['witness_pd'], 1)
        with self.assertRaises(ValueError):
            verify_enum_upper_bound(cert, {})


class CandidateTests(unittest.TestCase):
    def test_jones_hypotheses_do_not_become_witnesses(self):
        self.assertNotEqual(candidate_status('K', 'sum-jones', [1, 1], 0), 'WITNESS')
        self.assertEqual(candidate_status('K', 'sum-jones', [2, 3], 0), 'excluded:u>=2(jones)')
        self.assertEqual(candidate_status('K', 'unknown-method', [2, 3], 0), 'candidate:unverified')
        self.assertEqual(candidate_status(None, None, None, None), 'candidate:unidentified')
        self.assertEqual(candidate_status(None, None, None, 4), 'excluded:sigma')
        self.assertEqual(candidate_status('3_1', 'code', [1, 1], 2), 'WITNESS')

    def test_one_nontrivial_block_is_not_composite(self):
        known = {'3_1': [1, 1], '7_1': [3, 3]}
        self.assertEqual(connected_sum_interval(['0_1', '3_1'], known), [1, 1])
        self.assertEqual(connected_sum_interval(['3_1', '3_1'], known), [2, 2])
        # No additive lower bound is assumed for opposite-sign summands.
        self.assertEqual(connected_sum_interval(['7_1', '7_1'], known, 0), [2, 6])

    def test_mccoy_positive_negative_and_unknown(self):
        self.assertEqual(mccoy.check_diagram(knotinfo.pd_code('3_1'))['verdict'], 'u = 1')
        self.assertEqual(mccoy.check_diagram(knotinfo.pd_code('5_1'))['verdict'], 'u >= 2 (McCoy)')
        with patch.object(mccoy, 'greedy_reduce', side_effect=RuntimeError('injected failure')):
            result = mccoy.check_diagram(knotinfo.pd_code('5_1'))
        self.assertEqual(result['verdict'], 'undecided')
        self.assertEqual(len(result['undecided']), 5)
        with self.assertRaises(ValueError):
            mccoy.check_diagram(knotinfo.pd_code('8_19'))


if __name__ == '__main__':
    unittest.main()
