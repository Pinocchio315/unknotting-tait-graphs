"""Reject manuscript inputs that would conceal stale or incomplete results.

These tests exercise the boundary between authenticated mathematical records
and LaTeX inclusion. They do not rerun an obstruction or certify its proof.
"""
from pathlib import Path
import tempfile
import unittest

from consolidate_results import appendix_tag, manuscript_tex, reconstruct
from verify_manuscript import load_manuscript, verify


class ManuscriptInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        consolidated, table = reconstruct()
        cls.generated = manuscript_tex(consolidated, table)

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.inputs = self.root / 'code' / 'generated' / 'paper_v1_2'
        self.inputs.mkdir(parents=True)
        for name, text in self.generated.items():
            (self.inputs / name).write_text(text)
        self.manuscript = self.root / 'main_v1.2.tex'

    def include(self, *names):
        self.manuscript.write_text('\n'.join(
            rf'\input{{code/generated/paper_v1_2/{name}}}' for name in names))

    def test_authentic_inputs_include_both_greene_knots(self):
        self.include('counts.tex', 'rows_u2.tex')
        source, included = load_manuscript(self.manuscript)
        self.assertEqual(set(included), {'counts.tex', 'rows_u2.tex'})
        self.assertIn(r'12n491$^{\mathrm{G}}$', source)
        self.assertIn(r'13n3370$^{\mathrm{G}}$', source)

    def test_stale_count_or_missing_greene_case_is_rejected(self):
        self.include('counts.tex', 'rows_u2.tex')
        changes = {
            'counts.tex': lambda text: text.replace('1{,}227', '1{,}225'),
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
        self.manuscript.write_text(r'\input{main_v1.2.tex}')
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
            verify(source)

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


if __name__ == '__main__':
    unittest.main()
