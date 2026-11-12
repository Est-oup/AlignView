"""Tree geometry, Newick edge cases and integration with reports and projects."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_core import Seq
from av_tree import TREE_FIGURE, parse_newick, tree_svg
from av_report_v6 import build_report_v6
from alignview import Session, parse_args
from av_project import save_project, resume_project
from av_ui import UI

NS = {'s': 'http://www.w3.org/2000/svg'}
NEWICK = "((alpha:0.1,beta:0.2)0.97:0.3,gamma:0.4,delta:0.5);"


class TreeTests(unittest.TestCase):
    def test_lengths_supports_and_quoted_labels(self):
        tree = parse_newick("[&U] (('alpha,one':1e-2,'O''Brien':0)0.97:0.3,gamma[comment]:0.4);")
        clade = tree.children[0]
        self.assertEqual(clade.label, '0.97')
        self.assertEqual(clade.length, 0.3)
        self.assertEqual(clade.children[0].label, 'alpha,one')
        self.assertEqual(clade.children[0].length, 0.01)
        self.assertEqual(clade.children[1].label, "O'Brien")
        self.assertEqual(clade.children[1].length, 0)

    def test_invalid_trees_do_not_produce_figures(self):
        for source in ('', '(a,b)', '(a,b;', '(a,,b);', '(a,b));', '(a,b);(c,d);',
                       "('unclosed,b);", '(a[unclosed,b);', '(a:-0.1,b:1);',
                       '(a:NaN,b:1);', '(a:inf,b:1);', '(a:bad,b:1);', '(a:,b);',
                       '(a:1e308,(b:1e308,c:1):1e308);', 'a;'):
            with self.subTest(source=source), self.assertRaises(ValueError):
                tree_svg(source)

    def test_svg_preserves_distances_tip_names_and_support(self):
        svg = ET.fromstring(tree_svg(NEWICK))
        self.assertEqual(svg.attrib['data-layout'], 'phylogram')
        tips = {node.text: node for node in svg.findall('.//s:text[@class="tree-tip"]', NS)}
        self.assertEqual(set(tips), {'alpha', 'beta', 'gamma', 'delta'})
        # Root-to-tip distances: alpha=0.4, beta=0.5, gamma=0.4, delta=0.5.
        self.assertEqual(tips['alpha'].attrib['x'], tips['gamma'].attrib['x'])
        self.assertEqual(tips['beta'].attrib['x'], tips['delta'].attrib['x'])
        delta = float(tips['beta'].attrib['x']) - float(tips['alpha'].attrib['x'])
        self.assertAlmostEqual(delta, 720 * 0.1 / 0.5)
        self.assertEqual(svg.find('.//s:text[@class="tree-support"]', NS).text, '0.97')
        self.assertIsNotNone(svg.find('.//s:path[@class="tree-scale"]', NS))
        self.assertIn('without biological rooting', ''.join(svg.itertext()))

    def test_missing_and_zero_lengths_have_no_distance_scale(self):
        for source in ('((a,b),c);', '(a:0,b:0,c:0);', '(a:0.1,b,c:0.3);'):
            with self.subTest(source=source):
                svg = ET.fromstring(tree_svg(source))
                self.assertEqual(svg.attrib['data-layout'], 'cladogram')
                self.assertIsNone(svg.find('.//s:path[@class="tree-scale"]', NS))
                self.assertIn('distances are not represented', ''.join(svg.itertext()))

    def test_labels_are_xml_text_and_never_executable_markup(self):
        name = '<script>alert("x")</script>&_' + 'long_id_' * 20
        svg = ET.fromstring(tree_svg("('" + name + "':0.1,b:0.2);"))
        self.assertEqual(svg.find('.//s:text[@class="tree-tip"]', NS).text, name)
        self.assertEqual(svg.findall('.//s:script', NS), [])

    def test_large_unbalanced_tree_does_not_hit_recursion_limit(self):
        source = 'tip_0:0.1'
        for i in range(1, 1200):
            source = f'({source},tip_{i}:0.1):0.1'
        svg = ET.fromstring(tree_svg(source + ';'))
        self.assertEqual(svg.attrib['data-leaves'], '1200')
        self.assertEqual(len(svg.findall('.//s:text[@class="tree-tip"]', NS)), 1200)


class TreeIntegrationTests(unittest.TestCase):
    def test_report_catalog_and_stale_figure_removal(self):
        records = [Seq('alpha', 'alpha', 'ACGT'), Seq('beta', 'beta', 'ACGA'), Seq('gamma', 'gamma', 'ACGG')]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            tree = out / 'tree_fasttree.nwk'
            tree.write_text('(alpha:0.1,beta:0.2,gamma:0.3);')

            def report():
                return build_report_v6(out, records[0], records[1:], [], 'msa', {}, [], records).read_text()

            page = report()
            figure = out / 'figures' / (TREE_FIGURE + '.svg')
            self.assertTrue(figure.is_file())
            self.assertIn('data-group="Phylogeny" data-figure="' + TREE_FIGURE + '"', page)
            self.assertIn(figure.read_text(), page)
            self.assertIn(TREE_FIGURE + '.svg', (out / 'figures_catalog.tsv').read_text())
            self.assertIn('Download SVG', page)
            self.assertEqual(page.count('class="panel figcard"'), len(list((out / 'figures').glob('*.svg'))))
            tree.write_text('(alpha,,gamma);')
            page = report()
            self.assertIn('Figure d’arbre unavailable', page)
            self.assertFalse(figure.exists())
            tree.unlink()
            page = report()
            self.assertNotIn('data-figure="' + TREE_FIGURE + '"', page)
            self.assertNotIn(TREE_FIGURE, (out / 'figures_catalog.tsv').read_text())
            self.assertNotIn('Figure d’arbre unavailable', page)

    def test_failed_newick_preserves_previous_tree_and_svg(self):
        session = Session(parse_args([str(ROOT / 'examples/example_msa.fasta'), '--msa', 'view', '--no-menu']), UI(color=False))
        try:
            session.load_msa()
            tree = session.work / 'tree_fasttree.nwk'
            tree.write_text(NEWICK)
            (session.work / 'figures').mkdir()
            figure = session.work / 'figures' / (TREE_FIGURE + '.svg')
            figure.write_text('previous SVG')

            def invalid(alignment, output, *args, **kwargs):
                output.write_text('(a,,b);')

            with patch.object(session.tools, 'fasttree', side_effect=invalid), self.assertRaises(ValueError):
                session.make_tree()
            self.assertEqual(tree.read_text(), NEWICK)
            self.assertEqual(figure.read_text(), 'previous SVG')
        finally:
            session.close()

    def test_cli_exports_project_resume_and_alignment_invalidation(self):
        env = os.environ.copy()
        env['PATH'] = str(ROOT / 'tests/mockextra') + os.pathsep + env['PATH']
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'tree_project'
            result = subprocess.run([sys.executable, str(ROOT / 'alignview.py'),
                str(ROOT / 'examples/example_msa.fasta'), '--msa', 'view', '--tree', '--report',
                '--export-profile', 'publication', '--no-menu', '--module-policy', 'never',
                '--save', '--out', str(out)], env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            figure = Path('figures') / (TREE_FIGURE + '.svg')
            self.assertTrue((out / figure).is_file())
            self.assertEqual((out / figure).read_bytes(), (out / 'exports' / figure).read_bytes())
            self.assertIn(TREE_FIGURE, (out / 'exports' / 'report.html').read_text())
            session = resume_project(out, parse_args(['--resume', str(out), '--no-menu']), UI(color=False))
            try:
                # Old projects with Newick only gain a figure without rerunning FastTree.
                (session.work / figure).unlink()
                history = list(session.tools.history)
                page = session.make_report(open_preview=False).read_text()
                self.assertIn(TREE_FIGURE, page)
                self.assertEqual(session.tools.history, history)
                saved = save_project(session, Path(tmp) / 'saved_again')
                self.assertTrue((saved / figure).is_file())
                session.invalidate_derived()
                self.assertFalse((session.work / figure).exists())
                self.assertFalse((session.work / 'tree_fasttree.nwk').exists())
                page = session.make_report(open_preview=False).read_text()
                self.assertNotIn('data-figure="' + TREE_FIGURE + '"', page)
            finally:
                session.close()


if __name__ == '__main__':
    unittest.main()
