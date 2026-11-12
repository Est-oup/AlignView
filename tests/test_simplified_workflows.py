"""Regression checks for result reuse, explicit exports and Python 3.10 docs."""
from contextlib import redirect_stdout
import argparse
import csv
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,parse_args
from av_navigation import Workbench
from av_parameters import PARAMETERS
from av_tree import tree_text
from av_ui import UI
from scripts.generate_cli_docs import reference


class SimplifiedWorkflowsTests(unittest.TestCase):
    def setUp(self):
        self.session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--msa','view','--no-menu']),UI(False))
        self.addCleanup(self.session.close)

    def test_reading_and_exporting_preserve_clean_state_and_input(self):
        session=self.session;original=session.query_path.read_bytes()
        session.load_msa()
        self.assertFalse(session.dirty)
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'alignment.fa'
            session.export_artifact('msa',target)
            self.assertIn('>seq_1',target.read_text())
            with self.assertRaises(FileExistsError):session.export_artifact('msa',target)
            with patch.object(session.tools,'msa') as aligner,redirect_stdout(io.StringIO()):
                session.export_summary()
                aligner.assert_not_called()
        self.assertEqual(session.query_path.read_bytes(),original)
        self.assertFalse(session.dirty)
        self.assertTrue(Workbench(session).confirm_exit())
        self.assertEqual(session.exit_decision,'clean')

    def test_svg_export_uses_existing_tree_without_running_engine(self):
        session=self.session
        (session.work/'tree_fasttree.nwk').write_text('(seq_1:0.1,(seq_2:0.2,seq_3:0.3)0.95:0.1);')
        with tempfile.TemporaryDirectory() as tmp,patch.object(session,'make_tree') as calculate:
            target=session.export_artifact('tree',Path(tmp)/'tree.svg')
            self.assertIn('seq_3',target.read_text())
            self.assertIn('<svg',target.read_text())
            calculate.assert_not_called()
        self.assertFalse(session.dirty)

    def test_zero_matches_exports_a_valid_header(self):
        self.session.search_completed=True
        with tempfile.TemporaryDirectory() as tmp:
            target=self.session.export_artifact('hits',Path(tmp)/'matches.tsv')
            with target.open() as stream:
                reader=csv.DictReader(stream,delimiter='\t')
                self.assertIn('query_coverage',reader.fieldnames)
                self.assertEqual(list(reader),[])

    def test_tree_terminal_view_draws_topology_and_explains_scale(self):
        text=tree_text('(A:0.1,(B:0.2,C:0.3)0.95:0.1);')
        self.assertIn('├─ A: 0.1',text)
        self.assertIn('└─ C: 0.3',text)
        self.assertIn('not to scale',text)
        self.assertIn('Arbitrary drawing root',text)

    def test_cli_reference_ignores_argparse_version_specific_title(self):
        from alignview import build_parser
        expected=reference()
        def legacy_parser():
            parser=build_parser();parser._optionals.title='optional arguments'
            return parser
        with patch('scripts.generate_cli_docs.build_parser',side_effect=legacy_parser):
            self.assertEqual(reference(),expected)

    def test_boolean_help_uses_registry_text_on_every_python(self):
        from alignview import build_parser
        parser=build_parser()
        by_dest={action.dest:action for action in parser._actions}
        for parameter in PARAMETERS:
            if isinstance(by_dest[parameter.name],argparse.BooleanOptionalAction):
                self.assertEqual(by_dest[parameter.name].help,parameter.label+': '+parameter.effect)


if __name__=='__main__':unittest.main()
