"""Selections must preserve their biological rows and original MSA coordinates."""
from __future__ import annotations

import json
from contextlib import redirect_stderr, redirect_stdout
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from alignview import Session, parse_args, main
from av_core import Seq, read_fasta, sha256
from av_msa_tui import MSAViewer
from av_ui import UI
from test_release import Terminal

ROOT = Path(__file__).resolve().parents[1]


class ViewerScopeTests(unittest.TestCase):
    def viewer(self):
        return MSAViewer([
            Seq('alpha', 'alpha first', 'A-CGTAC'),
            Seq('beta', 'beta retained', 'ATCGTAT'),
            Seq('gamma', 'gamma retained', 'G-CGTAT'),
        ], color=False)

    def test_filtered_checked_region_exports_original_noncontiguous_columns(self):
        viewer = self.viewer()
        viewer.only_variable = True
        viewer.row_filter = 'retained'
        viewer.checked = {'alpha', 'gamma'}
        viewer.region = (1, 6)
        records, columns = viewer.export_selection()
        self.assertEqual(columns, [7])
        self.assertEqual([(r.id, r.seq) for r in records], [('gamma', 'T')])
        # A view must leave the original source sequences intact.
        self.assertEqual([r.seq for r in viewer.seqs], ['A-CGTAC', 'ATCGTAT', 'G-CGTAT'])

    def test_trace_distinguishes_gaps_and_native_residue_positions(self):
        viewer = self.viewer()
        viewer._handle_key('r')
        self.assertEqual(viewer.trace_id, 'alpha')
        viewer.col = 1
        self.assertTrue(any('Trace: alpha · residue gap' in line for line in viewer._frame()))
        viewer.col = 2
        self.assertTrue(any('Trace: alpha · residue 2' in line for line in viewer._frame()))
        viewer._handle_key('down')
        self.assertEqual(viewer.trace_id, 'alpha')
        viewer._handle_key('r')
        self.assertEqual(viewer.trace_id, 'beta')

    def test_no_matching_rows_cannot_select_a_hidden_sequence(self):
        viewer = self.viewer()
        viewer.row_filter = 'no such record'
        viewer._handle_key(' ')
        viewer._handle_key('r')
        self.assertEqual(viewer.checked, set())
        self.assertIsNone(viewer.trace_id)
        self.assertEqual(viewer.export_selection()[0], [])

    def test_saved_selection_records_source_digest_and_column_map(self):
        session = Session(parse_args([str(ROOT/'examples/example_msa.fasta'), '--msa', 'view', '--no-menu']), UI(False))
        self.addCleanup(session.close)
        session.load_msa()
        viewer = MSAViewer(session.msa_records, False)
        viewer.checked = {viewer.seqs[1].id}
        viewer.region = (1, 4)
        with tempfile.TemporaryDirectory() as folder:
            destination = Path(folder)/'exported scope'
            with patch.object(session.ui, 'ask', return_value=str(destination)):
                session.export_alignment_selection(viewer)
            self.assertEqual([r.id for r in read_fasta(destination/'alignment.fasta')], [viewer.seqs[1].id])
            self.assertEqual((destination/'columns.tsv').read_text(), 'output_column\tmsa_column\n1\t2\n2\t3\n3\t4\n4\t5\n')
            metadata = json.loads((destination/'selection.json').read_text())
            self.assertEqual(metadata['source_alignment_sha256'], sha256(session.msa_file))
            before = (destination/'alignment.fasta').read_bytes()
            with patch.object(session.ui, 'ask', return_value=str(destination)):
                with self.assertRaises(FileExistsError):
                    session.export_alignment_selection(viewer)
            self.assertEqual((destination/'alignment.fasta').read_bytes(), before)

    def test_ctrl_s_still_works_after_a_line_input_prompt(self):
        code = "import sys;sys.path.insert(0,sys.argv[1]);from av_core import Seq;from av_msa_tui import MSAViewer\nv=MSAViewer([Seq('a','a','ACGT'),Seq('b','b','ACGA')],False);v.run();print('RESULT='+str(v.result),flush=True)"
        terminal = Terminal(code)
        try:
            terminal.expect('AlignView / Alignment')
            terminal.send('g')
            terminal.expect('Go to MSA column')
            terminal.send('2\n')
            terminal.expect('MSA column 2')
            terminal.send('\x13')
            terminal.expect('RESULT=save')
            self.assertEqual(terminal.process.wait(timeout=3), 0)
        finally:
            terminal.close()

    def test_explicit_batch_tasks_reject_inputs_that_cannot_complete_the_task(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            unequal = root/'unequal.fasta'
            unequal.write_text('>a\nACGT\n>b\nACG\n')
            single = root/'single.fasta'
            single.write_text('>a\nACGT\n')
            for task, source in [('view', unequal), ('report', unequal), ('align', single)]:
                with self.subTest(task=task):
                    destination = root/task
                    html = root/(task+'.html')
                    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()), patch('av_tools.Toolchain.run') as run:
                        status = main([task, str(source), '--no-menu', '--save', '--out', str(destination), '--html-out', str(html)])
                    self.assertNotEqual(status, 0)
                    self.assertFalse(destination.exists())
                    self.assertFalse(html.exists())
                    run.assert_not_called()


if __name__ == '__main__':
    unittest.main()
