"""Loading a FASTA explains the file; readers and engines require a choice."""
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from av_cli import parse_args
from av_core import Seq,write_fasta
from av_session import Session
from av_workspace import Workspace
from av_ui import UI,menu_item
from test_release import Terminal

ROOT=Path(__file__).resolve().parents[1]


class FileMenuTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory(prefix='alignview_file_menu_')
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.source=self.folder/'my_sequences.fna'
        write_fasta([Seq('one','one','ACGT'*20),Seq('two','two','ACGA'*20),Seq('three','three','ACGG'*20)],self.source)

    def session(self):
        session=Session(parse_args([str(self.source)]),UI(False));self.addCleanup(session.close)
        return session

    def terminal(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.source,),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        self.addCleanup(terminal.close)
        return terminal

    def test_initial_frame_identifies_the_file_and_explains_all_five_actions(self):
        session=self.session();workspace=Workspace(session)
        self.assertEqual(workspace.tab,'Home')
        with patch.object(session.tools,'run') as engine:
            frame='\n'.join(workspace.render(80,24))
        for label in ('my_sequences.fna','3 DNA/RNA sequences','No analysis has been run',
                      'Look at the sequences','Align sequences','Find matching sequences',
                      'Other analyses','Files and previous results'):
            self.assertIn(label,frame)
        self.assertNotIn('Query:',frame);self.assertNotIn('Coord. units',frame)
        self.assertIsNone(session.msa_records);engine.assert_not_called()

    def test_narrow_menu_still_identifies_the_data_and_fits(self):
        workspace=Workspace(self.session())
        for width,height in ((40,14),(60,18),(80,24)):
            frame=workspace.render(width,height)
            self.assertLessEqual(len(frame),height-1)
            self.assertTrue(any('DNA/RNA' in line for line in frame))
            self.assertTrue(any('my_sequences.fna' in line for line in frame))

    def test_disabled_alignment_explains_why_without_hiding_other_actions(self):
        write_fasta([Seq('one','one','ACGT'*20)],self.source)
        workspace=Workspace(self.session())
        _,choices,_,_=workspace.landing()
        align=next(choice for choice in choices if choice.action=='align')
        self.assertFalse(align.enabled)
        self.assertIn('two sequences',align.reason)
        self.assertTrue(next(choice for choice in choices if choice.action=='compare').enabled)

    def test_actual_load_does_not_open_a_viewer_calculate_or_ask_to_save(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.expect('Files and previous results')
        self.assertNotIn(b'AlignView / Alignment',terminal.output)
        self.assertNotIn(b'Create an alignment',terminal.output)
        self.assertNotIn(b'Project directory',terminal.output)
        terminal.send('q');terminal.expect('Session closed.')
        self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_canceling_preparation_returns_to_same_file_menu_without_an_engine(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.send('2\r')
        terminal.expect('Create an alignment');terminal.send('q')
        terminal.expect('File loaded — choose an action');terminal.expect('No analysis has been run')
        terminal.send('q');terminal.expect('Session closed.')
        self.assertNotIn(b'Alignment ready',terminal.output)
        self.assertNotIn(b'Project directory',terminal.output)
        self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_browsing_sequences_is_explicit_and_escape_returns_to_the_file(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.send('1\r')
        terminal.expect('Look at the loaded sequences');terminal.send('1\r')
        terminal.expect('/ Sequences');terminal.send('\x1b')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Session closed.');self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_other_analyses_are_grouped_and_do_not_run_on_group_selection(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.send('4\r')
        terminal.expect('Other analyses');terminal.send('2\r')
        terminal.expect('Translation and ORFs');terminal.send('q')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Session closed.')
        self.assertNotIn(b'SIX-FRAME TRANSLATION',terminal.output)
        self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_selecting_a_single_record_is_explicit_before_a_sequence_analysis(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.send('4\r')
        terminal.expect('Other analyses');terminal.send('2\r')
        terminal.expect('Translation and ORFs');terminal.send('1\r')
        terminal.expect('Choose the sequence to analyse');terminal.send('q')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Session closed.')
        self.assertNotIn(b'SIX-FRAME TRANSLATION',terminal.output)
        self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_sequence_checks_also_wait_for_an_explicit_result_open(self):
        terminal=self.terminal()
        terminal.expect('File loaded — choose an action');terminal.send('4\r')
        terminal.expect('Other analyses');terminal.send('1\r')
        terminal.expect('Sequence checks, motifs and regions');terminal.send('1\r')
        terminal.expect('Sequence checks ready');terminal.expect('3 sequences inspected')
        self.assertNotIn(b'FASTA QUALITY CHECK',terminal.output)
        terminal.send('1\r');terminal.expect('ungapped')
        terminal.send('q');terminal.expect('Sequence checks ready');terminal.send('3\r')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Close workspace?');terminal.send('2\r')
        terminal.expect('without saving');self.assertEqual(terminal.process.wait(timeout=5),0)


if __name__=='__main__':unittest.main()
