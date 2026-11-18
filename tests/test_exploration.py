"""Free exploration: independent attempts, Back, and optional persistence."""
from contextlib import redirect_stdout, redirect_stderr
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from alignview import Session, parse_args
from av_app import _register_inputs
from av_core import Seq, write_fasta, read_fasta
from av_library import entries
from av_navigation import Workbench
from av_project import save_project, resume_project, read_project, local_path
from av_tools import CommandError
from av_ui import UI
from av_workspace import Workspace, VIEWS
from test_release import Terminal

ROOT=Path(__file__).resolve().parents[1]


class ExplorationTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory(prefix='alignview_exploration_')
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.source=self.folder/'sequences.fna'
        write_fasta([Seq('a','a','ACGT'*30),Seq('b','b','ACGA'*30),Seq('c','c','ACGG'*30)],self.source)

    def session(self,*args):
        session=Session(parse_args([*map(str,args or [self.source]),'--no-menu','--module-policy','never',
            '--tools-dir',str(ROOT/'tests/mockbin'),'--tools-dir',str(ROOT/'tests/mockextra')]),UI(False))
        self.addCleanup(session.close)
        session.library_root=self.folder/('library_'+str(len(list(self.folder.glob('library_*')))))
        with redirect_stdout(io.StringIO()):_register_inputs(session)
        return session

    def calculate(self,session,method,*args):
        with redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()),patch.object(session,'show_msa_window'):
            return getattr(session,method)(*args)

    def test_two_alignments_and_tree_are_independent_after_an_unrelated_analysis(self):
        session=self.session()
        self.calculate(session,'msa_run','all')
        first=entries(session,'result')[0]
        root,state=read_project(first['path']);original=local_path(root,state['msa']).read_bytes()
        self.calculate(session,'make_tree')
        tree=entries(session,'result')[-1]
        self.assertEqual(tree['inputs'],[first['id']])
        self.calculate(session,'translate_frames')
        session.a.aligner='muscle'
        self.calculate(session,'msa_run','all')
        second=entries(session,'result')[-1]
        self.assertNotEqual(first['id'],second['id'])
        self.assertEqual(local_path(root,state['msa']).read_bytes(),original)
        self.assertTrue((tree['path']/'tree_fasttree.nwk').is_file())
        self.assertEqual([row['analysis'] for row in entries(session,'result')],['msa','tree','translation','msa'])
        self.assertFalse(any(self.folder.glob('saved*')))

    def test_selected_records_are_frozen_in_each_alignment(self):
        session=self.session();session.sequence_selection=['b','c']
        self.calculate(session,'msa_run','selection')
        first=entries(session,'result')[0]
        root,data=read_project(first['path'])
        self.assertEqual([record.id for record in read_fasta(local_path(root,data['msa']))],['b','c'])
        session.sequence_selection=['a','b'];self.calculate(session,'msa_run','selection')
        self.assertEqual(data['calculation_runs']['msa']['preparation']['record_ids'],['b','c'])
        self.assertEqual([record.id for record in read_fasta(local_path(root,data['msa']))],['b','c'])

    def test_failure_and_cancel_keep_successful_attempts(self):
        session=self.session();self.calculate(session,'msa_run','all')
        previous=list(session.history);original=session.msa_file.read_bytes()
        with patch.object(session.tools,'msa',side_effect=CommandError('failed')):
            with self.assertRaises(CommandError):self.calculate(session,'msa_run','all')
        self.assertEqual(session.history,previous)
        self.assertEqual(session.msa_file.read_bytes(),original)
        bench=Workbench(session)
        with patch.object(bench,'choose',return_value=None),patch.object(session.tools,'msa') as engine:
            bench.align();engine.assert_not_called()
        self.assertEqual(session.history,previous)

    def test_filtering_search_never_removes_the_stored_alignment_and_tree(self):
        refs=self.folder/'refs.fna';write_fasta([Seq('r1','r1','ACGT'*30),Seq('r2','r2','ACGA'*30)],refs)
        session=self.session(self.source,refs,'--query-id','a','--mode','blastn')
        self.calculate(session,'search');self.calculate(session,'msa_run','hits');self.calculate(session,'make_tree')
        alignment=next(row for row in entries(session) if row.get('analysis')=='msa')
        tree=entries(session,'result')[-1]
        session.a.identity=100;session.a.coverage=100
        self.calculate(session,'apply_filters')
        root,data=read_project(alignment['path'])
        self.assertTrue(local_path(root,data['msa']).is_file())
        self.assertTrue((tree['path']/'tree_fasttree.nwk').is_file())

    def test_save_and_reopen_preserve_catalog_and_inputs_without_engines(self):
        session=self.session();self.calculate(session,'msa_run','all');self.calculate(session,'make_tree')
        project=save_project(session,self.folder/'saved')
        self.source.unlink()
        with patch('av_tools.Toolchain.run') as engine,redirect_stdout(io.StringIO()):
            restored=resume_project(project,parse_args(['resume',str(project),'--no-menu']),UI(False))
        self.addCleanup(restored.close);engine.assert_not_called()
        self.assertEqual(len(entries(restored,'data')),1)
        self.assertEqual(len(entries(restored,'result')),2)
        self.assertEqual(entries(restored,'result')[-1]['inputs'],[entries(restored,'result')[0]['id']])
        for row in entries(restored):read_project(row['path'])

    def test_back_restores_cursor_filter_and_viewer_position(self):
        session=self.session();workspace=Workspace(session);workspace.switch('Sequences')
        workspace.active=2;workspace.query='c';workspace.switch('Results');workspace.query='alignment'
        self.assertTrue(workspace.back())
        self.assertEqual((workspace.tab,workspace.active,workspace.query),('Sequences',2,'c'))
        workspace.switch('Results');self.assertEqual(workspace.query,'alignment')
        self.calculate(session,'msa_run','all')
        session.ui.tty=True;session.a.no_menu=False
        with patch('av_session.MSAViewer.run',return_value=None):session.show_msa_window()
        session._viewer.col=70;session._viewer.row=1;session._viewer.checked={'b'}
        workspace.remember();session._viewer=None
        with patch('av_session.MSAViewer.run',return_value=None):session.show_msa_window()
        self.assertEqual((session._viewer.col,session._viewer.row,session._viewer.checked),(70,1,{'b'}))

    def test_html_preview_never_asks_for_a_destination(self):
        session=self.session();session.a.no_menu=False;session.ui.tty=True
        with patch.object(session.ui,'ask',side_effect=AssertionError('unexpected destination prompt')),\
             patch.object(session,'_open_local_html') as browser,redirect_stdout(io.StringIO()):
            session.open_report()
        self.assertTrue(browser.call_args.kwargs['temporary'])
        self.assertTrue(browser.call_args.args[0].is_file())

    def test_new_analysis_uses_highlighted_sequence(self):
        session=self.session();bench=Workbench(session);bench.workspace=Workspace(session);bench.workspace.switch('Sequences')
        bench.workspace.active=1
        with patch.object(bench,'choose',return_value='motif'):bench.new_analysis()
        self.assertEqual(session.pending,('query_analysis',('b','motif')))

    def test_reusing_translation_keeps_its_source_link_without_export(self):
        from av_app import _apply_start
        session=self.session();self.calculate(session,'translate_frames')
        origin=session.workspace_entry['id']
        bench=Workbench(session);bench.workspace=Workspace(session);bench.reuse_result()
        request=session.pending[1]
        args=parse_args([]);args.no_menu=True;_apply_start(args,request,self.folder)
        reused=Session(args,UI(False));self.addCleanup(reused.close)
        reused.library_root=session.library_root;reused.history=session.history
        _register_inputs(reused)
        self.assertEqual(reused.workspace_entry['inputs'],[origin])
        self.assertEqual(len(reused.query_list),6)
        self.assertFalse(any(self.folder.glob('export*')))

    def test_small_terminal_views_do_not_execute_or_save(self):
        session=self.session();workspace=Workspace(session)
        with patch.object(session.tools,'run') as engine,patch('av_project.save_project') as save:
            for width,height in ((40,14),(60,18),(80,24),(110,32)):
                for view in VIEWS:
                    workspace.switch(view);frame=workspace.render(width,height)
                    self.assertLessEqual(len(frame),height-1)
                    self.assertFalse(any('Next [n]' in line or 'Choose a goal' in line for line in frame))
            engine.assert_not_called();save.assert_not_called()

    def test_terminal_explores_reopens_and_closes_without_a_save_destination(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.source,),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        try:
            terminal.expect('File loaded — choose an action');terminal.send('2\r')
            terminal.expect('Create an alignment');terminal.send('1\r')
            terminal.expect('Alignment ready');terminal.send('1\r')
            terminal.expect('AlignView / Alignment');terminal.send('g');terminal.expect('Go to MSA column')
            terminal.send('70\r');terminal.expect('MSA column 70');terminal.send('q')
            terminal.expect('Alignment ready');terminal.send('4\r')
            terminal.expect('Choose a previous result');terminal.send('1\r')
            terminal.expect('Alignment ready');terminal.send('1\r')
            terminal.expect('cursor 70');terminal.send('q');terminal.expect('Alignment ready')
            terminal.send('3\r');terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Close workspace?');terminal.send('2\r')
            terminal.expect('without saving');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertNotIn(b'Destination',terminal.output)
            self.assertNotIn(b'Choose your goal',terminal.output)
        finally:terminal.close()

    def test_browsing_saved_results_does_not_create_another_save_requirement(self):
        session=self.session();self.calculate(session,'msa_run','all')
        project=save_project(session,self.folder/'saved')
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['resume',sys.argv[2],'--no-color']))"
        terminal=Terminal(code,args=(project,))
        try:
            terminal.expect('Alignment ready');terminal.send('4\r')
            terminal.expect('Choose a previous result');terminal.send('1\r')
            terminal.expect('Alignment ready');terminal.send('3\r')
            terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Project retained:');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertNotIn(b'Close workspace?',terminal.output)
            self.assertNotIn(b'Destination',terminal.output)
        finally:terminal.close()

    def test_added_dataset_and_analysis_cancel_keep_both_files_accessible(self):
        second=self.folder/'proteins.faa'
        write_fasta([Seq('p1','p1','MKWVTFISLL'),Seq('p2','p2','MKWVTFISLA')],second)
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color']))"
        terminal=Terminal(code,args=(self.source,))
        try:
            terminal.expect('File loaded — choose an action');terminal.send('o')
            terminal.expect('Choose a file');terminal.send('p');terminal.expect('FASTA path')
            terminal.send(str(second)+'\r');terminal.expect('File loaded — choose an action')
            terminal.expect('2 protein sequences');terminal.send('d')
            terminal.expect('Choose a file');terminal.expect('sequences.fna');terminal.expect('proteins.faa')
            terminal.send('1\r');terminal.expect('File loaded — choose an action');terminal.send('4\r')
            terminal.expect('Other analyses');terminal.send('\x1bOQ')
            terminal.expect('File loaded — choose an action');terminal.send('d')
            terminal.expect('Choose a file');terminal.expect('proteins.faa');terminal.send('q')
            terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Session closed.');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertNotIn(b'Project directory',terminal.output)
        finally:terminal.close()


if __name__=='__main__':unittest.main()
