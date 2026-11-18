"""Selected FASTA queries keep independent evidence and portable grouped results."""
from contextlib import redirect_stdout,redirect_stderr
import io
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from alignview import Session,parse_args
from av_app import _register_inputs
from av_core import Seq,write_fasta,read_fasta
from av_library import entries
from av_project import read_project,resume_project,save_project
from av_report import build_comparison_summary
from av_start import choose_sequences
from av_tools import Toolchain,CommandError
from av_ui import UI
from test_release import Terminal

ROOT=Path(__file__).resolve().parents[1]


class SelectedComparisonsTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory(prefix='alignview_selected_queries_')
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.query=self.folder/'queries.fna';self.refs=self.folder/'references.fna'
        self.records=[Seq('query_a','query_a short','ACGT'*20),Seq('query_b','query_b medium','ACGT'*30),Seq('query_c','query_c long','ACGT'*40)]
        write_fasta(self.records,self.query)
        write_fasta([Seq('ref_a','ref_a','ACGT'*45),Seq('ref_b','ref_b','ACGA'*45)],self.refs)

    def session(self):
        session=Session(parse_args([str(self.query),str(self.refs),'--query-id','query_a','--mode','blastn',
            '--no-menu','--module-policy','never','--tools-dir',str(ROOT/'tests/mockbin'),
            '--tools-dir',str(ROOT/'tests/mockextra')]),UI(False))
        self.addCleanup(session.close);session.library_root=self.folder/'catalog'
        _register_inputs(session)
        return session

    def compare(self,session,identifiers):
        with redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()):
            return session.compare_many(identifiers,'blastn')

    def test_only_checked_queries_are_run_with_separate_coverage_denominators(self):
        session=self.session();original={path:path.read_bytes() for path in (self.query,self.refs)}
        group=self.compare(session,['query_a','query_c'])
        _,state=read_project(group)
        self.assertEqual(state['comparison_selection'],['query_a','query_c'])
        self.assertEqual([run['query_id'] for run in state['comparison_batch']['runs']],['query_a','query_c'])
        self.assertFalse(state['hsps']);self.assertFalse(state['search_completed'])
        leaves=[row for row in entries(session,'result') if row['analysis']=='search']
        self.assertEqual(len(leaves),2)
        restored=[]
        for row,identifier in zip(leaves,['query_a','query_c']):
            _,leaf=read_project(row['path'])
            self.assertEqual(leaf['active_query_id'],identifier)
            self.assertEqual({hsp['query'] for hsp in leaf['hsps']},{identifier})
            with redirect_stdout(io.StringIO()):
                result=resume_project(row['path'],parse_args(['resume',str(row['path']),'--no-menu']),UI(False))
            self.addCleanup(result.close);restored.append(result)
        self.assertGreater(restored[0].coverage_rows()[0]['query_coverage'],restored[1].coverage_rows()[0]['query_coverage'])
        self.assertEqual(set(state['workspace_entry']['members']),{row['id'] for row in leaves})
        self.assertEqual({path:path.read_bytes() for path in original},original)

    def test_a_failed_query_does_not_hide_successful_results(self):
        session=self.session();original=Toolchain.search
        def search(toolchain,mode,path,*args,**kwargs):
            if read_fasta(path)[0].id=='query_c':raise CommandError('synthetic failure for query_c')
            return original(toolchain,mode,path,*args,**kwargs)
        with patch.object(Toolchain,'search',new=search):group=self.compare(session,['query_a','query_c'])
        _,state=read_project(group);runs=state['comparison_batch']['runs']
        self.assertEqual([run['status'] for run in runs],['completed','failed'])
        self.assertIn('synthetic failure',runs[1]['error'])
        self.assertTrue((group/'attempts/0002/commands.json').is_file())
        self.assertEqual(len(state['workspace_entry']['members']),1)

    def test_interrupt_keeps_completed_query_and_marks_remaining_queries(self):
        session=self.session();original=Toolchain.search
        def search(toolchain,mode,path,*args,**kwargs):
            if read_fasta(path)[0].id=='query_c':raise KeyboardInterrupt
            return original(toolchain,mode,path,*args,**kwargs)
        with patch.object(Toolchain,'search',new=search):group=self.compare(session,['query_a','query_c','query_b'])
        _,state=read_project(group)
        self.assertEqual([run['status'] for run in state['comparison_batch']['runs']],['completed','canceled','not_run'])
        self.assertEqual(len(state['workspace_entry']['members']),1)

    def test_group_and_member_results_survive_save_reopen_without_engines(self):
        session=self.session();group=self.compare(session,['query_a','query_c'])
        with redirect_stdout(io.StringIO()):
            current=resume_project(group,parse_args(['resume',str(group),'--no-menu']),UI(False))
        self.addCleanup(current.close);current.history=list(session.history)
        project=save_project(current,self.folder/'saved')
        self.query.unlink();self.refs.unlink()
        with patch('av_tools.Toolchain.run') as engine,redirect_stdout(io.StringIO()):
            restored=resume_project(project,parse_args(['resume',str(project),'--no-menu']),UI(False))
        self.addCleanup(restored.close);engine.assert_not_called()
        self.assertEqual(restored.comparison_selection,['query_a','query_c'])
        member_ids=set(restored.workspace_entry['members'])
        self.assertEqual(len([row for row in entries(restored) if row['id'] in member_ids]),2)
        for row in entries(restored):read_project(row['path'])
        with patch('av_tools.Toolchain.run') as engine,redirect_stdout(io.StringIO()):report=restored.make_report(open_preview=False)
        engine.assert_not_called();self.assertIn('not pooled',report.read_text())

    def test_rerunning_a_selection_keeps_the_previous_group(self):
        session=self.session();first=self.compare(session,['query_a','query_c'])
        original=(first/'comparison_batch.json').read_bytes()
        second=self.compare(session,['query_b','query_c'])
        self.assertNotEqual(first,second)
        self.assertEqual((first/'comparison_batch.json').read_bytes(),original)
        self.assertEqual(len([row for row in entries(session) if row.get('analysis')=='search_batch']),2)

    def test_invalid_scopes_are_rejected_before_running_any_engine(self):
        session=self.session()
        with patch.object(Toolchain,'search') as search:
            for selected in ([],['query_a'],['query_a','query_a'],['query_a','unknown']):
                with self.assertRaises(ValueError):self.compare(session,selected)
            search.assert_not_called()

    def test_simple_selector_supports_one_many_all_and_cancel(self):
        ui=UI(False)
        for answer,expected in [('2',['query_b']),('1,3',['query_a','query_c']),('1-3',['query_a','query_b','query_c']),('all',['query_a','query_b','query_c'])]:
            with patch.object(ui,'ask',return_value=answer),redirect_stdout(io.StringIO()):selected=choose_sequences(ui,self.records)
            self.assertEqual([record.id for record in selected],expected)
        with patch.object(ui,'ask',return_value='q'),redirect_stdout(io.StringIO()):
            self.assertIsNone(choose_sequences(ui,self.records,['query_b']))

    def test_group_html_preserves_identifiers_as_data(self):
        identifier='β_<script>alert(1)</script>'
        batch={'runs':[{'query_id':identifier,'status':'completed','matches':0}]}
        report=build_comparison_summary(self.folder/'html',batch,{})
        source=report.read_text();self.assertNotIn('<script>',source)
        self.assertIn('β_&lt;script&gt;',source)

    def test_filter_keeps_checked_queries_that_are_no_longer_visible(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_start import choose_sequences;from av_core import Seq;from av_ui import UI\nrows=[Seq('query_a','query_a','ACGT'),Seq('query_b','query_b','ACGA'),Seq('query_c','query_c','ACGG')];chosen=choose_sequences(UI(False),rows,['query_a']);print('PICKED='+','.join(r.id for r in chosen) if chosen else 'CANCELED',flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('Choose sequences to compare');terminal.send('/query_c\r \r')
            terminal.expect('PICKED=query_a,query_c');self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()

    def test_canceling_checkbox_selection_runs_nothing(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color']))"
        terminal=Terminal(code,args=(self.query,))
        try:
            terminal.expect('File loaded — choose an action');terminal.send('3\r')
            terminal.expect('Choose sequences to compare');terminal.send('a\x1b')
            terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Session closed.');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertNotIn(b'Choose the target dataset',terminal.output)
            self.assertNotIn(b'SEARCH / EXECUTION',terminal.output)
        finally:terminal.close()

    def test_interactive_compare_command_offers_multiple_queries_before_running(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['compare',sys.argv[2],sys.argv[3],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.query,self.refs))
        try:
            terminal.expect('Choose sequences to compare');terminal.send('a\r')
            terminal.expect('Find matching sequences');terminal.expect('3 selected sequences');terminal.send('q')
            terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Session closed.');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertNotIn(b'SEARCH / EXECUTION',terminal.output)
        finally:terminal.close()

    def test_terminal_selects_two_queries_and_opens_their_independent_results(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.query,),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        try:
            terminal.expect('File loaded — choose an action');terminal.send('3\r')
            terminal.expect('Choose sequences to compare');terminal.send(' \x1b[B\x1b[B \r')
            terminal.expect('Choose the target dataset');terminal.send('\x1b[F\r')
            terminal.expect('Choose a file');terminal.send('p');terminal.expect('FASTA path')
            terminal.send(str(self.refs)+'\r');terminal.expect('Find matching sequences')
            terminal.expect('2 selected sequences')
            self.assertNotIn(b'SEARCH / EXECUTION',terminal.output)
            terminal.send('1\r');terminal.expect('Sequence comparisons ready');terminal.expect('2/2 queries completed')
            terminal.send('1\r');terminal.expect('Results by query sequence')
            terminal.expect('query_a');terminal.expect('query_c');terminal.send('2\r')
            terminal.expect('Sequence matches ready');terminal.expect('Query: query_c')
            terminal.send('1\r');terminal.expect('/ Matches');terminal.send('\x1b')
            terminal.expect('Sequence matches ready');terminal.send('3\r')
            terminal.expect('File loaded — choose an action');terminal.send('q')
            terminal.expect('Close workspace?');terminal.send('2\r');terminal.expect('without saving')
            self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()


if __name__=='__main__':unittest.main()
