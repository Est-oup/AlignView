"""Users can identify and change the local search target before running."""
from contextlib import redirect_stdout
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_app import _register_inputs
from av_cli import parse_args
from av_core import Seq,write_fasta
from av_navigation import Workbench
from av_project import read_project,resume_project
from av_session import Session
from av_ui import UI,menu_item
from av_workspace import Workspace,comparison_context
from test_release import Terminal


class ComparisonOrientationTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.query=self.folder/'my_queries.fna'
        self.target=self.folder/'first_targets.fna'
        self.other=self.folder/'other_targets.fna'
        write_fasta([Seq('q1','q1','ACGT'*30),Seq('q2','q2','ACGA'*30)],self.query)
        write_fasta([Seq('ref1','ref1','ACGT'*30),Seq('ref2','ref2','ACGA'*30)],self.target)
        write_fasta([Seq('other1','other1','ACGG'*30),Seq('other2','other2','ACGA'*30)],self.other)

    def session(self):
        result=Session(parse_args([str(self.query),str(self.target),'--query-id','q1','--mode','blastn','--no-menu']),UI(False))
        self.addCleanup(result.close)
        result.library_root=self.folder/'catalog'
        with redirect_stdout(io.StringIO()):_register_inputs(result)
        return result

    def test_source_and_target_remain_visible_in_small_preparation_and_result_views(self):
        session=self.session();session.comparison_selection=['q1','q2']
        bench=Workbench(session)
        with session.ui.section('Comparison'):bench.refresh_context()
        self.assertIn('2 selected sequences',session.ui.context[0])
        for width,height in ((80,24),(60,18),(40,14)):
            frame=session.ui._menu_frame('Find matching sequences',[menu_item('1','Run BLASTN search')],dimensions=(width,height))
            self.assertLessEqual(len(frame.lines),height-1)
            text='\n'.join(frame.lines)
            self.assertIn('my_queries.fna',text)
            self.assertIn('Target: first_targets.fna',text)
            view=Workspace(session);view.switch('Matches')
            result='\n'.join(view.render(width,height))
            self.assertIn('Target: first_targets.fna',result)
        self.assertEqual(comparison_context(session)[1],'Target: first_targets.fna · 2 sequences')

    def test_canceling_a_query_draft_does_not_relabel_existing_evidence(self):
        session=self.session()
        with patch.object(session.tools,'search',return_value=[]),redirect_stdout(io.StringIO()):session.search()
        bench=Workbench(session)
        with patch.object(bench,'choose',side_effect=['2',None]),patch('av_navigation.choose_sequences',return_value=session.query_list):
            self.assertFalse(bench.search())
        self.assertEqual(session.comparison_selection,['q1','q2'])
        self.assertEqual(comparison_context(session)[0],'Query: q1 · my_queries.fna')
        bench.refresh_context()
        self.assertEqual(session.ui.context[0],'Query: q1 · my_queries.fna')
        frame='\n'.join(Workspace(session).render(80,24))
        self.assertIn('Query: q1',frame)
        self.assertNotIn('2 selected sequences',frame)

    def test_review_inputs_shows_full_paths_and_selected_ids_without_an_engine(self):
        session=self.session();session.comparison_selection=['q2']
        bench=Workbench(session)
        with patch.object(session.ui,'page') as page,patch.object(session.tools,'run') as engine:
            bench.comparison_inputs('blastn')
        content=page.call_args.args[0]
        self.assertIn(str(self.query),content);self.assertIn(str(self.target),content)
        self.assertIn('q2 · 120 characters',content);self.assertNotIn('q1 ·',content)
        self.assertIn('all 2 target sequences',content);engine.assert_not_called()
        self.assertIn('ref1 · 120 characters · DNA/RNA',content)

    def test_reopening_an_earlier_result_shows_its_original_target(self):
        session=self.session()
        with patch.object(session.tools,'search',return_value=[]),redirect_stdout(io.StringIO()):
            session.search()
        snapshot=session.history[-1]
        session.source_refs_path=self.other
        with redirect_stdout(io.StringIO()):
            restored=resume_project(snapshot,parse_args(['resume',str(snapshot),'--no-menu']),UI(False))
        self.addCleanup(restored.close)
        self.assertIn('first_targets.fna',comparison_context(restored)[1])
        self.assertNotIn('other_targets',comparison_context(restored)[1])

    def test_canceling_target_picker_keeps_previous_target_and_query_selection(self):
        session=self.session();session.comparison_selection=['q1','q2']
        bench=Workbench(session);bench.workspace=Workspace(session)
        with patch.object(bench,'choose',side_effect=['6',None,None]),patch.object(session.tools,'run') as engine:
            self.assertFalse(bench.search())
        self.assertEqual(session.source_refs_path,self.target)
        self.assertEqual(session.comparison_selection,['q1','q2'])
        self.assertIsNone(session.pending);engine.assert_not_called()

    def test_terminal_changes_target_preserves_queries_and_saves_the_actual_search(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['compare',sys.argv[2],sys.argv[3],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.query,self.target),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        self.addCleanup(terminal.close)
        terminal.expect('Choose sequences to compare');terminal.send('a\r')
        terminal.expect('Find matching sequences');terminal.expect('Target: first_targets.fna')
        terminal.send('6\r');terminal.expect('Choose the target dataset');terminal.send('\x1b[F\r')
        terminal.expect('Choose a file — Choose target FASTA');terminal.send('p')
        terminal.expect('FASTA path');terminal.send(str(self.other)+'\r')
        terminal.expect('Find matching sequences');terminal.expect('2 selected sequences')
        terminal.expect('Target: other_targets.fna')
        self.assertNotIn(b'SEARCH / EXECUTION',terminal.output)
        terminal.send('1\r');terminal.expect('Sequence comparisons ready')
        terminal.expect('Target: other_targets.fna');terminal.expect('2/2 queries completed')
        terminal.send('1\r');terminal.expect('Results by query sequence');terminal.send('2\r')
        terminal.expect('Sequence matches ready');terminal.expect('Query: q2')
        terminal.expect('Target: other_targets.fna')
        terminal.send('3\r');terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Close workspace?');terminal.send('3\r')
        terminal.expect('Project directory');saved=self.folder/'saved';terminal.send(str(saved)+'\r')
        terminal.expect('Project saved:');self.assertEqual(terminal.process.wait(timeout=5),0)
        _,state=read_project(saved)
        self.assertEqual(state['inputs']['references']['original_path'],str(self.other))
        self.assertEqual(state['active_query_id'],'q2')
        self.assertEqual(state['mode'],'blastn')
        groups=[]
        for archive in (saved/'history').iterdir():
            _,item=read_project(archive)
            if item.get('comparison_batch'):
                groups.append(item)
        self.assertTrue(groups)
        self.assertEqual(groups[-1]['comparison_selection'],['q1','q2'])
        self.assertEqual(groups[-1]['inputs']['references']['original_path'],str(self.other))


if __name__=='__main__':unittest.main()
