"""New-user journeys distinguish results, filters and readable aligned regions."""
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
from av_cli import parse_args
from av_core import HSP,Seq,write_fasta
from av_navigation import Workbench
from av_session import Session
from av_ui import UI
from av_workspace import Workspace
from test_release import Terminal


class UserJourneyTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.query=self.folder/'query.fna';self.target=self.folder/'targets.fna'
        write_fasta([Seq('query_one','query_one','ACGT'*30)],self.query)
        write_fasta([Seq('target_one','target_one','ACGT'*30),Seq('target_two','target_two','ACGA'*30)],self.target)

    def session(self):
        result=Session(parse_args([str(self.query),str(self.target),'--mode','blastn','--no-menu']),UI(False))
        self.addCleanup(result.close)
        return result

    def test_zero_hit_result_is_distinct_from_an_unrun_search_and_old_vsearch(self):
        session=self.session();view=Workspace(session);view.switch('Matches')
        self.assertIn('prepare a comparison','\n'.join(view._body([],80,10)))
        session.search_completed=True
        session.vsearch_results={'search':{'summary':{'matches':4}}}
        text='\n'.join(view._body([],80,10))
        self.assertIn('Search completed: no local aligned segments',text)
        self.assertNotIn('Native VSEARCH results',text)
        self.assertEqual(view.primary()[0],'compare')

    def test_filtered_out_segments_are_not_reported_as_a_failed_or_empty_engine_search(self):
        session=self.session()
        session.hsps=[HSP('query_one','target_one',75,4,1,4,1,4,1e-8,40,0,0,'ACGT','ACGA')]
        session.search_completed=True;session.a.identity=90
        with redirect_stdout(io.StringIO()):session.apply_filters()
        view=Workspace(session);view.switch('Matches')
        text='\n'.join(view._body([],80,10))
        self.assertIn('No matches pass the current filters',text)
        self.assertIn('1 local segment was found',text)
        self.assertNotIn('no local aligned segments found',text)

    def test_segment_inspection_preserves_provenance_reverse_direction_and_native_units(self):
        session=self.session()
        session.mode='blastx'
        session.hsps=[HSP('query_one','target_one',75,4,12,1,5,8,1e-8,40,-1,0,'MMMM','MMMA')]
        session.search_completed=True
        with redirect_stdout(io.StringIO()):session.apply_filters()
        bench=Workbench(session);bench.workspace=Workspace(session);bench.workspace.switch('Matches')
        with patch.object(session.ui,'page') as page,patch.object(session.tools,'run') as engine:
            bench.inspect_reference()
        content,title=page.call_args.args
        self.assertEqual(title,'Aligned regions · target_one')
        for expected in (str(self.query),str(self.target),'Query 12 → 1 nt','Target 5 → 8 aa',
                         'aligned rows: aa','1-based inclusive','Q = query, S = target'):
            self.assertIn(expected,content)
        engine.assert_not_called()
        self.assertEqual(bench.workspace.primary(),('inspect','View aligned regions'))

    def test_native_global_search_has_an_enabled_honest_export_choice(self):
        session=self.session()
        session.vsearch_results={'search':{'summary':{'matches':0}}}
        bench=Workbench(session)
        with patch.object(bench,'choose',return_value=None) as choose:
            bench.export()
        choice=next(option for option in choose.call_args.args[1] if option.key=='hits')
        self.assertTrue(choice.enabled)
        self.assertEqual(choice.label,'TSV global matches')
        self.assertIn('Native VSEARCH',choice.description)
        session.search_completed=True
        with patch.object(bench,'choose',return_value=None) as choose:bench.export()
        choice=next(option for option in choose.call_args.args[1] if option.key=='hits')
        self.assertEqual(choice.label,'TSV matches')

    def test_actual_search_opens_aligned_regions_and_back_restores_matches(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['compare',sys.argv[2],sys.argv[3],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.query,self.target),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        self.addCleanup(terminal.close)
        terminal.expect('Find matching sequences');terminal.expect('Target: targets.fna')
        terminal.send('1\r');terminal.expect('Sequence matches ready')
        terminal.expect('matching targets');terminal.send('1\r')
        terminal.expect('/ Matches');terminal.expect('View aligned regions');terminal.send('\r')
        terminal.expect('Aligned regions · target_one');terminal.expect('Query: query_one')
        terminal.expect('Target sequence: target_one')
        terminal.send('q');terminal.expect('/ Matches');terminal.send('q')
        terminal.expect('Sequence matches ready');terminal.send('3\r')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Close workspace?');terminal.send('2\r')
        terminal.expect('without saving');self.assertEqual(terminal.process.wait(timeout=5),0)
        self.assertEqual(terminal.output.count(b'SEARCH / EXECUTION'),1)


if __name__=='__main__':unittest.main()
