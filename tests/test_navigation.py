"""Acceptance tests for the simplified workflows, prerequisites and 100-hit limit."""
from __future__ import annotations

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
from alignview import Session,parse_args,_apply_start
from av_core import HSP,Seq,read_fasta,select_hits
from av_navigation import Workbench,WorkspaceState,HOME_ITEMS,item
from av_project import save_project,resume_project
from av_start import StartRequest
from av_ui import UI,MenuItem
from test_release import Terminal


class DefaultsTests(unittest.TestCase):
    def test_default_and_all_sensitivity_profiles_keep_one_hundred_hits(self):
        for profile in ('standard','sensitive','strict'):
            args=parse_args(['--demo','--no-menu','--preset',profile])
            self.assertEqual(args.hits,100)
            explicit=parse_args(['--demo','--no-menu','--preset',profile,'--hits','7'])
            self.assertEqual(explicit.hits,7)

    def test_hit_selection_really_keeps_one_hundred_references(self):
        hsps=[HSP('q',f'ref_{i:03d}',95,4,1,4,1,4,1e-20,200-i,0,0,'ACGT','ACGT') for i in range(120)]
        selected=select_hits(hsps,4)
        self.assertEqual(len(selected),100)
        self.assertEqual(list(selected)[-1],'ref_099')
        self.assertEqual(len(select_hits(hsps,4,max_hits=7)),7)

    def test_session_search_and_exports_keep_one_hundred(self):
        with tempfile.TemporaryDirectory() as tmp:
            query=Path(tmp)/'query.fa';query.write_text('>q\nACGT\n')
            refs=Path(tmp)/'refs.fa';refs.write_text(''.join(f'>ref_{i:03d}\nACGT\n' for i in range(120)))
            session=Session(parse_args([str(query),str(refs),'--mode','blastn','--no-menu']),UI(False))
            try:
                hsps=[HSP('q',f'ref_{i:03d}',95,4,1,4,1,4,1e-20,200-i,0,0,'ACGT','ACGT') for i in range(120)]
                with patch.object(session.tools,'search',return_value=hsps),redirect_stdout(io.StringIO()):session.search()
                self.assertEqual(len(session.selected),100)
                self.assertEqual(len((session.work/'hsps_selected.tsv').read_text().splitlines()),101)
                self.assertTrue(WorkspaceState.from_session(session).searched)
            finally:session.close()

    def test_changing_a_profile_does_not_change_a_custom_limit(self):
        session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--hits','73','--no-menu']),UI(False))
        try:
            workbench=Workbench(session)
            with patch.object(workbench,'choose',return_value='strict'),redirect_stdout(io.StringIO()):workbench.choose_preset()
            self.assertEqual(session.a.hits,73)
            self.assertEqual(session.a.evalue,1e-10)
        finally:session.close()


class WorkflowStateTests(unittest.TestCase):
    def session(self,*extra):
        return Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--no-menu',*extra]),UI(False))

    def test_single_file_active_sequence_does_not_reduce_alignment_input(self):
        session=self.session('--query-id','seq_3')
        try:
            self.assertEqual(session.query.id,'seq_3')
            self.assertEqual(len(read_fasta(session.query_search_path)),4)
        finally:session.close()

    def test_viewer_starts_on_the_active_sequence(self):
        session=self.session('--msa','view','--query-id','seq_3')
        try:
            session.load_msa();session.ui.tty=True
            with patch('av_msa_tui.MSAViewer.run',return_value=None):session.show_msa_window()
            self.assertEqual(session._viewer.row,2)
        finally:session.close()

    def test_view_intent_opens_alignment_without_recalculation(self):
        args=parse_args([])
        _apply_start(args,StartRequest('view',ROOT/'examples/example_msa.fasta'),Path('/tmp'))
        ui=UI(False);ui.tty=True
        session=Session(args,ui)
        try:
            with patch.object(session,'show_msa_window') as viewer,patch.object(session.tools,'msa') as aligner:
                session.initial()
                self.assertIsNotNone(session.msa_records)
                viewer.assert_called_once();aligner.assert_not_called()
        finally:session.close()

    def test_guided_comparison_asks_which_query_before_calculation(self):
        args=parse_args([])
        _apply_start(args,StartRequest('compare',ROOT/'examples/example_msa.fasta',ROOT/'examples/example_msa.fasta'),Path('/tmp'))
        ui=UI(False);ui.tty=True
        with patch.object(ui,'menu',return_value='3') as menu:
            session=Session(args,ui)
            try:
                with patch.object(Workbench,'compare') as compare,patch.object(session.tools,'search') as search:
                    session.initial();search.assert_not_called();compare.assert_called_once()
                self.assertEqual(session.query.id,'seq_3')
                self.assertEqual(session.a.hits,100)
                menu.assert_called_once()
                self.assertEqual(menu.call_args.args[0],'Choose the active query')
            finally:session.close()

    def test_unaligned_file_opened_as_alignment_stays_available_for_alignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'raw.fa';source.write_text('>a\nACGT\n>b\nACGTA\n')
            args=parse_args([]);_apply_start(args,StartRequest('view',source),Path(tmp))
            ui=UI(False);ui.tty=True
            session=Session(args,ui)
            try:
                session.initial()
                self.assertIsNone(session.msa_records)
                self.assertTrue(WorkspaceState.from_session(session).can_align)
                self.assertIn('Different lengths',ui.notice)
            finally:session.close()

    def test_zero_hit_search_is_still_completed_after_project_resume(self):
        args=parse_args([str(ROOT/'examples/example_query.fna'),str(ROOT/'examples/example_refs.faa'),
                         '--mode','blastx','--no-menu'])
        session=Session(args,UI(False))
        try:
            with patch.object(session.tools,'search',return_value=[]),redirect_stdout(io.StringIO()):session.search()
            with tempfile.TemporaryDirectory() as tmp:
                project=save_project(session,Path(tmp)/'project')
                with redirect_stdout(io.StringIO()):resumed=resume_project(project,parse_args(['--resume',str(project),'--no-menu']),UI(False))
                try:
                    self.assertTrue(WorkspaceState.from_session(resumed).searched)
                    self.assertEqual(len(resumed.selected),0)
                finally:resumed.close()
        finally:session.close()

    def test_existing_saved_hit_limit_is_preserved(self):
        session=self.session('--hits','10')
        try:
            with tempfile.TemporaryDirectory() as tmp:
                project=save_project(session,Path(tmp)/'project')
                with redirect_stdout(io.StringIO()):resumed=resume_project(project,parse_args(['--resume',str(project),'--no-menu']),UI(False))
                try:self.assertEqual(resumed.a.hits,10)
                finally:resumed.close()
        finally:session.close()

    def test_search_settings_changes_are_distinct_from_filter_changes(self):
        session=self.session()
        try:
            session.search_completed=True
            session.search_parameters={'mode':session.mode,'evalue':session.a.evalue,'sensitivity':session.a.sensitivity,'task':session.a.task}
            session.a.identity=90
            self.assertFalse(WorkspaceState.from_session(session).search_stale)
            session.a.evalue=0.01
            self.assertTrue(WorkspaceState.from_session(session).search_stale)
            self.assertEqual(session.settings()['search_run']['evalue'],1e-5)
        finally:session.close()

    def test_home_six_actions_fit_small_terminal_with_complete_labels(self):
        ui=UI(False);ui.context=('Demo file · 4 sequences','Alignment ready · unsaved changes')
        for columns,rows in ((80,24),(60,24),(110,32)):
            with patch('shutil.get_terminal_size',return_value=os.terminal_size((columns,rows))):
                frame=ui._menu_frame('Project',[item(*entry) for entry in HOME_ITEMS],hint='Open Results to read the alignment.')
                self.assertLessEqual(len(frame.lines),rows-1)
                for _,label,_ in HOME_ITEMS:self.assertTrue(any(label in line for line in frame.lines))


class NavigationTerminalTests(unittest.TestCase):
    demo="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['--demo','--no-color']))"

    def discard(self,terminal,dirty=False):
        if not dirty:
            terminal.send('q');terminal.expect('Session closed.')
            self.assertEqual(terminal.process.wait(timeout=5),0);return
        terminal.send('q');terminal.expect('Quitter AlignView ?');terminal.send('2\r')
        terminal.expect('without saving');self.assertEqual(terminal.process.wait(timeout=5),0)

    def test_demo_to_viewer_and_f2_home(self):
        terminal=Terminal(self.demo)
        try:
            terminal.expect('AlignView / Alignment');terminal.send('q')
            terminal.expect('Results ready');terminal.expect('Advanced')
            terminal.send('1\r');terminal.expect('Read the alignment')
            terminal.send('1\r');terminal.expect('AlignView / Alignment')
            terminal.send('\x1bOQ');terminal.expect('Back to Home')
            self.discard(terminal)
        finally:terminal.close()

    def test_f2_leaves_nested_settings_without_quitting_project(self):
        terminal=Terminal(self.demo)
        try:
            terminal.expect('AlignView / Alignment');terminal.send('q')
            terminal.expect('Results ready');terminal.send('6\r')
            terminal.expect('Additional analyses');terminal.send('3\r')
            terminal.expect('Calculation settings');terminal.send('1\r')
            terminal.expect('Search settings');terminal.send('5\r')
            terminal.expect('Filter matches')
            terminal.send('\x1b[12~');terminal.expect('Back to Home')
            self.discard(terminal)
        finally:terminal.close()

    def test_comparison_without_references_opens_reference_picker(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[1]+'/examples/example_query.fna','--no-color']))"
        terminal=Terminal(code)
        try:
            terminal.expect('Choose Compare');terminal.send('2\r')
            terminal.expect('Choose a file — Choose references')
            terminal.send('\x1bOQ');terminal.expect('Back to Home')
            self.discard(terminal)
        finally:terminal.close()

    def test_add_references_and_change_query_before_comparison(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[1]+'/examples/example_msa.fasta','--module-policy','never','--no-color']))"
        terminal=Terminal(code)
        try:
            terminal.expect('Choose Compare');terminal.send('2\r')
            terminal.expect('Choose a file — Choose references');terminal.send('p')
            terminal.expect('FASTA path');terminal.send(str(ROOT/'examples/example_msa.fasta')+'\r')
            terminal.expect('Choose the active query');terminal.send('3\r')
            terminal.expect('Run comparison');terminal.expect('seq_3')
            terminal.send('2\r');terminal.expect('Choose the active query');terminal.send('2\r')
            terminal.expect('Selected query: seq_2');terminal.expect('Run comparison')
            self.assertNotIn(b'Search termin',terminal.output)
            terminal.send('q');terminal.expect('Choose Compare')
            self.discard(terminal,dirty=True)
        finally:terminal.close()

    def test_nested_details_page_restores_the_menu_and_keyboard(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_ui import UI\nu=UI(False);result=u.menu('Long list',[(str(i),'Sequence '+str(i),'Sequence description','') for i in range(250)]);print('SELECTED',result,flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('Long list');terminal.send('\x1bOP')
            terminal.expect('Details');terminal.send('q')
            terminal.expect('Long list');terminal.send('\x1b[F\r')
            terminal.expect('SELECTED 249');self.assertEqual(terminal.process.wait(timeout=5),0)
            self.assertEqual(terminal.output.count(b'\x1b[?1049h'),1)
            self.assertEqual(terminal.output.count(b'\x1b[?1049l'),1)
        finally:terminal.close()

    def test_guided_alignment_has_one_preview_then_opens_result(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;from unittest.mock import patch\nfrom av_start import StartRequest\nwith patch('alignview.startup',return_value=StartRequest('align',__import__('pathlib').Path(sys.argv[1])/'examples/example_msa.fasta')):\n raise SystemExit(main(['--module-policy','never','--no-color']))"
        terminal=Terminal(code,env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
        try:
            terminal.expect('Create an alignment');terminal.send('1\r')
            terminal.expect('AlignView / Alignment');terminal.send('q')
            terminal.expect('Results ready')
            self.discard(terminal,dirty=True)
        finally:terminal.close()

    def test_file_picker_can_return_to_start_with_f2(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['--no-color']))"
        terminal=Terminal(code)
        try:
            terminal.expect('What would you like to do?');terminal.send('1\r')
            terminal.expect('Choose a file');terminal.send('\x1bOQ')
            terminal.expect('What would you like to do?');terminal.send('q')
            self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()

    def test_reopening_viewer_does_not_repeat_a_previous_report_request(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_core import Seq;from av_msa_tui import MSAViewer\nv=MSAViewer([Seq('a','a','ACGT'),Seq('b','b','ACGA')],False);print('FIRST',v.run(),flush=True);print('SECOND',v.run(),flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('AlignView / Alignment');terminal.send('w')
            terminal.expect('FIRST web');terminal.expect('AlignView / Alignment');terminal.send('q')
            terminal.expect('SECOND None');self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()

    def test_f2_after_calculation_still_requires_saving_new_results(self):
        with tempfile.TemporaryDirectory() as tmp:
            session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--no-menu']),UI(False))
            try:project=save_project(session,Path(tmp)/'saved')
            finally:session.close()
            code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['--resume',sys.argv[2],'--module-policy','never','--no-color']))"
            terminal=Terminal(code,args=(project,),env={'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']})
            try:
                terminal.expect('Choose Compare');terminal.send('3\r')
                terminal.expect('Create an alignment');terminal.send('1\r')
                terminal.expect('AlignView / Alignment');terminal.send('\x1bOQ')
                terminal.expect('unsaved changes')
                self.discard(terminal,dirty=True)
                self.assertTrue(project.is_dir())
            finally:terminal.close()


if __name__=='__main__':unittest.main()
