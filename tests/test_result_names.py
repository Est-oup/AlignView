"""Result names are editable workspace metadata, independent of evidence."""
from contextlib import redirect_stdout
import io
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
from av_library import entries,rename_result
from av_navigation import Workbench
from av_project import read_project,resume_project,save_project
from av_session import Session
from av_ui import UI
from av_workspace import Workspace


class ResultNameTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.folder=Path(temporary.name)
        source=self.folder/'aligned.fa'
        write_fasta([Seq('a','a','ACGT'),Seq('b','b','ACGA')],source)
        self.session=Session(parse_args(['view',str(source),'--no-menu']),UI(False))
        self.addCleanup(self.session.close)
        self.session.library_root=self.folder/'catalog'
        with redirect_stdout(io.StringIO()):
            _register_inputs(self.session);self.session.load_msa()
        self.entry=entries(self.session,'result')[0]

    def test_rename_refreshes_cached_catalog_without_changing_evidence(self):
        s=self.session;workspace=Workspace(s)
        workspace.catalog('result')
        before={p.relative_to(self.entry['path']):p.read_bytes() for p in self.entry['path'].rglob('*') if p.is_file()}
        s.dirty=False
        rename_result(s,self.entry['id'],'Alignment with conserved region')
        workspace.remember()
        self.assertEqual(workspace.catalog('result')[0]['label'],'Alignment with conserved region')
        self.assertEqual(workspace.landing()[3],'Alignment with conserved region')
        self.assertTrue(s.dirty)
        after={p.relative_to(self.entry['path']):p.read_bytes() for p in self.entry['path'].rglob('*') if p.is_file()}
        self.assertEqual(before,after)
        read_project(self.entry['path'])
        self.assertEqual(s.result_sources['msa'],self.entry['id'])

    def test_name_survives_save_resume_and_opening_older_result(self):
        s=self.session;name='Alignment for publication'
        rename_result(s,self.entry['id'],name)
        project=save_project(s,self.folder/'saved')
        with redirect_stdout(io.StringIO()):
            restored=resume_project(project,parse_args(['resume',str(project),'--no-menu']),UI(False))
        self.addCleanup(restored.close)
        self.assertEqual(entries(restored,'result')[0]['label'],name)
        row=entries(restored,'result')[0]
        with redirect_stdout(io.StringIO()):
            child=resume_project(row['path'],parse_args(['resume',str(row['path']),'--no-menu']),UI(False))
        self.addCleanup(child.close)
        child.navigation_state=restored.navigation_state
        self.assertEqual(Workspace(child).landing()[3],name)

    def test_invalid_unknown_and_unchanged_names_leave_state_intact(self):
        s=self.session;s.dirty=False
        original=dict(s.navigation_state)
        for label in ('','   ','a\nb','\x1b[2J','x'*121,None):
            with self.subTest(label=label),self.assertRaises(ValueError):rename_result(s,self.entry['id'],label)
        with self.assertRaises(ValueError):rename_result(s,'missing','Name')
        with self.assertRaises(ValueError):rename_result(s,entries(s,'data')[0]['id'],'Name')
        rename_result(s,self.entry['id'],self.entry['label'])
        self.assertEqual(s.navigation_state,original)
        self.assertFalse(s.dirty)

    def test_canceling_rename_leaves_label_and_save_state_unchanged(self):
        s=self.session;s.dirty=False
        with patch.object(s.ui,'ask',return_value='q'):
            Workbench(s).rename_result()
        self.assertEqual(entries(s,'result')[0]['label'],self.entry['label'])
        self.assertFalse(s.dirty)

    def test_terminal_rename_is_visible_and_saved(self):
        from test_release import Terminal
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main([sys.argv[2],'--no-color','--module-policy','never']))"
        terminal=Terminal(code,args=(self.session.source_query_path,))
        self.addCleanup(terminal.close)
        terminal.expect('File loaded — choose an action');terminal.send('1\r')
        terminal.expect('Look at the loaded sequences');terminal.send('2\r')
        terminal.expect('Alignment');terminal.send('q')
        terminal.expect('Alignment ready');terminal.send('6\r')
        terminal.expect('Result name (q to cancel)');terminal.send('My reference alignment\r')
        terminal.expect('My reference alignment');terminal.send('3\r')
        terminal.expect('File loaded — choose an action');terminal.send('q')
        terminal.expect('Close workspace?');terminal.send('3\r')
        terminal.expect('Project directory')
        saved=self.folder/'terminal_saved';terminal.send(str(saved)+'\r')
        terminal.expect('Project saved:')
        self.assertEqual(terminal.process.wait(timeout=5),0)
        with redirect_stdout(io.StringIO()):
            restored=resume_project(saved,parse_args(['resume',str(saved),'--no-menu']),UI(False))
        self.addCleanup(restored.close)
        self.assertEqual(entries(restored,'result')[0]['label'],'My reference alignment')


if __name__=='__main__':unittest.main()
