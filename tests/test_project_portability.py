"""Reopened workspaces own their result catalog independently of saved sources."""
from contextlib import redirect_stdout
import io
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_app import _register_inputs
from av_cli import parse_args
from av_core import Seq,write_fasta
from av_library import entries
from av_project import read_project,resume_project,save_project
from av_session import Session
from av_ui import UI


class ProjectPortabilityTests(unittest.TestCase):
    def test_resumed_catalog_survives_removing_original_project_and_resaving(self):
        with tempfile.TemporaryDirectory() as temporary,redirect_stdout(io.StringIO()):
            folder=Path(temporary)
            source=folder/'aligned.fa'
            write_fasta([Seq('a','a','ACGT'),Seq('b','b','ACGA')],source)
            original=Session(parse_args(['view',str(source),'--no-menu']),UI(False))
            try:
                original.library_root=folder/'catalog'
                _register_inputs(original)
                original.load_msa()
                original.navigation_state={original.workspace_entry['id']:{'tab':'Alignment'}}
                expected=[(row['id'],row['kind'],row['label']) for row in entries(original)]
                project=save_project(original,folder/'saved')
            finally:original.close()
            # Reopening and browsing saved evidence never needs an engine.
            with patch('av_tools.Toolchain.run',side_effect=AssertionError('unexpected engine execution')):
                reopened=resume_project(project,parse_args(['resume',str(project),'--no-menu']),UI(False))
                try:
                    owned=list(reopened.history)
                    self.assertTrue(owned)
                    self.assertTrue(all(path.is_relative_to(Path(reopened.temporary.name)) for path in owned))
                    shutil.rmtree(project)
                    shutil.rmtree(folder/'catalog')
                    source.unlink()
                    self.assertEqual([(row['id'],row['kind'],row['label']) for row in entries(reopened)],expected)
                    for row in entries(reopened):
                        child=resume_project(row['path'],parse_args(['resume',str(row['path']),'--no-menu']),UI(False))
                        child.close()
                    saved_again=save_project(reopened,folder/'saved_again')
                finally:reopened.close()
                self.assertTrue(all(not path.exists() for path in owned))
                _,state=read_project(saved_again)
                self.assertEqual(state['navigation_state'],original.navigation_state)
                restored=resume_project(saved_again,parse_args(['resume',str(saved_again),'--no-menu']),UI(False))
                try:
                    self.assertEqual([(row['id'],row['kind'],row['label']) for row in entries(restored)],expected)
                    for row in entries(restored):read_project(row['path'])
                finally:restored.close()


if __name__=='__main__':unittest.main()
