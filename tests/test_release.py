"""Regression coverage for the 0.1 workflows and real terminal escape sequences."""
from __future__ import annotations

import argparse
import contextlib
import gzip
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,main,parse_args
from av_core import Seq,HSP,read_fasta
from av_msa_tui import MSAViewer
from av_project import read_project,resume_project,save_project
from av_setup import diagnostic,install_plan,setup_tools
from av_start import choose_file,is_fasta_name
from av_terminal import fit
from av_tools import Toolchain,CommandError
from av_ui import UI


class Terminal:
    """A real PTY, using only the standard library (no optional test dependency)."""
    def __init__(self,code,dimensions=(24,80),args=(),env=None):
        import fcntl,pty,struct,termios
        self.master,slave=pty.openpty()
        fcntl.ioctl(slave,termios.TIOCSWINSZ,struct.pack('HHHH',*dimensions,0,0))
        self.process=subprocess.Popen([sys.executable,'-c',code,str(ROOT),*map(str,args)],
            stdin=slave,stdout=slave,stderr=slave,env={**os.environ,'TERM':'xterm-256color',**(env or {})})
        os.close(slave);self.output=b'';self.offset=0

    def expect(self,needle,timeout=10):
        import select
        deadline=time.monotonic()+timeout
        needle=needle.encode()
        while time.monotonic()<deadline:
            found=self.output.find(needle,self.offset)
            if found>=0:self.offset=found+len(needle);return
            if select.select([self.master],[],[],.1)[0]:
                try:chunk=os.read(self.master,65536)
                except OSError:break
                self.output+=chunk
        raise AssertionError(f'Missing {needle!r} in PTY output: {self.output[-4000:]!r}')

    def send(self,value):os.write(self.master,value.encode())

    def resize(self,rows,columns):
        import fcntl,struct,termios
        fcntl.ioctl(self.master,termios.TIOCSWINSZ,struct.pack('HHHH',rows,columns,0,0))

    def close(self):
        if self.process.poll() is None:
            self.process.terminate()
            try:self.process.wait(timeout=3)
            except subprocess.TimeoutExpired:self.process.kill();self.process.wait()
        os.close(self.master)


class ViewerReleaseTests(unittest.TestCase):
    def viewer(self):return MSAViewer([Seq(f's{i}',f's{i}','ACGT'*300) for i in range(100)],False)

    def test_frame_fits_common_terminal_sizes(self):
        viewer=self.viewer()
        for columns,rows in ((80,24),(110,32),(60,18),(40,12),(160,50)):
            with patch('shutil.get_terminal_size',return_value=os.terminal_size((columns,rows))):
                frame=viewer._frame()
                self.assertLessEqual(len(frame),rows-1)
                for line in frame:
                    self.assertLessEqual(len(fit(line,columns-1)),columns-1)

    def test_fast_horizontal_and_vertical_navigation(self):
        viewer=self.viewer()
        with patch('shutil.get_terminal_size',return_value=os.terminal_size((80,24))):
            viewer._handle_key('ctrl-right')
            self.assertGreater(viewer.xoff,20)
            viewer._handle_key('pagedown')
            self.assertEqual(viewer.row,8)
            viewer._handle_key('end');viewer._handle_key('ctrl-end')
            self.assertEqual((viewer.col,viewer.row),(1199,99))
            viewer._handle_key('home');viewer._handle_key('ctrl-home')
            self.assertEqual((viewer.col,viewer.row),(0,0))

    def test_mouse_wheel_and_overview_jump(self):
        viewer=self.viewer()
        viewer._handle_key('mouse:65;30;10M')
        self.assertEqual(viewer.row,3)
        viewer._handle_key('mouse:69;30;10M')
        self.assertGreater(viewer.col,1)
        viewer._handle_key('mouse:0;75;3M')
        self.assertGreater(viewer.col,600)

    def test_search_preserves_exact_invariant_coordinate(self):
        viewer=MSAViewer([Seq('a','a','ATGCC'),Seq('b','b','ATGCT')],False)
        viewer.only_variable=True
        viewer._search('ATG')
        self.assertFalse(viewer.only_variable)
        self.assertEqual(viewer.cols[viewer.col],0)

    def test_invalid_alignment_has_actionable_error(self):
        for records in ([],[Seq('a','a','')],[Seq('a','a','AA'),Seq('b','b','A')]):
            with self.assertRaisesRegex(ValueError,'aligned'):MSAViewer(records)

    def test_utf8_and_modified_arrows_are_read_as_whole_keys(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_terminal import screen,read_key\nwith screen():\n print('READY',flush=True)\n for i in range(3):\n  key=read_key(timeout=5);print('KEY='+key,flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('READY');terminal.send('\x1b[1;5Cλ\x1b[6~')
            terminal.expect('KEY=ctrl-right');terminal.expect('KEY=λ');terminal.expect('KEY=pagedown')
            self.assertEqual(terminal.process.wait(timeout=3),0)
        finally:terminal.close()

    def test_real_msa_terminal_pages_search_and_resize(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_core import Seq;from av_msa_tui import MSAViewer\nv=MSAViewer([Seq('s'+str(i),'s'+str(i),'ACGT'*300) for i in range(40)],False);v.run();print('FINAL',v.row,v.col,flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('AlignView / Alignment')
            terminal.send('\x1b[6~\x1b[1;5C')
            terminal.expect('Row 9/40')
            terminal.resize(32,100)
            terminal.send('g');terminal.expect('Go to MSA column');terminal.send('75%\n')
            terminal.expect('Colonne MSA 900')
            terminal.send('q');terminal.expect('FINAL 8 899')
            self.assertEqual(terminal.process.wait(timeout=3),0)
        finally:terminal.close()

    def test_long_menu_can_reach_last_item(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_ui import UI\nu=UI(False);r=u.menu('LONG MENU',[(str(i),'Sequence '+str(i),'Description','DATA') for i in range(250)]);print('SELECTED',r,flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('LONG MENU');terminal.send('\x1b[F\r')
            terminal.expect('SELECTED 249')
            self.assertEqual(terminal.process.wait(timeout=3),0)
        finally:terminal.close()

    def test_pager_reaches_bottom_and_scrolls_sideways(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from av_ui import UI\nUI(False).page('\\n'.join('row'+str(i)+' '+('X'*200) for i in range(300)),'PAGER');print('DONE',flush=True)"
        terminal=Terminal(code)
        try:
            terminal.expect('PAGER');terminal.send('\x1b[F\x1b[C')
            terminal.expect('300/300');terminal.expect('horizontal offset 8')
            terminal.send('q');terminal.expect('DONE')
        finally:terminal.close()

    def test_file_browser_lists_all_files(self):
        class Browser(UI):
            def menu(self,title,items,*args):
                self.count=len(items)
                return next(key for key,label,_,_ in items if 'seq_149.fasta' in label)
        with tempfile.TemporaryDirectory() as tmp:
            for i in range(150):(Path(tmp)/f'seq_{i:03d}.fasta').write_text('>a\nACGT\n')
            ui=Browser(False)
            with contextlib.redirect_stdout(io.StringIO()):chosen=choose_file(ui,'test',tmp)
            self.assertEqual(chosen.name,'seq_149.fasta');self.assertEqual(ui.count,152)


class ProjectReleaseTests(unittest.TestCase):
    def session(self):
        session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--msa','view','--no-menu','--module-policy','never']),UI(False))
        session.load_msa()
        return session

    def test_project_roundtrip_detaches_inputs_and_msa(self):
        session=self.session()
        with tempfile.TemporaryDirectory() as tmp:
            try:
                project=save_project(session,Path(tmp)/'saved')
                resumed=resume_project(project,parse_args(['--resume',str(project),'--no-menu']),UI(False))
                shutil.rmtree(project)
                try:
                    self.assertEqual(len(read_fasta(resumed.query_path)),4)
                    self.assertEqual(len(read_fasta(resumed.query_search_path)),4)
                    self.assertEqual(resumed.msa_records,session.msa_records)
                    final=save_project(resumed,Path(tmp)/'again')
                    self.assertTrue((final/'session.json').is_file())
                finally:resumed.close()
            finally:session.close()

    def test_session_keeps_original_inputs_after_source_is_deleted(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'original.fasta';source.write_text('>a\nACGT\n>b\nACGA\n')
            session=Session(parse_args([str(source),'--no-menu']),UI(False))
            try:
                source.unlink()
                project=save_project(session,Path(tmp)/'saved')
                self.assertEqual(len(read_fasta(project/'inputs/query.fasta')),2)
            finally:session.close()

    def test_real_workbench_can_switch_query_and_restore_history(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'project'
            code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main\nraise SystemExit(main([sys.argv[1]+'/examples/example_msa.fasta',sys.argv[1]+'/examples/example_msa.fasta','--mode','global-dna','--query-id','seq_1','--module-policy','never']))"
            terminal=Terminal(code,dimensions=(32,110))
            try:
                terminal.expect('Choose Compare');terminal.send('c')
                terminal.expect('Files and query');terminal.send('1\r')
                terminal.expect('Choose the active query');terminal.send('2\r')
                terminal.expect('Selected query: seq_2')
                terminal.expect('Run comparison');terminal.send('q')
                terminal.expect('Choose Compare');terminal.send('c')
                terminal.expect('Files and query');terminal.send('3\r')
                terminal.expect('Previous analyses');terminal.send('1\r')
                terminal.expect('Project resumed: seq_1')
                terminal.expect('Choose Compare');terminal.send('s')
                terminal.expect('Project directory');terminal.send(str(target)+'\n')
                terminal.expect('Project saved:')
                terminal.expect('Choose Compare');terminal.send('x')
                terminal.expect('Project retained:');self.assertEqual(terminal.process.wait(timeout=5),0)
                _,state=read_project(target)
                self.assertEqual(state['active_query_id'],'seq_1')
                self.assertGreaterEqual(len(list((target/'history').iterdir())),2)
            finally:terminal.close()

    def test_changing_filters_invalidates_stale_alignment_and_tree(self):
        session=self.session()
        try:
            session.single=False
            session.hsps=[HSP('q','r',90,4,1,4,1,4,1e-8,30,0,0,'ACGT','ACGT')]
            session.selected={'r':session.hsps}
            (session.work/'tree_fasttree.nwk').write_text('(a,b);')
            (session.work/'report.html').write_text('old report')
            session.a.identity=99
            with patch.object(session,'export_hsps'),patch.object(session,'save_settings'):
                session.apply_filters()
            self.assertIsNone(session.msa_records)
            self.assertFalse((session.work/'tree_fasttree.nwk').exists())
            self.assertFalse((session.work/'report.html').exists())
        finally:session.close()

    def test_checksum_rejects_changed_input(self):
        session=self.session()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                folder=save_project(session,Path(tmp)/'saved')
                (folder/'inputs/query.fasta').write_text('>tampered\nAAA\n')
                with self.assertRaisesRegex(ValueError,'changed'):read_project(folder)
        finally:session.close()

    def test_project_paths_cannot_escape_project(self):
        session=self.session()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                folder=save_project(session,Path(tmp)/'saved')
                state=json.loads((folder/'session.json').read_text())
                state['inputs']['query']['path']='../../outside.fasta'
                (folder/'session.json').write_text(json.dumps(state))
                with self.assertRaisesRegex(ValueError,'outside the project'):read_project(folder)
        finally:session.close()

    def test_save_never_overwrites_existing_folder(self):
        session=self.session()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                target=Path(tmp)/'saved';target.mkdir();(target/'keep').write_text('keep')
                with self.assertRaises(ValueError):save_project(session,target)
                self.assertEqual((target/'keep').read_text(),'keep')
        finally:session.close()

    def test_failed_copy_does_not_publish_half_a_project(self):
        session=self.session()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                target=Path(tmp)/'saved'
                with patch('av_project.shutil.copyfile',side_effect=OSError('disk full')):
                    with self.assertRaises(OSError):save_project(session,target)
                self.assertFalse(target.exists());self.assertEqual(list(Path(tmp).iterdir()),[])
        finally:session.close()

    def test_nested_history_is_not_recursively_copied(self):
        session=self.session()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                archive=save_project(session,Path(tmp)/'previous',include_history=False)
                session.history=[archive]
                final=save_project(session,Path(tmp)/'final')
                self.assertTrue((final/'history/001/session.json').is_file())
                self.assertFalse((final/'history/001/history').exists())
        finally:session.close()

    def test_cli_resume_does_not_require_a_search_engine(self):
        with tempfile.TemporaryDirectory() as tmp:
            initial=Path(tmp)/'initial';resumed=Path(tmp)/'resumed'
            env={**os.environ,'PATH':str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']}
            run=subprocess.run([sys.executable,str(ROOT/'alignview.py'),str(ROOT/'examples/example_query.fna'),
                str(ROOT/'examples/example_refs.faa'),'--mode','blastx','--msa','hits','--save','--no-menu','--module-policy','never','--out',str(initial)],capture_output=True,text=True,env=env,timeout=30)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            run=subprocess.run([sys.executable,str(ROOT/'alignview.py'),'--resume',str(initial),'--save','--no-menu',
                '--module-policy','never','--out',str(resumed)],capture_output=True,text=True,timeout=30)
            self.assertEqual(run.returncode,0,run.stdout+run.stderr)
            a=json.loads((initial/'session.json').read_text());b=json.loads((resumed/'session.json').read_text())
            self.assertEqual(a['hsps'],b['hsps']);self.assertEqual(a['commands'],b['commands'])
            self.assertEqual((initial/'alignment_active.fasta').read_text(),(resumed/'alignment_active.fasta').read_text())

    def test_gzip_fasta_reaches_an_external_aligner_as_plain_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'input.fasta.gz'
            with gzip.open(source,'wt') as stream:stream.write('>a\nACGT\n>b\nACG\n')
            self.assertTrue(is_fasta_name(source))
            args=parse_args([str(source),'--no-menu','--module-policy','never','--tools-dir',str(ROOT/'tests/mockbin')])
            session=Session(args,UI(False))
            try:
                with contextlib.redirect_stdout(io.StringIO()):session.msa_run('all')
                self.assertEqual(len(session.msa_records),2)
            finally:session.close()


class SetupReleaseTests(unittest.TestCase):
    def test_plan_uses_strict_channels_and_dedicated_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            plan=install_plan('muscle',Path(tmp)/'tools',manager='/bin/conda')
            self.assertIn('create',plan);self.assertIn('--override-channels',plan)
            self.assertLess(plan.index('conda-forge'),plan.index('bioconda'))
            self.assertEqual(plan[-1],'muscle>=5')

    def test_existing_conda_prefix_is_extended(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp)/'conda-meta').mkdir()
            self.assertEqual(install_plan('mafft',tmp,manager='/bin/conda')[1],'install')

    def test_unrelated_existing_folder_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'without a Conda environment'):install_plan('mafft',tmp,manager='/bin/conda')

    def test_plan_never_installs_without_explicit_flag(self):
        with patch('av_setup.install_plan',return_value=['conda','create','mafft']),patch('av_setup.subprocess.run') as run:
            with contextlib.redirect_stdout(io.StringIO()):setup_tools(UI(False),'mafft',install=False)
            run.assert_not_called()

    def test_missing_tool_suggests_its_specific_installation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tools=Toolchain(Path(tmp),auto_modules=False)
            with patch.object(tools,'available',return_value=False):
                with self.assertRaisesRegex(CommandError,'--setup blast --install-tools'):tools.require('tblastn')

    def test_module_alias_failure_does_not_hide_other_alias(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'bin').mkdir()
            program=root/'bin/muscle5';program.write_text('#!/bin/sh\nexit 0\n');program.chmod(0o755)
            with patch.dict(os.environ,{'PATH':'/usr/bin'}):tools=Toolchain(root/'work')
            with patch.object(tools,'_module_env',return_value={'PATH':str(root/'bin')}):
                self.assertFalse(tools.available('muscle'))
                self.assertTrue(tools.available('muscle5'))

    def test_managed_binary_is_found_without_shell_activation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'bin').mkdir()
            executable=root/'bin/mmseqs';executable.write_text('#!/bin/sh\necho test\n');executable.chmod(0o755)
            with patch.dict(os.environ,{'ALIGNVIEW_TOOLS_PREFIX':tmp,'PATH':'/usr/bin'}):
                tools=Toolchain(root/'work',auto_modules=False)
                self.assertEqual(tools.executable('mmseqs'),str(executable))
                self.assertEqual(tools.source('mmseqs'),'AlignView environment')

    def test_nan_and_infinite_parameters_are_rejected(self):
        for flag,value in (('--evalue','nan'),('--evalue','inf'),('--bitscore','nan'),('--gapopen','-1'),('--ref','0')):
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):parse_args(['--demo','--no-menu',flag,value])


if __name__=='__main__':unittest.main()
