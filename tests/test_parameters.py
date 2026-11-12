"""Scientific controls must reach execution, survive portability, and describe results honestly."""
from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
import os
from pathlib import Path
import shlex
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,parse_args,main
from av_core import HSP,Seq,write_fasta
from av_navigation import Workbench,WorkspaceState
from av_parameters import relevant,load_profile,write_profile,snapshot
from av_project import save_project,resume_project
from av_tools import Toolchain,CommandError
from av_ui import UI
from test_release import Terminal


class ParameterTests(unittest.TestCase):
    def args(self,*extra):return parse_args(['--demo','--no-menu',*extra])

    def tools(self,path):return Toolchain(Path(path),auto_modules=False)

    def test_blast_scoring_masking_strand_are_used_without_emboss_leakage(self):
        a=self.args('--blast-word-size','7','--blast-reward','2','--blast-penalty','-3',
            '--blast-gapopen','5','--blast-gapextend','2','--blast-masking','no',
            '--blast-soft-masking','true','--blast-strand','minus','--blast-max-hsps','3',
            '--gapopen','17','--matrix','EBLOSUM80')
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp)
            argv=tools.search_commands('blastn',Path('q'),Path('r'),options=a)[-1].argv
            for flag,value in (('-word_size','7'),('-reward','2'),('-penalty','-3'),('-gapopen','5'),
                ('-gapextend','2'),('-dust','no'),('-soft_masking','true'),('-strand','minus'),('-max_hsps','3')):
                self.assertEqual(argv[argv.index(flag)+1],value)
            self.assertNotIn('17',argv);self.assertNotIn('-matrix',argv)
            argv=tools.search_commands('tblastx',Path('q'),Path('r'),options=a)[-1].argv
            self.assertNotIn('-gapopen',argv);self.assertNotIn('-gapextend',argv)
            self.assertNotIn('-reward',argv)

    def test_diamond_and_mmseqs_have_independent_sensitivity_and_native_options(self):
        a=self.args('--sensitivity','ultra-sensitive','--mmseqs-sensitivity','6.5',
            '--mmseqs-mask','0','--mmseqs-comp-bias','0','--diamond-masking','seg',
            '--diamond-max-hsps','0','--diamond-min-orf','20','--diamond-matrix','BLOSUM45')
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp)
            diamond=tools.search_commands('diamond-blastx',Path('q'),Path('r'),options=a,sensitivity=a.sensitivity)[-1].argv
            self.assertIn('--ultra-sensitive',diamond)
            self.assertEqual(diamond[diamond.index('--min-orf')+1],'20')
            self.assertEqual(diamond[diamond.index('--max-hsps')+1],'0')
            self.assertEqual(diamond[diamond.index('--masking')+1],'seg')
            mmseqs=tools.search_commands('mmseqs-blastp',Path('q'),Path('r'),max_target=731,options=a)[0].argv
            self.assertNotIn('--ultra-sensitive',mmseqs)
            self.assertIn('-a',mmseqs)
            for flag,value in (('-s','6.5'),('--max-seqs','731'),('--mask','0'),('--comp-bias-corr','0')):
                self.assertEqual(mmseqs[mmseqs.index(flag)+1],value)

    def test_native_defaults_are_omitted_and_task_incompatibility_fails_before_execution(self):
        a=self.args()
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp)
            argv=tools.search_commands('blastp',Path('q'),Path('r'),options=a)[-1].argv
            for flag in ('-word_size','-gapopen','-seg','-soft_masking','-matrix'):self.assertNotIn(flag,argv)
            with self.assertRaisesRegex(ValueError,'incompatible'):
                tools.search_commands('blastp',Path('q'),Path('r'),task='megablast',options=self.args('--task','megablast'))
            self.assertEqual(tools.history,[])

    def test_mafft_auto_rejects_ambiguous_iterations_and_other_engines_keep_own_knobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp)
            a=self.args('--mafft-maxiterate','50')
            with self.assertRaisesRegex(ValueError,'global, local or fast'):
                tools.msa_command(Path('q'),Path('o'),options=a)
            a.mafft='local';a.mafft_op=2.5;a.mafft_ep=.1
            argv=tools.msa_command(Path('q'),Path('o'),strategy=a.mafft,options=a).argv
            for flag,value in (('--maxiterate','50'),('--op','2.5'),('--ep','0.1')):self.assertEqual(argv[argv.index(flag)+1],value)
            a.aligner='muscle';a.muscle_method='super5';a.muscle_perm='abc';a.muscle_perturb=9
            argv=tools.msa_command(Path('q'),Path('o'),engine='muscle',options=a).argv
            self.assertIn('-super5',argv);self.assertNotIn('--op',argv)
            self.assertEqual(argv[argv.index('-perm')+1],'abc')
            a.aligner='clustalo';a.clustalo_iterations=2;a.clustalo_full=True;a.clustalo_full_iter=True
            argv=tools.msa_command(Path('q'),Path('o'),engine='clustalo',options=a).argv
            self.assertIn('--full',argv);self.assertIn('--full-iter',argv)
            self.assertEqual(argv[argv.index('--iterations')+1],'2')

    def test_tree_model_checks_alphabet_and_support_is_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp);a=self.args('--fasttree-model','GTR','--fasttree-gamma','--fasttree-boot','0','--fasttree-seed','42')
            with self.assertRaisesRegex(ValueError,'incompatible'):tools.fasttree_command(Path('q'),Path('o'),False,a)
            argv=tools.fasttree_command(Path('q'),Path('o'),True,a).argv
            for flag in ('-nt','-gtr','-gamma','-nosupport'):self.assertIn(flag,argv)
            self.assertNotIn('-boot',argv)
            self.assertEqual(argv[argv.index('-seed')+1],'42')

    def test_pairwise_end_penalties_are_separate_and_water_ignores_ends(self):
        a=self.args('--needle-endopen','3','--needle-endextend','.2')
        with tempfile.TemporaryDirectory() as tmp:
            tools=self.tools(tmp);q=Seq('q','q','MPEPTIDE');r=Seq('r','r','MPEPIDE')
            preview=tools.emboss('needle',q,r,Path(tmp)/'n',options=a,preview=True)
            self.assertIn('-gapopen 10',preview);self.assertIn('-endopen 3.0',preview)
            self.assertIn('-endextend 0.2',preview)
            preview=tools.emboss('water',q,r,Path(tmp)/'w',options=a,preview=True)
            self.assertNotIn('-endopen',preview);self.assertEqual(tools.history,[])

    def test_profiles_validate_atomically_and_cli_wins_including_negative_boolean(self):
        with tempfile.TemporaryDirectory() as tmp:
            filename=Path(tmp)/'profile.json'
            a=self.args('--hits','73','--mafft','local','--mafft-maxiterate','40','--clustalo-full','--blast-masking','no')
            write_profile(a,filename)
            loaded=self.args('--load-profile',str(filename),'--hits','100','--no-clustalo-full')
            self.assertEqual(loaded.mafft,'local');self.assertEqual(loaded.mafft_maxiterate,40)
            self.assertEqual(loaded.blast_masking,'no');self.assertEqual(loaded.hits,100)
            self.assertFalse(loaded.clustalo_full)
            self.assertEqual(self.args('--load-profile',str(filename)).hits,73)
            data=json.loads(filename.read_text());data['parameters']['query']='/tmp/unrelated'
            filename.write_text(json.dumps(data));before=vars(loaded).copy()
            with self.assertRaisesRegex(ValueError,'Unknown'):load_profile(loaded,filename)
            self.assertEqual(vars(loaded),before)
            data['parameters'].pop('query');data['parameters']['mmseqs_sensitivity']=float('nan')
            filename.write_text(json.dumps(data))
            with self.assertRaises(ValueError):load_profile(loaded,filename)
            self.assertEqual(vars(loaded),before)

    def test_profile_writer_does_not_run_tools_or_overwrite_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            file=Path(tmp)/'p.json'
            with patch.object(Toolchain,'run') as run,redirect_stdout(io.StringIO()):
                self.assertEqual(main(['--write-profile',str(file),'--no-menu']),0)
                run.assert_not_called()
            original=file.read_bytes()
            with redirect_stdout(io.StringIO()),patch('sys.stderr',new=io.StringIO()):
                self.assertEqual(main(['--write-profile',str(file),'--hits','7']),1)
            self.assertEqual(original,file.read_bytes())

    def test_contextual_fields_do_not_offer_unsupported_flags(self):
        a=self.args()
        names=lambda context,mode='',strategy='':{p.name for p in relevant(context,a,mode,strategy=strategy)}
        self.assertIn('mmseqs_sensitivity',names('search','mmseqs-blastp'))
        self.assertNotIn('sensitivity',names('search','mmseqs-blastp'))
        self.assertNotIn('blast_gapopen',names('search','tblastx'))
        self.assertNotIn('free_ends',names('water','blastp'))
        a.aligner='muscle'
        self.assertIn('mafft_op',names('msa',strategy='fragments'))
        self.assertNotIn('muscle_perm',names('msa',strategy='fragments'))


class SessionParametersTests(unittest.TestCase):
    def session(self,*extra,paired=False):
        paths=[str(ROOT/'examples/example_query.fna'),str(ROOT/'examples/example_refs.faa'),'--mode','blastx'] if paired else [str(ROOT/'examples/example_msa.fasta')]
        return Session(parse_args([*paths,'--no-menu','--module-policy','never',*extra]),UI(False))

    def test_search_preview_matches_execution_and_engine_changes_are_stale(self):
        session=self.session('--blast-word-size','3',paired=True)
        try:
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with redirect_stdout(io.StringIO()):
                preview=session.search_preview();self.assertEqual(session.tools.history,[])
                session.search()
            for row in session.tools.history:self.assertIn(row['command'],preview)
            self.assertEqual(session.search_parameters['blast_word_size'],3)
            session.a.identity=70
            self.assertFalse(WorkspaceState.from_session(session).search_stale)
            session.a.blast_masking='no'
            self.assertTrue(WorkspaceState.from_session(session).search_stale)
            self.assertEqual(session.calculation_runs['search']['parameters']['blast_masking'],'native')
        finally:session.close()

    def test_alignment_preview_is_not_a_calculation_and_provenance_survives_resume(self):
        session=self.session('--mafft','local','--mafft-op','2.2','--mafft-maxiterate','7')
        try:
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with redirect_stdout(io.StringIO()):preview=session.msa_run('all',preview=True)
            self.assertEqual(session.tools.history,[]);self.assertIsNone(session.msa_records)
            self.assertEqual(session.calculation_runs,{})
            with redirect_stdout(io.StringIO()):session.msa_run('all')
            self.assertIn(session.tools.history[0]['command'],preview)
            self.assertEqual(session.calculation_runs['msa']['parameters']['mafft_op'],2.2)
            session.a.mafft_op=3.1;session.save_settings()
            with tempfile.TemporaryDirectory() as tmp,redirect_stdout(io.StringIO()):
                project=save_project(session,Path(tmp)/'project')
                resumed=resume_project(project,parse_args(['--resume',str(project),'--no-menu']),UI(False))
                try:
                    self.assertEqual(resumed.a.mafft_op,3.1)
                    self.assertEqual(resumed.calculation_runs['msa']['parameters']['mafft_op'],2.2)
                    self.assertEqual(len(resumed.calculation_runs['msa']['commands']),1)
                finally:resumed.close()
        finally:session.close()

    def test_fragment_plan_and_execution_use_mafft_in_both_stages_with_tuned_reference_alignment(self):
        session=self.session('--aligner','muscle','--mafft','local','--mafft-maxiterate','3','--mafft-op','2.1',paired=True)
        try:
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with redirect_stdout(io.StringIO()):
                preview=session.msa_run('fragments',preview=True)
                session.msa_run('fragments')
            self.assertEqual(len(session.tools.history),2)
            for row in session.tools.history:self.assertIn(row['command'],preview)
            self.assertIn('--localpair',session.tools.history[0]['command'])
            self.assertIn('--maxiterate 3',session.tools.history[0]['command'])
            self.assertIn('--addfragments',session.tools.history[1]['command'])
            self.assertNotIn('--maxiterate',session.tools.history[1]['command'])
            self.assertEqual(session.calculation_runs['msa']['parameters']['aligner'],'mafft')
        finally:session.close()

    def test_hsp_overlap_and_strategy_actually_change_selected_alignment_input(self):
        session=self.session('--hsp-overlap','.8','--mafft','local',paired=True)
        try:
            anchor=HSP('q','ref_HSP_1',95,20,1,60,1,20,1e-20,200,1,0,'M'*20,'M'*20)
            other=HSP('q','ref_HSP_2',95,20,25,84,1,20,1e-20,100,1,0,'M'*20,'M'*20)
            session.hsps=[anchor,other];session.apply_filters()
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with redirect_stdout(io.StringIO()):session.msa_run('hsp')
            self.assertEqual(len(session.msa_records),2)
            self.assertIn('--localpair',session.tools.history[-1]['command'])
            session.a.hsp_overlap=.5
            with redirect_stdout(io.StringIO()):session.msa_run('hsp')
            self.assertEqual(len(session.msa_records),3)
        finally:session.close()

    def test_slurm_preserves_new_native_controls(self):
        session=self.session('--mmseqs-sensitivity','6.5','--mafft','local','--mafft-op','2.2',
            '--fasttree-model','LG','--fasttree-gamma','--needle-endopen','4')
        try:
            with patch('alignview.show_file'),redirect_stdout(io.StringIO()):session.slurm_script()
            command=shlex.split((session.work/'alignview_job.sbatch').read_text().splitlines()[-1])
            for flag,value in (('--mmseqs-sensitivity','6.5'),('--mafft-op','2.2'),('--fasttree-model','LG'),('--needle-endopen','4.0')):
                self.assertEqual(command[command.index(flag)+1],value)
            self.assertIn('--fasttree-gamma',command)
        finally:session.close()

    def test_failed_alignment_preserves_active_result_and_previous_logs(self):
        session=self.session()
        try:
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with redirect_stdout(io.StringIO()):session.msa_run('all');session.msa_run('all')
            logs=[row['stderr_log'] for row in session.tools.history]
            self.assertEqual(len(set(logs)),2)
            self.assertTrue(all((session.work/log).is_file() for log in logs))
            original=session.msa_file.read_bytes();run=session.calculation_runs['msa']
            def fail(input,output,**kwargs):
                output.write_text('>broken\nAAA\n>other\nAA\n')
                return output
            with patch.object(session.tools,'msa',side_effect=fail),redirect_stdout(io.StringIO()):
                with self.assertRaises(CommandError):session.msa_run('all')
            self.assertEqual(session.msa_file.read_bytes(),original)
            self.assertEqual(session.calculation_runs['msa'],run)
            self.assertEqual(session.a.msa,'all')
        finally:session.close()

    def test_profile_on_resume_changes_future_settings_and_preserves_executed_parameters(self):
        session=self.session('--mafft-op','2.2')
        try:
            session.tools.env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+os.environ['PATH']
            with tempfile.TemporaryDirectory() as tmp,redirect_stdout(io.StringIO()):
                session.msa_run('all');project=save_project(session,Path(tmp)/'project')
                filename=Path(tmp)/'profile.json'
                write_profile(parse_args(['--mafft-op','4.1','--hits','9']),filename)
                args=parse_args(['--resume',str(project),'--load-profile',str(filename),'--mafft-op','3.2','--hits','100','--no-menu'])
                resumed=resume_project(project,args,UI(False))
                try:
                    self.assertEqual(resumed.a.mafft_op,3.2);self.assertEqual(resumed.a.hits,100)
                    self.assertEqual(resumed.calculation_runs['msa']['parameters']['mafft_op'],2.2)
                    self.assertIsNotNone(resumed.msa_records);self.assertTrue(resumed.dirty)
                finally:resumed.close()
        finally:session.close()

    def test_terminal_reaches_actual_aligner_parameters_before_launch(self):
        code="import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(['--demo','--no-color']))"
        terminal=Terminal(code)
        try:
            terminal.expect('AlignView / Alignment');terminal.send('q')
            terminal.expect('Results ready');terminal.send('3\r')
            terminal.expect('Create an alignment');terminal.send('2\r')
            terminal.expect('Advanced alignment options');terminal.send('3\r')
            terminal.expect('Calculation parameters · MAFFT');terminal.expect('MAFFT strategy')
            terminal.send('2\r');terminal.expect('L-INS-i');terminal.send('3\r')
            terminal.expect('applies to the next calculation')
            terminal.send('\x1bOQ');terminal.expect('Back to Home')
            terminal.send('q');terminal.expect('Quitter AlignView ?');terminal.send('2\r')
            terminal.expect('without saving');self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()


if __name__=='__main__':unittest.main()
