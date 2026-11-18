"""Task-command acceptance through the public executable and real result files.

Scientific executables are isolated, deterministic fixtures. These checks exercise
command routing, saved evidence and output safety, not engine accuracy.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]


class TaskCommandTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='alignview_cli_v1_')
        self.addCleanup(self.temporary.cleanup)
        self.folder=Path(self.temporary.name)
        self.launch=self.folder/'launch';self.launch.mkdir()
        self.inputs=self.folder/'inputs';self.inputs.mkdir()
        self.runtime=self.folder/'runtime';self.runtime.mkdir()
        (self.runtime/'python3').symlink_to(sys.executable)
        self.env={**os.environ,'PATH':os.pathsep.join([str(ROOT/'tests/mockbin'),str(ROOT/'tests/mockextra'),str(self.runtime)]),
                  'ALIGNVIEW_TOOLS_PREFIX':str(self.folder/'no-managed-tools'),'ALIGNVIEW_AUTO_MODULES':'0','NO_COLOR':'1'}
        self.aligned=self.fasta('aligned.fasta',[('alpha','ACGTACGT'),('beta','ACGTACGA'),('gamma','ACGTTCGT')])
        self.raw=self.fasta('unaligned.fasta',[('alpha','ACGTACGTACGT'),('beta','ACGTACGA'),('gamma','ACGTTC')])
        self.dna=self.fasta('query.fna',[('query','ACGT'*50)])
        self.dna_refs=self.fasta('references.fna',[('ref_A','ACGT'*45),('ref_B','ACGT'*40)])
        self.protein=self.fasta('query.faa',[('protein_query','MPLDIGDAVKRHITATDGKVLWNIA')])
        self.protein_refs=self.fasta('references.faa',[('protein_A','MPLDIGDAVKRHITATDGKVLWNIA'),('protein_B','MPLDIADAVKRHITATDGRVLWKIA')])

    def fasta(self,name,records):
        path=self.inputs/name
        path.write_text(''.join('>'+identifier+'\n'+sequence+'\n' for identifier,sequence in records),encoding='utf-8')
        return path

    def cli(self,*arguments,env=None):
        return subprocess.run([sys.executable,str(ROOT/'alignview.py'),*map(str,arguments)],cwd=self.launch,
            env=self.env if env is None else env,capture_output=True,text=True,timeout=40)

    def succeeded(self,*arguments):
        run=self.cli(*arguments)
        self.assertEqual(run.returncode,0,run.stdout+'\n'+run.stderr)
        return run

    def project(self,name):
        return json.loads((self.launch/name/'session.json').read_text(encoding='utf-8'))

    def engines(self,project):
        return [Path((row.get('argv') or shlex.split(row['command']))[0]).name.lower() for row in project['commands']]

    def blocked_engine(self):
        marker=self.folder/'engine_was_run'
        executable=self.runtime/'blocked-engine'
        executable.write_text('#!'+sys.executable+'\nfrom pathlib import Path\nimport sys\nPath('+repr(str(marker))+').write_text("executed")\nraise SystemExit(91)\n')
        executable.chmod(0o755)
        arguments=[]
        for engine in ['mafft','FastTree','makeblastdb','blastn','blastx']:
            arguments+=['--executable',engine+'='+str(executable)]
        return marker,arguments

    def saved_alignment(self,*,tree=False,report=False):
        arguments=['align',self.raw,'--no-menu','--save','--out','project','--module-policy','never']
        if tree:arguments.append('--tree')
        if report:arguments.append('--report')
        self.succeeded(*arguments)
        return self.launch/'project'

    def test_help_separates_daily_tasks_from_native_parameters(self):
        concise=self.succeeded('--help').stdout
        complete=self.succeeded('--help-all').stdout
        for task in ['view FILE','align FILE','compare QUERY','report FILE_OR_PROJECT','resume PROJECT','tools']:
            self.assertIn(task,concise)
        self.assertNotIn('--blast-word-size',concise)
        self.assertIn('--blast-word-size',complete)
        self.assertIn('--mafft-maxiterate',complete)
        self.assertLess(len(concise),len(complete)//2)
        comparison=self.succeeded('compare','--help').stdout
        self.assertIn('--query-id',comparison)
        self.assertIn('--mode',comparison)
        self.assertIn('--hits',comparison)

    def test_view_reads_alignment_without_running_an_engine(self):
        original=self.aligned.read_bytes()
        self.succeeded('view',self.aligned,'--no-menu','--save','--out','view-project')
        data=self.project('view-project')
        self.assertEqual(data['mode'],'msa')
        self.assertEqual(self.engines(data),[])
        self.assertTrue(data['msa'])
        self.assertEqual(self.aligned.read_bytes(),original)

    def test_align_single_file_retains_input_records_and_runs_alignment(self):
        self.succeeded('align',self.raw,'--no-menu','--save','--out','alignment','--module-policy','never')
        data=self.project('alignment')
        self.assertEqual(self.engines(data),['mafft'])
        from av_core import read_fasta
        rows=read_fasta(self.launch/'alignment'/data['msa'])
        self.assertEqual([row.id for row in rows],['alpha','beta','gamma'])
        self.assertEqual(len({len(row.seq) for row in rows}),1)
        self.assertEqual([row.seq.replace('-','') for row in rows],['ACGTACGTACGT','ACGTACGA','ACGTTC'])

    def test_two_file_align_never_runs_a_similarity_search(self):
        self.succeeded('align',self.dna,self.dna_refs,'--no-menu','--save','--out','two-files','--module-policy','never')
        data=self.project('two-files')
        self.assertEqual(data['mode'],'global-dna')
        self.assertEqual(self.engines(data),['mafft'])
        self.assertFalse(data['search_completed'])
        self.assertEqual(data['hsps'],[])
        from av_core import read_fasta
        self.assertEqual([row.id for row in read_fasta(self.launch/'two-files'/data['msa'])],['query','ref_A','ref_B'])

    def test_compare_infers_native_and_translated_methods_from_inputs(self):
        cases=[(self.dna,self.dna_refs,'blastn'),(self.dna,self.protein_refs,'blastx'),
               (self.protein,self.protein_refs,'blastp'),(self.protein,self.dna_refs,'tblastn')]
        for query,references,mode in cases:
            with self.subTest(mode=mode):
                self.succeeded('compare',query,references,'--no-menu','--save','--out',mode,'--module-policy','never')
                data=self.project(mode)
                self.assertEqual(data['mode'],mode)
                self.assertTrue(data['search_completed'])
                self.assertTrue(data['hsps'])
                self.assertIn(mode,self.engines(data))
                self.assertFalse(data['msa'])

    def test_batch_comparison_requires_then_honors_explicit_query(self):
        queries=self.fasta('multiple.fna',[('first','ACGT'*50),('chosen','ACGA'*50)])
        rejected=self.cli('compare',queries,self.dna_refs,'--no-menu','--save','--out','rejected','--module-policy','never')
        self.assertNotEqual(rejected.returncode,0)
        self.assertIn('--query-id',rejected.stderr)
        self.assertFalse((self.launch/'rejected').exists())
        self.succeeded('compare',queries,self.dna_refs,'--query-id','chosen','--no-menu','--save','--out','chosen','--module-policy','never')
        data=self.project('chosen')
        self.assertEqual(data['active_query_id'],'chosen')
        self.assertTrue(all(row['query']=='chosen' for row in data['hsps']))

    def test_batch_comparison_requires_reference_input_before_creating_a_project(self):
        original=self.dna.read_bytes()
        marker,overrides=self.blocked_engine()
        rejected=self.cli('compare',self.dna,'--no-menu','--save','--out','missing-references',*overrides)
        self.assertNotEqual(rejected.returncode,0,rejected.stdout+'\n'+rejected.stderr)
        self.assertIn('reference',rejected.stderr.lower())
        self.assertFalse((self.launch/'missing-references').exists())
        self.assertFalse(marker.exists())
        self.assertEqual(self.dna.read_bytes(),original)

    def test_legacy_positional_flags_keep_working(self):
        self.succeeded(self.aligned,'--msa','view','--html-out','legacy.html','--no-menu')
        self.assertTrue((self.launch/'legacy.html').is_file())
        self.succeeded(self.dna,self.dna_refs,'--mode','blastn','--no-menu','--save','--out','legacy-comparison','--module-policy','never')
        self.assertEqual(self.project('legacy-comparison')['mode'],'blastn')

    def test_report_default_is_durable_and_relative_to_launch_directory(self):
        original=self.aligned.read_bytes()
        self.succeeded('report',self.aligned)
        target=self.launch/'aligned_alignview.html'
        self.assertTrue(target.is_file())
        self.assertIn('<html lang="en">',target.read_text())
        self.assertIn('id="reportData"',target.read_text())
        self.assertEqual(self.aligned.read_bytes(),original)
        self.assertFalse((self.inputs/'aligned_alignview.html').exists())
        self.assertFalse((self.launch/'session.json').exists())
        previous=target.read_bytes()
        duplicate=self.cli('report',self.aligned)
        self.assertNotEqual(duplicate.returncode,0)
        self.assertEqual(target.read_bytes(),previous)

    def test_shell_launcher_keeps_caller_relative_inputs_and_outputs(self):
        environment={**self.env,'PATH':self.env['PATH']+os.pathsep+'/usr/bin'+os.pathsep+'/bin'}
        run=subprocess.run(['/bin/bash',str(ROOT/'alignview.sh'),'report','../inputs/aligned.fasta'],
            cwd=self.launch,env=environment,capture_output=True,text=True,timeout=40)
        self.assertEqual(run.returncode,0,run.stdout+'\n'+run.stderr)
        self.assertTrue((self.launch/'aligned_alignview.html').is_file())
        self.assertFalse((self.inputs/'aligned_alignview.html').exists())

    def test_report_saved_project_uses_existing_results_without_engines(self):
        project=self.saved_alignment(tree=True)
        original=(project/'session.json').read_bytes();data=self.project('project')
        marker,overrides=self.blocked_engine()
        self.succeeded('report',project,'--html-out','reopened.html','--module-policy','never',*overrides)
        report=self.launch/'reopened.html'
        payload=json.loads(re.search(r'<script id="reportData" type="application/json">(.*?)</script>',report.read_text(),re.S)[1])
        self.assertEqual(payload['history'],data['commands'])
        self.assertIn('45_fasttree_tree',report.read_text())
        self.assertFalse(marker.exists())
        self.assertEqual((project/'session.json').read_bytes(),original)

    def test_report_rejects_two_inputs_for_files_and_saved_projects(self):
        project=self.saved_alignment()
        original=(project/'session.json').read_bytes()
        marker,overrides=self.blocked_engine()
        for index,(first,second) in enumerate([(self.aligned,self.dna_refs),(self.aligned,project),
                                                (project,self.dna_refs),(project,project)]):
            with self.subTest(first=first.name,second=second.name):
                output=self.launch/f'rejected-report-{index}.html'
                rejected=self.cli('report',first,second,'--html-out',output,*overrides)
                self.assertNotEqual(rejected.returncode,0,rejected.stdout+'\n'+rejected.stderr)
                self.assertFalse(output.exists())
        self.assertFalse(marker.exists())
        self.assertEqual((project/'session.json').read_bytes(),original)

    def test_resume_without_actions_does_not_rerun_or_rewrite_project(self):
        project=self.saved_alignment(tree=True,report=True)
        before=(project/'session.json').read_bytes()
        marker,overrides=self.blocked_engine()
        run=self.succeeded('resume',project,'--no-menu','--module-policy','never',*overrides)
        self.assertIn('Project retained:',run.stdout)
        self.assertNotIn('Report not retained',run.stdout+run.stderr)
        self.assertNotIn('temporary results removed',run.stdout)
        self.assertFalse(marker.exists())
        self.assertEqual((project/'session.json').read_bytes(),before)

    def test_tools_needs_no_dataset_and_json_stays_machine_readable(self):
        run=self.succeeded('tools','--diagnose-json','--module-policy','never')
        rows=json.loads(run.stdout)
        self.assertIsInstance(rows,list)
        self.assertTrue(any(row['command']=='mafft' and row['available'] for row in rows))
        self.assertTrue(all({'command','available','path','version'}<=set(row) for row in rows))
        self.assertFalse((self.launch/'session.json').exists())

    def test_tools_rejects_an_input_instead_of_silently_ignoring_it(self):
        original=self.dna.read_bytes()
        rejected=self.cli('tools',self.dna,'--diagnose-json','--module-policy','never')
        self.assertNotEqual(rejected.returncode,0,rejected.stdout+'\n'+rejected.stderr)
        self.assertFalse((self.launch/'session.json').exists())
        self.assertEqual(self.dna.read_bytes(),original)

    def test_output_paths_are_data_and_existing_symlinks_are_not_overwritten(self):
        # This filename would execute a command if routed through a shell.
        destination='report $(touch SHELL_EXECUTED).html'
        self.succeeded('report',self.aligned,'--html-out',destination)
        self.assertTrue((self.launch/destination).is_file())
        self.assertFalse((self.launch/'SHELL_EXECUTED').exists())
        original=self.aligned.read_bytes()
        link=self.launch/'input-link.html';link.symlink_to(self.aligned)
        rejected=self.cli('report',self.aligned,'--html-out',link)
        self.assertNotEqual(rejected.returncode,0)
        self.assertTrue(link.is_symlink())
        self.assertEqual(self.aligned.read_bytes(),original)
        directory=self.launch/'occupied';directory.mkdir();sentinel=directory/'keep';sentinel.write_text('original')
        rejected=self.cli('view',self.aligned,'--no-menu','--save','--out',directory)
        self.assertNotEqual(rejected.returncode,0)
        self.assertEqual(sentinel.read_text(),'original')
        self.assertEqual(list(directory.iterdir()),[sentinel])


if __name__=='__main__':unittest.main()
