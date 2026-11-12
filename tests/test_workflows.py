"""Storage, translation codes, VSEARCH and language integration regressions."""
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,parse_args
from av_core import Seq,translate_nt,orfs,write_fasta,read_fasta
from av_genetics import TABLES,code_table,context_dependent

from av_parameters import snapshot
from av_project import save_project,resume_project
from av_tools import Toolchain,CommandError
from av_ui import UI
from av_vsearch import commands,summarize,publish_directory,FIELDS
from av_navigation import Workbench
from av_tree import tree_svg
from av_setup import probe_version
from test_release import Terminal


class WorkflowTests(unittest.TestCase):

    def session(self,*extra):
        return Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--msa','view','--no-menu',*extra]),UI(False))

    def cli(self,*args,cwd):
        return subprocess.run([sys.executable,str(ROOT/'alignview.py'),*map(str,args)],cwd=cwd,
                              capture_output=True,text=True,timeout=30)

    def test_standalone_html_uses_launch_directory_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            write_fasta([Seq('a','a','ACGT'),Seq('b','b','ACGA')],tmp/'input.fasta')
            command=[str(ROOT/'alignview.sh'),'input.fasta','--msa','view','--html-out','reports/view.html','--no-menu']
            result=subprocess.run(command,cwd=tmp,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            report=tmp/'reports/view.html';original=report.read_bytes()
            self.assertIn(b'<!doctype html>',original)
            result=subprocess.run(command,cwd=tmp,capture_output=True,text=True,timeout=30)
            self.assertEqual(result.returncode,1)
            self.assertEqual(report.read_bytes(),original)
            self.assertFalse((tmp/'session.json').exists())

    def test_resume_rebuilds_report_and_honors_explicit_export_profile(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp);session=self.session()
            try:
                session.load_msa();project=save_project(session,tmp/'project')
            finally:session.close()
            result=self.cli('--resume',project,'--report','--export-profile','publication',
                            '--save','--out','updated','--no-menu',cwd=tmp)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertTrue((tmp/'updated/report.html').is_file())
            self.assertTrue((tmp/'updated/exports/msa_alignment.nex').is_file())
            previous=json.loads((project/'session.json').read_text())
            current=json.loads((tmp/'updated/session.json').read_text())
            self.assertEqual(previous['commands'],current['commands'])
            self.assertEqual(previous['msa'],current['msa'])

    def test_resume_explicit_settings_do_not_mutate_historical_results(self):
        session=self.session('--genetic-code','2')
        try:
            session.calculation_runs['translation']={'parameters':{'genetic_code':2}}
            with tempfile.TemporaryDirectory() as tmp:
                project=save_project(session,Path(tmp)/'project')
                resumed=resume_project(project,parse_args(['--resume',str(project),'--genetic-code','4','--no-menu']),UI(False))
                try:
                    self.assertEqual(resumed.a.genetic_code,4)
                    self.assertEqual(resumed.calculation_runs['translation']['parameters']['genetic_code'],2)
                    with patch.object(resumed.tools,'run') as run:resumed.resume_actions();run.assert_not_called()
                finally:resumed.close()
        finally:session.close()

    def test_legacy_preset_profile_loads_with_english_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            profile=Path(tmp)/'legacy.json'
            profile.write_text(json.dumps({'schema':1,'parameters':{'preset':'sensible','hits':42}}))
            args=parse_args(['--demo','--load-profile',str(profile),'--no-menu'])
            self.assertEqual(args.preset,'sensitive')
            self.assertEqual(args.hits,42)

    def test_english_is_default_and_has_no_language_selector(self):
        from alignview import build_parser
        parser=build_parser()
        self.assertNotIn('--lang',parser._option_string_actions)
        with tempfile.TemporaryDirectory() as tmp,patch.dict(os.environ,{'ALIGNVIEW_LANG':'fr'}):
            result=self.cli('--help',cwd=tmp)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('Input and resuming',result.stdout)
            self.assertIn('Command-line',self.cli('--docs','cli',cwd=tmp).stdout)
            self.assertNotIn('Langue',result.stdout)

    def test_english_output_preserves_identifiers(self):
        record=Seq('Δ_gene_42','Δ_gene_42','ACGT')
        self.assertEqual(record.id,'Δ_gene_42')
        svg=tree_svg("('Δ_gene_42':0.1,b:0.2);")
        self.assertIn('Exploratory tree',svg)
        self.assertIn('Δ_gene_42',svg)
        self.assertIn('substitutions per site',svg)

    def test_english_help_and_guides_are_available_without_data(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=self.cli('--help',cwd=tmp)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('Input and resuming',result.stdout)
            result=self.cli('--docs','quickstart','--no-menu',cwd=tmp)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertTrue(result.stdout.startswith('# Getting started'))
            self.assertIn('--html-out',result.stdout)
            result=self.cli('q.fa','r.fa','--msa','view','--no-menu',cwd=tmp)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('requires a single aligned FASTA',result.stderr)

    def test_guided_english_query_selector_reaches_filtered_identifier(self):
        code="from alignview import Session,parse_args;from av_ui import UI;import sys;a=parse_args([sys.argv[1]+'/examples/example_msa.fasta',sys.argv[1]+'/examples/example_msa.fasta',]);a._guided=True;s=Session(a,UI(False));print('ACTIVE='+s.query.id,flush=True);s.close()"
        terminal=Terminal(code)
        try:
            terminal.expect('Choose the active query');terminal.send('/seq_3\r')
            terminal.expect('ACTIVE=seq_3');self.assertEqual(terminal.process.wait(timeout=5),0)
        finally:terminal.close()


class GeneticCodeTests(unittest.TestCase):

    def test_alternative_codes_and_ambiguity(self):
        self.assertEqual(translate_nt('ATATGAAGA',table=1),'I*R')
        self.assertEqual(translate_nt('ATATGAAGA',table=2),'MW*')
        self.assertEqual(translate_nt('ATATGAAGA',table=4),'IWR')
        self.assertEqual(translate_nt('GCN',table=1),'A')
        self.assertEqual(translate_nt('UGA',table=2),'W')
        self.assertEqual(translate_nt('TCATAT',frame=-1,table=2),'MW')
        for identifier in TABLES:self.assertEqual(len(code_table(identifier)),64)
        with self.assertRaises(ValueError):translate_nt('ATG',table=7)

    def test_orf_starts_and_stop_codons_follow_table(self):
        record=Seq('x','x','GTGAAATGA')
        self.assertEqual(orfs(record,min_aa=1,table=11,starts='atg'),[])
        found=orfs(record,min_aa=1,table=11,starts='table')
        self.assertEqual(found[0]['protein'],'MK')
        self.assertEqual(orfs(record,min_aa=1,table=4,starts='table'),[])
        for identifier in (27,28,31):
            self.assertTrue(context_dependent(identifier))
            with self.assertRaises(ValueError):orfs(record,min_aa=1,table=identifier)

    def test_native_search_commands_receive_query_and_reference_codes(self):
        options=parse_args(['--genetic-code','2','--reference-genetic-code','4'])
        with tempfile.TemporaryDirectory() as tmp:
            tool=Toolchain(Path(tmp),auto_modules=False)
            command=tool.search_commands('tblastx',Path('q'),Path('r'),options=options)[-1].argv
            self.assertEqual(command[command.index('-query_gencode')+1],'2')
            self.assertEqual(command[command.index('-db_gencode')+1],'4')
            diamond=tool.search_commands('diamond-blastx',Path('q'),Path('r'),options=options)[-1].argv
            self.assertEqual(diamond[diamond.index('--query-gencode')+1],'2')
            self.assertEqual(snapshot('search',options,'tblastx')['genetic_code'],2)

    def test_internal_translation_keeps_its_own_provenance(self):
        session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--genetic-code','2','--no-menu']),UI(False))
        try:
            with contextlib.redirect_stdout(io.StringIO()):session.translate_frames()
            session.a.genetic_code=4
            self.assertEqual(session.calculation_runs['translation']['parameters']['genetic_code'],2)
            self.assertIn('table=2',(session.work/'six_frames.faa').read_text())
        finally:session.close()


class VsearchTests(unittest.TestCase):
    def options(self,*extra):return parse_args(['--vsearch','cluster',*extra])

    def test_commands_have_native_metrics_and_explicit_heuristics(self):
        options=self.options('--vsearch-id','.91','--hits','23','--vsearch-maxrejects','0','--vsearch-minwordmatches','0')
        search=commands('search',Path('q'),Path('r'),Path('out'),options)[0].argv
        self.assertIn('--usearch_global',search)
        self.assertEqual(search[search.index('--maxhits')+1],'23')
        self.assertEqual(search[search.index('--id')+1],'0.91')
        fields=search[search.index('--userfields')+1]
        self.assertNotIn('evalue',fields);self.assertNotIn('bits',fields)
        self.assertIn('qrow',fields)
        self.assertEqual(len(commands('chimeras',Path('q'),None,Path('out'),options)),2)
        with self.assertRaises(ValueError):commands('search',Path('q'),None,Path('out'),options)
        with self.assertRaises(ValueError):commands('cluster',Path('q'),None,Path('out'),self.options('--vsearch-min-length','20','--vsearch-max-length','10'))

    def test_output_table_retains_native_fields_and_validates_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            values=['q','ref','99.0','100','1','100','1','100','100','100','+','A'*100,'A'*100]
            (folder/'matches_raw.tsv').write_text('\t'.join(values)+'\n')
            summary,preview=summarize('search',folder,1,1)
            self.assertEqual(summary['matches'],1)
            self.assertEqual(preview[0]['id'],99.0)
            self.assertEqual((folder/'matches.tsv').read_text().splitlines()[0],'\t'.join(FIELDS))
            values[2]='nan';(folder/'matches_raw.tsv').write_text('\t'.join(values)+'\n')
            with self.assertRaises(ValueError):summarize('search',folder,1,1)

    def test_failed_vsearch_run_keeps_previous_successful_results(self):
        session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--no-menu']),UI(False))
        try:
            previous=session.work/'vsearch/cluster';previous.mkdir(parents=True)
            (previous/'centroids.fasta').write_text('original')
            with patch.object(session.tools,'execute',side_effect=CommandError('failed')),self.assertRaises(CommandError):
                session.vsearch_run('cluster')
            self.assertEqual((previous/'centroids.fasta').read_text(),'original')
        finally:session.close()

    def test_protein_input_is_rejected_before_engine_execution(self):
        session=Session(parse_args([str(ROOT/'examples/example_refs.faa'),'--no-menu']),UI(False))
        try:
            with patch.object(session.tools,'execute') as execute,self.assertRaises(ValueError):session.vsearch_run('cluster')
            execute.assert_not_called()
        finally:session.close()

    def test_engine_version_is_not_mistaken_for_the_citation(self):
        response=subprocess.CompletedProcess([],0,'Rognes T et al. (2016)\n','vsearch v2.32.0_linux_x86_64\n')
        with patch('av_setup.subprocess.run',return_value=response):
            self.assertEqual(probe_version('/tools/vsearch',{}),'vsearch v2.32.0_linux_x86_64')

    def test_failed_directory_publication_restores_previous_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            previous=Path(tmp)/'cluster';previous.mkdir();(previous/'result').write_text('previous')
            with self.assertRaises(FileNotFoundError):publish_directory(Path(tmp)/'missing',previous)
            self.assertEqual((previous/'result').read_text(),'previous')

    def test_home_opens_explore_when_vsearch_results_are_ready(self):
        session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--no-menu']),UI(False))
        try:
            session.vsearch_results={'cluster':{'summary':{'representatives':2}}}
            workbench=Workbench(session)
            with patch.object(workbench,'choose',return_value=None) as choose,patch.object(workbench,'confirm_exit',return_value=True):
                workbench.run()
            self.assertEqual(choose.call_args.kwargs['default'],'1')
        finally:session.close()


if __name__=='__main__':unittest.main()
