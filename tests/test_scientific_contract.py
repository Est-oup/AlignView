"""Scientific evidence and complete-snapshot integrity regressions for v1."""
from contextlib import redirect_stdout
from dataclasses import replace
import io
import json
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import (HSP,Seq,coverage_data,select_hits,normalize_search_record,parse_blast_tab,
                     validate_search_results,validate_alignment_results,write_fasta,read_fasta,sequence_stats)
from av_units import coord_audit
from av_project import read_project,save_project,resume_project,INTEGRITY_FILE
from alignview import Session,parse_args
from av_ui import UI
from av_tools import Toolchain,CommandError
from av_tree import parse_newick,relabel_newick,tree_svg


class EvidenceTests(unittest.TestCase):
    def hsp(self,**changes):
        row=HSP('q','ref',80,4,1,4,1,4,1e-5,20,0,0,'ACGT','ACGA')
        return replace(row,**changes)

    def test_reference_ranking_does_not_depend_on_engine_order(self):
        weak=self.hsp(identity=100,bitscore=10,evalue=1e-3)
        strong=self.hsp(identity=80,bitscore=30,evalue=1e-8,qstart=3,qend=6,sstart=3,send=6)
        other=self.hsp(subject='other',bitscore=25)
        selected=select_hits([weak,other,strong],10,max_hits=1)
        self.assertEqual(list(selected),['ref'])
        self.assertEqual(selected['ref'][0],strong)
        row=coverage_data([weak,other,strong],10,{'ref':8,'other':4})[0]
        self.assertEqual(row['subject'],'ref')
        self.assertEqual(row['best_identity'],80)
        self.assertEqual(row['max_identity'],100)
        self.assertEqual(row['best_bitscore'],30)
        self.assertEqual(row['query_aligned_length'],6)
        self.assertEqual(row['query_coverage'],60)

    def test_malformed_metrics_and_rows_cannot_enter_project_evidence(self):
        for changes in ({'identity':float('nan')},{'identity':101},{'evalue':-1},
                        {'bitscore':float('inf')},{'qstart':0},{'align_length':0},
                        {'qseq':'ACG'},{'sseq':''},{'qframe':4}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):self.hsp(**changes)
        row=self.hsp(qseq='a.GT',sseq='A-TA')
        self.assertEqual((row.matches,row.mismatches,row.gap_columns),(1,2,1))

    def test_search_validation_checks_identity_bounds_and_native_units(self):
        query=Seq('q','q','ACGTACGT');refs=[Seq('ref','ref','ACGTACGT')]
        self.assertEqual(validate_search_results([],query,refs,'blastn'),[])
        self.assertEqual(validate_search_results([self.hsp()],query,refs,'blastn'),[self.hsp()])
        for changes in ({'query':'foreign'},{'subject':'unknown'},
                        {'qstart':7,'qend':10},{'qend':5},
                        {'qseq':'','sseq':''}):
            with self.subTest(changes=changes),self.assertRaises(ValueError):
                validate_search_results([self.hsp(**changes)],query,refs,'blastn')
        translated=self.hsp(qend=12,qframe=1)
        self.assertEqual(validate_search_results([translated],Seq('q','q','ATG'*4),refs,'blastx'),[translated])
        with self.assertRaises(ValueError):validate_search_results([replace(translated,qframe=0)],Seq('q','q','ATG'*4),refs,'blastx')
        with self.assertRaisesRegex(ValueError,'reading frame disagrees'):
            validate_search_results([replace(translated,qframe=2)],Seq('q','q','ATG'*4),refs,'blastx')

    def test_untranslated_audit_does_not_tolerate_two_extra_positions(self):
        row=coord_audit([self.hsp(qend=6)],'blastn')[0]
        self.assertFalse(row['query_consistent'])
        self.assertTrue(row['reference_consistent'])

    def test_fractional_engine_frame_is_not_silently_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'bad.tsv'
            target.write_text('q\tref\t100\t4\t1\t4\t1\t4\t1e-5\t20\t100\t1.5\t0\tACGT\tACGT\n')
            with self.assertRaisesRegex(ValueError,'frames must be integers'):parse_blast_tab(target)

    def test_normalization_preserves_unknown_positions_and_identifiers(self):
        record=Seq('rna','rna original description','ac-u.?g')
        normalized=normalize_search_record(record)
        self.assertEqual(normalized.seq,'ACTNG')
        self.assertEqual(normalized.id,record.id)
        self.assertEqual(record.seq,'ac-u.?g')
        with self.assertRaises(ValueError):normalize_search_record(Seq('empty','empty','---'))
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'records.fa'
            write_fasta([Seq('stable_id','a description without the identifier','ACGT')],target)
            self.assertEqual(read_fasta(target)[0].id,'stable_id')
        self.assertEqual(sequence_stats(Seq('lower','lower','gcgc'))['gc_pct'],100)

    def test_alignment_validation_restores_order_and_preserves_residues(self):
        expected=[Seq('q','q','ACG'),Seq('ref','ref','ATCG')]
        aligned=[Seq('ref','ref','ATCG'),Seq('q','q','A-CG')]
        self.assertEqual([r.id for r in validate_alignment_results(aligned,expected)],['q','ref'])
        for invalid in ([aligned[0]],aligned+[aligned[0]],
                        [aligned[0],Seq('q','q','A-TG')],
                        [aligned[0],Seq('q','q','ACG')]):
            with self.subTest(invalid=invalid),self.assertRaises(ValueError):validate_alignment_results(invalid,expected)
        rna=[Seq('rna','rna','AU')];dna=[Seq('rna','rna','AT')]
        self.assertEqual(validate_alignment_results(dna,rna,nucleotide=True),dna)
        with self.assertRaises(ValueError):validate_alignment_results(dna,rna)


class CompleteProjectTests(unittest.TestCase):
    def session(self):
        result=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--msa','view','--no-menu']),UI(False))
        result.load_msa();self.addCleanup(result.close)
        return result

    def test_alignment_tree_metadata_and_logs_are_verified(self):
        for name in ('alignment_active.fasta','tree_fasttree.nwk','logs/engine.log','session.json'):
            with self.subTest(name=name),tempfile.TemporaryDirectory() as tmp:
                session=self.session()
                (session.work/'tree_fasttree.nwk').write_text('(seq_1:0.1,seq_2:0.2);')
                (session.work/'logs').mkdir(exist_ok=True)
                (session.work/'logs/engine.log').write_text('recorded execution')
                target=save_project(session,Path(tmp)/'saved')
                self.assertTrue((target/INTEGRITY_FILE).is_file())
                read_project(target)
                file=target/name
                if name=='session.json':
                    state=json.loads(file.read_text());state['started_utc']='edited';file.write_text(json.dumps(state))
                else:file.write_text('changed')
                with self.assertRaisesRegex(ValueError,'changed'):read_project(target)

    def test_missing_manifest_rejected_but_legacy_project_still_loads(self):
        with tempfile.TemporaryDirectory() as tmp:
            target=save_project(self.session(),Path(tmp)/'saved')
            (target/INTEGRITY_FILE).unlink()
            with self.assertRaises(ValueError):read_project(target)
            state=json.loads((target/'session.json').read_text());state.pop('integrity')
            (target/'session.json').write_text(json.dumps(state))
            read_project(target)

    def test_reference_selection_roundtrips_without_original_paths(self):
        args=parse_args([str(ROOT/'examples/example_query.fna'),str(ROOT/'examples/example_refs.faa'),'--mode','blastx','--no-menu'])
        session=Session(args,UI(False));self.addCleanup(session.close)
        session.reference_selection=['ref_HSP_2','ref_HSP_1']
        with tempfile.TemporaryDirectory() as tmp:
            target=save_project(session,Path(tmp)/'saved')
            with redirect_stdout(io.StringIO()):resumed=resume_project(target,parse_args(['--resume',str(target),'--no-menu']),UI(False))
            self.addCleanup(resumed.close)
            self.assertEqual(resumed.reference_selection,session.reference_selection)

    def test_invalid_reference_selection_cannot_publish_a_project(self):
        session=self.session();session.reference_selection=['not_loaded']
        with tempfile.TemporaryDirectory() as tmp:
            target=Path(tmp)/'invalid'
            with self.assertRaisesRegex(ValueError,'reference selection'):save_project(session,target)
            self.assertFalse(target.exists())


class TreeIdentityTests(unittest.TestCase):
    def test_alias_restoration_preserves_comments_supports_and_lengths(self):
        source="(TREE_1[&note='a:b']:1e-05,(TREE_2:0.2,TREE_3:0.3)0.97:0.4)[&R];"
        mapping={'TREE_1':"query:O'Brien_β",'TREE_2':'ref:gene,2','TREE_3':'Δ_isoform(3)'}
        restored=relabel_newick(source,mapping)
        self.assertIn("'query:O''Brien_β'",restored)
        self.assertIn("[&note='a:b']:1e-05",restored)
        self.assertIn(')0.97:0.4)[&R];',restored)
        tree=parse_newick(restored)
        self.assertEqual(tree.children[0].label,mapping['TREE_1'])
        self.assertEqual(tree.children[0].length,1e-5)
        self.assertEqual(tree.children[1].label,'0.97')
        self.assertEqual(tree.children[1].length,0.4)
        svg=tree_svg(restored)
        self.assertIn('Δ_isoform(3)',svg)
        self.assertIn('query:O&#x27;Brien_β',svg)

    def test_tip_set_changes_and_ambiguous_names_are_rejected(self):
        mapping={'TREE_1':'query','TREE_2':'reference'}
        for source in ('(TREE_1:0.1,TREE_1:0.2);','(TREE_1:0.1,OTHER:0.2);','TREE_1;'):
            with self.subTest(source=source),self.assertRaises(ValueError):relabel_newick(source,mapping)
        with self.assertRaises(ValueError):relabel_newick('(TREE_1,TREE_2);',{'TREE_1':'same','TREE_2':'same'})
        with self.assertRaises(ValueError):relabel_newick('(TREE_1,TREE_2);',{'TREE_1':'bad\nlabel','TREE_2':'okay'})


class FreshOutputTests(unittest.TestCase):
    def test_interrupt_and_broken_terminal_stop_child_before_output_rollback(self):
        class SilentUI:
            tty=True
            def stage(self,*args):pass
            def info(self,*args):pass
            def success(self,*args):pass
        class BrokenTerminal:
            def __init__(self,error):self.error=error
            def write(self,*args):raise self.error
            def flush(self):pass
        for error in (KeyboardInterrupt(),OSError('terminal disconnected'),subprocess.TimeoutExpired('fixture',1)):
            with self.subTest(error=type(error).__name__),tempfile.TemporaryDirectory() as tmp:
                folder=Path(tmp);tools=Toolchain(folder,SilentUI(),auto_modules=False)
                output=folder/'result.txt';output.write_text('previous successful result')
                tools.version_cache[sys.executable]='test Python'
                processes=[];popen=subprocess.Popen
                def capture(*args,**kwargs):
                    process=popen(*args,**kwargs);processes.append(process);return process
                expected=KeyboardInterrupt if isinstance(error,KeyboardInterrupt) else CommandError
                with self.assertRaises(expected):
                    with tools.fresh_output(output),patch('av_tools.subprocess.Popen',side_effect=capture),patch('av_tools.sys.stdout',BrokenTerminal(error)):
                        tools.run([sys.executable,'-c','import time;print("new result",flush=True);time.sleep(30)'],
                                  'child fixture','child.log',stdout_path=output)
                self.assertEqual(len(processes),1)
                self.assertIsNotNone(processes[0].poll())
                self.assertEqual(output.read_text(),'previous successful result')
                self.assertTrue(tools.history[-1]['failed'])

    def test_success_without_new_search_file_cannot_reuse_previous_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);tools=Toolchain(folder,auto_modules=False)
            raw=folder/'raw_hsps.tsv';raw.write_text('previous successful evidence')
            with patch.object(tools,'execute'),self.assertRaises(CommandError):
                tools.search('blastn',folder/'query.fa',folder/'refs.fa')
            self.assertEqual(raw.read_text(),'previous successful evidence')

    def test_malformed_search_restores_old_file_but_empty_new_search_is_valid(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);tools=Toolchain(folder,auto_modules=False)
            raw=folder/'raw_hsps.tsv';raw.write_text('previous successful evidence')
            with patch.object(tools,'execute',side_effect=lambda _:raw.write_text('not a valid result row')),self.assertRaises(ValueError):
                tools.search('blastn',folder/'query.fa',folder/'refs.fa')
            self.assertEqual(raw.read_text(),'previous successful evidence')
            with patch.object(tools,'execute',side_effect=lambda _:raw.write_text('')):
                self.assertEqual(tools.search('blastn',folder/'query.fa',folder/'refs.fa'),[])
            self.assertEqual(raw.read_text(),'')

    def test_missing_new_msa_restores_previous_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);tools=Toolchain(folder,auto_modules=False)
            source=folder/'input.fa';source.write_text('>q\nACGT\n')
            result=folder/'alignment.fa';result.write_text('previous successful alignment')
            with patch.object(tools,'execute'),self.assertRaises(CommandError):tools.msa(source,result)
            self.assertEqual(result.read_text(),'previous successful alignment')

    def test_partial_pairwise_outputs_restore_both_previous_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp);tools=Toolchain(folder,auto_modules=False)
            text=folder/'comparison.txt';aligned=folder/'comparison.fasta'
            text.write_text('previous text');aligned.write_text('previous alignment')
            with patch.object(tools,'execute',side_effect=lambda _:text.write_text('partial new text')),self.assertRaises(CommandError):
                tools.emboss('needle',Seq('q','q','ACGT'),Seq('r','r','ACGA'),folder/'comparison',mode='blastn')
            self.assertEqual(text.read_text(),'previous text')
            self.assertEqual(aligned.read_text(),'previous alignment')


if __name__=='__main__':unittest.main()
