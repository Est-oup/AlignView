"""Session integration contracts; mock engines do not assess biological accuracy."""
from __future__ import annotations

from contextlib import redirect_stdout,redirect_stderr
import copy
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_cli import parse_args
from av_core import HSP,Seq,read_fasta,write_fasta
from av_project import INTEGRITY_FILE,read_project
from av_navigation import Workbench
from av_session import Session
from av_tools import CommandError
from av_tree import TREE_FIGURE,parse_newick
from av_ui import UI
from av_workspace import Workspace,evidence_rows


class WorkspaceSessionTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory(prefix='alignview_workspace_contract_')
        self.addCleanup(self.temporary.cleanup)
        self.folder=Path(self.temporary.name)

    def session(self,query,references=None,*,mode=None,extra=()):
        source=self.folder/('input_'+str(len(list(self.folder.glob('input_*')))))
        source.mkdir()
        qfile=source/'query.fa';write_fasta(query,qfile)
        argv=[str(qfile)]
        if references is not None:
            rfile=source/'references.fa';write_fasta(references,rfile);argv.append(str(rfile))
        if mode:argv+=['--mode',mode]
        argv+=['--no-menu','--no-open','--module-policy','never',
               '--tools-dir',str(ROOT/'tests/mockbin'),'--tools-dir',str(ROOT/'tests/mockextra'),*extra]
        result=Session(parse_args(argv),UI(False));self.addCleanup(result.close)
        return result

    def comparison(self):
        query=[Seq('query','query synthetic fixture','ACGTACGTTACGTA')]
        refs=[Seq('reference_a','reference_a synthetic fixture','ACGTACGTTACGTA'),
              Seq('reference_b','reference_b synthetic fixture','ACGTTCGTTACGTA'),
              Seq('reference_c','reference_c synthetic fixture','ACGTACGTTACGTT')]
        return self.session(query,refs,mode='blastn')

    def compute_baseline(self,session):
        with redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()),patch.object(session,'show_msa_window'):
            session.search();session.msa_run('hits');session.make_tree()
        self.assertTrue(session.selected)
        self.assertTrue(session.msa_file.is_file())
        self.assertTrue((session.work/'tree_fasttree.nwk').is_file())
        session.dirty=False

    @staticmethod
    def active_state(session):
        files={name:path.read_bytes() for name,path in (
            ('msa',session.msa_file),('tree',session.work/'tree_fasttree.nwk'),
            ('tree_figure',session.work/'figures'/(TREE_FIGURE+'.svg')),
            ('all_hsps',session.work/'hsps_all.tsv'),('selected_hsps',session.work/'hsps_selected.tsv'))}
        return (copy.deepcopy(session.hsps),copy.deepcopy(session.selected),
                copy.deepcopy(session.msa_records),session.msa_file,
                copy.deepcopy(session.calculation_runs),files)

    def assert_active_state(self,session,previous):
        self.assertEqual(self.active_state(session),previous)
        self.assertFalse(session.dirty)

    def test_selected_references_keep_biological_ids_and_safe_tree_labels(self):
        shared="shared:O'Brien_β"
        query=[Seq(shared,shared+' synthetic query','ACGTACGTTACGTA')]
        refs=[Seq(shared,shared+' synthetic reference','ACGTTCGTTACGTA'),
              Seq('other(2)','other(2) synthetic reference','ACGTACGTTACGTAAC'),
              Seq('excluded','excluded synthetic reference','TTTTACGTTACGTA')]
        session=self.session(query,refs,mode='global-dna')
        inputs={p:p.read_bytes() for p in (session.source_query_path,session.source_refs_path)}
        session.reference_selection=['other(2)',shared]
        with redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()),patch.object(session,'show_msa_window'):
            session.msa_run('selection');session.make_tree()
        expected=['query:'+shared,'other(2)','ref:'+shared]
        self.assertEqual([r.id for r in session.msa_records],expected)
        self.assertEqual([r.id for r in read_fasta(session.msa_file)],expected)
        self.assertEqual(session.calculation_runs['msa']['preparation']['reference_ids'],['other(2)',shared])
        self.assertEqual(session.calculation_runs['msa']['preparation']['reference_alignment_ids'],
                         {'other(2)':'other(2)',shared:'ref:'+shared})
        self.assertNotIn('excluded',[r.id for r in session.msa_records])
        self.assertEqual(session.msa_records[0].seq.replace('-',''),query[0].seq)
        self.assertEqual(session.msa_records[1].seq.replace('-',''),refs[1].seq)
        self.assertEqual(session.msa_records[2].seq.replace('-',''),refs[0].seq)
        newick=(session.work/'tree_fasttree.nwk').read_text()
        tree=parse_newick(newick);nodes=[tree];tips=[]
        while nodes:
            node=nodes.pop()
            if node.children:nodes.extend(node.children)
            else:tips.append(node.label)
        self.assertEqual(set(tips),set(expected))
        self.assertIn("O''Brien",newick)
        self.assertTrue((session.work/'figures'/(TREE_FIGURE+'.svg')).is_file())
        tree_aliases=session.calculation_runs['tree']['sequence_aliases']
        self.assertEqual(set(tree_aliases.values()),set(expected))
        self.assertTrue(all(key.startswith('TREE_') for key in tree_aliases))
        for path,content in inputs.items():self.assertEqual(path.read_bytes(),content)

    def test_translated_checked_reference_alignment_scopes_hsps_and_keeps_safe_aliases(self):
        coding='ATGCCGGAACCGATTACCGATGAAAAATTA'
        for mode in ('tblastn','tblastx'):
            with self.subTest(mode=mode):
                query=Seq('active_query','active_query synthetic fixture',
                          'MPEPITDEKL' if mode=='tblastn' else coding)
                refs=[Seq(name,name+' synthetic fixture',coding)
                      for name in ('excluded','QUERY_HSP','selected:β')]
                session=self.session([query],refs,mode=mode)
                qend,qframe=(4,0) if mode=='tblastn' else (12,1)
                hsps=[HSP(query.id,ref.id,100,4,1,qend,1,12,1e-20,100-i*10,
                          qframe,1,'MPEP','MPEP') for i,ref in enumerate(refs)]
                with patch.object(session.tools,'search',return_value=hsps),redirect_stdout(io.StringIO()):
                    session.search()
                self.assertEqual(set(session.selected),{r.id for r in refs})
                session.reference_selection=['selected:β','QUERY_HSP']
                commands={command.key:command for command in Workspace(session).commands()}
                self.assertTrue(commands['selection'].enabled)
                self.assertTrue(commands['translation_table'].enabled)
                workbench=Workbench(session)
                with patch.object(workbench,'choose',return_value='1'),patch.object(session,'show_msa_window'),redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()):
                    self.assertTrue(workbench.align('selection'))
                self.assertEqual([record.id for record in session.msa_records],
                                 ['query:QUERY_HSP','selected:β','QUERY_HSP'])
                self.assertTrue(all(record.seq.replace('-','')=='MPEP' for record in session.msa_records))
                preparation=session.calculation_runs['msa']['preparation']
                self.assertEqual(preparation['anchor_subject'],'selected:β')
                self.assertEqual(preparation['checked_reference_scope'],session.reference_selection)
                self.assertEqual(preparation['reference_ids'],session.reference_selection)
                self.assertEqual(preparation['reference_alignment_ids'],
                                 {'selected:β':'selected:β','QUERY_HSP':'QUERY_HSP'})
                self.assertEqual([record.id for record in read_fasta(session.work/'msa_hsp_input.fasta')],
                                 ['HSP_1','HSP_2','HSP_3'])
                self.assertNotIn('excluded',[record.id for record in session.msa_records])

    def test_command_palette_does_not_offer_alignment_for_one_sequence(self):
        session=self.session([Seq('only','only synthetic fixture','ACGT')])
        commands={command.key:command for command in Workspace(session).commands()}
        self.assertFalse(commands['align'].enabled)
        self.assertIn('At least two sequences',commands['align'].reason)
        self.assertFalse(commands['selection'].enabled)
        self.assertFalse(commands['open_alignment'].enabled)

    def test_malformed_search_does_not_replace_existing_evidence_msa_or_tree(self):
        session=self.comparison();self.compute_baseline(session)
        previous=self.active_state(session)
        invalid=HSP('different_query','reference_a',100,4,1,4,1,4,1e-8,50,0,0,'ACGT','ACGT')
        with patch.object(session.tools,'search',return_value=[invalid]),redirect_stdout(io.StringIO()),self.assertRaisesRegex(ValueError,'unexpected query'):
            session.search()
        self.assert_active_state(session,previous)

    def test_vsearch_replaces_active_local_evidence_and_keeps_archived_attempt(self):
        from av_library import entries
        from av_project import save_project,resume_project
        session=self.comparison()
        session.library_root=self.folder/'catalog'
        session.history=[]
        with redirect_stdout(io.StringIO()):session.search()
        local=entries(session,'result')[0]
        for matched in (True,False):
            def execute(commands):
                output=session.work/'pending_vsearch_search/matches_raw.tsv'
                output.write_text('query\treference_b\t95\t14\t1\t14\t1\t14\t100\t100\t+\tACGT\tACGT\n' if matched else '')
            with patch.object(session.tools,'execute',side_effect=execute),redirect_stdout(io.StringIO()):
                session.vsearch_run('search')
            rows,kind=evidence_rows(session)
            self.assertEqual(kind,'VSEARCH global search')
            self.assertEqual([r['subject'] for r in rows],['reference_b'] if matched else [])
            self.assertEqual(session.mode,'vsearch')
            self.assertFalse(session.search_completed)
            self.assertFalse(session.hsps)
            self.assertFalse(session.selected)
            self.assertEqual(read_project(local['path'])[1]['mode'],'blastn')
        saved=save_project(session,self.folder/'vsearch_saved')
        with redirect_stdout(io.StringIO()):
            resumed=resume_project(saved,parse_args(['--resume',str(saved),'--no-menu']),UI(False))
        self.addCleanup(resumed.close)
        self.assertEqual(evidence_rows(resumed),([],'VSEARCH global search'))

    def test_failed_vsearch_preserves_local_comparison(self):
        session=self.comparison();self.compute_baseline(session)
        previous=self.active_state(session)
        with patch.object(session.tools,'execute',side_effect=CommandError('failed')),redirect_stdout(io.StringIO()):
            with self.assertRaises(CommandError):session.vsearch_run('search')
        self.assertEqual(session.mode,'blastn')
        self.assertTrue(session.search_completed)
        self.assert_active_state(session,previous)

    def test_malformed_msa_preserves_previous_alignment_tree_and_provenance(self):
        session=self.comparison();self.compute_baseline(session)
        previous=self.active_state(session)
        for fault in ('missing_sequence','changed_residues'):
            def bad_alignment(source,output,**kwargs):
                records=read_fasta(source)
                if fault=='missing_sequence':records=records[:-1]
                else:
                    record=records[0]
                    first='C' if record.seq[0]!='C' else 'A'
                    records[0]=Seq(record.id,record.description,first+record.seq[1:])
                length=max(len(r.seq) for r in records)
                write_fasta([Seq(r.id,r.description,r.seq.ljust(length,'-')) for r in records],output)
                return output
            with self.subTest(fault=fault),patch.object(session.tools,'msa',side_effect=bad_alignment),redirect_stdout(io.StringIO()),self.assertRaises((ValueError,CommandError)):
                session.msa_run('all')
            self.assert_active_state(session,previous)

    def test_native_coverage_uses_prepared_ungapped_lengths(self):
        query=Seq('gapped_query','gapped_query synthetic fixture','A-C.G?TU')
        reference=Seq('gapped_reference','gapped_reference synthetic fixture','AC.G-TTU')
        session=self.session([query],[reference],mode='blastn')
        source_bytes={p:p.read_bytes() for p in (session.source_query_path,session.source_refs_path)}
        self.assertEqual(session.query.seq,'ACGNTT')
        self.assertEqual(session.references[0].seq,'ACGTTT')
        self.assertEqual(session.query_list[0].seq,query.seq)
        self.assertEqual(read_fasta(session.query_search_path)[0].seq,'ACGNTT')
        with redirect_stdout(io.StringIO()),redirect_stderr(io.StringIO()):session.search()
        row=session.coverage_rows()[0]
        self.assertEqual(row['query_aligned_length'],6)
        self.assertEqual(row['reference_length'],6)
        self.assertEqual(row['query_coverage'],100)
        self.assertEqual(row['reference_coverage'],100)
        preparation=session.calculation_runs['search']['preparation']
        self.assertEqual(preparation['query_gap_characters_removed'],2)
        self.assertEqual(preparation['reference_gap_characters_removed'],2)
        for path,content in source_bytes.items():self.assertEqual(path.read_bytes(),content)

    def test_save_shortcut_creates_new_snapshot_only_after_changes(self):
        records=[Seq('one','one synthetic fixture','ACGT'),Seq('two','two synthetic fixture','ACGA')]
        target=self.folder/'saved_project'
        session=self.session(records,extra=('--msa','view','--out',str(target)))
        session.load_msa()
        with patch.object(session.ui,'ask',return_value=str(target)),redirect_stdout(io.StringIO()):
            first=session.save_checkpoint()
        original_manifest=(first/INTEGRITY_FILE).read_bytes()
        original_state=(first/'session.json').read_bytes()
        with patch.object(session,'export_summary') as export,patch.object(session.ui,'ask') as ask,redirect_stdout(io.StringIO()):
            self.assertEqual(session.save_checkpoint(),first)
            export.assert_not_called();ask.assert_not_called()
        self.assertEqual(list(self.folder.glob('saved_project*')),[first])
        session.a.hits=3;session.save_settings()
        with patch.object(session.ui,'ask') as ask,redirect_stdout(io.StringIO()):second=session.save_checkpoint()
        ask.assert_not_called()
        self.assertNotEqual(second,first)
        self.assertEqual(second.parent,first.parent)
        self.assertEqual((first/INTEGRITY_FILE).read_bytes(),original_manifest)
        self.assertEqual((first/'session.json').read_bytes(),original_state)
        _,original=read_project(first);_,updated=read_project(second)
        self.assertEqual(original['args']['hits'],100)
        self.assertEqual(updated['args']['hits'],3)
        self.assertEqual(session.last_saved,second)
        self.assertFalse(session.dirty)
        with patch.object(session,'export_summary') as export,redirect_stdout(io.StringIO()):
            self.assertEqual(session.save_checkpoint(),second)
            export.assert_not_called()
        self.assertEqual(len(list(self.folder.glob('saved_project*'))),2)


if __name__=='__main__':unittest.main()
