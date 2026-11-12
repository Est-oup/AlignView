"""Tests without external data; mock tools support integration tests."""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import (Seq,HSP,merge_intervals,missing_intervals,interval_length,select_hits,
                     reverse_complement,translate_nt,find_motifs,orfs,msa_metrics,distance_matrix,
                     dotplot,parse_blast_tab)
from av_visual import build_report

class CoreTests(unittest.TestCase):
    def test_intervals(self):
        self.assertEqual(merge_intervals([(10,20),(18,30),(5,8),(55,50)],100),[(5,8),(10,30),(50,55)])
        self.assertEqual(missing_intervals([(2,6),(5,10),(18,20)],20),[(1,1),(11,17)])
        self.assertEqual(interval_length(merge_intervals([(2,10),(6,20)])),19)

    def test_hsps_multiple_union(self):
        a=HSP('q','ref',98,10,1,30,1,10,1e-10,120,3,0,'ABCDEFGHIJ','ABCDEFGHIJ')
        b=HSP('q','ref',88,10,20,49,4,13,1e-9,115,3,0,'ABCDEFGHIJ','ABCDEFGHIK')
        c=HSP('q','other',75,5,75,89,1,5,1e-3,80,3,0,'ABCDE','ABCDE')
        s=select_hits([a,b,c],100,identity=80,coverage=45,max_hits=10)
        self.assertEqual(list(s),['ref'])
        self.assertEqual(len(s['ref']),2)
        self.assertEqual(s['ref'][0].matches,10)
        self.assertEqual(s['ref'][1].mismatches,1)
        self.assertEqual(interval_length(merge_intervals([x.qspan for x in s['ref']])),49)

    def test_translation(self):
        self.assertEqual(translate_nt('ATGGGATAA',1),'MG*')
        self.assertEqual(translate_nt('CAT',-1),'M')
        self.assertEqual(reverse_complement('ATGCRY'),'RYGCAT')

    def test_motif_orfs(self):
        s=Seq('X','X','ATGAAATAA')
        self.assertEqual(len(orfs(s,min_aa=1)),1)
        self.assertEqual(orfs(s,min_aa=1)[0]['protein'],'MK')
        hits=find_motifs(s,'ATG',True)
        self.assertEqual(hits[0]['start'],1)
        self.assertEqual(hits[0]['strand'],'+')
        self.assertEqual(len(find_motifs(Seq('s','s','AAAA'),'AA',False)),3)

    def test_msa(self):
        seqs=[Seq('a','a','ATGC-'),Seq('b','b','ATGT-'),Seq('c','c','ATGTG'),Seq('d','d','ATGCG')]
        m=msa_metrics(seqs)
        self.assertEqual(m['n'],4)
        self.assertEqual(m['variable'],1)
        self.assertEqual(m['parsimony_informative'],1)
        self.assertEqual(m['consensus'][:3],'ATG')
        self.assertEqual(distance_matrix(seqs)[0][1],.25)

    def test_dotplot(self):
        d=dotplot('ATGACATG','ATGCAATG',k=3)
        self.assertTrue(d['points'])
        self.assertEqual(d['xlen'],8)

    def test_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            q=Seq('q<script>','q<script>','ATGGCTATGGCT')
            r=[Seq('ref','ref','ATGGCT')]
            h=[HSP(q.id,'ref',100,4,1,4,1,4,1e-5,50,0,0,'ATGG','ATGG')]
            path=build_report(Path(tmp),q,r,h,'blastn',{'evalue':1e-5},[])
            text=path.read_text()
            self.assertIn('<html',text)
            self.assertNotIn('"q<script>"',text)
            self.assertTrue((Path(tmp)/'coverage.svg').exists())

class IntegrationTests(unittest.TestCase):
    def run_cli(self,*argv,env=None):
        e=os.environ.copy();e['PATH']=str(ROOT/'tests'/'mockbin')+os.pathsep+e['PATH']
        if env:e.update(env)
        return subprocess.run(['bash',str(ROOT/'alignview.sh'),*map(str,argv)],
                              cwd=self.temp,env=e,capture_output=True,text=True,timeout=90)

    def setUp(self):self.temp=Path(tempfile.mkdtemp(prefix='alignview_unit_'))
    def tearDown(self):shutil.rmtree(self.temp)

    def test_save_single(self):
        result=self.run_cli(ROOT/'examples'/'example_msa.fasta','--msa','view','--report','--save','--no-menu','--out','output')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.temp/'output'/'report.html').is_file())
        self.assertTrue((self.temp/'output'/'msa_consensus.fasta').is_file())

    def test_no_save(self):
        result=self.run_cli(ROOT/'examples'/'example_query.fna','--no-menu','--out','should_not_exist')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertFalse((self.temp/'should_not_exist').exists())

    def test_blastx(self):
        result=self.run_cli(ROOT/'examples'/'example_query.fna',ROOT/'examples'/'example_refs.faa',
                            '--mode','blastx','--msa','full','--report','--save','--no-menu','--out','output')
        self.assertEqual(result.returncode,0,result.stderr)
        out=self.temp/'output'
        self.assertTrue((out/'alignment_hits.fasta').is_file())
        self.assertTrue((out/'coverage.svg').is_file())
        self.assertTrue((out/'hsp_query.bed').is_file())
        self.assertEqual(len((out/'hsps_selected.tsv').read_text().splitlines()),4)

    def test_diamond(self):
        result=self.run_cli(ROOT/'examples'/'example_query.fna',ROOT/'examples'/'example_refs.faa',
                            '--mode','diamond-blastx','--save','--no-menu','--out','output')
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.temp/'output'/'hsps_all.tsv').is_file())
        command=(self.temp/'output'/'commands.json').read_text()
        self.assertIn('qseq_translated',command)
        self.assertIn('--min-orf 1',command)

    def test_global_and_pairwise(self):
        result=self.run_cli(ROOT/'examples'/'example_msa.fasta',ROOT/'examples'/'example_msa.fasta',
                            '--mode','global-dna','--msa','all','--pairwise','global',
                            '--no-menu','--save','--out','output')
        self.assertEqual(result.returncode,1) # Multi-sequence QUERY deliberately rejected
        self.assertIn('--query-id',result.stderr)

    def test_invalid_sequence(self):
        f=self.temp/'bad.fasta';f.write_text('hello\n')
        result=self.run_cli(f,'--no-menu')
        self.assertNotEqual(result.returncode,0)
        self.assertIn('FASTA',result.stderr)

if __name__=='__main__':unittest.main()
