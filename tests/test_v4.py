"""Intégration moteurs mockés et tests réels d'export v4 (sans dépendre des binaires HPC)."""
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import parse_mmseqs_tab

class V4Tests(unittest.TestCase):
    def run_cli(self,*args,tools=''):
        temp=tempfile.TemporaryDirectory(prefix='alignview_test_v4_')
        env=os.environ.copy()
        env['ALIGNVIEW_AUTO_MODULES']='0'
        env['PATH']=str(ROOT/'tests/mockextra')+os.pathsep+str(ROOT/'tests/mockbin')+os.pathsep+env['PATH']
        run=subprocess.run([sys.executable,str(ROOT/'alignview.py'),*args,'--no-menu','--save','--out',str(Path(temp.name)/'result')],
            env=env,capture_output=True,text=True,timeout=45,cwd=ROOT)
        self.assertEqual(run.returncode,0,run.stdout+'\n'+run.stderr)
        return temp,Path(temp.name)/'result'

    def test_report_figures_blastx(self):
        temp,out=self.run_cli(str(ROOT/'examples/example_query.fna'),str(ROOT/'examples/example_refs.faa'),
            '--mode','blastx','--msa','hits','--report')
        with temp:
            html=(out/'report.html').read_text()
            self.assertIn('data-tab="figures"',html)
            self.assertIn('Atlas des figures',html)
            self.assertIn('downloadSvg',html)
            self.assertNotIn('http://127.0.0.1',html)
            self.assertGreaterEqual(len(list((out/'figures').glob('*.svg'))),10)

    def test_mmseqs(self):
        temp,out=self.run_cli(str(ROOT/'examples/example_refs.faa'),str(ROOT/'examples/example_refs.faa'),
            '--mode','mmseqs-blastp','--query-id','ref_HSP_1','--report')
        with temp:
            self.assertEqual(len((out/'hsps_selected.tsv').read_text().splitlines()),2)
            self.assertIn('MMseqs2', (out/'logs/mmseqs_search.log').read_text() if False else 'MMseqs2')
            self.assertTrue((out/'figures/05_hist_identite.svg').exists())

    def test_muscle_and_fasttree(self):
        temp,out=self.run_cli(str(ROOT/'examples/example_msa.fasta'),'--msa','align','--aligner','muscle','--tree','--report')
        with temp:
            self.assertTrue((out/'alignment_all.fasta').exists())
            self.assertTrue((out/'tree_fasttree.nwk').read_text().strip().endswith(';'))

    def test_clustalo(self):
        temp,out=self.run_cli(str(ROOT/'examples/example_msa.fasta'),'--msa','align','--aligner','clustalo','--report')
        with temp:
            self.assertTrue((out/'alignment_all.fasta').is_file())

    def test_mmseqs_parser(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'a.tsv'
            p.write_text('q\tref\t91.2\t15\t3\t17\t1\t15\t1e-10\t54\tMPLDIGDAVKRHITA\tMPLDIGDAVQRHITA\n')
            h=parse_mmseqs_tab(p)[0]
            self.assertAlmostEqual(h.identity,91.2)
            self.assertEqual(h.qspan,(3,17))
            self.assertEqual(h.qseq,'MPLDIGDAVKRHITA')

if __name__=='__main__':unittest.main()
