"""v3 workflows: persistent HTML, startup wizard and query selection."""
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_start import suggest_modes, is_fasta_name

class V3Tests(unittest.TestCase):
    def test_methods(self):
        self.assertEqual(suggest_modes('DNA/RNA','protein')[0],'blastx')
        self.assertEqual(suggest_modes('protein','DNA/RNA')[0],'tblastn')
        self.assertEqual(suggest_modes('DNA/RNA','DNA/RNA')[1],'global-dna')
        self.assertTrue(is_fasta_name('reads.fna'))
        self.assertFalse(is_fasta_name('config.txt'))

    def test_diagnose_without_fasta(self):
        env=os.environ.copy();env['ALIGNVIEW_AUTO_MODULES']='0'
        run=subprocess.run([sys.executable,str(ROOT/'alignview.py'),'--diagnose'],env=env,capture_output=True,text=True,timeout=15)
        self.assertEqual(run.returncode,0,run.stderr)
        self.assertIn('DIAGNOSTIC',run.stdout)

    def test_multiple_query_selection(self):
        with tempfile.TemporaryDirectory() as temp:
            env=os.environ.copy();env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+env['PATH']
            q=ROOT/'examples/example_msa.fasta'
            result=subprocess.run([sys.executable,str(ROOT/'alignview.py'),str(q),str(q),
                '--mode','global-dna','--query-id','seq_2','--msa','all','--save','--no-menu','--out','out'],
                cwd=temp,env=env,capture_output=True,text=True,timeout=40)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertTrue((Path(temp)/'out'/'alignment_all.fasta').is_file())
            settings=(Path(temp)/'out'/'parameters.json').read_text()
            self.assertIn('seq_2',settings)

    def test_web_report_has_no_network_dependency(self):
        with tempfile.TemporaryDirectory() as temp:
            result=subprocess.run([sys.executable,str(ROOT/'alignview.py'),
                str(ROOT/'examples/example_msa.fasta'),'--report','--save','--no-menu','--out','out'],
                cwd=temp,capture_output=True,text=True,timeout=40)
            self.assertEqual(result.returncode,0,result.stderr)
            content=(Path(temp)/'out'/'report.html').read_text()
            self.assertIn('<!doctype html>',content)
            self.assertIn('const D=',content)
            self.assertNotIn('<script src=',content)
            self.assertNotIn('http://127.0.0.1',content)
            self.assertNotIn('fetch(',content)

if __name__=='__main__':unittest.main()
