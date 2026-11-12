"""MSA viewer tests, including a simulated interactive terminal."""
from pathlib import Path
import os
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq
from av_msa_tui import MSAViewer, column_summary, msa_static


class MSAViewerTests(unittest.TestCase):
    def setUp(self):
        self.seqs=[Seq('s1','s1','ATGC-ATACG'),Seq('s2','s2','ATGT-ATACG'),
                   Seq('s3','s3','ATGTGATACG'),Seq('s4','s4','ATGC-ATATG')]

    def test_conservation_and_variability(self):
        con,var,gaps,scores=column_summary(self.seqs)
        self.assertEqual(len(con),10)
        self.assertEqual([i+1 for i,x in enumerate(var) if x],[4,9])
        self.assertAlmostEqual(gaps[4],.75)
        self.assertEqual(con[:3],'ATG')
        self.assertEqual(scores[0],1.0)

    def test_variable_navigation_and_search(self):
        v=MSAViewer(self.seqs,False)
        v._search('c:ATG')
        self.assertEqual(v.search_hits[0],0)
        v.only_variable=True
        self.assertEqual(len(v.cols),2)
        v._goto_genomic(8)
        self.assertEqual(v.cols[v.col],8)
        v._goto_genomic(0)
        self.assertEqual(v.cols[v.col],3)

    def test_static_display(self):
        output=msa_static(self.seqs,width=20)
        self.assertIn('CONSENSUS',output)
        self.assertIn('VARIABLE SITES',output)
        self.assertIn('s1',output)

    def test_real_terminal_navigation(self):
        try:import pexpect
        except ImportError:self.skipTest('pexpect unavailable')
        script=("import sys;sys.path.insert(0,sys.argv[1]);"
                "from av_core import Seq;from av_msa_tui import MSAViewer;"
                "s=[Seq('q','q','ATGC-ATACG'),Seq('r','r','ATGTGATATG')];"
                "print('RESULT:',MSAViewer(s,False).run())")
        p=pexpect.spawn(sys.executable,['-c',script,str(ROOT)],encoding='utf-8',timeout=8,env={**os.environ,'TERM':'xterm-256color'},dimensions=(32,110))
        try:
            p.expect('AlignView / Alignment')
            p.send('\x1b[C')
            p.send('v')
            p.send('g')
            p.expect('Go to MSA column')
            p.sendline('9')
            p.send('w')
            p.expect('RESULT: web')
            p.expect(pexpect.EOF)
            p.close()
            self.assertEqual(p.exitstatus,0)
        finally:
            if p.isalive():p.close(force=True)

if __name__=='__main__':unittest.main()
