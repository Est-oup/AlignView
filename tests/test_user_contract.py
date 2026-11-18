"""Native segment reading and genuinely noninteractive CLI behavior."""
from contextlib import redirect_stdout
import io
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_cli import parse_args
from av_core import HSP,Seq,write_fasta
from av_session import Session
from av_ui import UI
from test_release import Terminal


class UserContractTests(unittest.TestCase):
    def setUp(self):
        temporary=tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup);self.folder=Path(temporary.name)
        self.input=self.folder/'aligned.fa'
        write_fasta([Seq('q','q','ACGTACGT'),Seq('s','s','ACGTTCGT')],self.input)
        self.session=Session(parse_args([str(self.input),'--no-menu']),UI(False))
        self.addCleanup(self.session.close)

    def alignment(self,mode,hsp,width):
        self.session.mode=mode
        output=io.StringIO()
        with redirect_stdout(output):self.session.print_alignment(hsp,width)
        text=output.getvalue()
        rows=re.findall(r'^  ([QS])\s+(\d+|—)\s+(\S+)\s+(\d+|—)$',text,re.M)
        return text,rows

    def test_native_reverse_segment_coordinates_handle_dot_gaps_and_lowercase(self):
        h=HSP('q','s',50,8,10,15,50,44,1e-5,30,1,-1,'Ac.-GTAC','ATGG-TAC')
        text,rows=self.alignment('blastn',h,4)
        self.assertEqual(rows,[('Q','10','Ac.-','11'),('S','50','ATGG','47'),
                               ('Q','12','GTAC','15'),('S','46','-TAC','44')])
        self.assertIn('|.  ',text)
        self.assertIn('1-based inclusive',text)
        self.assertNotIn('Coordinates are local',text)
        self.assertNotIn('identiques',text)

    def test_translated_reverse_query_and_target_use_native_codon_ranges(self):
        cases=[('blastx',HSP('q','s',50,6,100,89,50,54,1e-5,30,-2,0,'MP.-LK','MPA-LK'),
                [('Q','100','MP.','95'),('S','50','MPA','52'),('Q','94','-LK','89'),('S','53','-LK','54')]),
               ('tblastn',HSP('q','s',50,6,11,15,48,37,1e-5,30,0,-1,'MPA-LK','MP.-LK'),
                [('Q','11','MPA','13'),('S','48','MP.','43'),('Q','14','-LK','15'),('S','42','-LK','37')])]
        for mode,hsp,expected in cases:
            with self.subTest(mode=mode):
                text,rows=self.alignment(mode,hsp,3)
                self.assertEqual(rows,expected)
                self.assertIn('alignment columns (aa)',text)

    def test_all_gap_block_never_fabricates_a_native_residue_range(self):
        _,rows=self.alignment('blastn',HSP('q','s',50,4,10,11,20,23,1e-5,30,1,1,'-.AC','ACAC'),2)
        self.assertEqual(rows[0],('Q','—','-.','—'))
        self.assertEqual(rows[2],('Q','10','AC','11'))

    def test_native_global_matches_export_preserves_fields_and_empty_tables(self):
        s=self.session;native=s.work/'vsearch/search/matches.tsv';native.parent.mkdir(parents=True)
        header='query\ttarget\tid\tqcov\ttcov\tqstrand\n'
        s.vsearch_results={'search':{'summary':{'matches':1}}}
        for i,content in enumerate((header+'q\ts\t95\t80\t60\t-\n',header)):
            native.write_text(content)
            target=s.export_artifact('hits',self.folder/f'matches_{i}.tsv')
            self.assertEqual(target.read_text(),content)
            with self.assertRaises(FileExistsError):s.export_artifact('hits',target)
        s.search_completed=True
        target=s.export_artifact('hits',self.folder/'local.tsv')
        self.assertIn('best_identity',target.read_text())
        self.assertNotIn('qstrand',target.read_text())

    def test_no_menu_and_report_commands_exit_on_a_real_terminal_without_saving_questions(self):
        for args in ([str(self.input),'--msa','view','--no-menu'],
                     ['--demo','--no-menu'],
                     ['align',str(self.input),'--no-menu','--tools-dir',str(ROOT/'tests/mockbin')],
                     [str(self.input),str(self.input),'--query-id','q','--mode','global-dna','--pairwise','global',
                      '--no-menu','--tools-dir',str(ROOT/'tests/mockbin')],
                     ['report',str(self.input),'--html-out',str(self.folder/'report.html')]):
            with self.subTest(args=args):
                code='import sys;sys.path.insert(0,sys.argv[1]);from alignview import main;raise SystemExit(main(sys.argv[2:]))'
                terminal=Terminal(code,args=args)
                try:
                    terminal.expect('Finished: temporary results removed')
                    self.assertEqual(terminal.process.wait(timeout=3),0)
                    self.assertNotIn(b'Keep this session',terminal.output)
                    self.assertNotIn(b'Project directory',terminal.output)
                    self.assertNotIn(b'\x1b[?1049h',terminal.output)
                finally:terminal.close()

    def test_unsaved_data_never_claims_that_a_project_was_saved(self):
        output=io.StringIO()
        self.session.dirty=False
        with redirect_stdout(output):self.session.show_dashboard()
        self.assertIn('temporary session · not saved',output.getvalue())
        self.assertNotIn('project saved',output.getvalue())

    def test_batch_save_never_prompts_to_change_an_existing_destination_even_on_tty(self):
        s=self.session;s.ui.tty=True;s.a.save=True;s.a.out=str(self.folder)
        with patch.object(s.ui,'ask',side_effect=AssertionError('batch requested interactive input')),redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(ValueError,'Output directory already exists'):s.save_and_exit()
        self.assertTrue(self.input.is_file())


if __name__=='__main__':unittest.main()
