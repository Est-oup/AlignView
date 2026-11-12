"""v6 regressions for units, figures and output modes."""
from __future__ import annotations
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq,HSP
from av_units import units_for,coord_audit
from av_report_v6 import tracks_svg,frame_tracks_svg,atlas_additions
from alignview import parse_args,Session
from av_ui import UI

class V6Tests(unittest.TestCase):
    def test_units_across_modes(self):
        expected={'blastn':('nt','nt','nt'),'blastp':('aa','aa','aa'),
                  'blastx':('nt','aa','aa'),'diamond-blastx':('nt','aa','aa'),
                  'tblastn':('aa','nt','aa'),'tblastx':('nt','nt','aa')}
        for mode,triplet in expected.items():
            u=units_for(mode)
            self.assertEqual((u.query,u.subject,u.aligned),triplet)

    def test_coordinate_translation_and_reverse(self):
        h=HSP('q','s',98,20,100,41,4,23,1e-15,99,-2,0,'M'*20,'M'*20)
        row=coord_audit([h],'blastx')[0]
        self.assertEqual(row['query_span'],60)
        self.assertEqual(row['reference_span'],20)
        self.assertTrue(row['query_consistent'])
        self.assertTrue(row['reference_consistent'])
        self.assertEqual(row['query_unit'],'nt')
        self.assertEqual(row['alignment_unit'],'aa')
        svg=tracks_svg([h],150,{'s':30},'blastx')
        self.assertIn('QUERY: nt · REFERENCE: aa',svg)
        self.assertIn('150 nt',svg)
        self.assertIn('30 aa',svg)
        self.assertNotIn('<line x1=',svg)
        self.assertIn('Cadre -2',frame_tracks_svg([h],150))

    def test_preset_cli_explicit_precedence(self):
        a=parse_args(['examples/example_query.fna','examples/example_refs.faa','--mode','blastx','--preset','strict','--no-menu'])
        self.assertEqual(a.evalue,1e-10)
        self.assertEqual(a.hits,100)
        b=parse_args(['examples/example_query.fna','examples/example_refs.faa','--mode','blastx','--preset','sensitive', '--evalue','0.5','--no-menu'])
        self.assertEqual(b.evalue,0.5)
        self.assertEqual(b.hits,100)

    def test_translated_dotplot_uses_amino_acids(self):
        a=parse_args(['examples/example_query.fna','examples/example_refs.faa','--mode','blastx','--no-menu'])
        session=Session(a,UI(color=False))
        try:
            ref=session.references[0]
            h=HSP(session.query.id,ref.id,100,15,3,47,1,15,1e-20,120,3,0,'A'*15,'A'*15)
            session.selected={ref.id:[h]}
            dot=session._compute_dotplot(ref,3)
            self.assertEqual((dot['xunit'],dot['yunit']),('aa','aa'))
            self.assertEqual(dot['qframe'],3)
            self.assertIn('translated',dot['note'])
        finally:session.close()

    def test_export_report_and_figures(self):
        env=os.environ.copy()
        env['PATH']=str(ROOT/'tests/mockbin')+os.pathsep+str(ROOT/'tests/mockextra')+os.pathsep+env['PATH']
        with tempfile.TemporaryDirectory() as tmp:
            destination=Path(tmp)/'results'
            p=subprocess.run([sys.executable,str(ROOT/'alignview.py'),str(ROOT/'examples/example_query.fna'),
              str(ROOT/'examples/example_refs.faa'),'--mode','blastx','--msa','hits','--preset','sensitive',
              '--export-profile','publication','--save','--out',str(destination),'--no-menu',
              '--module-policy','never'],env=env,capture_output=True,text=True,timeout=55)
            self.assertEqual(p.returncode,0,p.stdout+'\n'+p.stderr)
            report=destination/'report.html'
            self.assertTrue(report.is_file())
            page=report.read_text()
            figures=list((destination/'figures').glob('*.svg'))
            self.assertGreaterEqual(len(figures),30)
            self.assertEqual(page.count('class="panel figcard"'),len(figures))
            self.assertEqual(len(re.findall(r'data-figure="',page)),len(figures))
            self.assertTrue((destination/'exports'/'report.html').is_file())
            self.assertTrue((destination/'exports'/'coordinate_units.tsv').is_file())
            self.assertEqual(len((destination/'exports'/'figures_catalog.tsv').read_text().splitlines())-1,len(figures))
            self.assertIn('QUERY: nt · REFERENCE: aa',(destination/'figures'/'23_hsp_mapping.svg').read_text())
            self.assertIn('<strong>Units:</strong> Translated search',page)

if __name__=='__main__':unittest.main()

class UnitsReportTests(unittest.TestCase):
    def test_subject_nucleotide_frames(self):
        from av_report_v6 import atlas_additions
        prot=Seq('p','p','M'+'A'*54)
        genome=Seq('g','g','ATG'*110)
        h=HSP('p','g',100,20,5,24,180,121,1e-12,150,0,-1,'A'*20,'A'*20)
        figs,_=atlas_additions(prot,[genome],[h],None,'tblastn')
        self.assertIn('39_ratio_nt_aa_reference',figs)
        self.assertIn('37_reference_frames',figs)
        self.assertIn('QUERY: aa · REFERENCE: nt',figs['35_pistes_query_reference'] if '35_pistes_query_reference' in figs else tracks_svg([h],55,{'g':330},'tblastn'))
        self.assertNotIn('36_query_frame_coverage',figs)

    def test_both_nucleotide_translated(self):
        from av_report_v6 import atlas_additions
        dna=Seq('q','q','ATG'*40)
        genome=Seq('g','g','ATG'*100)
        h=HSP('q','g',95,20,60,1,180,121,1e-25,130,-1,-2,'M'*20,'M'*20)
        rows=coord_audit([h],'tblastx')
        self.assertTrue(rows[0]['query_consistent'])
        self.assertTrue(rows[0]['reference_consistent'])
        figs,_=atlas_additions(dna,[genome],[h],None,'tblastx')
        self.assertIn('36_query_frame_coverage',figs)
        self.assertIn('37_reference_frames',figs)
        self.assertIn('38_ratio_nt_aa_query',figs)
        self.assertIn('39_ratio_nt_aa_reference',figs)
