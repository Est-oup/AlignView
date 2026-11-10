"""V5 : exports interopérables, portabilité, rendu HTML et absence de sauvegarde implicite."""
from __future__ import annotations
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq,HSP
from av_exports import _phylip, _nexus
from av_report_v5 import mapping_svg

class V5Tests(unittest.TestCase):
    def run_app(self,*options,external=True,expect=0):
        env=os.environ.copy();env['ALIGNVIEW_AUTO_MODULES']='0'
        if external:env['PATH']=os.pathsep.join([str(ROOT/'tests/mockextra'),str(ROOT/'tests/mockbin'),env['PATH']])
        run=subprocess.run([sys.executable,str(ROOT/'alignview.py'),*options,'--no-menu','--module-policy','never'],cwd=ROOT,env=env,capture_output=True,text=True,timeout=40)
        self.assertEqual(run.returncode,expect,run.stdout+'\n'+run.stderr)
        return run

    def test_export_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            output=Path(tmp)/'saved'
            self.run_app('examples/example_query.fna','examples/example_refs.faa','--mode','blastx',
                '--msa','hits','--report','--export-profile','publication','--save','--out',str(output))
            export=output/'exports'
            for name in ('MANIFEST.tsv','hsps_selected.csv','hsps_selected.json','hsp_query.bed',
                         'hsp_reference.bed','hsp_features.gff3','msa_alignment.nex',
                         'msa_alignment.phy','msa_alignment.sto','msa_alignment.aln','msa_positions.csv','msa_pairwise_distances.tsv',
                         'provenance.json','README_exports.txt'):
                self.assertTrue((export/name).is_file(),name)
            with (export/'hsps_selected.csv').open(newline='') as stream:
                records=list(csv.DictReader(stream))
            self.assertGreater(len(records),0)
            bed=(export/'hsp_query.bed').read_text().splitlines()[1].split('\t')
            self.assertEqual(int(bed[1]),min(int(records[0]['qstart']),int(records[0]['qend']))-1)
            self.assertEqual(bed[2],str(max(int(records[0]['qstart']),int(records[0]['qend']))))
            self.assertEqual(json.loads((export/'provenance.json').read_text())['version'],'0.1.0')
            self.assertIn('data-tab="exports"',(output/'report.html').read_text())
            self.assertGreaterEqual(len(list((output/'figures').glob('*.svg'))),20)

    def test_export_complete_without_save(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'export'
            self.run_app('examples/example_query.fna','examples/example_refs.faa','--mode','blastx',
                '--export-profile','complete','--export-dir',str(folder))
            self.assertTrue((folder/'hsps_all.jsonl').is_file())
            self.assertTrue((folder/'query_without_hsp.bed').is_file())
            self.assertTrue((folder/'MANIFEST.tsv').is_file())
            self.assertFalse(any(p.name.startswith('alignview_') for p in Path(tmp).iterdir()))
            self.assertEqual(set(Path(tmp).iterdir()),{folder})

    def test_no_implicit_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.run_app('examples/example_query.fna','examples/example_refs.faa','--mode','blastx',
                '--report','--out',str(Path(tmp)/'results'))
            self.assertEqual(list(Path(tmp).iterdir()),[])

    def test_single_fasta_publication(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'out'
            self.run_app('examples/example_msa.fasta','--msa','view','--report',
                '--export-profile','publication','--save','--out',str(out))
            self.assertTrue((out/'exports/msa_alignment.fasta').exists())
            self.assertTrue((out/'exports/msa_positions.csv').exists())
            self.assertTrue((out/'exports/msa_alignment.nex').exists())
            self.assertFalse((out/'exports/hsp_query.bed').stat().st_size>200)

    def test_custom_tool_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'out'
            self.run_app('examples/example_query.fna','examples/example_refs.faa',
                '--mode','blastx','--tools-dir',str(ROOT/'tests/mockbin'),
                '--save','--out',str(out),external=False)
            self.assertTrue((out/'hsps_all.tsv').exists())

    def test_custom_executable_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'out'
            self.run_app('examples/example_query.fna','examples/example_refs.faa',
                '--mode','blastx',
                '--executable',f'blastx={ROOT}/tests/mockbin/blastx',
                '--executable',f'makeblastdb={ROOT}/tests/mockbin/makeblastdb',
                '--save','--out',str(out),external=False)
            self.assertTrue((out/'hsps_all.tsv').exists())

    def test_phylip_ids_collision(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.phy'
            _phylip([Seq('a/b','a/b','ACGT'),Seq('a:b','a:b','ACGT')],path)
            lines=path.read_text().splitlines()
            self.assertEqual(lines[1].split()[0],'a_b')
            self.assertEqual(lines[2].split()[0],'a_b_2')

    def test_nexus_protein(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'msa.nex'
            _nexus([Seq('REF-1','REF-1','MPLDAA'),Seq('REF-2','REF-2','MPLDA-')],path)
            self.assertIn('DATATYPE=PROTEIN',path.read_text())

    def test_mapping_svg_escaped(self):
        h=HSP('q','<script>',95.,10,1,30,1,10,1e-9,50.,3,0,'M'*10,'M'*10)
        svg=mapping_svg([h],90,{'<script>':20})
        self.assertNotIn('<script>',svg)
        self.assertIn('&lt;script&gt;',svg)

    def test_unknown_tool_path_fails_cleanly(self):
        self.run_app('examples/example_query.fna','examples/example_refs.faa','--mode','blastx',
            '--executable','blastx=/does/not/exist',external=False,expect=1)

if __name__=='__main__':unittest.main()
