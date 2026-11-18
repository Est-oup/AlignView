"""Results workspace data integrity, safe embedding and compatibility boundaries."""
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import HSP,Seq
from av_report import build_payload,build_report,render_report,build_comparison_summary,MAX_MSA_ROWS,MAX_HSPS
from av_report_v4 import build_report_v4
from av_report_v5 import build_report_v5
from av_report_v6 import build_report_v6


def payload_from_page(page):
    found=re.search(r'<script id="reportData" type="application/json">(.*?)</script>',page,re.S)
    return json.loads(found[1])


class ReportWorkspaceTests(unittest.TestCase):
    def records(self):
        return [Seq('query','query','ACGTACGT'),Seq('reference','reference','ACGTACGA'),Seq('third','third','ACGTTCGT')]

    def test_comparison_names_original_query_and_target_files_and_current_query(self):
        records=self.records()
        settings={'original_query':'/data/my_queries.fna','original_references':'/databases/target_genomes.fna',
                  'query':'/tmp/prepared/query.fasta','references':'/tmp/prepared/references.fasta'}
        data=build_payload(records[0],records[1:],[],'blastn',settings,[])
        context=data['comparisonContext']
        self.assertEqual(context['queryFile'],'my_queries.fna')
        self.assertEqual(context['targetFile'],'target_genomes.fna')
        self.assertEqual(context['targetPath'],'/databases/target_genomes.fna')
        self.assertEqual(context['queryIds'],['query'])
        self.assertEqual(context['targetCount'],2)
        self.assertEqual(context['method'],'BLASTN')
        self.assertIsNone(build_payload(records[0],records[1:],[],'msa',{},[])['comparisonContext'])

    def test_group_summary_names_both_datasets_and_preserves_paths_as_data(self):
        batch={'query_ids':['q1','q2'],'reference_count':8,'mode':'blastn',
               'runs':[{'query_id':'q1','status':'completed','matches':1}]}
        settings={'original_query':'/data/queries.fna','original_references':'/db/target_<script>.fna'}
        with tempfile.TemporaryDirectory() as tmp:
            page=build_comparison_summary(Path(tmp),batch,settings).read_text()
        self.assertIn('Query sequences:</strong> q1, q2',page)
        self.assertIn('From: queries.fna',page)
        self.assertIn('target_&lt;script&gt;.fna · 8 target sequences',page)
        self.assertNotIn('<script>',page)
        self.assertIn('BLASTN · supplied local FASTA targets',page)

    def test_user_identifiers_remain_data_in_json_svg_and_html(self):
        identifier='Δ_</script><script>window.pwned=true</script>'
        records=[Seq(identifier,identifier,'ACGT'),Seq('Ω_ref','Ω_ref','ACGA')]
        with tempfile.TemporaryDirectory() as tmp:
            report=build_report(Path(tmp),records[0],records[1:],[],'msa',{},[],records)
            source=report.read_text()
            data=payload_from_page(source)
            self.assertEqual(data['queryName'],identifier)
            self.assertEqual(data['msa']['seqs'][0]['id'],identifier)
            self.assertNotIn('<script>window.pwned=true</script>',source)
            self.assertIn('\\u003c/script\\u003e',source)
            self.assertNotIn('MutationObserver',source)
            self.assertNotRegex(source,r'<(?:script|link)[^>]+(?:src|href)="https?://')

    def test_preview_caps_are_explicit_and_do_not_change_full_metrics(self):
        records=[Seq(f'seq_{i}',f'seq_{i}','ACGT') for i in range(MAX_MSA_ROWS+1)]
        data=build_payload(records[0],records[1:],[],'msa',{},[],records)
        self.assertEqual(len(data['msa']['seqs']),MAX_MSA_ROWS)
        self.assertEqual(data['msa']['metrics']['n'],MAX_MSA_ROWS+1)
        self.assertEqual(data['limits']['msaTotal'],MAX_MSA_ROWS+1)
        self.assertEqual(data['limits']['msaIncluded'],MAX_MSA_ROWS)
        self.assertIn('full files remain in the saved project',render_report(data))

    def test_strongest_hsp_identity_and_native_units_are_preserved(self):
        query=Seq('dna','dna','ATG'*40);ref=Seq('protein','protein','M'*30)
        hsps=[HSP('dna','protein',100,4,1,12,1,4,1e-6,20,1,0,'MMMM','MMMM'),
              HSP('dna','protein',75,4,31,42,10,13,1e-9,80,1,0,'MMMM','MMMA')]
        data=build_payload(query,[ref],hsps,'blastx',{'search_run':{'evalue':1e-8}},[])
        self.assertEqual(data['coverage'][0]['best_identity'],75)
        self.assertEqual(data['units'],{'query':'nt','subject':'aa','aligned':'aa','msa':'aa'})
        self.assertEqual(data['settings']['search_run']['evalue'],1e-8)
        self.assertEqual(data['hsps'][0]['qstart'],1)

    def test_same_identifier_in_query_and_reference_keeps_both_sources(self):
        query=Seq('shared','shared','ACGT')
        reference=Seq('shared','shared','ACGTACGT')
        aligned=[Seq('query:shared','shared','ACGT----'),Seq('ref:shared','shared','ACGTACGT')]
        data=build_payload(query,[reference],[],'blastn',{},[],aligned)
        self.assertEqual([row['length'] for row in data['sequenceStats']],[4,8])
        self.assertEqual([row['role'] for row in data['sequenceStats']],['query','reference'])
        self.assertEqual(data['limits']['sequenceRowsTotal'],2)
        self.assertEqual(data['alignmentSources'][0]['msaId'],'query:shared')
        self.assertEqual(data['alignmentSources'][1]['msaId'],'ref:shared')

    def test_recorded_collision_safe_aliases_override_prefix_guessing(self):
        query=Seq('shared','shared','ACGT')
        references=[Seq('shared','shared','ACGTACGT'),Seq('query:shared','query:shared','ACGTACGA')]
        alignment=[Seq('query:query:shared','shared','ACGT----'),Seq('ref:shared','shared','ACGTACGT'),references[1]]
        settings={'calculation_runs':{'msa':{'preparation':{
            'sequence_aliases':{'QUERY':'query:query:shared'},
            'reference_alignment_ids':{'shared':'ref:shared','query:shared':'query:shared'}}}}}
        data=build_payload(query,references,[],'blastn',settings,[],alignment)
        self.assertEqual([source['msaId'] for source in data['alignmentSources']],['query:query:shared','ref:shared','query:shared'])

    def test_shared_metrics_and_distance_work_budget_avoid_repeated_scans(self):
        from av_core import msa_metrics
        from av_figures import build_figures
        rows=self.records();metrics=msa_metrics(rows)
        with patch('av_figures.msa_metrics',side_effect=AssertionError('unexpected repeated MSA scan')),patch('av_figures.MAX_DISTANCE_WORK',1),patch('av_figures.distance_matrix') as distance:
            figures,_=build_figures(rows[0],rows[1:],[],rows,'msa',metrics=metrics)
            self.assertIn('09_conservation',figures)
            self.assertNotIn('12_p_distance',figures)
            distance.assert_not_called()

    def test_legacy_builders_delegate_to_the_single_renderer(self):
        records=self.records()
        for builder in [build_report_v4,build_report_v5,build_report_v6]:
            with self.subTest(builder=builder.__name__),tempfile.TemporaryDirectory() as tmp:
                with patch('av_report.render_report',wraps=render_report) as renderer:
                    path=builder(Path(tmp),records[0],records[1:],[],'msa',{},[],records)
                    self.assertTrue(path.is_file())
                    renderer.assert_called_once()
                    self.assertIn('id="msaReference"',path.read_text())

    def test_catalog_tracks_only_current_generated_figures(self):
        records=self.records()
        hsp=HSP('query','reference',87.5,8,1,8,1,8,1e-8,80,0,0,'ACGTACGT','ACGTACGA')
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            build_report(folder,records[0],records[1:],[hsp],'blastn',{},[],records)
            previous=folder/'figures/05_identity_histogram.svg'
            self.assertTrue(previous.exists())
            custom=folder/'figures/custom.svg';custom.write_text('user attachment')
            report=build_report(folder,records[0],records[1:],[],'msa',{},[])
            self.assertFalse(previous.exists())
            self.assertEqual(custom.read_text(),'user attachment')
            data=payload_from_page(report.read_text())
            self.assertEqual(len(data['figureCatalog']),report.read_text().count('class="panel figcard"'))
            self.assertFalse(list(folder.glob('.report-*.html')))


if __name__=='__main__':unittest.main()
