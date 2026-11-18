#!/usr/bin/env python3
"""Exercise the real offline results workspace in Chrome, without a browser library."""
from pathlib import Path
import json
import re
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq,HSP
from av_report import build_report


def check(browser,report,script,marker,profile):
    source=report.read_text(encoding='utf-8')
    report.write_text(source.replace('</body>','<script>'+script+'</script></body>'),encoding='utf-8')
    run=subprocess.run([browser,'--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage',
        '--disable-background-networking','--no-first-run','--window-size=1440,1050','--virtual-time-budget=1800',
        '--user-data-dir='+str(profile),'--dump-dom',report.as_uri()],capture_output=True,text=True,timeout=45)
    found=re.search('<pre id="'+marker+'">([^<]+)</pre>',run.stdout)
    if run.returncode or not found:raise SystemExit('Browser check did not initialize: '+run.stderr[-1800:])
    print(found[1])
    if not found[1].startswith('PASS:'):raise SystemExit(1)


browser=shutil.which('google-chrome') or shutil.which('chromium') or shutil.which('chromium-browser')
if not browser:raise SystemExit('Chrome/Chromium is required for this optional check.')
with tempfile.TemporaryDirectory(prefix='alignview_workspace_browser_') as temporary:
    folder=Path(temporary)
    # Large enough to require both horizontal and sequence pages.
    rows=[Seq(f'seq_{i:03d}',f'seq_{i:03d}','ACGT'*300) for i in range(120)]
    rows[1].seq='G'+rows[1].seq[1:]
    rows[-1].seq=rows[-1].seq[:-1]+'A'
    report=build_report(folder/'alignment',rows[0],rows[1:],[],'msa',{},[
        {'command':'mafft --auto input.fasta','version':'1.2-test','duration_s':1,'exit_code':0,'stderr_log':'logs/test.log','stdout':'aligned.fasta'}],rows)
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',async()=>{
 const result=document.createElement('pre');result.id='check-alignment';document.body.appendChild(result);
 try {
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  assert(document.querySelector('.section.active').id==='msa','available alignment must open directly');
  assert(!document.querySelector('.metro-workflow'),'no imposed route');
  assert(document.getElementById('viewBack').disabled,'Back begins at the first viewed result');
  tab('sequences');assert(!document.getElementById('viewBack').disabled,'another view enables Back');
  document.getElementById('viewBack').click();
  await new Promise(resolve=>setTimeout(resolve,100));
  assert(document.querySelector('.section.active').id==='msa','Back restores the previous reading view');
  document.getElementById('viewForward').click();
  await new Promise(resolve=>setTimeout(resolve,100));
  assert(document.querySelector('.section.active').id==='sequences','Forward restores the next visited view');
  tab('msa');
  assert(document.querySelector('[data-tab="coverage"]').hidden,'no unused search view');
  assert(document.querySelector('[data-tab="dotplot"]').hidden,'no unused dotplot');
  assert(document.getElementById('history').textContent.includes('version 1.2-test'),'recorded tool version');
  assert(document.getElementById('history').textContent.includes('\n  stdout: aligned.fasta'),'stdout has a separate line');
  const board=document.getElementById('msaView');board.scrollTop=180;const top=board.scrollTop;
  assert(top>0,'fixture must scroll vertically');
  document.getElementById('msaNext').click();
  assert(Number(document.getElementById('msaStart').value)===81,'one horizontal page advances 80 columns');
  assert(board.scrollTop===top,'horizontal navigation preserves vertical position');
  board.dispatchEvent(new WheelEvent('wheel',{deltaX:100,shiftKey:true,cancelable:true}));
  assert(Number(document.getElementById('msaStart').value)===91,'trackpad moves columns');
  assert(board.scrollTop===top,'trackpad preserves vertical position');
  document.getElementById('msaLast').click();
  assert(Number(document.getElementById('msaStart').value)===1121,'last page stays within the alignment');
  document.getElementById('msaRowsNext').click();
  assert(board.querySelector('.msa-name').textContent==='POSITION','ruler stays visible');
  assert(board.querySelectorAll('.msa-name')[1].textContent==='seq_100','sequence pagination reaches row 101');
  const exportedIds=alignmentExport(true).split('\n').filter(line=>line.startsWith('>')).map(line=>line.slice(1).split(' ')[0]);
  const displayedIds=[...board.querySelectorAll('[data-trace]')].map(row=>row.dataset.trace);
  assert(JSON.stringify(exportedIds)===JSON.stringify(displayedIds),'view export contains exactly the current sequence page');
  assert(alignmentExport(false).split('\n').filter(line=>line.startsWith('>')).length===120,'full export retains all filtered rows');
  document.getElementById('msaRowsPrev').click();
  assert(alignmentExport(true).split('\n').filter(line=>line.startsWith('>')).length===100,'first page export contains 100 rows');
  document.getElementById('msaRowsNext').click();
  const trace=document.getElementById('msaReference');trace.value='r1';trace.onchange();
  tab('provenance');tab('msa');assert(trace.value==='r1','comparison trace persists across views');
  document.getElementById('msaFirst').click();
  const filter=document.getElementById('msaSeqFilter');filter.value='seq_119';filter.oninput();
  assert(board.querySelectorAll('.msa-name')[1].textContent==='seq_119','filter reaches later embedded sequences');
  document.getElementById('msaVar').checked=true;document.getElementById('msaVar').onchange();
  assert(board.querySelectorAll('[data-column]').length===2,'only original variable columns are shown');
  assert([...board.querySelectorAll('[data-column]')].map(x=>x.dataset.column).join(',')==='1,1200','masked columns keep native MSA coordinates');
  let blob=null,name='';const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
  URL.createObjectURL=value=>{blob=value;return create(value);};HTMLAnchorElement.prototype.click=function(){name=this.download;};
  document.querySelector('[data-export="msa_visible"]').click();
  assert(name==='alignview_alignment_view.fasta'&&blob,'scoped FASTA download');
  const fasta=await blob.text();assert(fasta==='>seq_119 alignview_msa_columns=1,1200\nAA\n','FASTA uses filtered rows and actual displayed residues');
  document.querySelector('[data-export="msa_columns"]').click();const map=await blob.text();
  assert(map==='export_column\tmsa_column\n1\t1\n2\t1200\n','column sidecar preserves original MSA positions');
  URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;
  document.getElementById('msaVar').checked=false;document.getElementById('msaVar').onchange();
  filter.value='';filter.oninput();document.getElementById('msaRowsPrev').click();document.getElementById('msaFirst').click();
  const cell=board.querySelector('[data-sequence="seq_000"][data-column="1"]');
  assert(cell.classList.contains('change'),'differences use the chosen sequence, not an implicit consensus');
  cell.click();assert(document.getElementById('msaCell').textContent.includes('ungapped position 1'),'residue inspection distinguishes coordinate systems');
  const motif=document.getElementById('msaMotif');motif.value='ACGT';motif.oninput();document.getElementById('msaFindNext').click();
  assert(Number(document.getElementById('msaStart').value)===5,'motif search follows the selected trace');
  result.textContent='PASS: direct alignment workspace, preserved view state, two-dimensional paging, trace comparison, coordinate inspection, exact motif search and scoped FASTA/column-map exports';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-alignment',folder/'profile-a')

    # A local comparison with Unicode and HTML-like identifiers, recorded score
    # ordering, a real tree and an alignment that can be linked to references.
    identifier='Δ_query<script>'
    query=Seq(identifier,identifier,'ACGT'*50)
    refs=[Seq('Ω_ref','Ω_ref','ACGT'*10),Seq('ref_B','ref_B','ACGT'*10)]
    aligned=[Seq(identifier,identifier,'ACGTACGT'),Seq('Ω_ref','Ω_ref','ACGTTCGA'),Seq('ref_B','ref_B','ACGTACGA')]
    hsps=[HSP(identifier,'Ω_ref',100,4,20,23,12,15,1e-6,20,0,0,'ACGT','ACGT'),
          HSP(identifier,'Ω_ref',75,8,1,8,1,8,1e-8,60,0,0,'ACGTACGT','ACGTTCGA'),
          HSP(identifier,'ref_B',87.5,8,1,8,1,8,1e-7,40,0,0,'ACGTACGT','ACGTACGA')]
    comparison=folder/'comparison';comparison.mkdir()
    (comparison/'tree_fasttree.nwk').write_text("(('"+identifier+"':0.1,'Ω_ref':0.2)0.97:0.3,ref_B:0.4);")
    report=build_report(comparison,query,refs,hsps,'blastn',{'search_completed':True,'evalue':.01,'search_run':{'evalue':1e-8},
        'original_query':'/samples/queries.fna','original_references':'/database/targets.fna'},[],aligned)
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',async()=>{
 const result=document.createElement('pre');result.id='check-evidence';document.body.appendChild(result);
 try {
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  assert(document.documentElement.lang==='en','English document');
  assert(!document.querySelector('[data-tab="coverage"]').hidden,'matches remain directly accessible');
  assert(!window.pwned,'identifiers must remain data');
  assert(document.getElementById('searchSummary').textContent.includes('1e-8'),'executed threshold retained');
  assert(!document.getElementById('comparisonContext').hidden,'comparison input roles remain visible outside result tabs');
  assert(document.getElementById('contextQueries').textContent==='Δ_query<script>','exact query identifier is visible safely');
  assert(document.getElementById('contextQueryFile').textContent==='From: queries.fna','original query filename is visible');
  assert(document.getElementById('contextTargetFile').textContent==='targets.fna','target filename is visible');
  assert(document.getElementById('contextTargetCount').textContent==='2 target sequences','target count is explicit');
  assert(document.getElementById('contextPaths').textContent.includes('/database/targets.fna'),'full paths distinguish identically named target files');
  tab('coverage');
  const first=document.querySelector('#hitRows [data-subject]');assert(first.dataset.subject==='Ω_ref','strongest bitscore orders references');
  assert(first.children[4].textContent==='75.0%','identity comes from the strongest segment, not maximum identity');
  first.querySelector('[data-select]').click();
  document.getElementById('matchesToMsa').click();
  assert(document.querySelector('.section.active').id==='msa','selected reference links to alignment');
  assert(document.getElementById('msaPage').textContent.includes('of 2 filtered'),'linked alignment includes query and selected reference');
  document.getElementById('openTree').click();
  assert(document.querySelector('.section.active').id==='tree','tree has a dedicated reading view');
  assert(document.querySelectorAll('#treeCanvas .tree-tip').length===3,'all labels preserved');
  assert(document.querySelector('#treeCanvas .tree-tip').textContent==='Δ_query<script>','Unicode/HTML-like labels remain literal data');
  assert(document.querySelector('#treeCanvas .tree-support').textContent==='0.97','local support is visible');
  assert(document.querySelector('#treeCanvas .tree-scale'),'branch length scale');
  const tip=[...document.querySelectorAll('#treeCanvas .tree-tip')].find(x=>x.textContent==='Ω_ref');tip.dispatchEvent(new MouseEvent('click',{bubbles:true}));
  assert(document.getElementById('msaReference').selectedOptions[0].textContent==='Ω_ref','tree label links to comparison trace');
  let blob=null,name='';const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
  URL.createObjectURL=value=>{blob=value;return create(value);};HTMLAnchorElement.prototype.click=function(){name=this.download;};
  document.querySelector('#tree .downloadSvg').click();
  assert(name==='45_fasttree_tree.svg'&&blob,'actual SVG download');
  const svg=await blob.text();assert(svg.includes('tree-tip')&&svg.includes('0.97'),'download contains labels and support');
  URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;
  tab('figures');document.getElementById('figFilter').value='conservation';document.getElementById('figFilter').oninput();
  const cards=[...document.querySelectorAll('.figcard')];assert(cards.some(x=>x.hidden)&&cards.some(x=>!x.hidden),'supplementary figures are searchable');
  const details=document.querySelector('.figure-view');assert(!details.open,'figures start collapsed');
  window.dispatchEvent(new Event('beforeprint'));assert(details.open,'print opens figure content');window.dispatchEvent(new Event('afterprint'));assert(!details.open,'print restores disclosure state');
  tab('provenance');assert(document.getElementById('advancedNav').closest('details').open,'active secondary view remains visible');
  result.textContent='PASS: ranked local evidence, honest identity semantics, linked selected references/tree/trace, safe Unicode labels, real SVG download, figure filtering, print and methods';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-evidence',folder/'profile-b')

    native=folder/'vsearch';(native/'vsearch').mkdir(parents=True)
    (native/'vsearch/analyses.json').write_text(json.dumps({'search':{'summary':{'matches':1,'matched_references':1},'settings':{'vsearch_id':.97},'commands':[],
        'preview':[{'query':identifier,'target':'Ω_ref','id':100,'qcov':100,'tcov':100,'qstrand':'+'}]}},ensure_ascii=False))
    report=build_report(native,query,refs,[],'vsearch',{'hits':100},[])
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',()=>{
 const result=document.createElement('pre');result.id='check-vsearch';document.body.appendChild(result);
 try{
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  assert(document.querySelector('.section.active').id==='vsearch','native results open directly');
  assert(document.querySelector('[data-tab="coverage"]').hidden,'VSEARCH is not a local-HSP view');
  assert(document.querySelector('[data-tab="msa"]').hidden,'no alignment entry without an alignment');
  assert(document.querySelector('#vsearch td').textContent==='Δ_query<script>','query identifier preserved');
  assert(document.querySelector('#vsearch td:nth-child(2)').textContent==='Ω_ref','reference identifier preserved');
  assert(document.querySelector('[data-export="msa_visible"]').disabled,'unavailable exports explain state through disabled controls');
  result.textContent='PASS: native VSEARCH workspace, correct empty states and preserved identifiers';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-vsearch',folder/'profile-c')

    shared=folder/'shared-id'
    report=build_report(shared,Seq('shared','shared','ACGT'),[Seq('shared','shared','ACGTACGT')],
        [HSP('shared','shared',100,4,1,4,5,8,1e-8,40,0,0,'ACGT','ACGT')],
        'blastn',{'search_completed':True},[],
        [Seq('query:shared','shared','ACGT----'),Seq('ref:shared','shared','ACGTACGT')])
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',()=>{
 const result=document.createElement('pre');result.id='check-shared';document.body.appendChild(result);
 try{
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  tab('coverage');document.getElementById('hitToMsa').click();
  assert(document.getElementById('msaReference').selectedOptions[0].textContent==='ref:shared','reference evidence resolves to the reference-qualified alignment row');
  assert(document.querySelector('[data-sequence="query:shared"][data-column="5"]').classList.contains('change'),'gap differences are highlighted against the chosen trace');
  const query=document.querySelector('#seqRows [data-role="query"]');query.click();
  assert(document.getElementById('msaReference').selectedOptions[0].textContent==='query:shared','sequence data resolves query role independently');
  assert(document.querySelectorAll('#seqRows tr').length===2,'both inputs survive a shared identifier');
  result.textContent='PASS: identical query/reference identifiers retain source roles and link to distinct alignment rows';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-shared',folder/'profile-d')

    gapped_trace='T'*90+'AC-G.TAC-G.T?A'+'T'*100
    aligned=[Seq('gapped_trace','gapped_trace',gapped_trace),Seq('same_trace','same_trace',gapped_trace)]
    report=build_report(folder/'gapped-motifs',aligned[0],aligned[1:],[],'msa',{},[],aligned)
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',()=>{
 const result=document.createElement('pre');result.id='check-gapped-motifs';document.body.appendChild(result);
 try{
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  const motif=document.getElementById('msaMotif'),find=document.getElementById('msaFindNext'),status=document.getElementById('msaMotifStatus');
  document.getElementById('msaReference').value='r0';document.getElementById('msaReference').onchange();
  document.getElementById('msaVar').checked=true;document.getElementById('msaVar').onchange();
  motif.value='acgt';motif.oninput();find.click();
  assert(Number(document.getElementById('msaStart').value)===91,'ungapped search locates the original first alignment column');
  assert(status.textContent.includes('ungapped residues 91–94 · MSA columns 91–96'),'hit includes biological and alignment ranges');
  assert(!document.getElementById('msaVar').checked&&status.textContent.includes('Showing all columns'),'a masked motif is revealed explicitly');
  find.click();assert(status.textContent.includes('ungapped residues 95–98 · MSA columns 97–102'),'next hit skips gaps and advances in biological coordinates');
  find.click();assert(status.textContent.startsWith('Wrapped to start.'),'search wraps and explains its position');
  motif.value='?A';motif.oninput();find.click();
  assert(status.textContent.includes('ungapped residues 99–100 · MSA columns 103–104'),'unknown residues retain a biological position');
  document.querySelector('[data-sequence="gapped_trace"][data-column="103"]').click();
  assert(document.getElementById('msaCell').textContent.includes('ungapped position 99'),'residue inspection uses the same unknown-residue coordinate mapping');
  motif.value='AC-G.T';motif.oninput();find.click();
  assert(status.textContent.includes('Enter residues only'),'gap symbols explain the selected search mode');
  document.getElementById('msaMotifUngapped').checked=false;document.getElementById('msaMotifUngapped').onchange();find.click();
  assert(status.textContent.includes('MSA columns 91–96')&&!status.textContent.includes('ungapped residues'),'literal alignment search remains available');
  motif.value='GGGG';motif.oninput();find.click();assert(status.textContent.startsWith('Motif not found'),'empty search result is explicit');
  result.textContent='PASS: gap-aware motif search, original alignment and ungapped coordinates, next/wrap navigation and literal search';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-gapped-motifs',folder/'profile-e')

    # The segment reader must map the displayed window, not merely repeat
    # whole-HSP bounds. Reverse rows and translated codons use native units.
    qrow='A'*70+'-'+'C'*79
    srow='A'*69+'G'+'T'*80
    query=Seq('reverse_query','reverse_query','A'*600)
    target=Seq('target_seq','target_seq','A'*400)
    reverse=HSP(query.id,target.id,46,150,500,352,20,169,1e-15,120,-1,1,qrow,srow)
    second=HSP(query.id,target.id,100,4,10,13,300,303,1e-5,20,1,1,'ACGT','ACGT')
    report=build_report(folder/'segment-windows',query,[target],[reverse,second],'blastn',{'search_completed':True},[])
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',async()=>{
 const result=document.createElement('pre');result.id='check-segment-window';document.body.appendChild(result);
 try{
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  const info=document.getElementById('hspInfo'),windowInfo=document.getElementById('hspWindow');
  assert(info.textContent.includes('Q 500–352 nt (reverse strand)'),'whole segment states native reverse orientation');
  assert(windowInfo.textContent.includes('Displayed columns 1–140/150 (nt) · Q 500–362 nt · T 20–159 nt'),'first window advances native coordinates only for non-gap letters');
  const focused=[...document.querySelectorAll('#selectedTrack .hsp-focus')];
  assert(focused.length===2&&focused[0].title==='Selected segment: 352–500','chosen segment is marked independently of all covered regions');
  document.getElementById('hspNext').click();
  assert(windowInfo.textContent.includes('Displayed columns 141–150/150 (nt) · Q 361–352 nt · T 160–169 nt'),'later page maps reverse and forward native spans correctly');
  assert(document.getElementById('hspNext').disabled&&!document.getElementById('hspPrev').disabled,'paging bounds reflect available segment columns');
  let blob=null,name='';const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;
  URL.createObjectURL=value=>{blob=value;return create(value);};HTMLAnchorElement.prototype.click=function(){name=this.download;};
  document.getElementById('hspDownload').click();const fasta=await blob.text();
  assert(name==='alignview_selected_segment.fasta','selected segment is directly downloadable');
  assert(fasta.includes('>query:reverse_query native_range=500-352 native_unit=nt')&&fasta.includes('>target:target_seq native_range=20-169 native_unit=nt'),'download records exact original role, coordinates and units');
  assert(fasta.split('\n')[1].length===150&&fasta.split('\n')[1].includes('-'),'download contains the whole selected aligned segment rather than only its displayed page');
  document.getElementById('hspPicker').value='1';document.getElementById('hspPicker').onchange();
  assert(windowInfo.textContent.includes('columns 1–4/4')&&info.textContent.includes('Q 10–13'),'choosing another segment resets paging');
  assert(document.querySelector('#selectedTrack .hsp-focus').title==='Selected segment: 10–13','coverage focus follows the chosen segment');
  URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;
  document.getElementById('hitSearch').value='absent';document.getElementById('hitSearch').oninput();
  assert(document.getElementById('hspDownload').disabled&&document.getElementById('hspNext').disabled,'no-result state disables stale segment controls');
  result.textContent='PASS: native HSP window coordinates, reverse strands, selected coverage focus, complete segment download and empty controls';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-segment-window',folder/'profile-f')

    amino_query='M'*70+'-'+'K'*79
    amino_target='M'*50+'.'+'R'*99
    query=Seq('translated_query','translated_query','ATG'*300)
    target=Seq('translated_target','translated_target','ATG'*300)
    translated=HSP(query.id,target.id,30,150,700,254,101,547,1e-15,120,-3,2,amino_query,amino_target)
    report=build_report(folder/'translated-segment',query,[target],[translated],'tblastx',{'search_completed':True},[])
    check(browser,report,r'''
window.addEventListener('DOMContentLoaded',()=>{
 const result=document.createElement('pre');result.id='check-translated-window';document.body.appendChild(result);
 try{
  const assert=(condition,message)=>{if(!condition)throw Error(message);};
  const info=document.getElementById('hspInfo'),windowInfo=document.getElementById('hspWindow');
  assert(info.textContent.includes('reverse strand · frame -3')&&info.textContent.includes('forward strand · frame +2'),'translated strands and reading frames are explicit');
  assert(windowInfo.textContent.includes('columns 1–140/150 (aa) · Q 700–284 nt · T 101–517 nt'),'amino-acid windows map to three native nucleotides per non-gap letter');
  document.getElementById('hspNext').click();
  assert(windowInfo.textContent.includes('Q 283–254 nt · T 518–547 nt'),'translated second page retains inclusive codon endpoints');
  const h=hspCurrent();h.qend=253;showHsp();
  assert(windowInfo.textContent.includes('Q unavailable (inconsistent segment coordinates)'),'legacy inconsistent spans are not mapped to invented residue coordinates');
  h.qseq='-'.repeat(140)+'M'.repeat(10);h.qstart=30;h.qend=1;hspOffset=0;showHsp();
  assert(windowInfo.textContent.includes('Q gap-only nt'),'a gap-only page invents no native query coordinates');
  document.getElementById('hspNext').click();
  assert(windowInfo.textContent.includes('Q 30–1 nt'),'residues after a gap-only page start at the original native coordinate');

  result.textContent='PASS: translated HSP amino-acid windows retain native nucleotide coordinates, reading frames and dot gaps';
 }catch(error){result.textContent='FAIL: '+error.message;}
});''','check-translated-window',folder/'profile-g')
