#!/usr/bin/env python3
"""Exercise the offline MSA controls in headless Chrome, without a browser library."""
from pathlib import Path
import json
import re
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq,msa_metrics
from av_visual import report_html
from av_report_v6 import build_report_v6

records=[Seq(f'seq_{i:03d}',f'seq_{i:03d}','ACGT'*300) for i in range(120)]
payload=dict(units={'query':'nt','subject':'nt','aligned':'nt','msa':'nt'},queryName='example',queryLen=1200,
    queryStats={'length':1200,'kind':'ADN/ARN'},referenceCount=0,mode='msa',
    settings={'evalue':.01,'search_run':{'evalue':1e-8}},
    history=[{'command':'mafft --auto input.fasta','version':'1.2-test','duration_s':1,'exit_code':0,
        'stderr_log':'logs/test.log','stdout':'aligned.fasta'}],
    coverage=[],combined=[],coverageSVG='',profileSVG='',hsps=[],sequenceStats=[],dotplot=None,
    msa={'seqs':[{'id':r.id,'seq':r.seq} for r in records],'metrics':msa_metrics(records),'svg':''})
checks=r'''<script>window.addEventListener('DOMContentLoaded',()=>{
const result=document.createElement('pre');result.id='alignview-browser-check';document.body.appendChild(result);
try {
 const assert=(ok,why)=>{if(!ok)throw Error(why);};
 assert(document.getElementById('history').textContent.includes('version 1.2-test'),'executed version must remain visible');
 assert(document.getElementById('history').textContent.includes('\n  stdout : aligned.fasta'),'stdout destination must appear on its own line');
 assert(document.querySelector('[data-tab="coverage"]').hidden,'unused search tab must be hidden');
 assert(document.querySelector('[data-tab="dotplot"]').hidden,'unused dotplot tab must be hidden');
 assert(document.querySelector('[data-open-tab="msa"]'),'reading guide must link to the available alignment');
 const figure=document.querySelector('.figure-view');
 assert(figure&&!figure.open,'figures must be collapsed initially');
 window.dispatchEvent(new Event('beforeprint'));assert(figure.open,'printing must expand figures');
 window.dispatchEvent(new Event('afterprint'));assert(!figure.open,'printing must restore collapsed figures');
 tab('msa');
 assert(location.hash==='#msa','navigation must update browser history');
 assert(!document.getElementById('backToSummary').hidden,'summary return must be visible in an alignment');
 document.getElementById('backToSummary').click();
 assert(document.getElementById('overview').classList.contains('active'),'summary return must work');
 tab('msa');
 const board=document.getElementById('msaView'),position=document.getElementById('msaStart');
 board.scrollTop=180;const top=board.scrollTop;
 assert(top>0,'fixture must scroll vertically');
 document.getElementById('msaNext').click();
 assert(Number(position.value)===81,'horizontal page must advance by 80');
 assert(board.scrollTop===top,'horizontal page must preserve vertical scroll');
 board.dispatchEvent(new WheelEvent('wheel',{deltaX:100,shiftKey:true,cancelable:true}));
 assert(Number(position.value)>81,'trackpad must advance horizontally');
 assert(board.scrollTop===top,'trackpad must preserve vertical scroll');
 document.getElementById('msaLast').click();
 assert(Number(position.value)===1121,'last page must end at the last column');
 document.getElementById('msaRowsNext').click();
 assert(board.querySelector('.msa-name').textContent==='POSITION','ruler must remain visible');
 assert(board.querySelectorAll('.msa-name')[1].textContent==='seq_100','sequence page must start at row 101');
 const filter=document.getElementById('msaSeqFilter');filter.value='seq_119';filter.oninput();
 assert(board.querySelectorAll('.msa-name')[1].textContent==='seq_119','filter must reach a sequence beyond the first page');
 const variables=document.getElementById('msaVar');variables.checked=true;variables.onchange();
 assert(board.querySelectorAll('[data-column]').length===0,'invariant fixture must have zero variable columns');
 result.textContent='PASS: guided navigation, empty tabs, expandable figures, print, summary return, horizontal pages, trackpad, vertical position, bounds, row pages, sequence filter, empty variables';
} catch(error) {result.textContent='FAIL: '+error.message;}
});</script>'''

browser=shutil.which('google-chrome') or shutil.which('chromium') or shutil.which('chromium-browser')
if not browser:raise SystemExit('Chrome/Chromium is required for this optional browser check.')
with tempfile.TemporaryDirectory(prefix='alignview_browser_') as tmp:
    tmp=Path(tmp);page=tmp/'report.html'
    markup=report_html(payload).replace('<footer class="footer">',
        '<section id="figures" class="section"><article class="panel figcard"><h2>Example figure</h2><div class="figure-svg"><svg viewBox="0 0 10 10"><rect width="10" height="10"/></svg></div></article></section><footer class="footer">',1)
    page.write_text(markup.replace('</body>',checks+'</body>'),encoding='utf-8')
    command=[browser,'--headless','--no-sandbox','--disable-gpu','--disable-dev-shm-usage',
        '--disable-background-networking','--no-first-run','--window-size=1200,900',
        '--user-data-dir='+str(tmp/'profile'),'--dump-dom',page.as_uri()]
    run=subprocess.run(command,capture_output=True,text=True,timeout=30)
    marker=re.search(r'<pre id="alignview-browser-check">([^<]+)</pre>',run.stdout)
    if run.returncode or not marker:raise SystemExit('Browser check could not complete: '+run.stderr[-1500:])
    print(marker[1])
    if not marker[1].startswith('PASS:'):raise SystemExit(1)

    # Also exercise the actual report builder chain, including its atlas,
    # provenance and export sections, rather than only the synthetic viewport.
    full=tmp/'full';full.mkdir()
    tiny=[Seq('a','a','ACGTACGT'),Seq('b','b','ACGTACGA'),Seq('c','c','ACGTTCGT')]
    (full/'tree_fasttree.nwk').write_text('((a:0.1,b:0.2)0.97:0.3,c:0.4);')
    report=build_report_v6(full,tiny[0],tiny[1:],[],'blastn',
        {'hits':100,'search_completed':True,'evalue':.01,'search_run':{'evalue':1e-8}},[],msa_seqs=tiny)
    full_checks=r'''<script>window.addEventListener('DOMContentLoaded',async ()=>{
const result=document.createElement('pre');result.id='alignview-full-report-check';document.body.appendChild(result);
try {
 const assert=(ok,why)=>{if(!ok)throw Error(why);};
 assert(document.getElementById('searchSummary').textContent.includes('1e-8'),'search summary must show the executed threshold');
 assert(document.getElementById('readingGuide'),'reading guide must exist');
 assert(document.querySelectorAll('.figure-view').length===document.querySelectorAll('.figcard').length,'all atlas figures must be expandable');
 assert(document.querySelector('[data-tab="figures"]').textContent==='Figures','atlas entry must have a plain label');
 tab('figures');assert(document.getElementById('pageTitle').textContent==='Atlas des figures','atlas title must survive report composition');
 document.getElementById('figFilter').value='conservation';document.getElementById('figFilter').dispatchEvent(new Event('input'));
 const cards=[...document.querySelectorAll('.figcard')];assert(cards.some(x=>x.hidden)&&cards.some(x=>!x.hidden),'atlas filtering must remain usable');
 document.getElementById('openTree').click();
 const tree=document.querySelector('.figcard[data-figure="45_arbre_fasttree"]');
 assert(tree&&!tree.hidden&&tree.querySelector('details').open,'tree shortcut must reveal and expand the figure');
 assert(cards.filter(x=>!x.hidden).length===1,'phylogeny filter must isolate the tree');
 assert(tree.querySelectorAll('.tree-tip').length===3,'all tree labels must be rendered');
 assert(tree.querySelector('.tree-support').textContent==='0.97','internal support must remain visible');
 assert(tree.querySelector('.tree-scale'),'tree must include a branch length scale');
 const create=URL.createObjectURL,click=HTMLAnchorElement.prototype.click;let downloaded=null,filename='';
 URL.createObjectURL=blob=>{downloaded=blob;return create(blob);};
 HTMLAnchorElement.prototype.click=function(){filename=this.download;};
 tree.querySelector('.downloadSvg').click();
 URL.createObjectURL=create;HTMLAnchorElement.prototype.click=click;
 assert(filename==='45_arbre_fasttree.svg'&&downloaded,'tree SVG download must be available');
 const svg=await downloaded.text();assert(svg.includes('tree-tip')&&svg.includes('0.97'),'download must contain the actual tree and support');
 tab('provenance');assert(document.getElementById('advancedNav').closest('details').open,'advanced active navigation must be visible');
 tab('exports');assert(document.getElementById('exports').classList.contains('active'),'exports must remain accessible');
 result.textContent='PASS: complete report, real atlas, filtering, tree shortcut, supports, SVG download, parameters and exports';
} catch(error) {result.textContent='FAIL: '+error.message;}
});</script>'''
    report.write_text(report.read_text().replace('</body>',full_checks+'</body>'),encoding='utf-8')
    run=subprocess.run([*command[:-1],report.as_uri()],capture_output=True,text=True,timeout=30)
    marker=re.search(r'<pre id="alignview-full-report-check">([^<]+)</pre>',run.stdout)
    if run.returncode or not marker:raise SystemExit('Full report check could not complete: '+run.stderr[-1500:])
    print(marker[1])
    raise SystemExit(0 if marker[1].startswith('PASS:') else 1)
