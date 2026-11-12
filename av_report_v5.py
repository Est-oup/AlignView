"""Figure atlas and export center for the standalone v5 HTML report."""
from __future__ import annotations
import collections
import html
import math
from pathlib import Path
from av_display import figure_title
from av_core import Seq,HSP,coverage_data,msa_metrics,sequence_stats,merge_intervals
from av_report_v4 import build_report_v4, bars_svg, hist_svg, plot_svg, scatter_svg


def mapping_svg(hsps:list[HSP],query_length:int,ref_lengths:dict):
    """Map HSP segments using each sequence's native coordinates."""
    limit=70;W=1050;L=100;T=80;R=55;B=55;H=520
    shown=hsps[:limit]
    max_s=max((ref_lengths.get(h.subject) or max(h.sstart,h.send) for h in shown),default=1)
    max_s=max(1,max_s); qn=max(1,query_length)
    x=lambda p:L+(p-1)/qn*(W-L-R)
    y=lambda p:T+(p-1)/max_s*(H-T-B)
    a=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="HSP correspondence map">',
       '<rect width="1050" height="520" fill="#102039" rx="14"/>',
       '<text x="28" y="33" fill="#cfe0f2" font-size="17" font-family="Arial">QUERY ↔ reference matches (native coordinates)</text>',
       '<text x="28" y="54" fill="#91a8c3" font-size="12" font-family="Arial">Each line represents an HSP, not global collinearity. References have different lengths.</text>',
       f'<rect x="{L}" y="{T}" width="{W-L-R}" height="{H-T-B}" fill="#172743" stroke="#435574"/>']
    for i in range(6):
        xx=L+i/5*(W-L-R); yy=T+i/5*(H-T-B)
        a.append(f'<line x1="{xx:.2f}" y1="{T}" x2="{xx:.2f}" y2="{H-B}" stroke="#344762" stroke-dasharray="2,5"/>')
        a.append(f'<line x1="{L}" y1="{yy:.2f}" x2="{W-R}" y2="{yy:.2f}" stroke="#344762" stroke-dasharray="2,5"/>')
        a.append(f'<text x="{xx:.1f}" y="{H-B+21}" text-anchor="middle" fill="#9bb4ce" font-size="12">{round(i/5*qn)}</text>')
        a.append(f'<text x="{L-7}" y="{yy+4:.1f}" text-anchor="end" fill="#9bb4ce" font-size="12">{round(i/5*max_s)}</text>')
    for h in shown:
        hue='#53d6b7' if h.qstart<=h.qend and h.sstart<=h.send else '#f9bd7b'
        label=html.escape(h.subject,quote=True)
        a.append(f'<line x1="{x(h.qstart):.2f}" y1="{y(h.sstart):.2f}" x2="{x(h.qend):.2f}" y2="{y(h.send):.2f}" stroke="{hue}" stroke-width="2.8" opacity=".82"><title>{label}: Q {h.qstart}-{h.qend}; S {h.sstart}-{h.send}; id {h.identity:.2f}%</title></line>')
    a.append(f'<text x="{W/2}" y="{H-7}" fill="#acc1d7" text-anchor="middle" font-size="13">QUERY — positions ({"nucleotides" if any(h.qframe for h in hsps) else "residues"})</text>')
    a.append(f'<text x="{W-220}" y="{H-7}" fill="#9eafc3" font-size="11">{len(shown)} / {len(hsps)} HSP</text></svg>')
    return ''.join(a)


def density_svg(hsps,qlength):
    if not hsps or qlength<=0:return None
    perref=collections.defaultdict(list)
    for h in hsps:perref[h.subject].append(h.qspan)
    groups=[merge_intervals(xs,qlength) for xs in perref.values()]
    n=min(250,qlength); data=[]
    for i in range(n):
        p=max(1,min(qlength,int((i+.5)*qlength/n)+1))
        num=sum(any(a<=p<=b for a,b in spans) for spans in groups)
        data.append(dict(center=p,value=num))
    return plot_svg(data,'References with an HSP along query positions','references')


def build_report_v5(out:Path,query:Seq,ref_seqs:list[Seq],hsps:list[HSP],mode:str,settings:dict,
                    history:list,msa_seqs:list[Seq]|None=None,dot:dict|None=None,sha_query='',sha_refs='',
                    motifs:list|None=None,orf_rows:list|None=None):
    report=build_report_v4(out,query,ref_seqs,hsps,mode,settings,history,msa_seqs,dot,sha_query,sha_refs,motifs,orf_rows)
    figures={}
    explanations={}
    refs_len={r.id:len(r.seq) for r in ref_seqs}
    if hsps:
        coverage=coverage_data(hsps,len(query.seq),refs_len)
        figures['16_hsp_count_per_reference']=bars_svg(collections.Counter(h.subject for h in hsps).most_common(35),'HSPs per reference')
        explanations['16_hsp_count_per_reference']='HSPs can overlap: this count does not represent independent regions.'
        figures['17_query_coverage_per_reference']=bars_svg([(r['subject'],round(r['query_coverage'],2)) for r in coverage[:35]],'Query coverage, union of HSPs (%)')
        explanations['17_query_coverage_per_reference']='Percentage of query positions covered by at least one HSP from the reference.'
        with_lengths=[(r['subject'],round(r['reference_coverage'],2)) for r in coverage if r['reference_coverage'] is not None]
        if with_lengths:
            figures['18_reference_coverage']=bars_svg(with_lengths[:35],'Reference coverage (%)')
            explanations['18_reference_coverage']='Coverage uses each reference\'s own coordinate system.'
        figures['19_bitscore_distribution']=hist_svg([h.bitscore for h in hsps],'HSP bitscores','Bitscore')
        figures['20_significance_strength']=hist_svg([-math.log10(max(h.evalue,1e-300)) for h in hsps],'Significance −log₁₀(E-value)','−log₁₀(E-value)')
        figures['21_hsp_gap_rate']=hist_svg([100*h.gap_columns/max(1,h.align_length) for h in hsps],'Columns with gaps per HSP (%)','Gaps / columns (%)')
        figures['22_identity_vs_bitscore']=scatter_svg([(h.identity,h.bitscore,h.subject) for h in hsps],'Local identity × bitscore','Local identity (%)','Bitscore',100,max(1,max(h.bitscore for h in hsps)))
        figures['23_hsp_mapping']=mapping_svg(hsps,len(query.seq),refs_len)
        d=density_svg(hsps,len(query.seq))
        if d:figures['24_reference_density_on_query']=d
        explanations['23_hsp_mapping']='The Y axis depends on each reference: compare segments without inferring a complete genome structure.'
        explanations['24_reference_density_on_query']='References with at least one HSP at the sampled position; not an abundance measure.'
    if msa_seqs:
        met=msa_metrics(msa_seqs)
        figures['25_gap_fraction_per_sequence']=bars_svg([(r.id,round(100*sum(x in '-.?' for x in r.seq)/max(1,len(r.seq)),2)) for r in msa_seqs[:35]],'Gaps per MSA sequence (%)')
        figures['26_msa_entropy_distribution']=hist_svg([p['entropy'] for p in met['positions']],'Distribution of MSA column entropy','Shannon (bits)')
        explanations['25_gap_fraction_per_sequence']='Gaps can reflect biology, different lengths or alignment artifacts.'
        explanations['26_msa_entropy_distribution']='Calculated from observed residues, excluding gaps and missing characters.'
    if not figures:
        figures['16_sequence_ambiguity']=hist_svg([s['n_pct'] for s in map(sequence_stats,[query,*ref_seqs])],'% N per sequence','N (%)')
    figure_dir=out/'figures';figure_dir.mkdir(exist_ok=True)

    for name,svg in list(figures.items()):
        figures[name]=str(svg)
        (figure_dir/f'{name}.svg').write_text(figures[name],encoding='utf-8')
    extra_cards=''.join(f'<article class="panel figcard" data-figure="{html.escape(k,quote=True)}"><div class="row"><h2>{html.escape(figure_title(figures[k],k))}</h2><button class="downloadSvg right" type="button" data-fig="{html.escape(k,quote=True)}">Download SVG</button></div><div class="figure-svg">{svg}</div><p class="hint">{html.escape(explanations.get(k,"Exploratory figure; interpretation depends on inputs and parameters."))}</p></article>' for k,svg in figures.items())
    page=report.read_text(encoding='utf-8')
    original_figures=page.count('class="panel figcard"')
    page=page.replace('</div></section><footer class="footer">',extra_cards+'</div></section><footer class="footer">',1)
    import re
    page=re.sub(r'\d+ vector SVG figures',f'{original_figures+len(figures)} vector SVG figures',page,count=1)
    # The v4 atlas immediately precedes the footer; rendering requires no network.
    page=page.replace('</nav>', '<button data-tab="exports">08 · Export data</button><button data-tab="guide">09 · Interpret</button></nav>',1)
    page=page.replace("figures:'Figure atlas'}", "figures:'Figure atlas',exports:'Export results',guide:'Methods and limitations'}",1)
    guide=f'''<section class="section" id="guide"><div class="panel"><h2>Understand the calculations</h2>
    <p>A local result is not whole-sequence identity. Two sequences can share a highly similar domain while their unaligned regions are long or divergent.</p>
    <details open><summary>BLAST / DIAMOND / MMseqs2</summary><p>Find local HSPs. Scores, E-values and identities describe the detected segments. E-values depend on the search space.</p></details>
    <details><summary>Coverage, overlap and undetected regions</summary><p>Query coverage is the union of intervals covered by each reference. Uncovered intervals have no detected HSP; they are not proof of non-homology.</p></details>
    <details><summary>Global and multiple alignment</summary><p>Needle performs global pairwise alignment. MAFFT, MUSCLE and Clustal Omega build MSAs and can produce different results. Alignment is not a homology search.</p></details>
    <details><summary>MSA statistics and phylogeny</summary><p>Variable sites, entropy and p-distances describe the calculated MSA, which must be biologically homogeneous. FastTree produces an exploratory tree; it does not replace phylogenetic validation.</p></details>
    <details><summary>BLASTX and coordinates</summary><p>BLASTX queries use nucleotides and references use amino acids, so HSP positions have different units. BED intervals always use the coordinates of their respective sequences.</p></details>
    <details><summary>Reproducibility</summary><p>Export parameters, checksums and logs. Figures alone cannot reproduce the search.</p></details>
    </div></section>'''
    exports='''<section id="exports" class="section"><div class="panel"><h2>Browser exports</h2>
    <p>The downloads below use data embedded in the report, which may be <strong>limited</strong> to the first 2,000 HSPs, 301 sequence-quality records and 1,000 MSA sequences. Use application exports for complete results.</p>
    <div class="row"><button data-export="hsp_csv">HSP · CSV</button><button data-export="hsp_json">HSP · JSON</button><button data-export="cov_tsv">Coverage · TSV</button><button data-export="qc_csv">Sequences · CSV</button><button data-export="msa_fasta">MSA visible · FASTA</button><button data-export="msa_positions">MSA sites · TSV</button><button data-export="provenance">Parameters · JSON</button><button id="exportPdf">Print / PDF</button></div>
    <p class="hint" id="exportStatus" aria-live="polite">All files are created locally in your browser without a server.</p>
    </div><div class="panel"><h2>Recommended formats by use</h2><table><thead><tr><th>Purpose</th><th>Format</th><th>Why</th></tr></thead><tbody>
    <tr><td>Figures for a paper</td><td>SVG</td><td>Vector graphics, independently exportable from the atlas</td></tr>
    <tr><td>Statistics for R / Python</td><td>TSV or CSV</td><td>Readable tabular data</td></tr>
    <tr><td>Aligned / uncovered regions</td><td>BED / GFF3</td><td>Explicitly documented coordinates</td></tr>
    <tr><td>MSA for other software</td><td>FASTA / NEXUS / PHYLIP</td><td>Interoperable alignment exchange</td></tr>
    <tr><td>Report to share</td><td>Standalone HTML</td><td>One page, without network or server</td></tr>
    </tbody></table></div></section>'''
    page=page.replace('<footer class="footer">AlignView',exports+guide+'<footer class="footer">AlignView',1)
    css='''<style>#guide details{padding:13px 16px;margin:10px 0;background:var(--side);border:1px solid var(--edge);border-radius:9px}#guide summary{cursor:pointer;font-weight:750}#guide details p{color:var(--muted);margin:13px 0 4px}#exports button{white-space:nowrap}button:focus-visible{outline:2px solid var(--mint)}@media print{#exports{display:none!important}}</style>'''
    script=r'''<script>
(function(){
const q=(x)=>String(x===null||x===undefined?'':x);
const csvCell=(x)=>'"'+q(x).replace(/"/g,'""')+'"';
const table=(rows,delimiter)=>{if(!rows.length)return '';const cols=Object.keys(rows[0]);return [cols.map(csvCell).join(delimiter),...rows.map(r=>cols.map(k=>csvCell(r[k])).join(delimiter))].join('\n')+'\n';};
function dl(content,name,type){let blob=new Blob([content],{type:type+';charset=utf-8'});let url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),2000);}
const actions={
 hsp_csv:()=>dl(table(D.hsps||[],','),'alignview_hsp_visible.csv','text/csv'),
 hsp_json:()=>dl(JSON.stringify(D.hsps||[],null,2),'alignview_hsp_visible.json','application/json'),
 cov_tsv:()=>dl(table(D.coverage||[],'\t'),'alignview_couvertures_visibles.tsv','text/tab-separated-values'),
 qc_csv:()=>dl(table(D.sequenceStats||[],','),'alignview_sequences_visibles.csv','text/csv'),
 msa_fasta:()=>{let msa=D.msa?.seqs||[];if(!msa.length)throw Error('No MSA available');dl(msa.map(x=>'>'+x.id+'\n'+x.seq+'\n').join(''),'alignview_msa_visible.fasta','text/plain');},
 msa_positions:()=>{let pos=D.msa?.metrics?.positions||[];if(!pos.length)throw Error('No MSA available');dl(table(pos,'\t'),'alignview_msa_positions.tsv','text/tab-separated-values');},
 provenance:()=>dl(JSON.stringify({settings:D.settings,history:D.history,shaQuery:D.shaQuery,shaRefs:D.shaRefs},null,2),'alignview_provenance.json','application/json')
};
document.querySelectorAll('[data-export]').forEach(b=>b.addEventListener('click',()=>{let status=document.getElementById('exportStatus');try{actions[b.dataset.export]();status.textContent='Export started: '+b.textContent;}catch(err){status.textContent=err.message;}}));
document.getElementById('exportPdf').onclick=()=>window.print();
})();
</script>'''
    page=page.replace('</head>',css+'</head>',1).replace('</body>',script+'</body>',1)
    page=page.replace('Sequence explorer · v4.0','Sequence explorer · v5.0')
    report.write_text(page,encoding='utf-8')
    return report
