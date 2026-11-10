"""Atlas de figures et centre d'exports intégrés au rapport HTML v5 autonome."""
from __future__ import annotations
import collections
import html
import math
from pathlib import Path
from av_core import Seq,HSP,coverage_data,msa_metrics,sequence_stats,merge_intervals
from av_report_v4 import build_report_v4, bars_svg, hist_svg, plot_svg, scatter_svg


def mapping_svg(hsps:list[HSP],query_length:int,ref_lengths:dict):
    """Cartographie des segments HSP en coordonnées propres à chaque séquence."""
    limit=70;W=1050;L=100;T=80;R=55;B=55;H=520
    shown=hsps[:limit]
    max_s=max((ref_lengths.get(h.subject) or max(h.sstart,h.send) for h in shown),default=1)
    max_s=max(1,max_s); qn=max(1,query_length)
    x=lambda p:L+(p-1)/qn*(W-L-R)
    y=lambda p:T+(p-1)/max_s*(H-T-B)
    a=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Carte de correspondance HSP">',
       '<rect width="1050" height="520" fill="#102039" rx="14"/>',
       '<text x="28" y="33" fill="#cfe0f2" font-size="17" font-family="Arial">Correspondances QUERY ↔ références (coordonnées originales)</text>',
       '<text x="28" y="54" fill="#91a8c3" font-size="12" font-family="Arial">Chaque trait représente un HSP, pas une colinéarité globale. Les références ont des tailles différentes.</text>',
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
    a.append(f'<text x="{W/2}" y="{H-7}" fill="#acc1d7" text-anchor="middle" font-size="13">QUERY — positions ({"nucléotides" if any(h.qframe for h in hsps) else "résidus"})</text>')
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
    return plot_svg(data,'Nombre de références avec HSP selon la position QUERY','références')


def build_report_v5(out:Path,query:Seq,ref_seqs:list[Seq],hsps:list[HSP],mode:str,settings:dict,
                    history:list,msa_seqs:list[Seq]|None=None,dot:dict|None=None,sha_query='',sha_refs='',
                    motifs:list|None=None,orf_rows:list|None=None):
    report=build_report_v4(out,query,ref_seqs,hsps,mode,settings,history,msa_seqs,dot,sha_query,sha_refs,motifs,orf_rows)
    figures={}
    explanations={}
    refs_len={r.id:len(r.seq) for r in ref_seqs}
    if hsps:
        coverage=coverage_data(hsps,len(query.seq),refs_len)
        figures['16_nombre_hsp_par_reference']=bars_svg(collections.Counter(h.subject for h in hsps).most_common(35),'HSP par référence')
        explanations['16_nombre_hsp_par_reference']='Plusieurs HSP peuvent se chevaucher : ce nombre ne représente pas des régions indépendantes.'
        figures['17_couverture_query_par_reference']=bars_svg([(r['subject'],round(r['query_coverage'],2)) for r in coverage[:35]],'Couverture QUERY, union des HSP (%)')
        explanations['17_couverture_query_par_reference']='Pourcentage de positions QUERY touchées par au moins un HSP de la référence.'
        with_lengths=[(r['subject'],round(r['reference_coverage'],2)) for r in coverage if r['reference_coverage'] is not None]
        if with_lengths:
            figures['18_couverture_references']=bars_svg(with_lengths[:35],'Couverture dans les références (%)')
            explanations['18_couverture_references']='Couverture calculée dans les coordonnées propres à chaque référence.'
        figures['19_distribution_bitscore']=hist_svg([h.bitscore for h in hsps],'Bitscore des HSP','Bitscore')
        figures['20_force_significativite']=hist_svg([-math.log10(max(h.evalue,1e-300)) for h in hsps],'Significativité −log₁₀(E-value)','−log₁₀(E-value)')
        figures['21_taux_gaps_hsp']=hist_svg([100*h.gap_columns/max(1,h.align_length) for h in hsps],'% colonnes avec gap par HSP','Gaps / colonnes (%)')
        figures['22_identites_vs_bitscore']=scatter_svg([(h.identity,h.bitscore,h.subject) for h in hsps],'Identité locale × bitscore','Identité locale (%)','Bitscore',100,max(1,max(h.bitscore for h in hsps)))
        figures['23_cartographie_hsp']=mapping_svg(hsps,len(query.seq),refs_len)
        d=density_svg(hsps,len(query.seq))
        if d:figures['24_densite_references_query']=d
        explanations['23_cartographie_hsp']='Le repère Y dépend de chaque référence : comparer les segments, pas inférer une structure génomique complète.'
        explanations['24_densite_references_query']='Références ayant au moins un HSP à la position échantillonnée; ne mesure pas l’abondance.'
    if msa_seqs:
        met=msa_metrics(msa_seqs)
        figures['25_fraction_gaps_par_sequence']=bars_svg([(r.id,round(100*sum(x in '-.?' for x in r.seq)/max(1,len(r.seq)),2)) for r in msa_seqs[:35]],'Gaps par séquence du MSA (%)')
        figures['26_distribution_entropie_msa']=hist_svg([p['entropy'] for p in met['positions']],'Distribution de l’entropie des colonnes du MSA','Shannon (bits)')
        explanations['25_fraction_gaps_par_sequence']='Les gaps peuvent refléter de la biologie, des longueurs différentes ou des artefacts d’alignement.'
        explanations['26_distribution_entropie_msa']='Calcul sur résidus observés, gaps et caractères manquants exclus.'
    if not figures:
        figures['16_ambiguite_sequences']=hist_svg([s['n_pct'] for s in map(sequence_stats,[query,*ref_seqs])],'% N par séquence','N (%)')
    figure_dir=out/'figures';figure_dir.mkdir(exist_ok=True)
    for name,svg in figures.items():(figure_dir/f'{name}.svg').write_text(svg,encoding='utf-8')
    extra_cards=''.join(f'<article class="panel figcard" data-figure="{html.escape(k,quote=True)}"><div class="row"><h2>{html.escape(k[3:].replace("_"," ").capitalize())}</h2><button class="downloadSvg right" type="button" data-fig="{html.escape(k,quote=True)}">Télécharger SVG</button></div><div class="figure-svg">{svg}</div><p class="hint">{html.escape(explanations.get(k,"Figure exploratoire ; interprétation conditionnée par les entrées et paramètres."))}</p></article>' for k,svg in figures.items())
    page=report.read_text(encoding='utf-8')
    original_figures=page.count('class="panel figcard"')
    page=page.replace('</div></section><footer class="footer">',extra_cards+'</div></section><footer class="footer">',1)
    import re
    page=re.sub(r'\d+ figures SVG vectorielles',f'{original_figures+len(figures)} figures SVG vectorielles',page,count=1)
    # L'atlas v4 est immédiatement suivi du footer ; construction sans réseau.
    page=page.replace('</nav>', '<button data-tab="exports">08 · Exporter les données</button><button data-tab="guide">09 · Interpréter</button></nav>',1)
    page=page.replace("figures:'Atlas des figures'}", "figures:'Atlas des figures',exports:'Exporter les résultats',guide:'Méthodes & limites'}",1)
    guide=f'''<section class="section" id="guide"><div class="panel"><h2>Comprendre ce qui est calculé</h2>
    <p>Un résultat local n'est pas une identité globale. Deux séquences peuvent partager un domaine très similaire alors que leurs régions non alignées sont longues ou divergentes.</p>
    <details open><summary>BLAST / DIAMOND / MMseqs2</summary><p>Localisent des HSP. Les scores, E-values et identités concernent les segments détectés. Les E-values dépendent de l’espace de recherche.</p></details>
    <details><summary>Couverture, chevauchement et zones non détectées</summary><p>La couverture QUERY est calculée par union des intervalles couverts par référence. Les intervalles non couverts représentent des positions sans HSP détecté, pas une preuve de non-homologie.</p></details>
    <details><summary>Alignement global et multiple</summary><p>Needle réalise un alignement global pairwise. MAFFT, MUSCLE et Clustal Omega construisent des MSA et peuvent différer. L’alignement n’est pas une recherche d’homologues.</p></details>
    <details><summary>Statistiques MSA et phylogénie</summary><p>Sites variables, entropie et distances p décrivent le MSA calculé, qui doit être biologiquement homogène. FastTree produit un arbre exploratoire : il ne remplace pas une validation phylogénétique.</p></details>
    <details><summary>BLASTX et coordonnées</summary><p>QUERY est nucléotidique et les références sont protéiques ; les positions HSP sont donc exprimées dans des unités distinctes. Les BED concernent toujours les coordonnées de leur séquence respective.</p></details>
    <details><summary>Reproductibilité</summary><p>Exporter les paramètres, les checksums et les logs. Les figures seules ne suffisent pas à reproduire la recherche.</p></details>
    </div></section>'''
    exports='''<section id="exports" class="section"><div class="panel"><h2>Exports depuis le navigateur</h2>
    <p>Les exports ci-dessous utilisent les données intégrées au rapport : elles peuvent être <strong>limitées</strong> aux 2 000 premiers HSP, 301 statistiques de séquences et 1 000 séquences de MSA. Pour les résultats complets, utiliser le centre d’export du logiciel.</p>
    <div class="row"><button data-export="hsp_csv">HSP · CSV</button><button data-export="hsp_json">HSP · JSON</button><button data-export="cov_tsv">Couverture · TSV</button><button data-export="qc_csv">Séquences · CSV</button><button data-export="msa_fasta">MSA visible · FASTA</button><button data-export="msa_positions">Sites du MSA · TSV</button><button data-export="provenance">Paramètres · JSON</button><button id="exportPdf">Imprimer / PDF</button></div>
    <p class="hint" id="exportStatus" aria-live="polite">Tous les fichiers sont créés localement dans ton navigateur, sans serveur.</p>
    </div><div class="panel"><h2>Format recommandé selon l’usage</h2><table><thead><tr><th>Besoin</th><th>Format</th><th>Pourquoi</th></tr></thead><tbody>
    <tr><td>Figures pour article</td><td>SVG</td><td>Vectoriel, exportable indépendamment depuis l'atlas</td></tr>
    <tr><td>Statistiques pour R / Python</td><td>TSV ou CSV</td><td>Données tabulaires lisibles</td></tr>
    <tr><td>Régions alignées / non couvertes</td><td>BED / GFF3</td><td>Coordonnées clairement documentées</td></tr>
    <tr><td>MSA vers d’autres logiciels</td><td>FASTA / NEXUS / PHYLIP</td><td>Compatibilité interopérable après alignement</td></tr>
    <tr><td>Rapport à partager</td><td>HTML autonome</td><td>Une seule page, sans réseau ni serveur</td></tr>
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
 msa_fasta:()=>{let msa=D.msa?.seqs||[];if(!msa.length)throw Error('Aucun MSA disponible');dl(msa.map(x=>'>'+x.id+'\n'+x.seq+'\n').join(''),'alignview_msa_visible.fasta','text/plain');},
 msa_positions:()=>{let pos=D.msa?.metrics?.positions||[];if(!pos.length)throw Error('Aucun MSA disponible');dl(table(pos,'\t'),'alignview_msa_positions.tsv','text/tab-separated-values');},
 provenance:()=>dl(JSON.stringify({settings:D.settings,history:D.history,shaQuery:D.shaQuery,shaRefs:D.shaRefs},null,2),'alignview_provenance.json','application/json')
};
document.querySelectorAll('[data-export]').forEach(b=>b.addEventListener('click',()=>{let status=document.getElementById('exportStatus');try{actions[b.dataset.export]();status.textContent='Export lancé : '+b.textContent;}catch(err){status.textContent=err.message;}}));
document.getElementById('exportPdf').onclick=()=>window.print();
})();
</script>'''
    page=page.replace('</head>',css+'</head>',1).replace('</body>',script+'</body>',1)
    page=page.replace('Sequence explorer · v4.0','Sequence explorer · v5.0')
    report.write_text(page,encoding='utf-8')
    return report
