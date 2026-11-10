"""Rapport enrichi v4 : figures SVG autonomes sans JavaScript externe.
Les données de recherche sont déjà filtrées par le moteur ; les figures indiquent
explicitement quand une statistique est locale (HSP) et non globale.
"""
from __future__ import annotations
import collections
import html
import math
import statistics
from pathlib import Path
from av_core import Seq,HSP,sequence_stats,msa_metrics,distance_matrix,sliding_scores,coverage_data
from av_units import units_for
from av_visual import build_report,plot_svg

INK='#cfe0f2'; MUTED='#8da5c1'; BLUE='#6d9ef6'; TEAL='#54d2b5'; ORANGE='#ffc47f'; RED='#f48d9b'

def txt(s):return html.escape(str(s),quote=True)

def shell(title,xmax=780,height=285):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {xmax} {height}" role="img" aria-label="{txt(title)}">',
            f'<rect width="{xmax}" height="{height}" rx="14" fill="#102039"/>',
            f'<text x="24" y="29" fill="{INK}" font-family="Arial,sans-serif" font-size="15" font-weight="bold">{txt(title)}</text>']

def finish(parts):parts.append('</svg>');return ''.join(parts)

def hist_svg(values,title,x_label,bins=14,percent=False):
    values=[float(x) for x in values if x is not None and math.isfinite(float(x))]
    parts=shell(title)
    if not values:
        parts.append(f'<text x="32" y="145" fill="{MUTED}" font-size="15">Aucune observation</text>');return finish(parts)
    minimum=min(values);maximum=max(values)
    if maximum==minimum:maximum=minimum+1
    bins=max(1,min(bins,30)); counts=[0]*bins
    for v in values:counts[min(bins-1,int((v-minimum)/(maximum-minimum)*bins))]+=1
    chartx=55;charty=60;w=690;h=165;peak=max(counts)
    for g in range(5):
        y=charty+h-g*h/4
        parts.append(f'<line x1="{chartx}" y1="{y:.1f}" x2="{chartx+w}" y2="{y:.1f}" stroke="#334761" stroke-width=".7"/>')
        parts.append(f'<text x="{chartx-9}" y="{y+4:.1f}" text-anchor="end" font-family="sans-serif" font-size="10" fill="{MUTED}">{round(g*peak/4)}</text>')
    for i,n in enumerate(counts):
        bw=w/bins-3;x=chartx+i*w/bins;bh=(n/max(1,peak))*h
        parts.append(f'<rect x="{x:.2f}" y="{charty+h-bh:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="2" fill="{TEAL}"><title>{n} observations · {minimum+i*(maximum-minimum)/bins:.2f} à {minimum+(i+1)*(maximum-minimum)/bins:.2f}</title></rect>')
    parts.extend([f'<text x="55" y="248" fill="{MUTED}" font-size="11" font-family="sans-serif">{minimum:,.2f}</text>',
                  f'<text x="745" y="248" text-anchor="end" fill="{MUTED}" font-family="sans-serif" font-size="11">{max(values):,.2f}</text>',
                  f'<text x="390" y="269" fill="{INK}" text-anchor="middle" font-family="sans-serif" font-size="11">{txt(x_label)} · n={len(values)}</text>'])
    return finish(parts)

def bars_svg(items,title,nmax=18):
    items=sorted(items,key=lambda v:v[1],reverse=True)[:nmax]
    height=max(180,len(items)*31+65);parts=shell(title,height=height)
    if not items:parts.append(f'<text x="30" y="85" fill="{MUTED}">Aucune donnée</text>');return finish(parts)
    peak=max(float(k) for _,k in items) or 1
    for i,(name,val) in enumerate(items):
        y=46+i*31
        parts.append(f'<text x="20" y="{y+13}" fill="{INK}" font-family="monospace" font-size="11">{txt(str(name)[:20])}</text>')
        parts.append(f'<rect x="185" y="{y}" width="{max(2,float(val)/peak*500):.2f}" height="18" rx="3" fill="{BLUE}"><title>{txt(name)} : {val}</title></rect>')
        parts.append(f'<text x="700" y="{y+14}" fill="{MUTED}" font-size="12" font-family="monospace">{float(val):,.1f}</text>')
    return finish(parts)

def scatter_svg(points,title,xlabel,ylabel,xmax=100,ymax=100):
    parts=shell(title);l=62;t=48;w=670;h=174
    for i in range(6):
        x=l+i*w/5;y=t+h-i*h/5
        parts.append(f'<path d="M{x:.1f} {t}V{t+h} M{l} {y:.1f}H{l+w}" stroke="#32455f" fill="none"/>')
        parts.append(f'<text x="{x:.1f}" y="{t+h+16}" fill="{MUTED}" font-size="10" text-anchor="middle">{i*xmax/5:.0f}</text>')
        parts.append(f'<text x="{l-10}" y="{y+4:.1f}" fill="{MUTED}" font-size="10" text-anchor="end">{i*ymax/5:.0f}</text>')
    for x,y,label in points[:1500]:
        if x is None or y is None:continue
        px=l+max(0,min(1,float(x)/xmax))*w;py=t+h-max(0,min(1,float(y)/ymax))*h
        parts.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="{TEAL}" fill-opacity=".73" stroke="#0c1c2f"><title>{txt(label)} · {x:.1f} / {y:.1f}</title></circle>')
    parts.extend([f'<text x="390" y="270" text-anchor="middle" fill="{INK}" font-size="11">{txt(xlabel)}</text>',
                  f'<text transform="translate(15,150) rotate(-90)" fill="{INK}" font-size="11">{txt(ylabel)}</text>'])
    return finish(parts)

def heatmap_svg(names,mat,title):
    n=min(len(names),45)
    if n<2:return bars_svg([],title)
    s=max(10,min(26,560/n)); margin=145; side=n*s; width=margin+side+95;height=margin+side+35
    parts=shell(title,xmax=width,height=height)
    for i in range(n):
        y=margin+i*s;parts.append(f'<text x="{margin-5}" y="{y+s*.7:.1f}" text-anchor="end" fill="{INK}" font-size="{max(8,min(12,s*.65))}" font-family="monospace">{txt(names[i][:18])}</text>')
    for j in range(n):
        x=margin+j*s;parts.append(f'<text transform="translate({x+s*.6:.1f},{margin-7}) rotate(-60)" fill="{MUTED}" font-size="{max(8,min(12,s*.65))}" font-family="monospace">{txt(names[j][:15])}</text>')
    for i in range(n):
        for j in range(n):
            v=mat[i][j]
            frac=max(0,min(1,float(v))) if v is not None else None
            fill=(f'rgb({int(34+frac*202)},{int(181-frac*99)},{int(153-frac*83)})' if frac is not None else '#28354b')
            parts.append(f'<rect x="{margin+j*s:.1f}" y="{margin+i*s:.1f}" width="{s:.1f}" height="{s:.1f}" fill="{fill}" stroke="#13263d" stroke-width=".4"><title>{txt(names[i])} / {txt(names[j])}: {"NA" if v is None else f"{v:.4f}"}</title></rect>')
    parts.append(f'<text x="{margin}" y="{height-10}" fill="{MUTED}" font-size="11">Distance p · gaps exclus · {n} séquences max</text>')
    return finish(parts)

def annotation_svg(length,covered,title='Régions non couvertes (QUERY)'):
    parts=shell(title,height=173)
    if length<=0:return finish(parts)
    parts.append(f'<rect x="45" y="66" width="690" height="26" fill="#3b4356" rx="5"/>')
    for a,b in covered:
        x=45+(a-1)/length*690;w=(b-a+1)/length*690
        parts.append(f'<rect x="{x:.2f}" y="66" width="{w:.2f}" height="26" fill="{TEAL}" rx="2"><title>HSP couvert {a}–{b}</title></rect>')
    parts.append(f'<text x="45" y="120" fill="{MUTED}" font-size="12">1</text><text x="735" y="120" text-anchor="end" fill="{MUTED}" font-size="12">{length:,}</text>')
    parts.append(f'<text x="45" y="146" fill="{MUTED}" font-size="12">Vert : au moins un HSP · Gris : aucune région détectée par la recherche</text>')
    return finish(parts)

def build_report_v4(out:Path,query:Seq,ref_seqs:list[Seq],hsps:list[HSP],mode:str,settings:dict,
                    history:list,msa_seqs:list[Seq]|None=None,dot:dict|None=None,sha_query='',sha_refs='',
                    motifs:list|None=None,orf_rows:list|None=None):
    report=build_report(out,query,ref_seqs,hsps,mode,settings,history,msa_seqs,dot,sha_query,sha_refs)
    figures={}
    allseq=[query,*ref_seqs]
    stats=[sequence_stats(s) for s in allseq]
    summary=collections.Counter(query.seq.replace('-','').replace('.',''))
    figures['01_composition_query']=bars_svg(summary.most_common(18),'Composition de la séquence QUERY')
    u=units_for(mode,sequence_stats(query)['kind'])
    if u.query==u.subject:
        figures['02_longueurs_sequences']=hist_svg([x['ungapped'] for x in stats],
            'Distribution des longueurs',f'Longueur ({u.query}, sans gaps)')
    else:
        # Les longueurs ADN/protéine ne sont pas comparables sur un axe commun.
        if ref_seqs:
            figures['02b_longueurs_references']=hist_svg([sequence_stats(x)['ungapped'] for x in ref_seqs],
                'Distribution des longueurs des références',f'RÉFÉRENCE ({u.subject}, sans gaps)')
    nt_stats=[x for x in stats if x['kind']=='ADN/ARN' and x['gc_pct'] is not None]
    if nt_stats:
        figures['03_qualite_gc']=scatter_svg([(x['ungapped'],x['gc_pct'],x['id']) for x in nt_stats],
            'Longueur nucléotidique / GC (%)','Longueur (nt)','GC (%)',max(1,max(x['ungapped'] for x in nt_stats)),100)
    if hsps:
        figures['04_identite_couverture']=scatter_svg([(h.identity,100*(h.qspan[1]-h.qspan[0]+1)/max(1,len(query.seq)),h.subject) for h in hsps],
            'Identité locale / étendue HSP sur QUERY','Identité HSP (%)','Étendue QUERY (%)')
        figures['05_hist_identite']=hist_svg([h.identity for h in hsps],'Identités des HSP (%)','Identité locale (%)')
        figures['06_hist_hsp']=hist_svg([h.align_length for h in hsps],'Longueurs alignées (HSP)',f'Colonnes alignées ({u.aligned})')
        frames=collections.Counter(str(h.qframe) for h in hsps if h.qframe)
        if frames:figures['07_cadres']=bars_svg(sorted(frames.items()),'Cadres de lecture QUERY (HSP)')
        covered=[]
        from av_core import merge_intervals
        covered=merge_intervals([h.qspan for h in hsps],len(query.seq))
        figures['08_segments']=annotation_svg(len(query.seq),covered,f'Régions sans HSP (QUERY, {u.query})')
    if msa_seqs:
        m=msa_metrics(msa_seqs)
        for suffix,title,key in [('09_conservation','Conservation MSA (%)','conservation'),('10_entropie','Entropie de Shannon / colonne','entropy'),('11_gaps','Colonnes avec gaps (%)','gap_pct')]:
            figures[suffix]=plot_svg([dict(center=x['position'],value=x[key]) for x in m['positions']],title,key)
        if 2<=len(msa_seqs)<=45:
            mat=distance_matrix(msa_seqs)
            figures['12_distance_p']=heatmap_svg([s.id for s in msa_seqs],mat,'Distances p entre séquences (MSA)')
        counts=collections.Counter(str(round(p['gap_pct']//10)*10)+'–'+str(min(100,round(p['gap_pct']//10)*10+10))+' %' for p in m['positions'])
        figures['13_gaps_distribution']=bars_svg(counts.items(),'Distribution des taux de gaps par colonne')
    if motifs:figures['14_motifs']=bars_svg(collections.Counter(r.get('seqid','?') for r in motifs).items(),'Occurrences de motifs par séquence')
    if orf_rows:figures['15_orfs']=hist_svg([r['length_aa'] for r in orf_rows],'Longueurs des ORF détectées','Acides aminés')
    directory=out/'figures';directory.mkdir(exist_ok=True)
    for name,svg in figures.items():(directory/(name+'.svg')).write_text(svg,encoding='utf-8')
    cards=''.join(f'<article class="panel figcard" data-figure="{txt(k)}"><div class="row"><h2>{txt(k[3:].replace("_"," ").capitalize())}</h2><button type="button" class="downloadSvg right" data-fig="{txt(k)}" aria-label="Exporter SVG">Télécharger SVG</button></div><div class="figure-svg">{svg}</div></article>' for k,svg in figures.items())
    section=f'''<section id="figures" class="section"><div class="panel"><h2>Atlas des figures</h2>
<p>Les figures sont calculées à partir des données de la session. Les statistiques des HSP représentent des <strong>alignements locaux</strong> et non l’identité globale des génomes ou protéines. Les données MSA sont disponibles après alignement.</p>
<div class="row"><label>Filtrer les figures <input id="figFilter" type="search" placeholder="Ex. identité, conservation, gaps, GC" style="min-width:280px"></label><span class="pill">{len(figures)} figures SVG vectorielles</span><button id="printFigures" type="button">Imprimer / PDF</button></div></div>
<div class="atlas">{cards}</div></section>'''
    page=report.read_text(encoding='utf-8')
    page=page.replace('</nav>', '<button data-tab="figures">07 · Atlas des figures</button></nav>',1)
    page=page.replace('<footer class="footer">',section+'<footer class="footer">',1)
    page=page.replace("provenance:'Reproductibilité'}", "provenance:'Reproductibilité',figures:'Atlas des figures'}",1)
    css='''<style>.atlas{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,510px),1fr));gap:16px}.figcard{margin-bottom:0}.figcard svg{max-height:480px;display:block}.figcard h2{font-size:15px;margin:0}.figure-svg{width:100%;overflow:auto}@media print{.sidebar,.top button,.downloadSvg,#figFilter,#printFigures{display:none!important}.layout{display:block}.section{display:none!important}#figures{display:block!important}.figcard{break-inside:avoid;min-height:200px;page-break-inside:avoid}.atlas{grid-template-columns:1fr 1fr}}</style>'''
    js='''<script>document.querySelectorAll('.downloadSvg').forEach(btn=>btn.addEventListener('click',()=>{const svg=btn.closest('.figcard').querySelector('svg');const clone=svg.cloneNode(true);clone.setAttribute('xmlns','http://www.w3.org/2000/svg');const text=new XMLSerializer().serializeToString(clone);const url=URL.createObjectURL(new Blob([text],{type:'image/svg+xml;charset=utf-8'}));const a=document.createElement('a');a.href=url;a.download=btn.dataset.fig+'.svg';a.click();setTimeout(()=>URL.revokeObjectURL(url),5000);}));document.getElementById('figFilter').oninput=e=>{const needle=e.target.value.trim().toLowerCase();document.querySelectorAll('.figcard').forEach(card=>card.hidden=!card.dataset.figure.toLowerCase().includes(needle));};document.getElementById('printFigures').onclick=()=>window.print();</script>'''
    page=page.replace('</head>',css+'</head>',1).replace('</body>',js+'</body>',1)
    report.write_text(page,encoding='utf-8')
    return report
