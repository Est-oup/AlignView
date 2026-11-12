"""Atlas v6: vues complémentaires conditionnelles, unités natives et navigation ciblée.

Les graphiques sont descriptifs ; pas de statistiques inférentielles sur HSP non indépendants.
"""
from __future__ import annotations
import collections
import html
import json
import math
import re
from pathlib import Path
from av_i18n import figure_title
from av_core import Seq,HSP,coverage_data,merge_intervals,interval_length,msa_metrics,sequence_stats
from av_units import units_for,coord_audit
from av_report_v5 import build_report_v5
from av_tree import TREE_FIGURE, tree_svg
from av_report_v4 import bars_svg,hist_svg,scatter_svg,plot_svg,shell,finish,txt,INK,MUTED,BLUE,TEAL


def tracks_svg(hsps:list[HSP],qlen:int,rlen:dict,mode:str,max_refs=20):
    """Pistes indépendantes : une échelle native pour QUERY, une pour chaque référence.

    Ne dessine aucune diagonale nt/aa qui laisserait penser à une échelle commune.
    """
    u=units_for(mode);groups=collections.defaultdict(list)
    for h in hsps:groups[h.subject].append(h)
    groups=dict(list(groups.items())[:max_refs]);width=1060
    height=max(215,136+len(groups)*78)
    left=205;right=45;usable=width-left-right
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="Pistes coordonnées HSP">',
         f'<rect width="{width}" height="{height}" rx="12" fill="#102039"/>',
         '<text x="25" y="30" fill="#d5e5f9" font-size="17">Couverture par référence · deux coordonnées indépendantes</text>',
         f'<text x="25" y="51" fill="#9fb4ca" font-size="12">QUERY : {u.query} · RÉFÉRENCE : {u.subject} · alignement HSP : {u.aligned} · intervalles 1-based inclusifs</text>']
    for i,(subject,rows) in enumerate(groups.items()):
        top=82+i*78; slen=max(1,rlen.get(subject, max(max(h.sspan) for h in rows)))
        out.append(f'<text x="20" y="{top+11}" fill="#d5e5f9" font-size="12">{txt(subject[:26])}</text>')
        for side,length,spans,y,color in [
            ('Q',max(1,qlen),merge_intervals([h.qspan for h in rows],qlen),top,TEAL),
            ('S',slen,merge_intervals([h.sspan for h in rows],slen),top+28,BLUE)]:
            out.extend([f'<text x="{left-22}" y="{y+10}" text-anchor="end" fill="#9fb4ca" font-size="11">{side}</text>',
                f'<rect x="{left}" y="{y}" width="{usable}" height="13" fill="#33455e" rx="4"/>'])
            for a,b in spans:
                x=left+(a-1)/length*usable; w=(b-a+1)/length*usable
                out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(1,w):.1f}" height="13" fill="{color}" rx="3"><title>{side} {a}–{b} / {length} {u.query if side=="Q" else u.subject}</title></rect>')
            out.append(f'<text x="{width-11}" y="{y+10}" text-anchor="end" fill="#e1eaf7" font-size="11">{length:,} {u.query if side=="Q" else u.subject}</text>')
    out.append(f'<text x="25" y="{height-12}" fill="#9fb4ca" font-size="11">Vert = positions QUERY couvertes · bleu = positions référence couvertes · gris = sans HSP détecté</text>')
    out.append('</svg>');return ''.join(out)


def frame_tracks_svg(hsps:list[HSP],qlen:int):
    frames=sorted({h.qframe for h in hsps if h.qframe},key=lambda v:(v<0,abs(v)))
    if not frames:return None
    w=1000;height=90+len(frames)*49;l=115;usable=810
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {height}">',f'<rect width="{w}" height="{height}" rx="12" fill="#102039"/>',
        '<text x="22" y="30" font-size="16" fill="#dae8f5">Cadres QUERY · pistes de HSP (positions nt)</text>']
    for i,fr in enumerate(frames):
        y=55+i*49
        out.append(f'<text x="20" y="{y+10}" fill="#c5d4e6" font-size="12">Cadre {fr:+d}</text>')
        out.append(f'<rect x="{l}" y="{y}" width="{usable}" height="19" fill="#304058" rx="4"/>')
        for h in (h for h in hsps if h.qframe==fr):
            a,b=h.qspan;x=l+(a-1)/max(1,qlen)*usable;length=(b-a+1)/max(1,qlen)*usable
            out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(2,length):.1f}" height="19" fill="{"#64d9bd" if fr>0 else "#f6ad70"}" opacity=".75"><title>{txt(h.subject)} · {a}–{b} nt · cadre {fr:+d}</title></rect>')
        out.append(f'<text x="{l}" y="{y+34}" fill="#91a8c3" font-size="10">1</text><text x="{l+usable}" y="{y+34}" text-anchor="end" fill="#91a8c3" font-size="10">{qlen:,} nt</text>')
    out.append('</svg>');return ''.join(out)


def atlas_additions(query:Seq,ref_seqs:list[Seq],hsps:list[HSP],msa_seqs:list[Seq]|None,mode:str):
    u=units_for(mode,sequence_stats(query)['kind']); figures={};meta={}
    def put(key,svg,group,description,kind='HSP'):
        if svg:
            figures[key]=svg;meta[key]={'group':group,'description':description,'source':kind}
    ref_lengths={r.id:len(r.seq) for r in ref_seqs}
    if hsps:
        spans_q=[h.qspan[1]-h.qspan[0]+1 for h in hsps]
        spans_s=[h.sspan[1]-h.sspan[0]+1 for h in hsps]
        put('27_etendue_query_hsp',hist_svg(spans_q,'Étendue de chaque HSP sur QUERY',f'Positions QUERY ({u.query})'),'HSP et coordonnées','Coordonnées sur QUERY dans son alphabet d’origine, jamais des colonnes alignées.')
        put('28_etendue_reference_hsp',hist_svg(spans_s,'Étendue de chaque HSP sur sa référence',f'Positions RÉFÉRENCE ({u.subject})'),'HSP et coordonnées','La longueur du segment est exprimée dans l’unité propre aux références.')
        put('29_pourcentage_hsp_query',hist_svg([100*n/max(1,len(query.seq)) for n in spans_q], 'Étendue de chaque HSP / QUERY (%)','Étendue relative QUERY (%)'),'HSP et coordonnées','Chaque HSP est mesuré séparément, les chevauchements ne sont pas soustraits.')
        put('30_gaps_nombres',hist_svg([h.gap_columns for h in hsps],'Colonnes avec gap par HSP','Colonnes alignées ('+u.aligned+')'),'Qualité des HSP','Nombre de colonnes contenant au moins un gap ; en BLASTX, les colonnes représentent des acides aminés.')
        put('31_substitutions_hsp',hist_svg([h.mismatches for h in hsps],'Substitutions par HSP','Substitutions (colonnes '+u.aligned+')'),'Qualité des HSP','Comparaison des séquences alignées ; les gaps ne sont pas comptés comme substitutions.')
        put('32_substitutions_vs_gaps',scatter_svg([(h.mismatches,h.gap_columns,h.subject) for h in hsps], 'Substitutions et gaps des HSP', 'Substitutions (colonnes '+u.aligned+')','Colonnes avec gap',max(1,max(h.mismatches for h in hsps)),max(1,max(h.gap_columns for h in hsps))),'Qualité des HSP','Les valeurs dépendent de la longueur des alignements locaux ; comparaison descriptive.')
        cov=coverage_data(hsps,len(query.seq),ref_lengths)
        both=[r for r in cov if r['reference_coverage'] is not None]
        if both:put('33_couverture_relative_q_vs_s',scatter_svg([(r['query_coverage'],r['reference_coverage'],r['subject']) for r in both],'Couverture relative des deux séquences','Couverture QUERY (%)','Couverture RÉFÉRENCE (%)'),'Couvertures','Pourcentages calculés sur des dénominateurs distincts ; ne pas les convertir en longueurs comparables.')
        accum=[];allspans=[]
        for i,r in enumerate(cov,1):
            allspans.extend(r['query_segments']);current=merge_intervals(allspans,len(query.seq))
            accum.append(dict(center=i,value=100*interval_length(current)/max(1,len(query.seq))))
        if len(accum)>1:put('34_couverture_cumulative_hits',plot_svg(accum,'Couverture cumulée QUERY, ajout des références','QUERY couverte (%)'),'Couvertures','Union des positions couvertes après ajout successif des références classées par bitscore ; ce n’est pas une courbe de saturation d’échantillonnage.')
        if u.translated_query:
            put('36_cadres_query_couverture',frame_tracks_svg(hsps,len(query.seq)),'Traduction','Les cadres +/− sont représentés sur les coordonnées nucléotidiques originales de la requête.')
        if u.translated_subject:
            subject_frames=collections.Counter(str(h.sframe) for h in hsps if h.sframe)
            if subject_frames:put('37_cadres_reference',bars_svg(sorted(subject_frames.items()),'Cadres de lecture des RÉFÉRENCES'),'Traduction','Cadres définis dans les nucléotides de la référence, indépendamment des coordonnées de QUERY.')
        # Diagnostic de correspondance 3 nt ↔ 1 aa, sans considérer une colonne comme une base.
        if u.translated_query or u.translated_subject:
            audit=coord_audit(hsps,mode)
            def scale(r,which):
                a=r[f'{which}_span']; b=r[f'{which}_aligned_chars'];return round(a/b,3) if b else None
            if u.translated_query:put('38_ratio_nt_aa_query',hist_svg([scale(r,'query') for r in audit if r['query_aligned_chars']], 'Ratio étendue QUERY nt / résidus alignés','nt par aa aligné'),'Traduction','Contrôle qualité : ratio théorique de 3 nt/aa ; un écart peut indiquer un HSP atypique ou un problème de coordonnées.')
            if u.translated_subject:put('39_ratio_nt_aa_reference',hist_svg([scale(r,'reference') for r in audit if r['reference_aligned_chars']], 'Ratio étendue RÉFÉRENCE nt / résidus alignés','nt par aa aligné'),'Traduction','Contrôle qualité : ratio théorique de 3 nt/aa ; ne pas mélanger les deux types de positions.')
    if msa_seqs and len(msa_seqs)>=2:
        m=msa_metrics(msa_seqs);pos=m['positions'];length=m['length'];n=m['n']
        put('40_occupation_msa',plot_svg([dict(center=p['position'],value=100-p['gap_pct']) for p in pos], 'Occupancy du MSA','Séquences sans gap (%)'),'MSA','Position exprimée en colonnes du MSA ('+u.msa+'), pas en positions natives de QUERY.','MSA')
        put('41_entropie_vs_gaps',scatter_svg([(p['gap_pct'],p['entropy'],str(p['position'])) for p in pos],'Variabilité en fonction des gaps','Gaps (%)','Entropie (bits)',100,max(1,max(p['entropy'] for p in pos))),'MSA','Une entropie élevée sur des colonnes peu occupées doit être interprétée avec prudence.','MSA')
        put('42_consensus_composition',bars_svg(collections.Counter(m['consensus'].replace('-','')).most_common(25),'Composition du consensus MSA'),'MSA','Consensus calculé par colonne ; n’est pas une séquence biologique observée.','MSA')
        put('43_conservation_distribution',hist_svg([p['conservation'] for p in pos],'Distribution de la conservation du MSA','Conservation (%)'),'MSA','Conservation calculée sur les résidus non-gap, dépend des séquences retenues.','MSA')
        win=min(max(10,length//30),100);half=win//2
        smooth=[]
        if length>1:
            step=max(1,length//180)
            for i in range(0,length,step):
                rows=pos[max(0,i-half):min(length,i+half+1)]
                smooth.append(dict(center=i+1,value=sum(p['entropy'] for p in rows)/len(rows)))
            put('44_entropie_msa_lissee',plot_svg(smooth,f'Entropie locale (fenêtre ≈ {win} colonnes)','Entropie moyenne (bits)'),'MSA','Lissage descriptif de l’entropie autour des colonnes du MSA.','MSA')
    return figures,meta


def build_report_v6(out:Path,query:Seq,ref_seqs:list[Seq],hsps:list[HSP],mode:str,settings:dict,
                    history:list,msa_seqs:list[Seq]|None=None,dot:dict|None=None,sha_query='',sha_refs='',motifs=None,orf_rows=None):
    report=build_report_v5(out,query,ref_seqs,hsps,mode,settings,history,msa_seqs,dot,sha_query,sha_refs,motifs,orf_rows)
    figures,meta=atlas_additions(query,ref_seqs,hsps,msa_seqs,mode)
    u=units_for(mode,sequence_stats(query)['kind']); figdir=out/'figures';figdir.mkdir(exist_ok=True)
    tree=out/'tree_fasttree.nwk';tree_error=''
    # Never reuse an SVG belonging to an earlier tree (including resumed projects).
    (figdir/(TREE_FIGURE+'.svg')).unlink(missing_ok=True)
    if tree.is_file():
        try:
            figures[TREE_FIGURE]=tree_svg(tree.read_text(encoding='utf-8'))
            meta[TREE_FIGURE]={'group':'Phylogénie','source':'FastTree',
                'description':'Arbre FastTree issu du MSA ; longueurs en substitutions/site lorsqu’elles sont disponibles. '
                'Racine de dessin arbitraire, sans enracinement biologique. Supports locaux SH-like, pas des bootstraps.'}
        except (ValueError,OSError) as exc:
            tree_error='Figure d’arbre indisponible : '+str(exc)
    from av_i18n import localize_svg
    for name,svg in list(figures.items()):
        figures[name]=localize_svg(svg)
        (figdir/(name+'.svg')).write_text(figures[name],encoding='utf-8')
    page=report.read_text(encoding='utf-8')
    from av_vsearch import report_section
    vsearch_section=report_section(out)
    if vsearch_section:
        page=page.replace('</nav>','<button data-tab="vsearch">VSEARCH</button></nav>',1)
        page=page.replace('<footer class="footer">',vsearch_section+'<footer class="footer">',1)
        page=page.replace('</body>','''<script>
titles.vsearch='Analyses VSEARCH';
const openVsearch=document.createElement('button');openVsearch.textContent='Voir les analyses VSEARCH';
openVsearch.onclick=()=>tab('vsearch');document.getElementById('readingActions').appendChild(openVsearch);
</script></body>''',1)
    cards=''.join('<article class="panel figcard" data-group="'+html.escape(m['group'],quote=True)+'" data-figure="'+html.escape(name,quote=True)+'">'
        +'<div class="row"><h2>'+html.escape(figure_title(figures[name],name))+'</h2><button class="downloadSvg right" data-fig="'+html.escape(name,quote=True)+'">Télécharger SVG</button></div>'
        +'<div class="figure-svg">'+figures[name]+'</div><p class="hint"><b>'+html.escape(m['group'])+'.</b> '+html.escape(m['description'])+'</p></article>'
        for name,m in meta.items())
    page=page.replace('</div></section><section id="exports"',cards+'</div></section><section id="exports"',1)
    if tree_error:
        page=page.replace('<div class="atlas">','<p class="notice" role="status">'+html.escape(tree_error)+'</p><div class="atlas">',1)
    n=page.count('class="panel figcard"')
    page=re.sub(r'\d+ figures SVG vectorielles',f'{n} figures SVG vectorielles',page,count=1)
    # Replacement of the original misleading mixed-reference HSP mapping.
    replacement=localize_svg(tracks_svg(hsps,len(query.seq),{r.id:len(r.seq) for r in ref_seqs},mode)) if hsps else None
    if replacement:
        old=(figdir/'23_cartographie_hsp.svg')
        if old.exists():old.write_text(replacement,encoding='utf-8')
        page=re.sub(r'(<article class="panel figcard" data-figure="23_cartographie_hsp">.*?<div class="figure-svg">).*?(</div><p class="hint">)',lambda match:match.group(1)+replacement+match.group(2),page,flags=re.S,count=1)
        page=page.replace('Le repère Y dépend de chaque référence : comparer les segments, pas inférer une structure génomique complète.', 'Pistes indépendantes : les positions QUERY et RÉFÉRENCE utilisent leurs unités propres (nt ou aa).')
    unit_notice=('Recherche traduite : coordonnées QUERY en '+u.query+', références en '+u.subject+', alignement HSP en '+u.aligned+'. '
        'Ne pas convertir un intervalle en longueur protéique par simple soustraction, ni comparer des axes nt et aa comme s’ils avaient la même unité.'
        if u.translated_query or u.translated_subject else 'Coordonnées QUERY/RÉFÉRENCE en '+u.query+'/'+u.subject+' ; alignement en '+u.aligned+'.')
    page=page.replace('<section id="coverage" class="section">','<section id="coverage" class="section"><div class="notice" role="note"><strong>Unités :</strong> '+html.escape(unit_notice)+'</div>',1)
    page=page.replace('<div class="panel"><h2>Comprendre ce qui est calculé</h2>', '<div class="panel"><h2>Comprendre ce qui est calculé</h2><p><b>Unités de cette recherche :</b> '+html.escape(unit_notice)+'</p>',1)
    # Classer également les figures héritées de v4/v5 dans les domaines affichés.
    for i in range(1,27):
        group=('Séquences' if i in (1,2,3,14,15) else
               'MSA' if i in (9,10,11,12,13,25,26) else
               'Traduction' if i==7 else
               'Couvertures' if i in (8,17,18,23,24) else
               'Qualité des HSP' if i in (5,6,19,20,21,22) else 'HSP et coordonnées')
        page=re.sub(r'(<article class="panel figcard" )(data-figure="'+str(i).zfill(2)+r'[^"]+")',
                    lambda m:m.group(1)+'data-group="'+group+'" '+m.group(2),page)
    # Clean, useful filtering by scientific theme: robust to previous figures without data-group.
    filters='''<label style="margin-left:8px">Domaine <select id="figGroup"><option value="">Tous les domaines</option><option value="Séquences">Séquences</option><option value="HSP et coordonnées">HSP et coordonnées</option><option value="Qualité des HSP">Qualité des HSP</option><option value="Couvertures">Couvertures</option><option value="Traduction">Traduction</option><option value="MSA">MSA</option><option value="Phylogénie">Phylogénie</option></select></label><span id="figVisible" class="pill"></span>'''
    page=page.replace('<span class="pill">'+str(n)+' figures SVG vectorielles</span>','<span class="pill">'+str(n)+' figures SVG vectorielles</span>'+filters,1)
    page=page.replace('</head>','''<style>
.figcard[data-figure='45_arbre_fasttree']{grid-column:1/-1}
.figcard[data-figure='45_arbre_fasttree'] .figure-svg{max-height:650px;overflow:auto}
.figcard[data-figure='45_arbre_fasttree'] svg{min-width:1052px;max-height:none}
@media print{.figcard[data-figure='45_arbre_fasttree'] .figure-svg{max-height:none;overflow:visible}
.figcard[data-figure='45_arbre_fasttree'] svg{min-width:0}}
</style></head>''',1)
    page=page.replace('</body>','''<script>(function(){
  const text=document.getElementById('figFilter'),group=document.getElementById('figGroup'),count=document.getElementById('figVisible');
  if(!text||!group)return;
  function filter(){let n=0;document.querySelectorAll('.figcard').forEach(card=>{
    const hay=(card.dataset.figure+' '+card.textContent).toLowerCase();
    const ok=hay.includes(text.value.trim().toLowerCase())&&(!group.value||(card.dataset.group||'')===group.value);
    card.hidden=!ok;if(ok)n++;
  });count.textContent=n+' affichées';}
  text.addEventListener('input',filter);group.addEventListener('change',filter);filter();
  const tree=document.querySelector(".figcard[data-figure='45_arbre_fasttree']");
  if(tree){
    const button=document.createElement('button');button.id='openTree';button.textContent='Voir l’arbre';
    button.onclick=()=>{
      text.value='';group.value='Phylogénie';filter();tab('figures');
      const details=tree.querySelector('.figure-view');if(details)details.open=true;
      tree.scrollIntoView({block:'start'});
    };
    document.getElementById('readingActions').appendChild(button);
  }
})();</script></body>''',1)
    # Track sequence and graphics units for subsequent PDF / figure export.
    catalog=[]
    from av_i18n import translate
    for group,name in re.findall(r'<article class="panel figcard" data-group="([^"]+)" data-figure="([^"]+)"',page):
        m=meta.get(name,{})
        catalog.append(dict(file=f'figures/{name}.svg',category=translate(group),
            description=translate(m.get('description','Figure descriptive : consulter sa légende dans le rapport HTML')),
            coordinate_query='NA' if name==TREE_FIGURE else u.query,
            coordinate_subject='NA' if name==TREE_FIGURE else u.subject,
            alignment_columns='NA' if name==TREE_FIGURE else u.aligned))
    (out/'figures_catalog.tsv').write_text('file\tcategory\tdescription\tquery_unit\tsubject_unit\talignment_unit\n'+
        ''.join(f"{v['file']}\t{v['category']}\t{v['description']}\t{v['coordinate_query']}\t{v['coordinate_subject']}\t{v['alignment_columns']}\n" for v in catalog),encoding='utf-8')
    page=page.replace('Sequence explorer · v5.0','Sequence explorer · v6.0')
    from av_i18n import localize_html
    report.write_text(localize_html(page),encoding='utf-8');return report
