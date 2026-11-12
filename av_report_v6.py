"""v6 atlas: conditional supplementary views, native units and direct navigation.

Charts are descriptive; no inferential statistics are applied to dependent HSPs.
"""
from __future__ import annotations
import collections
import html
import json
import math
import re
from pathlib import Path
from av_display import figure_title
from av_core import Seq,HSP,coverage_data,merge_intervals,interval_length,msa_metrics,sequence_stats
from av_units import units_for,coord_audit
from av_report_v5 import build_report_v5
from av_tree import TREE_FIGURE, tree_svg
from av_report_v4 import bars_svg,hist_svg,scatter_svg,plot_svg,shell,finish,txt,INK,MUTED,BLUE,TEAL


def tracks_svg(hsps:list[HSP],qlen:int,rlen:dict,mode:str,max_refs=20):
    """Independent tracks: one native scale for QUERY and one for each reference.

    Avoid drawing a nucleotide/amino-acid diagonal that implies a shared scale.
    """
    u=units_for(mode);groups=collections.defaultdict(list)
    for h in hsps:groups[h.subject].append(h)
    groups=dict(list(groups.items())[:max_refs]);width=1060
    height=max(215,136+len(groups)*78)
    left=205;right=45;usable=width-left-right
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-label="HSP coordinate tracks">',
         f'<rect width="{width}" height="{height}" rx="12" fill="#102039"/>',
         '<text x="25" y="30" fill="#d5e5f9" font-size="17">Coverage per reference · two independent coordinate systems</text>',
         f'<text x="25" y="51" fill="#9fb4ca" font-size="12">QUERY: {u.query} · REFERENCE: {u.subject} · HSP alignment: {u.aligned} · 1-based inclusive intervals</text>']
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
    out.append(f'<text x="25" y="{height-12}" fill="#9fb4ca" font-size="11">Green = covered QUERY positions · blue = covered reference positions · gray = no detected HSP</text>')
    out.append('</svg>');return ''.join(out)


def frame_tracks_svg(hsps:list[HSP],qlen:int):
    frames=sorted({h.qframe for h in hsps if h.qframe},key=lambda v:(v<0,abs(v)))
    if not frames:return None
    w=1000;height=90+len(frames)*49;l=115;usable=810
    out=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {height}">',f'<rect width="{w}" height="{height}" rx="12" fill="#102039"/>',
        '<text x="22" y="30" font-size="16" fill="#dae8f5">QUERY frames · HSP tracks (nt positions)</text>']
    for i,fr in enumerate(frames):
        y=55+i*49
        out.append(f'<text x="20" y="{y+10}" fill="#c5d4e6" font-size="12">Cadre {fr:+d}</text>')
        out.append(f'<rect x="{l}" y="{y}" width="{usable}" height="19" fill="#304058" rx="4"/>')
        for h in (h for h in hsps if h.qframe==fr):
            a,b=h.qspan;x=l+(a-1)/max(1,qlen)*usable;length=(b-a+1)/max(1,qlen)*usable
            out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(2,length):.1f}" height="19" fill="{"#64d9bd" if fr>0 else "#f6ad70"}" opacity=".75"><title>{txt(h.subject)} · {a}–{b} nt · frame {fr:+d}</title></rect>')
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
        put('27_query_hsp_span',hist_svg(spans_q,'Each HSP\'s query span',f'Positions QUERY ({u.query})'),'HSPs and coordinates','Query coordinates use its native alphabet, never alignment columns.')
        put('28_reference_hsp_span',hist_svg(spans_s,'Each HSP\'s reference span',f'Positions REFERENCE ({u.subject})'),'HSPs and coordinates','Segment length uses the reference\'s native unit.')
        put('29_query_hsp_percentage',hist_svg([100*n/max(1,len(query.seq)) for n in spans_q], 'Each HSP span / query (%)','Relative query span (%)'),'HSPs and coordinates','Each HSP is measured separately; overlaps are not subtracted.')
        put('30_gap_counts',hist_svg([h.gap_columns for h in hsps],'Columns with gaps per HSP','Aligned columns ('+u.aligned+')'),'HSP quality','Columns containing at least one gap; in BLASTX, columns represent amino acids.')
        put('31_substitutions_hsp',hist_svg([h.mismatches for h in hsps],'Substitutions par HSP','Substitutions (columns '+u.aligned+')'),'HSP quality','Comparison of aligned sequences; gaps are not counted as mismatches.')
        put('32_substitutions_vs_gaps',scatter_svg([(h.mismatches,h.gap_columns,h.subject) for h in hsps], 'HSP substitutions and gaps', 'Substitutions (columns '+u.aligned+')','Columns with gaps',max(1,max(h.mismatches for h in hsps)),max(1,max(h.gap_columns for h in hsps))),'HSP quality','Values depend on local-alignment lengths; this comparison is descriptive.')
        cov=coverage_data(hsps,len(query.seq),ref_lengths)
        both=[r for r in cov if r['reference_coverage'] is not None]
        if both:put('33_relative_query_vs_reference_coverage',scatter_svg([(r['query_coverage'],r['reference_coverage'],r['subject']) for r in both],'Relative coverage of both sequences','Query coverage (%)','Reference coverage (%)'),'Coverage','Percentages use different denominators; do not turn them into directly comparable lengths.')
        accum=[];allspans=[]
        for i,r in enumerate(cov,1):
            allspans.extend(r['query_segments']);current=merge_intervals(allspans,len(query.seq))
            accum.append(dict(center=i,value=100*interval_length(current)/max(1,len(query.seq))))
        if len(accum)>1:put('34_cumulative_hit_coverage',plot_svg(accum,'Cumulative QUERY coverage as references are added','QUERY covered (%)'),'Coverage','Union of covered positions as references are added in bitscore order; this is not a sampling saturation curve.')
        if u.translated_query:
            put('36_query_frame_coverage',frame_tracks_svg(hsps,len(query.seq)),'Translation','Positive and negative frames use the query\'s original nucleotide coordinates.')
        if u.translated_subject:
            subject_frames=collections.Counter(str(h.sframe) for h in hsps if h.sframe)
            if subject_frames:put('37_reference_frames',bars_svg(sorted(subject_frames.items()),'REFERENCE reading frames'),'Translation','Frames use reference nucleotide coordinates, independently of QUERY coordinates.')
        # Check the 3 nt to 1 aa relationship without treating a column as a base.
        if u.translated_query or u.translated_subject:
            audit=coord_audit(hsps,mode)
            def scale(r,which):
                a=r[f'{which}_span']; b=r[f'{which}_aligned_chars'];return round(a/b,3) if b else None
            if u.translated_query:put('38_ratio_nt_aa_query',hist_svg([scale(r,'query') for r in audit if r['query_aligned_chars']], 'QUERY nucleotide span / aligned residues','nt per aligned aa'),'Translation','Quality check: expected ratio of 3 nt/aa; a difference can indicate an atypical HSP or a coordinate issue.')
            if u.translated_subject:put('39_ratio_nt_aa_reference',hist_svg([scale(r,'reference') for r in audit if r['reference_aligned_chars']], 'REFERENCE nucleotide span / aligned residues','nt per aligned aa'),'Translation','Quality check: expected ratio of 3 nt/aa; keep the two coordinate types separate.')
    if msa_seqs and len(msa_seqs)>=2:
        m=msa_metrics(msa_seqs);pos=m['positions'];length=m['length'];n=m['n']
        put('40_msa_occupancy',plot_svg([dict(center=p['position'],value=100-p['gap_pct']) for p in pos], 'MSA occupancy','Sequences without gaps (%)'),'MSA','Positions use MSA columns ('+u.msa+'), not native QUERY coordinates.','MSA')
        put('41_entropy_vs_gaps',scatter_svg([(p['gap_pct'],p['entropy'],str(p['position'])) for p in pos],'Variability versus gaps','Gaps (%)','Entropy (bits)',100,max(1,max(p['entropy'] for p in pos))),'MSA','High entropy in sparsely occupied columns requires cautious interpretation.','MSA')
        put('42_consensus_composition',bars_svg(collections.Counter(m['consensus'].replace('-','')).most_common(25),'MSA consensus composition'),'MSA','Column-wise consensus; not an observed biological sequence.','MSA')
        put('43_conservation_distribution',hist_svg([p['conservation'] for p in pos],'Distribution of MSA conservation','Conservation (%)'),'MSA','Conservation uses non-gap residues and depends on the selected sequences.','MSA')
        win=min(max(10,length//30),100);half=win//2
        smooth=[]
        if length>1:
            step=max(1,length//180)
            for i in range(0,length,step):
                rows=pos[max(0,i-half):min(length,i+half+1)]
                smooth.append(dict(center=i+1,value=sum(p['entropy'] for p in rows)/len(rows)))
            put('44_smoothed_msa_entropy',plot_svg(smooth,f'Local entropy (window ≈ {win} columns)','Mean entropy (bits)'),'MSA','Descriptive smoothing of entropy around MSA columns.','MSA')
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
            meta[TREE_FIGURE]={'group':'Phylogeny','source':'FastTree',
                'description':'FastTree tree inferred from the MSA; lengths use substitutions/site when available. '
                'Arbitrary drawing root, without biological rooting. Local SH-like support, not bootstrap values.'}
        except (ValueError,OSError) as exc:
            tree_error='Figure d’arbre unavailable: '+str(exc)

    for name,svg in list(figures.items()):
        figures[name]=str(svg)
        (figdir/(name+'.svg')).write_text(figures[name],encoding='utf-8')
    page=report.read_text(encoding='utf-8')
    from av_vsearch import report_section
    vsearch_section=report_section(out)
    if vsearch_section:
        page=page.replace('</nav>','<button data-tab="vsearch">VSEARCH</button></nav>',1)
        page=page.replace('<footer class="footer">',vsearch_section+'<footer class="footer">',1)
        page=page.replace('</body>','''<script>
titles.vsearch='VSEARCH analyses';
const openVsearch=document.createElement('button');openVsearch.textContent='View VSEARCH analyses';
openVsearch.onclick=()=>tab('vsearch');document.getElementById('readingActions').appendChild(openVsearch);
</script></body>''',1)
    cards=''.join('<article class="panel figcard" data-group="'+html.escape(m['group'],quote=True)+'" data-figure="'+html.escape(name,quote=True)+'">'
        +'<div class="row"><h2>'+html.escape(figure_title(figures[name],name))+'</h2><button class="downloadSvg right" data-fig="'+html.escape(name,quote=True)+'">Download SVG</button></div>'
        +'<div class="figure-svg">'+figures[name]+'</div><p class="hint"><b>'+html.escape(m['group'])+'.</b> '+html.escape(m['description'])+'</p></article>'
        for name,m in meta.items())
    page=page.replace('</div></section><section id="exports"',cards+'</div></section><section id="exports"',1)
    if tree_error:
        page=page.replace('<div class="atlas">','<p class="notice" role="status">'+html.escape(tree_error)+'</p><div class="atlas">',1)
    n=page.count('class="panel figcard"')
    page=re.sub(r'\d+ vector SVG figures',f'{n} vector SVG figures',page,count=1)
    # Replacement of the original misleading mixed-reference HSP mapping.
    replacement=str(tracks_svg(hsps,len(query.seq),{r.id:len(r.seq) for r in ref_seqs},mode)) if hsps else None
    if replacement:
        old=(figdir/'23_hsp_mapping.svg')
        if old.exists():old.write_text(replacement,encoding='utf-8')
        page=re.sub(r'(<article class="panel figcard" data-figure="23_hsp_mapping">.*?<div class="figure-svg">).*?(</div><p class="hint">)',lambda match:match.group(1)+replacement+match.group(2),page,flags=re.S,count=1)
        page=page.replace('The Y axis depends on each reference: compare segments without inferring a complete genome structure.', 'Independent tracks: query and reference positions use their own units (nt or aa).')
    unit_notice=('Translated search: QUERY coordinates in '+u.query+', references in '+u.subject+', HSP alignment in '+u.aligned+'. '
        'Do not obtain protein length by subtracting nucleotide coordinates or compare nucleotide and amino-acid axes as if they had the same units.'
        if u.translated_query or u.translated_subject else 'QUERY/REFERENCE coordinates in '+u.query+'/'+u.subject+' ; alignment in '+u.aligned+'.')
    page=page.replace('<section id="coverage" class="section">','<section id="coverage" class="section"><div class="notice" role="note"><strong>Units:</strong> '+html.escape(unit_notice)+'</div>',1)
    page=page.replace('<div class="panel"><h2>Understand what is calculated</h2>', '<div class="panel"><h2>Understand what is calculated</h2><p><b>Search units:</b> '+html.escape(unit_notice)+'</p>',1)
    # Group inherited v4/v5 figures into the displayed categories.
    for i in range(1,27):
        group=('Sequences' if i in (1,2,3,14,15) else
               'MSA' if i in (9,10,11,12,13,25,26) else
               'Translation' if i==7 else
               'Coverage' if i in (8,17,18,23,24) else
               'HSP quality' if i in (5,6,19,20,21,22) else 'HSPs and coordinates')
        page=re.sub(r'(<article class="panel figcard" )(data-figure="'+str(i).zfill(2)+r'[^"]+")',
                    lambda m:m.group(1)+'data-group="'+group+'" '+m.group(2),page)
    # Clean, useful filtering by scientific theme: robust to previous figures without data-group.
    filters='''<label style="margin-left:8px">Category <select id="figGroup"><option value="">All categories</option><option value="Sequences">Sequences</option><option value="HSPs and coordinates">HSPs and coordinates</option><option value="HSP quality">HSP quality</option><option value="Coverage">Coverage</option><option value="Translation">Translation</option><option value="MSA">MSA</option><option value="Phylogeny">Phylogeny</option></select></label><span id="figVisible" class="pill"></span>'''
    page=page.replace('<span class="pill">'+str(n)+' vector SVG figures</span>','<span class="pill">'+str(n)+' vector SVG figures</span>'+filters,1)
    page=page.replace('</head>','''<style>
.figcard[data-figure='45_fasttree_tree']{grid-column:1/-1}
.figcard[data-figure='45_fasttree_tree'] .figure-svg{max-height:650px;overflow:auto}
.figcard[data-figure='45_fasttree_tree'] svg{min-width:1052px;max-height:none}
@media print{.figcard[data-figure='45_fasttree_tree'] .figure-svg{max-height:none;overflow:visible}
.figcard[data-figure='45_fasttree_tree'] svg{min-width:0}}
</style></head>''',1)
    page=page.replace('</body>','''<script>(function(){
  const text=document.getElementById('figFilter'),group=document.getElementById('figGroup'),count=document.getElementById('figVisible');
  if(!text||!group)return;
  function filter(){let n=0;document.querySelectorAll('.figcard').forEach(card=>{
    const hay=(card.dataset.figure+' '+card.textContent).toLowerCase();
    const ok=hay.includes(text.value.trim().toLowerCase())&&(!group.value||(card.dataset.group||'')===group.value);
    card.hidden=!ok;if(ok)n++;
  });count.textContent=n+' shown';}
  text.addEventListener('input',filter);group.addEventListener('change',filter);filter();
  const tree=document.querySelector(".figcard[data-figure='45_fasttree_tree']");
  if(tree){
    const button=document.createElement('button');button.id='openTree';button.textContent='View the tree';
    button.onclick=()=>{
      text.value='';group.value='Phylogeny';filter();tab('figures');
      const details=tree.querySelector('.figure-view');if(details)details.open=true;
      tree.scrollIntoView({block:'start'});
    };
    document.getElementById('readingActions').appendChild(button);
  }
})();</script></body>''',1)
    # Track sequence and graphics units for subsequent PDF / figure export.
    catalog=[]

    for group,name in re.findall(r'<article class="panel figcard" data-group="([^"]+)" data-figure="([^"]+)"',page):
        m=meta.get(name,{})
        catalog.append(dict(file=f'figures/{name}.svg',category=str(group),
            description=str(m.get('description','Descriptive figure: see its caption in the HTML report')),
            coordinate_query='NA' if name==TREE_FIGURE else u.query,
            coordinate_subject='NA' if name==TREE_FIGURE else u.subject,
            alignment_columns='NA' if name==TREE_FIGURE else u.aligned))
    (out/'figures_catalog.tsv').write_text('file\tcategory\tdescription\tquery_unit\tsubject_unit\talignment_unit\n'+
        ''.join(f"{v['file']}\t{v['category']}\t{v['description']}\t{v['coordinate_query']}\t{v['coordinate_subject']}\t{v['alignment_columns']}\n" for v in catalog),encoding='utf-8')
    page=page.replace('Sequence explorer · v5.0','Sequence explorer · v6.0')

    report.write_text(str(page),encoding='utf-8');return report
