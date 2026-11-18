"""Pure scientific SVG plots and conditional figure composition.

This module has no HTML template mutation or report-file side effects.
"""
from __future__ import annotations
import collections
import html
import math
from av_core import Seq,HSP,sequence_stats,msa_metrics,distance_matrix,sliding_scores,coverage_data,merge_intervals,interval_length
from av_units import units_for,coord_audit
from av_visual import plot_svg

MAX_DISTANCE_WORK=20_000_000

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
        parts.append(f'<text x="32" y="145" fill="{MUTED}" font-size="15">No observations</text>');return finish(parts)
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
        parts.append(f'<rect x="{x:.2f}" y="{charty+h-bh:.2f}" width="{bw:.2f}" height="{bh:.2f}" rx="2" fill="{TEAL}"><title>{n} observations · {minimum+i*(maximum-minimum)/bins:.2f} to {minimum+(i+1)*(maximum-minimum)/bins:.2f}</title></rect>')
    parts.extend([f'<text x="55" y="248" fill="{MUTED}" font-size="11" font-family="sans-serif">{minimum:,.2f}</text>',
                  f'<text x="745" y="248" text-anchor="end" fill="{MUTED}" font-family="sans-serif" font-size="11">{max(values):,.2f}</text>',
                  f'<text x="390" y="269" fill="{INK}" text-anchor="middle" font-family="sans-serif" font-size="11">{txt(x_label)} · n={len(values)}</text>'])
    return finish(parts)

def bars_svg(items,title,nmax=18):
    items=sorted(items,key=lambda v:v[1],reverse=True)[:nmax]
    height=max(180,len(items)*31+65);parts=shell(title,height=height)
    if not items:parts.append(f'<text x="30" y="85" fill="{MUTED}">No data</text>');return finish(parts)
    peak=max(float(k) for _,k in items) or 1
    for i,(name,val) in enumerate(items):
        y=46+i*31
        parts.append(f'<text class="data-label" x="20" y="{y+13}" fill="{INK}" font-family="monospace" font-size="11">{txt(str(name)[:20])}</text>')
        parts.append(f'<rect x="185" y="{y}" width="{max(2,float(val)/peak*500):.2f}" height="18" rx="3" fill="{BLUE}"><title class="data-label">{txt(name)}: {val}</title></rect>')
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
    parts.append(f'<text x="{margin}" y="{height-10}" fill="{MUTED}" font-size="11">p-distance · gaps excluded · {n} sequences max</text>')
    return finish(parts)

def annotation_svg(length,covered,title='Uncovered regions (QUERY)'):
    parts=shell(title,height=173)
    if length<=0:return finish(parts)
    parts.append(f'<rect x="45" y="66" width="690" height="26" fill="#3b4356" rx="5"/>')
    for a,b in covered:
        x=45+(a-1)/length*690;w=(b-a+1)/length*690
        parts.append(f'<rect x="{x:.2f}" y="66" width="{w:.2f}" height="26" fill="{TEAL}" rx="2"><title>Covered HSP {a}–{b}</title></rect>')
    parts.append(f'<text x="45" y="120" fill="{MUTED}" font-size="12">1</text><text x="735" y="120" text-anchor="end" fill="{MUTED}" font-size="12">{length:,}</text>')
    parts.append(f'<text x="45" y="146" fill="{MUTED}" font-size="12">Green: at least one HSP · Gray: no detected region</text>')
    return finish(parts)

def mapping_svg(hsps:list[HSP],query_length:int,ref_lengths:dict):
    """Compatibility helper using independent native-coordinate tracks."""
    query_translated=any(h.qframe for h in hsps)
    reference_translated=any(h.sframe for h in hsps)
    mode='tblastx' if query_translated and reference_translated else 'blastx' if query_translated else 'tblastn' if reference_translated else 'blastn'
    return tracks_svg(hsps,query_length,ref_lengths,mode)

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
        out.append(f'<text x="20" y="{y+10}" fill="#c5d4e6" font-size="12">Frame {fr:+d}</text>')
        out.append(f'<rect x="{l}" y="{y}" width="{usable}" height="19" fill="#304058" rx="4"/>')
        for h in (h for h in hsps if h.qframe==fr):
            a,b=h.qspan;x=l+(a-1)/max(1,qlen)*usable;length=(b-a+1)/max(1,qlen)*usable
            out.append(f'<rect x="{x:.1f}" y="{y}" width="{max(2,length):.1f}" height="19" fill="{"#64d9bd" if fr>0 else "#f6ad70"}" opacity=".75"><title>{txt(h.subject)} · {a}–{b} nt · frame {fr:+d}</title></rect>')
        out.append(f'<text x="{l}" y="{y+34}" fill="#91a8c3" font-size="10">1</text><text x="{l+usable}" y="{y+34}" text-anchor="end" fill="#91a8c3" font-size="10">{qlen:,} nt</text>')
    out.append('</svg>');return ''.join(out)

def atlas_additions(query:Seq,ref_seqs:list[Seq],hsps:list[HSP],msa_seqs:list[Seq]|None,mode:str,metrics=None):
    u=units_for(mode,sequence_stats(query)['kind']); figures={};meta={}
    def put(key,svg,group,description,kind='HSP'):
        if svg:
            figures[key]=svg;meta[key]={'group':group,'description':description,'source':kind}
    ref_lengths={r.id:len(r.seq) for r in ref_seqs}
    if hsps:
        spans_q=[h.qspan[1]-h.qspan[0]+1 for h in hsps]
        spans_s=[h.sspan[1]-h.sspan[0]+1 for h in hsps]
        put('27_query_hsp_span',hist_svg(spans_q,'Each HSP\'s query span',f'QUERY positions ({u.query})'),'HSPs and coordinates','Query coordinates use its native alphabet, never alignment columns.')
        put('28_reference_hsp_span',hist_svg(spans_s,'Each HSP\'s reference span',f'REFERENCE positions ({u.subject})'),'HSPs and coordinates','Segment length uses the reference\'s native unit.')
        put('29_query_hsp_percentage',hist_svg([100*n/max(1,len(query.seq)) for n in spans_q], 'Each HSP span / query (%)','Relative query span (%)'),'HSPs and coordinates','Each HSP is measured separately; overlaps are not subtracted.')
        put('30_gap_counts',hist_svg([h.gap_columns for h in hsps],'Columns with gaps per HSP','Aligned columns ('+u.aligned+')'),'HSP quality','Columns containing at least one gap; in BLASTX, columns represent amino acids.')
        put('31_substitutions_hsp',hist_svg([h.mismatches for h in hsps],'Substitutions per HSP','Substitutions (columns '+u.aligned+')'),'HSP quality','Comparison of aligned sequences; gaps are not counted as mismatches.')
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
        m=metrics if metrics is not None else msa_metrics(msa_seqs);pos=m['positions'];length=m['length'];n=m['n']
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

def build_figures(query,ref_seqs,hsps,msa_seqs,mode,motifs=None,orf_rows=None,metrics=None):
    """Compose independent plots once, with an explicit catalog for every figure."""
    figures={}
    allseq=[query,*ref_seqs]
    stats=[sequence_stats(s) for s in allseq]
    summary=collections.Counter(query.seq.replace('-','').replace('.',''))
    figures['01_query_composition']=bars_svg(summary.most_common(18),'Query sequence composition')
    u=units_for(mode,sequence_stats(query)['kind'])
    if u.query==u.subject:
        figures['02_sequence_lengths']=hist_svg([x['ungapped'] for x in stats],
            'Length distribution',f'Length ({u.query}, without gaps)')
    else:
        # DNA and protein lengths cannot be compared on a shared axis.
        if ref_seqs:
            figures['02b_reference_lengths']=hist_svg([sequence_stats(x)['ungapped'] for x in ref_seqs],
                'Reference length distribution',f'REFERENCE ({u.subject}, without gaps)')
    nt_stats=[x for x in stats if x['kind']=='DNA/RNA' and x['gc_pct'] is not None]
    if nt_stats:
        figures['03_gc_quality']=scatter_svg([(x['ungapped'],x['gc_pct'],x['id']) for x in nt_stats],
            'Nucleotide length / GC (%)','Length (nt)','GC (%)',max(1,max(x['ungapped'] for x in nt_stats)),100)
    if hsps:
        figures['04_identity_coverage']=scatter_svg([(h.identity,100*(h.qspan[1]-h.qspan[0]+1)/max(1,len(query.seq)),h.subject) for h in hsps],
            'Local identity / HSP span on the query','HSP identity (%)','Query span (%)')
        figures['05_identity_histogram']=hist_svg([h.identity for h in hsps],'HSP identities (%)','Local identity (%)')
        figures['06_hsp_length_histogram']=hist_svg([h.align_length for h in hsps],'Aligned HSP lengths',f'Aligned columns ({u.aligned})')
        frames=collections.Counter(str(h.qframe) for h in hsps if h.qframe)
        if frames:figures['07_frames']=bars_svg(sorted(frames.items()),'Query reading frames (HSP)')
        covered=[]
        from av_core import merge_intervals
        covered=merge_intervals([h.qspan for h in hsps],len(query.seq))
        figures['08_segments']=annotation_svg(len(query.seq),covered,f'Regions without HSPs (QUERY, {u.query})')
    if msa_seqs:
        m=metrics if metrics is not None else msa_metrics(msa_seqs)
        for suffix,title,key in [('09_conservation','MSA conservation (%)','conservation'),('10_entropy','Shannon entropy per column','entropy'),('11_gaps','Columns with gaps (%)','gap_pct')]:
            figures[suffix]=plot_svg([dict(center=x['position'],value=x[key]) for x in m['positions']],title,key)
        if 2<=len(msa_seqs)<=45 and len(msa_seqs)**2*len(msa_seqs[0].seq)<=MAX_DISTANCE_WORK:
            mat=distance_matrix(msa_seqs)
            figures['12_p_distance']=heatmap_svg([s.id for s in msa_seqs],mat,'Pairwise p-distances (MSA)')
        counts=collections.Counter(str(round(p['gap_pct']//10)*10)+'–'+str(min(100,round(p['gap_pct']//10)*10+10))+' %' for p in m['positions'])
        figures['13_gaps_distribution']=bars_svg(counts.items(),'Distribution of per-column gap rates')
    if motifs:figures['14_motifs']=bars_svg(collections.Counter(r.get('seqid','?') for r in motifs).items(),'Motif occurrences per sequence')
    if orf_rows:figures['15_orfs']=hist_svg([r['length_aa'] for r in orf_rows],'Detected ORF lengths','Amino acids')
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
        figures['23_hsp_mapping']=tracks_svg(hsps,len(query.seq),refs_len,mode)
        d=density_svg(hsps,len(query.seq))
        if d:figures['24_reference_density_on_query']=d
        explanations['23_hsp_mapping']='Independent query and reference tracks use native coordinate units and their own scales.'
        explanations['24_reference_density_on_query']='References with at least one HSP at the sampled position; not an abundance measure.'
    if msa_seqs:
        met=metrics if metrics is not None else msa_metrics(msa_seqs)
        figures['25_gap_fraction_per_sequence']=bars_svg([(r.id,round(100*sum(x in '-.?' for x in r.seq)/max(1,len(r.seq)),2)) for r in msa_seqs[:35]],'Gaps per MSA sequence (%)')
        figures['26_msa_entropy_distribution']=hist_svg([p['entropy'] for p in met['positions']],'Distribution of MSA column entropy','Shannon (bits)')
        explanations['25_gap_fraction_per_sequence']='Gaps can reflect biology, different lengths or alignment artifacts.'
        explanations['26_msa_entropy_distribution']='Calculated from observed residues, excluding gaps and missing characters.'
    if not figures:
        figures['16_sequence_ambiguity']=hist_svg([s['n_pct'] for s in map(sequence_stats,[query,*ref_seqs])],'% N per sequence','N (%)')
    additions,metadata=atlas_additions(query,ref_seqs,hsps,msa_seqs,mode,metrics)
    figures.update(additions)
    for name in figures:
        if name not in metadata:
            number=int(name[:2])
            group=('Sequences' if number in (1,2,3,14,15,16) else
                   'MSA' if number in (9,10,11,12,13,25,26) else
                   'Translation' if number==7 else
                   'Coverage' if number in (8,17,18,23,24) else 'HSP quality')
            metadata[name]={'group':group,'source':'MSA' if group=='MSA' else 'HSP' if hsps else 'Sequence',
                'description':explanations.get(name,'Descriptive result calculated from the recorded inputs. Local HSP statistics are not whole-sequence identity.')}
    return figures,metadata
