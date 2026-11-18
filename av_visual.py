"""Standalone HTML/CSS/JS report and exportable SVG figures, without a CDN or external dependencies."""
from __future__ import annotations
import html
import json
import math
from pathlib import Path
from av_units import units_for
from av_core import (Seq, HSP, coverage_data, group_hsps, merge_intervals, missing_intervals,
                     interval_length, sequence_stats, sliding_scores, msa_metrics)


def coverage_svg(query_len,rows,combined,limit=100,unit="nt"):
    rows=rows[:limit]; W=1350; label=300; right=110; plot=W-label-right
    height=155+40*(len(rows)+1); ybase=110
    x=lambda k:label+(k-1)/max(1,query_len)*plot
    a=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {height}" role="img" aria-label="Coverage map">',
       f'<rect width="{W}" height="{height}" fill="#0b1220"/>',
       f'<text x="28" y="39" font-family="Arial" font-weight="700" font-size="22" fill="#eaf4ff">AlignView  /  HSP coverage</text>',
       f'<text x="28" y="64" font-family="Arial" font-size="13" fill="#a5bdd8">Query: {query_len} {unit} · no detected HSP does not prove non-homology</text>']
    for t in range(6):
        pos=1+round((query_len-1)*t/5);xx=x(pos)
        a.append(f'<line x1="{xx:.1f}" x2="{xx:.1f}" y1="91" y2="{height-25}" stroke="#26354b" stroke-dasharray="3,4"/>')
        a.append(f'<text x="{xx:.1f}" y="88" text-anchor="middle" font-family="Arial" font-size="12" fill="#9db5ce">{pos}</text>')
    def track(name,spans,index,color):
        y=ybase+index*40
        covered=interval_length(spans);cov=100*covered/max(1,query_len)
        a.append(f'<text x="28" y="{y+5}" font-family="Arial" font-size="13" fill="#e5efff">{html.escape(name[:33])}</text>')
        a.append(f'<rect x="{label}" y="{y-11}" width="{plot}" height="19" rx="4" fill="#2a3850"/>')
        for s,e in spans:
            left=x(s);w=max(1.6,plot*(e-s+1)/max(1,query_len))
            a.append(f'<rect x="{left:.1f}" y="{y-11}" width="{w:.1f}" height="19" rx="2" fill="{color}"><title>{html.escape(name)}: {s}–{e}</title></rect>')
        a.append(f'<text x="{W-96}" y="{y+4}" font-family="Arial" font-size="13" fill="#d8e9ff">{cov:.1f} %</text>')
    track('Reference union',combined,0,'#41d8ba')
    for i,row in enumerate(rows,1):track(row['subject'],row['query_segments'],i,'#70a6ff')
    a.append('</svg>')
    return '\n'.join(a)


def plot_svg(values:list[dict],title='Profile',ylabel='Value',color='#6da5ff'):
    W=1120;H=320;L=72;R=30;T=57;B=49
    valid=[r for r in values if r.get('value') is not None]
    if len(valid)>3500:
        stride=math.ceil(len(valid)/1600)
        sampled=[]
        for i in range(0,len(valid),stride):
            part=valid[i:i+stride]
            sampled.extend([min(part,key=lambda r:r['value']),max(part,key=lambda r:r['value'])])
        valid=sorted({(r['center'],r['value']):r for r in sampled}.values(),key=lambda r:r['center'])
    if not valid:
        return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 120"><text x="30" y="60">No data</text></svg>'
    ymin=min(r['value'] for r in valid);ymax=max(r['value'] for r in valid)
    if abs(ymax-ymin)<1e-12: ymin-=1;ymax+=1
    xmin=min(r['center'] for r in valid);xmax=max(r['center'] for r in valid)
    x=lambda p:L+(p-xmin)/max(1,xmax-xmin)*(W-L-R)
    y=lambda p:H-B-(p-ymin)/(ymax-ymin)*(H-B-T)
    pts=' '.join(f'{x(r["center"]):.1f},{y(r["value"]):.1f}' for r in valid)
    a=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">',
       f'<rect width="{W}" height="{H}" fill="#0d1628"/>',
       f'<text x="24" y="31" font-family="Arial" font-size="19" font-weight="700" fill="#edf5ff">{html.escape(title)}</text>']
    for i in range(5):
        yy=T+i*(H-T-B)/4;v=ymax-(ymax-ymin)*i/4
        a.append(f'<line x1="{L}" y1="{yy:.1f}" x2="{W-R}" y2="{yy:.1f}" stroke="#30415a"/>')
        a.append(f'<text x="{L-8}" y="{yy+4:.1f}" text-anchor="end" font-family="Arial" font-size="12" fill="#9fb4ce">{v:.1f}</text>')
    a.append(f'<polyline fill="none" stroke="{color}" stroke-width="3" stroke-linejoin="round" points="{pts}"/>')
    a.append(f'<text x="{W/2}" y="{H-11}" text-anchor="middle" font-family="Arial" fill="#a9bdd4" font-size="12">Position</text>')
    a.append('</svg>')
    return '\n'.join(a)


def dotplot_svg(dot:dict,x_name='QUERY',y_name='REFERENCE'):
    W=1100;H=820;L=96;T=72;R=34;B=77;Pw=W-L-R;Ph=H-T-B
    a=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}">',
       f'<rect width="{W}" height="{H}" fill="#0d1628"/>',
       f'<text x="30" y="36" fill="#e9f5ff" font-family="Arial" font-size="22" font-weight="bold">Dotplot k={dot["k"]}</text>',
       f'<rect x="{L}" y="{T}" width="{Pw}" height="{Ph}" fill="#17233a" stroke="#425471"/>']
    for i in range(6):
        xx=L+i*Pw/5; yy=T+i*Ph/5
        a.append(f'<line x1="{xx:.1f}" y1="{T}" x2="{xx:.1f}" y2="{H-B}" stroke="#304359"/>')
        a.append(f'<line x1="{L}" y1="{yy:.1f}" x2="{W-R}" y2="{yy:.1f}" stroke="#304359"/>')
        a.append(f'<text x="{xx:.1f}" y="{H-B+19}" text-anchor="middle" fill="#9fb5ca" font-family="Arial" font-size="12">{round(i*dot["xlen"]/5)}</text>')
        a.append(f'<text x="{L-11}" y="{yy+3:.1f}" text-anchor="end" fill="#9fb5ca" font-family="Arial" font-size="12">{round(i*dot["ylen"]/5)}</text>')
    for px,py in dot['points'][:30000]:
        a.append(f'<circle cx="{L+px/max(1,dot["xlen"])*Pw:.1f}" cy="{T+py/max(1,dot["ylen"])*Ph:.1f}" r="1.4" fill="#6cdbcd" opacity=".6"/>')
    a.append(f'<text x="{W/2}" y="{H-17}" fill="#d8e6fa" text-anchor="middle" font-family="Arial" font-size="14">{html.escape(x_name)}</text>')
    a.append(f'<text transform="translate(25 {H/2}) rotate(-90)" fill="#d8e6fa" text-anchor="middle" font-family="Arial" font-size="14">{html.escape(y_name)}</text>')
    a.append('</svg>');return '\n'.join(a)


def report_html(payload):
    """Compatibility renderer; new reports use the single workspace renderer."""
    from av_report import render_report
    return render_report(payload)


def build_report(*args,**kwargs):
    from av_report import build_report as build
    return build(*args,**kwargs)
