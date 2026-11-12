"""Rapport autonome HTML/CSS/JS, figures SVG exportables, sans CDN ni dépendances."""
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
       f'<text x="28" y="39" font-family="Arial" font-weight="700" font-size="22" fill="#eaf4ff">AlignView  /  couverture des HSP</text>',
       f'<text x="28" y="64" font-family="Arial" font-size="13" fill="#a5bdd8">Requête : {query_len} {unit} · absence de HSP ≠ absence d’homologie</text>']
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
            a.append(f'<rect x="{left:.1f}" y="{y-11}" width="{w:.1f}" height="19" rx="2" fill="{color}"><title>{html.escape(name)} : {s}–{e}</title></rect>')
        a.append(f'<text x="{W-96}" y="{y+4}" font-family="Arial" font-size="13" fill="#d8e9ff">{cov:.1f} %</text>')
    track('Union des références',combined,0,'#41d8ba')
    for i,row in enumerate(rows,1):track(row['subject'],row['query_segments'],i,'#70a6ff')
    a.append('</svg>')
    return '\n'.join(a)


def plot_svg(values:list[dict],title='Profil',ylabel='Valeur',color='#6da5ff'):
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
        return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 600 120"><text x="30" y="60">Aucune donnée</text></svg>'
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


def dotplot_svg(dot:dict,x_name='QUERY',y_name='RÉFÉRENCE'):
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
    data=json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    page=r'''<!doctype html>
<html lang="fr"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>AlignView · Sequence Explorer</title>
<style>
:root{color-scheme:dark;--bg:#090f1b;--side:#0d1727;--panel:#122038;--edge:#2c405d;--fg:#eff5ff;--muted:#9db2cd;--blue:#7ba9ff;--mint:#5bdbbb;--amber:#ffca7c;--red:#ff8fa1;--accent:#4f94ef;--track:#304259}
:root.light{color-scheme:light;--bg:#edf3fa;--side:#f9fcff;--panel:#fff;--edge:#bdcbde;--fg:#16283f;--muted:#576c86;--blue:#2062be;--mint:#137e68;--amber:#9b580e;--red:#b12c46;--accent:#306cc7;--track:#d1deed}
*{box-sizing:border-box}body{font:14px/1.53 ui-sans-serif,system-ui,-apple-system,Segoe UI,Arial,sans-serif;background:var(--bg);color:var(--fg);margin:0}button,input,select{font:inherit}
.layout{display:grid;grid-template-columns:242px minmax(0,1fr);min-height:100vh}.sidebar{padding:26px 18px;background:var(--side);border-right:1px solid var(--edge);position:sticky;top:0;height:100vh}
.brand{font-size:24px;font-weight:800;letter-spacing:-.9px}.brand b{color:var(--mint)}.version{font-size:11px;letter-spacing:1.3px;text-transform:uppercase;color:var(--muted);margin-top:0}
.nav{display:flex;flex-direction:column;gap:5px;margin-top:32px}.nav button{width:100%;text-align:left;background:transparent;color:var(--muted);border:0;border-radius:9px;padding:12px;cursor:pointer}
.nav button[hidden]{display:none}.nav-extra{margin-top:24px;color:var(--muted)}.nav-extra summary{cursor:pointer}.nav-extra .nav{margin-top:10px}.reading-actions{display:flex;gap:10px;flex-wrap:wrap;margin-top:16px}.figure-view summary{cursor:pointer;color:var(--mint);padding:10px 0}.nav button.active,.nav button:hover{background:var(--panel);color:var(--fg)}.nav button.active{border-left:3px solid var(--mint)}
.main{padding:28px 35px;max-width:1680px;width:100%;min-width:0}.top{display:flex;align-items:flex-start;justify-content:space-between;gap:16px;margin-bottom:22px}.eyebrow{text-transform:uppercase;letter-spacing:1.5px;color:var(--mint);font-size:11px;font-weight:800}
h1{font-size:30px;letter-spacing:-1.1px;margin:2px 0 5px}h2{font-size:19px;letter-spacing:-.3px;margin:0 0 14px}h3{font-size:15px;margin:0 0 8px}p{margin:8px 0 12px}.muted,.hint{color:var(--muted)}.hint{font-size:12px}.section{display:none}.section.active{display:block}
.cards{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:14px;margin-bottom:18px}.card,.panel{border-radius:14px;background:var(--panel);border:1px solid var(--edge);padding:20px;min-width:0}.panel{margin-bottom:17px}
.kicker{text-transform:uppercase;letter-spacing:1.2px;color:var(--muted);font-size:10px;font-weight:700}.value{font-size:27px;font-weight:780;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;margin-top:9px}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:17px}.accent{color:var(--mint)}.row{display:flex;gap:12px;align-items:center;flex-wrap:wrap;margin:10px 0}
button,select,input{background:var(--side);color:var(--fg);border:1px solid var(--edge);border-radius:8px;padding:8px 10px;min-width:0}button{cursor:pointer}button:hover{border-color:var(--accent)}.right{margin-left:auto}
label{color:var(--muted);font-size:12px}input[type=range]{accent-color:var(--mint);padding:0}.search{width:280px;max-width:100%}
.tablewrap{max-height:620px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border-bottom:1px solid var(--edge);padding:9px;text-align:left;vertical-align:middle}th{position:sticky;top:0;background:var(--panel);color:var(--muted);white-space:nowrap}tr[data-name]{cursor:pointer}tr[data-name]:hover{background:color-mix(in srgb,var(--accent) 13%,transparent)}
.mini{font:12px ui-monospace,monospace;color:var(--muted)}.track{height:17px;background:var(--track);position:relative;border-radius:5px;overflow:hidden;min-width:160px}.track span{position:absolute;top:0;bottom:0;background:var(--blue);border-radius:3px}
.pill{font:11px ui-monospace,monospace;display:inline-block;padding:3px 7px;border:1px solid var(--edge);border-radius:6px;color:var(--mint)}
pre,.code{background:var(--side);border:1px solid var(--edge);border-radius:10px;padding:15px;overflow:auto;font:12.5px/1.5 ui-monospace,SFMono-Regular,Consolas,monospace;white-space:pre}
.hspblock{margin:12px 0 20px}.hsp-header{padding:10px 13px;border:1px solid var(--edge);background:var(--side);border-radius:8px}.base-m{color:var(--mint)}.base-v{color:var(--amber)}.base-gap{color:var(--red)}
svg{width:100%;max-height:670px;height:auto}.notice{border-left:3px solid var(--amber);padding:10px 15px;background:color-mix(in srgb,var(--amber) 7%,transparent);font-size:12px;color:var(--muted);border-radius:4px}
.statbar{height:8px;border-radius:10px;background:var(--track);overflow:hidden}.statbar b{display:block;height:100%;background:var(--mint)}
.msa{font:12px/1.85 ui-monospace,SFMono-Regular,Consolas,monospace;overflow-x:auto}
.msa-board{font:13px/1.85 ui-monospace,SFMono-Regular,Consolas,monospace;overflow:auto;max-height:65vh;min-height:130px;background:var(--side);border:1px solid var(--edge);border-radius:10px;padding:12px;outline-offset:3px}
.msa-board .msa-line{white-space:nowrap;display:flex;width:max-content;min-width:100%;align-items:center}
.msa-board .msa-name{display:inline-block;position:sticky;left:0;background:var(--side);color:var(--muted);min-width:225px;width:225px;max-width:225px;overflow:hidden;text-overflow:ellipsis;z-index:1;padding-right:12px;border-right:1px solid var(--edge)}
.msa-board .resid{display:inline-block;width:1.15ch;min-width:1.15ch;text-align:center;cursor:crosshair}
.msa-board .resid:hover{background:var(--accent);color:#fff}
.msa-board .msa-ruler{color:var(--muted)}.msa-board .msa-different{background:color-mix(in srgb,var(--red) 25%,transparent);color:var(--red)}
.msa-line{white-space:pre;display:block}.msa-name{color:var(--muted)}.dna-a{color:#68ceaa}.dna-t{color:#ffb36e}.dna-g{color:#77aafd}.dna-c{color:#d695f0}.dna-gap{color:#f27895}.resid{font-weight:600}
#dotcanvas{display:block;width:100%;background:var(--side);border-radius:10px;border:1px solid var(--edge);max-height:600px;touch-action:none}.footer{color:var(--muted);font-size:11px;margin:30px 0}
@media(max-width:1100px){.cards{grid-template-columns:repeat(2,minmax(0,1fr))}.grid2{grid-template-columns:1fr}.layout{grid-template-columns:180px 1fr}.main{padding:22px}.sidebar{padding:24px 9px}}
@media(max-width:720px){.layout{display:block}.sidebar{height:auto;position:static;padding:12px}.nav{flex-direction:row;flex-wrap:wrap;margin:12px 0}.nav button{width:auto;padding:7px}.main{padding:15px}.cards{grid-template-columns:repeat(2,minmax(0,1fr))}h1{font-size:24px}}
</style></head><body><div class="layout"><aside class="sidebar"><div class="brand">AlignView</div><div class="version">Séquences et alignements · __ALIGNVIEW_VERSION__</div><nav class="nav">
<button data-tab="overview" class="active">Résumé</button><button data-tab="coverage">Correspondances</button>
<button data-tab="sequences">Séquences</button><button data-tab="msa">Alignement</button>
<button data-tab="dotplot">Dotplot</button><button data-tab="provenance">Paramètres</button></nav><details class="nav-extra"><summary>Paramètres et aide</summary><div id="advancedNav" class="nav"></div></details><div class="footer">Rapport autonome<br>Aucun appel réseau · données locales<br><br>Une absence de HSP ne signifie pas une absence d'homologie.</div></aside>
<main class="main"><header class="top"><div><div class="eyebrow">Exploration bioinformatique</div><h1 id="pageTitle">Vue générale</h1><p class="muted" id="subtitle"></p></div><div class="row"><button id="backToSummary" hidden>← Résumé</button><button onclick="document.documentElement.classList.toggle('light')">◐ Thème</button></div></header>
<section id="overview" class="section active"><div class="panel" id="readingGuide"><h2>Par où commencer ?</h2><p class="muted">Ce rapport présente les données enregistrées. Choisir une vue ci-dessous ; les nouveaux calculs se lancent dans l’application AlignView.</p><div id="readingActions" class="reading-actions"></div></div><div id="stats" class="cards"></div><div class="grid2"><article class="panel"><h2 id="summaryHeading">Résumé de la recherche</h2><div id="searchSummary"></div></article><article class="panel"><h2 id="queryHeading">Qualité de la requête</h2><div id="queryStats"></div></article></div>
<div class="panel" id="overviewCoverage"><h2>Couverture de la requête</h2><div id="coverageSvg"></div><p class="hint">Les intervalles sont calculés à partir de l'union des HSP. Les positions d'une requête BLASTX sont nucléotidiques.</p></div>
<div class="panel"><h2>Profils le long de la requête</h2><div id="profile"></div></div></section>
<section id="coverage" class="section"><div class="panel"><h2>Explorer les références</h2><div class="row"><input id="hitSearch" class="search" placeholder="Filtrer les IDs des références…"><label>Trier <select id="sort"><option value="bits">Bitscore</option><option value="coverage">Couverture</option><option value="identity">Identité</option><option value="evalue">E-value</option></select></label><span class="pill" id="hitCount"></span></div>
<div class="tablewrap"><table><thead><tr><th>#</th><th>Référence</th><th>HSP</th><th>Id. max</th><th>Couverture Q</th><th>E-value</th><th>Bits</th><th>Cadres</th></tr></thead><tbody id="hitRows"></tbody></table></div></div>
<div class="panel"><h2>Positions couvertes et non couvertes</h2><p id="selectedName" class="accent"></p><div id="selectedTrack"></div><div id="selectedMissing" class="hint"></div></div>
<div class="panel"><h2>Alignement local, colonne par colonne</h2><div class="row"><label>HSP <select id="hspPicker"></select></label><label><input type="checkbox" id="onlyDiff"> Montrer seulement les blocs contenant des différences</label></div><div id="hspInfo"></div><div id="hspAlign"></div></div><p class="notice">Attention : la couverture est calculée par union des coordonnées des HSP retenus, sans compter plusieurs fois les chevauchements. Identité = proportion dans le HSP, pas identité sur l'ensemble des deux séquences. Une région sans HSP peut néanmoins être homologue.</p></section>
<section id="sequences" class="section"><div class="panel"><h2>Statistiques FASTA</h2><div class="tablewrap"><table><thead><tr><th>ID</th><th>Type probable</th><th>Longueur</th><th>Sans gaps</th><th>% GC</th><th>% N</th><th>Entropie</th></tr></thead><tbody id="seqRows"></tbody></table></div></div><div class="panel"><h2>Profil GC / entropie</h2><p class="hint">Analyse exploratoire en fenêtres glissantes. Le GC ne s'applique qu'aux séquences nucléotidiques.</p><div id="profiles"></div></div></section>
<section id="msa" class="section"><div id="msaBody"></div></section>
<section id="dotplot" class="section"><div class="panel"><h2>Correspondances exactes par k-mers</h2><p class="muted">Un dotplot représente les positions partageant un même mot de longueur k. Les répétitions abondantes sont limitées pour préserver la lisibilité. Il ne remplace pas un alignement.</p><div class="row"><label>Afficher les points : <input type="range" min="1" max="100" value="100" id="dotRange"></label><span id="dotCount" class="pill"></span></div><canvas id="dotcanvas" width="1050" height="620"></canvas><p id="dotNote" class="hint"></p></div></section>
<section id="provenance" class="section"><div class="panel"><h2>Réglages actuels et paramètres des calculs effectués</h2><p>Les blocs search_run et calculation_runs conservent les réglages ayant produit les résultats. Les autres valeurs sont les réglages actuels pour le prochain calcul. null/native signifie : option laissée au moteur installé.</p><pre id="prov"></pre></div><div class="panel"><h2>Commandes externes exécutées</h2><pre id="history"></pre></div><p class="notice">Les fichiers complets TSV/JSON, alignements FASTA, cartes SVG et logs sont disponibles uniquement si la session a été enregistrée.</p></section>
<footer class="footer">AlignView · Rapport exploratoire reproductible · Tous les calculs sont locaux</footer></main></div>
<script>
const D=__DATA__;
const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const num=x=>Number(x||0).toLocaleString('fr-FR');
const percent=(v,d=1)=>Number(v||0).toFixed(d)+' %';
const nfmt=xs=>xs&&xs.length?xs.map(([a,b])=>`${num(a)}–${num(b)} (${num(b-a+1)})`).join(', '):'aucune';
const titles={overview:'Vue générale',coverage:'Couvertures & HSP',sequences:'Séquences & qualité',msa:'Alignement multiple',dotplot:'Dotplot',provenance:'Reproductibilité'};
let activeTab='overview';const tabScroll={};
function tab(name,record=true){
 if(!/^[a-z][a-z0-9_-]*$/i.test(String(name)))name='overview';
 const button=document.querySelector('[data-tab="'+name+'"]'),section=document.getElementById(name);
 if(!button||button.hidden||!section)name='overview';
 if(activeTab!==name)tabScroll[activeTab]=window.scrollY;
 document.querySelectorAll('[data-tab]').forEach(x=>{const active=x.dataset.tab===name;x.classList.toggle('active',active);if(active)x.setAttribute('aria-current','page');else x.removeAttribute('aria-current');});
 document.querySelectorAll('.section').forEach(x=>x.classList.toggle('active',x.id===name));
 const parentNav=document.querySelector('[data-tab="'+name+'"]').closest('details');if(parentNav)parentNav.open=true;
 document.getElementById('pageTitle').textContent=titles[name]||name;
 document.getElementById('backToSummary').hidden=name==='overview';
 activeTab=name;window.scrollTo(0,tabScroll[name]||0);
 if(record&&location.hash!=='#'+name){try{history.pushState(null,'','#'+name);}catch(error){/* Local report remains navigable without history access. */}}
 if(name==='dotplot')drawDot();
}
window.addEventListener('popstate',()=>tab(location.hash.slice(1)||'overview',false));
document.getElementById('backToSummary').onclick=()=>tab('overview');
document.querySelectorAll('[data-tab]').forEach(x=>x.onclick=()=>tab(x.dataset.tab));
const metric=(name,value)=>`<div class="card"><div class="kicker">${esc(name)}</div><div class="value">${esc(value)}</div></div>`;
const S=D.queryStats||{};const entries=D.coverage||[];
const cover=D.combined||[];const cov=cover.reduce((s,[a,b])=>s+b-a+1,0);
document.getElementById('subtitle').textContent=`${D.mode} · ${D.queryName} · ${entries.length} références retenues`;
document.getElementById('stats').innerHTML=metric('Requête',num(D.queryLen)+' '+D.units.query)+metric('Couverture union',percent(D.queryLen?100*cov/D.queryLen:0))+metric('Références',num(entries.length))+metric('HSP retenus',num(D.hsps.length));
document.getElementById('searchSummary').innerHTML=`<p><b>Mode :</b> ${esc(D.mode)}</p><p><b>Références au total :</b> ${num(D.referenceCount)}</p><p><b>E-value max :</b> ${esc((D.settings.search_run?.evalue??D.settings.evalue)||'—')}</p><p><b>Identité minimale :</b> ${esc(D.settings.identity||0)} %</p><p><b>Couverture minimale (union) :</b> ${esc(D.settings.coverage||0)} %</p>`;
document.getElementById('queryStats').innerHTML=`<p><b>Type :</b> ${esc(S.kind||'—')}</p><p><b>GC :</b> ${S.gc_pct===null?'non applicable':percent(S.gc_pct)}</p><p><b>Ambiguïtés :</b> ${esc(S.ambiguous??'—')}</p><p><b>Entropie :</b> ${esc(S.entropy??'—')} bits</p>`;
document.getElementById('coverageSvg').innerHTML=D.coverageSVG||'<p class="muted">Aucune recherche de similarité.</p>';
document.getElementById('profile').innerHTML=D.profileSVG||'<p class="muted">Aucune donnée.</p>';
const localSearch=!!(D.hsps?.length||D.settings?.search_completed||D.referenceCount&&D.mode!=='msa'&&D.mode!=='vsearch'&&!D.mode.startsWith('global-'));
const tabLabels={overview:'Résumé',coverage:'Correspondances',sequences:'Séquences',msa:'Alignement',dotplot:'Dotplot',figures:'Figures',exports:'Exporter',provenance:'Paramètres',guide:'Comprendre'};
const availableTabs={coverage:localSearch,msa:!!D.msa,dotplot:!!D.dotplot};
document.querySelectorAll('[data-tab]').forEach(button=>{
 const name=button.dataset.tab;button.textContent=tabLabels[name]||button.textContent.replace(/^\d+\s*·\s*/,'');
 if(availableTabs[name]===false)button.hidden=true;
 if(name==='provenance'||name==='guide')document.getElementById('advancedNav').appendChild(button);
});
const actions=[['sequences','Consulter les séquences']];
if(D.msa)actions.unshift(['msa','Lire l’alignement']);
if(localSearch)actions.push(['coverage','Examiner les correspondances']);
if(document.getElementById('exports'))actions.push(['exports','Exporter les résultats']);
document.getElementById('readingActions').innerHTML=actions.map(([name,label])=>`<button data-open-tab="${name}">${label}</button>`).join('');
document.querySelectorAll('[data-open-tab]').forEach(button=>button.onclick=()=>tab(button.dataset.openTab));
if(!localSearch){
 document.getElementById('overviewCoverage').hidden=true;
 document.getElementById('summaryHeading').textContent='Données et résultats disponibles';
 document.getElementById('queryHeading').textContent='Qualité de la première séquence';
 document.getElementById('searchSummary').innerHTML=`<p>${num(D.msa?.metrics?.n||D.sequenceStats.length)} séquences intégrées au rapport.</p><p>${D.msa?'Un alignement est disponible dans la visionneuse.':'Aucun alignement multiple n’a été calculé.'}</p>`;
 document.getElementById('stats').innerHTML=metric('Séquences',num(D.msa?.metrics?.n||D.sequenceStats.length))+metric('Alignement',D.msa?num(D.msa.metrics.length)+' colonnes':'À calculer')+(D.msa?metric('Colonnes variables',num(D.msa.metrics.variable)):metric('Type',S.kind||'—'))+metric('Fichier','Rapport autonome');
}
document.addEventListener('DOMContentLoaded',()=>{
 document.querySelectorAll('.figcard .figure-svg').forEach(figure=>{
  const details=document.createElement('details');details.className='figure-view';
  const summary=document.createElement('summary');summary.textContent='Afficher la figure';
  figure.before(details);details.append(summary,figure);
 });
 tab(location.hash.slice(1)||'overview',false);
});
let printedFigures=[];
window.addEventListener('beforeprint',()=>{printedFigures=[...document.querySelectorAll('.figure-view')].filter(x=>!x.open);printedFigures.forEach(x=>x.open=true);});
window.addEventListener('afterprint',()=>printedFigures.forEach(x=>x.open=false));

function bar(spans,len){return `<div class="track">${(spans||[]).map(([a,b])=>`<span title="${a}–${b}" style="left:${100*(a-1)/Math.max(1,len)}%;width:${Math.max(.15,100*(b-a+1)/Math.max(1,len))}%"></span>`).join('')}</div>`;}
let picked=entries[0]?.subject||'';
function renderHits(){const search=document.getElementById('hitSearch').value.toLowerCase(),sort=document.getElementById('sort').value;
let selected=entries.filter(e=>e.subject.toLowerCase().includes(search));selected.sort((a,b)=>sort==='coverage'?b.query_coverage-a.query_coverage:sort==='identity'?b.best_identity-a.best_identity:sort==='evalue'?a.best_evalue-b.best_evalue:b.best_bitscore-a.best_bitscore);
document.getElementById('hitCount').textContent=selected.length+' références';
document.getElementById('hitRows').innerHTML=selected.map((e,i)=>`<tr data-name="${esc(e.subject)}"><td>${i+1}</td><td><b>${esc(e.subject)}</b></td><td>${e.hsps}</td><td>${percent(e.best_identity)}</td><td>${bar(e.query_segments,D.queryLen)}<span class="mini">${percent(e.query_coverage)}</span></td><td class="mini">${Number(e.best_evalue).toExponential(2)}</td><td>${e.best_bitscore.toFixed(1)}</td><td>${esc(e.frames.join(',')||'—')}</td></tr>`).join('');
document.querySelectorAll('tr[data-name]').forEach(tr=>tr.onclick=()=>{picked=tr.dataset.name;showSelected();});}
document.getElementById('hitSearch').oninput=renderHits;document.getElementById('sort').onchange=renderHits;
function showSelected(){const e=entries.find(x=>x.subject===picked);if(!e){document.getElementById('selectedTrack').innerHTML='Aucune référence';return;}
document.getElementById('selectedName').textContent=e.subject;
document.getElementById('selectedTrack').innerHTML=`<p>QUERY · ${num(D.queryLen)} ${D.units.query} · ${percent(e.query_coverage)}</p>${bar(e.query_segments,D.queryLen)}<p>RÉFÉRENCE · ${num(e.reference_length)} ${D.units.subject} · ${e.reference_coverage===null?'inconnue':percent(e.reference_coverage)}</p>${bar(e.reference_segments,e.reference_length)}`;
document.getElementById('selectedMissing').innerHTML=`<p><b>QUERY sans HSP :</b> ${esc(nfmt(e.query_uncovered))}</p><p><b>RÉFÉRENCE sans HSP :</b> ${e.reference_length?esc(nfmt(e.reference_uncovered)):'longueur inconnue'}</p>`;
const hs=D.hsps.filter(h=>h.subject===picked);const picker=document.getElementById('hspPicker');picker.innerHTML=hs.map((h,i)=>`<option value="${i}">HSP ${i+1} · ${num(h.qstart)}–${num(h.qend)} (Q, ${D.units.query})</option>`).join('');picker.onchange=showHsp;showHsp();}
function colorLine(q,s,which){let v='';for(let i=0;i<q.length;i++){const c=which==='q'?q[i]:s[i];const cl=q[i]==='-'||s[i]==='-'?'base-gap':q[i]===s[i]?'base-m':'base-v';v+=`<span class="${cl}">${esc(c)}</span>`;}return v;}
function showHsp(){const hs=D.hsps.filter(h=>h.subject===picked),which=Number(document.getElementById('hspPicker').value)||0,h=hs[which];if(!h){document.getElementById('hspInfo').innerHTML='Aucun HSP';document.getElementById('hspAlign').innerHTML='';return;}
document.getElementById('hspInfo').innerHTML=`<p class="mini">Q ${num(h.qstart)}–${num(h.qend)} ${D.units.query} · S ${num(h.sstart)}–${num(h.send)} ${D.units.subject} · alignement ${D.units.aligned} · identité ${percent(h.identity)} · bits ${h.bitscore} · E ${Number(h.evalue).toExponential(2)} · cadres Q ${h.qframe} / S ${h.sframe}</p><p class="mini">${num(h.matches)} identités · ${num(h.mismatches)} substitutions · ${num(h.gap_columns)} colonnes avec gap</p>`;
const q=h.qseq||'',s=h.sseq||'';const diff=document.getElementById('onlyDiff').checked;let blocks='';const n=Math.min(q.length,s.length);
for(let i=0;i<n;i+=70){let qa=q.slice(i,i+70),sa=s.slice(i,i+70);if(diff && [...qa].every((c,j)=>c===sa[j]))continue;
let marks=[...qa].map((c,j)=>c===sa[j]&&c!=='-'?'|':c==='-'||sa[j]==='-'?' ':'.').join('');
blocks+=`<div class="hspblock"><pre>Q  ${colorLine(qa,sa,'q')}\n   ${esc(marks)}\nS  ${colorLine(qa,sa,'s')}</pre><span class="hint">Colonnes ${i+1}–${Math.min(n,i+70)} (HSP, non coordonnées génomiques)</span></div>`;}
document.getElementById('hspAlign').innerHTML=blocks||'<p class="muted">Aucune différence dans les blocs affichés.</p>';
}
document.getElementById('onlyDiff').onchange=showHsp;renderHits();showSelected();
const seqstats=D.sequenceStats||[];
document.getElementById('seqRows').innerHTML=seqstats.map(s=>`<tr><td>${esc(s.id)}</td><td>${esc(s.kind)}</td><td>${num(s.length)} ${s.kind==='ADN/ARN'?'nt':'aa'}</td><td>${num(s.ungapped)} ${s.kind==='ADN/ARN'?'nt':'aa'}</td><td>${s.gc_pct===null?'—':percent(s.gc_pct)}</td><td>${percent(s.n_pct)}</td><td>${s.entropy}</td></tr>`).join('');
document.getElementById('profiles').innerHTML=D.profileSVG||'<p>Aucune mesure.</p>';
function msa(){
 const m=D.msa;
 if(!m){document.getElementById('msaBody').innerHTML='<div class="panel"><h2>Pas de MSA dans cette session</h2><p>Lancer MAFFT ou charger un alignement pour explorer les colonnes et les sites variables.</p></div>';return;}
 const list=m.seqs||[], metrics=m.metrics||{}, L=metrics.length||0, positions=metrics.positions||[];
 const variants=positions.filter(x=>x.variable).map(x=>x.position-1);
 const body=document.getElementById('msaBody');
 body.innerHTML=`<div class="cards">${metric('Séquences',metrics.n)}${metric('Colonnes',num(L))}${metric('Sites variables',num(metrics.variable))}${metric('Sites informatifs',num(metrics.parsimony_informative))}</div>
 <div class="panel"><h2>Conservation des colonnes</h2>${m.svg||''}</div>
 <div class="panel"><h2>MSA interactif</h2>
 <p class="hint">Déplacement : flèches ← → pour avancer/reculer d'une colonne, Page précédente/suivante pour des blocs. Clique sur un résidu pour inspecter sa colonne. Les coordonnées sont celles du MSA, gaps inclus.</p>
 <div class="row">
 <button id="msaPrev" title="Page précédente">◀ −80</button><button id="msaPrevOne" title="Une colonne en arrière">−1</button>
 <label>Colonne du MSA <input type="number" min="1" max="${L}" value="1" id="msaStart" style="width:100px"></label>
 <button id="msaNextOne" title="Une colonne en avant">+1</button><button id="msaNext" title="Page suivante">+80 ▶</button>
 <button id="msaFirst" title="Début de l’alignement">Début</button><button id="msaLast" title="Fin de l’alignement">Fin</button>
 <label>Colonnes <select id="msaWidth"><option>80</option><option>160</option><option>320</option></select></label>
 <button id="msaNextVar" title="Prochain site variable">Variable suivante</button>
 </div>
 <div class="row"><input type="range" min="1" max="${L}" value="1" id="msaSlider" style="width:min(440px,100%)">
 <label><input type="checkbox" id="msaVar"> Sites variables uniquement</label>
 <label><input type="checkbox" id="msaDiff"> Différences au consensus</label>
 </div>
 <div class="row"><label>Filtrer les séquences <input type="search" id="msaSeqFilter" placeholder="Identifiant…" style="width:230px"></label>
 <button id="msaRowsPrev" title="Séquences précédentes">↑ 100</button><button id="msaRowsNext" title="Séquences suivantes">↓ 100</button>
 <span id="msaPage" class="pill"></span><span id="msaCell" class="pill"></span></div>
 <div id="msaView" class="msa-board" tabindex="0" aria-label="Alignement multiple navigable au clavier"></div>
 <p class="hint">Raccourcis dans la visionneuse : ← → = une colonne · PageUp/PageDown = un bloc · v = sites variables · d = différences. Le filtrage masque des colonnes mais ne modifie pas l'alignement.</p></div>`;
 let start=0;
 const startInput=document.getElementById('msaStart'),slider=document.getElementById('msaSlider');
 const variableCheckbox=document.getElementById('msaVar'),diffCheckbox=document.getElementById('msaDiff');
 const filterInput=document.getElementById('msaSeqFilter'),board=document.getElementById('msaView');
 let pageWidth=80,rowStart=0;
 function displayed(){return variableCheckbox.checked ? variants : null;}
 function anchorToIndex(k){
  const cols=displayed();if(!cols)return Math.max(0,Math.min(L-1,k));
  if(!cols.length)return 0;
  let left=0,right=cols.length;
  while(left<right){let mid=(left+right)>>1;if(cols[mid]<k)left=mid+1;else right=mid;}
  return Math.min(cols.length-1,left);
 }
 function indexToAnchor(k){const cols=displayed();return cols ? (cols[Math.max(0,Math.min(cols.length-1,k))]||0) : k;}
 function go(real){start=anchorToIndex(real);draw();}
 const clazz=c=>({'A':'dna-a','T':'dna-t','U':'dna-t','G':'dna-g','C':'dna-c','-':'dna-gap','*':'dna-gap'}[c]||'');
 function draw(){
  const savedTop=board.scrollTop,savedLeft=board.scrollLeft;
  const cols=displayed();const count=cols?cols.length:L;
  start=Math.max(0,Math.min(Math.max(0,count-pageWidth),start));
  const idx=Array.from({length:Math.min(pageWidth,Math.max(0,count-start))},(_,i)=>cols?cols[start+i]:start+i);
  const highlighted=diffCheckbox.checked;
  const paint=(name,seq)=>`<div class="msa-line"><span class="msa-name" title="${esc(name)}">${esc(name.slice(0,26))}</span>`+
   idx.map(j=>{const c=seq[j]||'-',dif=highlighted && c!=='-' && c!==(metrics.consensus[j]||'-');
   return `<span class="resid ${clazz(c)} ${dif?'msa-different':''}" data-column="${j}" title="${esc(name)} · colonne ${j+1} · ${esc(c)} · consensus ${esc(metrics.consensus[j])}">${esc(c)}</span>`;}).join('')+'</div>';
  const tick=(j)=>j%10===0?String(j+1).slice(-1):'·';
  const ticks=`<div class="msa-line"><span class="msa-name">POSITION</span>${idx.map(j=>`<span class="resid msa-ruler" title="Position ${j+1}">${tick(j)}</span>`).join('')}</div>`;
  const vard=`<div class="msa-line"><span class="msa-name">VARIABLES</span>${idx.map(j=>`<span class="resid msa-ruler" title="${j+1}: ${positions[j]?.variable?'variable':'non variable'}">${positions[j]?.variable?'▲':'·'}</span>`).join('')}</div>`;
  const query=filterInput.value.trim().toLowerCase();let visible=list.filter(s=>s.id.toLowerCase().includes(query));
  rowStart=Math.max(0,Math.min(Math.floor(Math.max(0,visible.length-1)/100)*100,rowStart));
  const rows=visible.slice(rowStart,rowStart+100);
  board.innerHTML=ticks+rows.map(s=>paint(s.id,s.seq)).join('')+paint('CONSENSUS',metrics.consensus)+vard;
  board.scrollTop=savedTop;board.scrollLeft=savedLeft;
  startInput.value=indexToAnchor(start)+1;slider.value=startInput.value;
  document.getElementById('msaPage').textContent=`${idx.length ? idx[0]+1 : '—'} – ${idx.length ? idx[idx.length-1]+1 : '—'} / ${num(L)} · séquences ${visible.length ? rowStart+1 : 0}–${Math.min(visible.length,rowStart+100)}/${visible.length} intégrées (total ${metrics.n})`;
  document.getElementById('msaCell').textContent=(cols?'Variables : ':'Toutes : ')+num(count)+' colonnes';
 }
 function step(n){const count=displayed()?.length||L;start=Math.max(0,Math.min(Math.max(0,count-1),start+n));draw();}
 document.getElementById('msaPrev').onclick=()=>step(-pageWidth);
 document.getElementById('msaNext').onclick=()=>step(pageWidth);
 document.getElementById('msaPrevOne').onclick=()=>step(-1);
 document.getElementById('msaNextOne').onclick=()=>step(1);
 document.getElementById('msaNextVar').onclick=()=>{const current=indexToAnchor(start);let next=variants.find(i=>i>current);go(next===undefined?(variants[0]||0):next);};
 startInput.onchange=()=>go(Math.max(0,Math.min(L-1,(Number(startInput.value)||1)-1)));
 slider.oninput=()=>go(Number(slider.value)-1);
 variableCheckbox.onchange=()=>go(Number(startInput.value)-1);
 diffCheckbox.onchange=draw;filterInput.oninput=()=>{rowStart=0;board.scrollTop=0;draw();};
 document.getElementById('msaWidth').onchange=e=>{pageWidth=Number(e.target.value);draw();};
 document.getElementById('msaFirst').onclick=()=>go(0);
 document.getElementById('msaLast').onclick=()=>go(L-1);
 document.getElementById('msaRowsPrev').onclick=()=>{rowStart-=100;board.scrollTop=0;draw();};
 document.getElementById('msaRowsNext').onclick=()=>{rowStart+=100;board.scrollTop=0;draw();};
 board.addEventListener('wheel',e=>{
  if(e.shiftKey||Math.abs(e.deltaX)>Math.abs(e.deltaY)){
   e.preventDefault();step(Math.sign(e.deltaX||e.deltaY)*Math.max(1,Math.round(pageWidth/8)));
  }
 },{passive:false});
 board.onclick=e=>{const el=e.target.closest('[data-column]');if(!el)return;const j=Number(el.dataset.column);let row=el.parentElement.querySelector('.msa-name')?.textContent||'';
  document.getElementById('msaCell').textContent=`${row} · colonne ${j+1} · ${el.textContent} · consensus ${metrics.consensus[j]} · conservation ${positions[j]?.conservation}% · gaps ${positions[j]?.gap_pct}%`;
 };
 board.onkeydown=e=>{
  if(e.key==='ArrowRight'){e.preventDefault();step(e.ctrlKey?pageWidth:e.shiftKey?10:1);}
  else if(e.key==='ArrowLeft'){e.preventDefault();step(e.ctrlKey?-pageWidth:e.shiftKey?-10:-1);}
  else if(e.key==='Home'){e.preventDefault();go(0);}
  else if(e.key==='End'){e.preventDefault();go(L-1);}
  else if(e.key==='PageDown'){e.preventDefault();step(pageWidth);}
  else if(e.key==='PageUp'){e.preventDefault();step(-pageWidth);}
  else if(e.key.toLowerCase()==='v'){variableCheckbox.checked=!variableCheckbox.checked;variableCheckbox.onchange();}
  else if(e.key.toLowerCase()==='d'){diffCheckbox.checked=!diffCheckbox.checked;draw();}
 };
 draw();
} msa();
const dot=D.dotplot;
function drawDot(){const canvas=document.getElementById('dotcanvas'),ctx=canvas.getContext('2d');const W=canvas.width,H=canvas.height;ctx.clearRect(0,0,W,H);ctx.fillStyle='#101d30';ctx.fillRect(0,0,W,H);ctx.strokeStyle='#3d506d';
for(let i=0;i<6;i++){const x=55+i*(W-85)/5,y=25+i*(H-70)/5;ctx.beginPath();ctx.moveTo(x,25);ctx.lineTo(x,H-45);ctx.moveTo(55,y);ctx.lineTo(W-30,y);ctx.stroke();}
if(!dot){document.getElementById('dotNote').textContent='Aucun dotplot calculé. Créer un dotplot depuis le terminal.';return;}
const cap=Math.floor(dot.points.length*Number(document.getElementById('dotRange').value)/100);ctx.fillStyle='#58d9bd';const step=Math.max(1,Math.ceil(dot.points.length/Math.max(1,cap)));
for(let i=0;i<dot.points.length;i+=step){const [a,b]=dot.points[i];ctx.fillRect(55+(a/dot.xlen)*(W-85),25+(b/dot.ylen)*(H-70),2,2);}
ctx.fillStyle='#c9dbf1';ctx.font='13px sans-serif';ctx.fillText('QUERY '+(dot.xunit||D.units.query)+' · '+dot.xlen+' positions',W/2-80,H-15);
ctx.save();ctx.translate(15,H/2+90);ctx.rotate(-Math.PI/2);ctx.fillText('RÉFÉRENCE '+(dot.yunit||D.units.subject)+' · '+dot.ylen+' positions',0,0);ctx.restore();
document.getElementById('dotCount').textContent=cap.toLocaleString('fr-FR')+' points';document.getElementById('dotNote').textContent=`k = ${dot.k} · ${dot.note||'Unités natives'} · données ${dot.truncated?'plafonnées':'complètes'} · ces correspondances exactes ne constituent pas des alignements.`;
}
document.getElementById('dotRange').oninput=drawDot;
document.getElementById('prov').textContent=JSON.stringify(D.settings,null,2)+'\n\nSHA256 QUERY : '+(D.shaQuery||'—')+'\nSHA256 REFS : '+(D.shaRefs||'—');
document.getElementById('history').textContent=(D.history||[]).map(x=>`$ ${x.command}\n  version ${x.version||"non disponible"} · durée ${x.duration_s}s · statut ${x.exit_code} · log ${x.stderr_log}${x.stdout?"\n  stdout : "+x.stdout:""}`).join('\n\n')||'Aucune commande externe exécutée.';
</script></body></html>'''
    return page.replace('__DATA__',data)


def build_report(out:Path,query:Seq,ref_seqs:list[Seq],hsps:list[HSP],mode:str,settings:dict,
                 history:list,msa_seqs:list[Seq]|None=None,dot:dict|None=None,
                 sha_query='',sha_refs=''):
    refs_len={r.id:len(r.seq) for r in ref_seqs}
    rows=coverage_data(hsps,len(query.seq),refs_len)
    combined=merge_intervals([h.qspan for h in hsps],len(query.seq))
    u=units_for(mode,sequence_stats(query)['kind'])
    svg=coverage_svg(len(query.seq),rows,combined,unit=u.query)
    (out/'coverage.svg').write_text(svg,encoding='utf-8')
    kind=sequence_stats(query)['kind']
    metric='gc' if kind=='ADN/ARN' else 'entropy'
    scores=sliding_scores(query.seq,window=min(101,len(query.seq)),step=max(1,min(101,len(query.seq))//4),kind=metric)
    profile=plot_svg(scores,title='Profil '+('GC' if metric=='gc' else 'entropie')+f' · positions QUERY ({u.query})',ylabel=metric)
    (out/'profile.svg').write_text(profile,encoding='utf-8')
    msa=None
    if msa_seqs:
        stats=msa_metrics(msa_seqs)
        vals=[dict(center=p['position'],value=p['conservation']) for p in stats['positions']]
        conservation=plot_svg(vals,'Conservation du MSA (%)')
        (out/'msa_conservation.svg').write_text(conservation,encoding='utf-8')
        msa=dict(seqs=[dict(id=s.id,seq=s.seq) for s in msa_seqs[:1000]],metrics=stats,svg=conservation)
    if dot: (out/'dotplot.svg').write_text(dotplot_svg(dot),encoding='utf-8')
    u=units_for(mode,kind)
    payload=dict(units=dict(query=u.query,subject=u.subject,aligned=u.aligned,msa=u.msa),queryName=query.id,queryLen=len(query.seq),queryStats=sequence_stats(query),
                 referenceCount=len(ref_seqs),mode=mode,settings=settings,history=history,
                 coverage=rows,combined=combined,coverageSVG=svg,profileSVG=profile,
                 hsps=[h.record() for h in hsps[:2000]],sequenceStats=[sequence_stats(query)]+[sequence_stats(r) for r in ref_seqs[:300]],
                 msa=msa,dotplot=dot,shaQuery=sha_query,shaRefs=sha_refs)
    from av_version import VERSION
    report=out/'report.html';report.write_text(report_html(payload).replace('__ALIGNVIEW_VERSION__',VERSION),encoding='utf-8')
    return report
