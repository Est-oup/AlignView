"""Exports scientifiques documentés : formats ouverts, coordonnées explicites, provenance.

Aucun export n'est persistant tant qu'une destination n'est pas explicitement demandée.
Les positions HSP entrantes sont 1-based inclusives; les BED sortants sont 0-based half-open.
"""
from __future__ import annotations
import csv
import datetime
import hashlib
import json
import math
import platform
import re
import shutil
import sys
from pathlib import Path
from urllib.parse import quote
from av_version import VERSION
from av_units import units_for,coord_audit
from av_core import (Seq, HSP, coverage_data, sequence_stats, msa_metrics, distance_matrix,
                     merge_intervals, missing_intervals, write_fasta, export_table, sha256)

PRESETS={
    'standard':'FASTA, QC, résultats par référence, HSP, couvertures, BED et provenance',
    'publication':'Standard + CSV, GFF3, MSA (NEXUS/PHYLIP), distances et figures',
    'complete':'Publication + données JSON/JSONL détaillées, régions sans HSP et inventaire',
}

def _json(path, payload):
    path.write_text(json.dumps(payload,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def _table(base:Path, rows:list[dict], kind='tsv', columns=None):
    if not rows and not columns: return None
    filename=base.with_suffix('.csv' if kind=='csv' else '.tsv')
    export_table(filename,rows,fieldnames=columns,sep=',' if kind=='csv' else '\t')
    return filename

def _sanitize(label): return re.sub(r'[^A-Za-z0-9_.-]+','_',label).strip('._') or 'seq'

def _gff(value): return quote(str(value),safe='-_.~:')

def _quality_report(records):
    types={r['kind'] for r in records}
    warnings=[]
    if len(types)>1:warnings.append('Types ADN/protéines mixtes ou ambigus : vérifier les fichiers FASTA')
    if len({r['id'] for r in records})!=len(records):warnings.append('Identifiants identiques dans les deux jeux FASTA : les exporter avec leur source')
    for r in records:
        if r.get('n_pct',0)>5:warnings.append(f"{r['id']}: plus de 5 % de N")
    return warnings

def _nexus(records:list[Seq],path:Path):
    if not records:return
    length=len(records[0].seq)
    if any(len(s.seq)!=length for s in records):raise ValueError('NEXUS : alignement requis')
    from av_core import sequence_type
    types={sequence_type(s.seq) for s in records}
    if types not in ({'ADN/ARN'},{'protéine'}):raise ValueError('NEXUS : types de séquences mélangés, export impossible')
    molecule='DNA' if types=={'ADN/ARN'} else 'PROTEIN'
    with path.open('w',encoding='utf-8') as out:
        out.write(f'#NEXUS\nBEGIN DATA;\n DIMENSIONS NTAX={len(records)} NCHAR={length};\n FORMAT DATATYPE={molecule} MISSING=? GAP=-;\n MATRIX\n')
        for r in records:out.write("'"+r.id.replace("'","''")+"' "+r.seq.replace('.','-')+'\n')
        out.write(' ;\nEND;\n')

def _phylip(records:list[Seq],path:Path):
    if not records:return
    n=len(records[0].seq)
    if any(len(s.seq)!=n for s in records):raise ValueError('PHYLIP : longueurs incohérentes')
    # Format PHYLIP relaxed, pas d'identifiants tronqués à 10 caractères.
    ids=[];seen=set()
    for s in records:
        key=_sanitize(s.id)
        base=key;idx=2
        while key in seen:key=f'{base}_{idx}';idx+=1
        ids.append(key);seen.add(key)
    with path.open('w',encoding='utf-8') as f:
        f.write(f'{len(records)} {n}\n')
        for key,s in zip(ids,records):f.write(f'{key} {s.seq.replace(".","-")}\n')
    export_table(path.with_name('msa_phylip_ids.tsv'),[dict(original=s.id,phylip=k) for k,s in zip(ids,records)])

def _stockholm(records:list[Seq],path:Path):
    """Stockholm 1.0 : recommandé pour des MSA de protéines et les outils HMMER."""
    if not records:return
    if len({len(s.seq) for s in records})!=1:raise ValueError('Stockholm : MSA de longueurs identiques requis')
    aliases=[];used=set()
    for r in records:
        key=_sanitize(r.id)
        base=key;n=2
        while key in used:key=f'{base}_{n}';n+=1
        used.add(key);aliases.append(key)
    with path.open('w',encoding='utf-8') as f:
        f.write('# STOCKHOLM 1.0\n')
        for name,r in zip(aliases,records):f.write(f'{name} {r.seq.replace(".","-")}\n')
        f.write('//\n')
    export_table(path.with_name('msa_stockholm_ids.tsv'),[dict(original=r.id,stockholm=name) for r,name in zip(records,aliases)])


def _clustal(records:list[Seq],path:Path):
    """Format CLUSTAL standard, sans calcul artificiel de consensus."""
    if not records:return
    length=len(records[0].seq)
    if any(len(r.seq)!=length for r in records):raise ValueError('CLUSTAL : séquences pas alignées')
    names=[];used=set()
    for s in records:
        base=_sanitize(s.id)[:37];name=base;k=2
        while name in used:name=(base[:33]+f'_{k}');k+=1
        used.add(name);names.append(name)
    with path.open('w',encoding='utf-8') as f:
        f.write('CLUSTAL AlignView 5.0 multiple sequence alignment\n\n')
        for start in range(0,length,70):
            for name,s in zip(names,records):f.write(f'{name:<40}{s.seq[start:start+70].replace(".","-")}\n')
            f.write('\n')
    export_table(path.with_name('msa_clustal_ids.tsv'),[dict(original=r.id,clustal=name) for r,name in zip(records,names)])


def _write_bed(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        f.write('# BED 0-based half-open ; source : HSP BLAST/DIAMOND/MMseqs2\n')
        for r in rows:f.write('\t'.join(str(x) for x in r)+'\n')

def generate_exports(session,profile='standard'):
    """Génère une arborescence temporaire. Aucune écriture permanente."""
    if profile not in PRESETS:raise ValueError('Profil export invalide : '+profile)
    folder=session.work/'exports'
    if folder.exists():shutil.rmtree(folder)
    folder.mkdir(parents=True)
    rich=profile in ('publication','complete')
    extended=profile=='complete'
    a=session
    query=a.query
    allseq=[*a.query_list,*a.references] if not a.single else list(a.query_list)
    qc=[]
    for origin,records in [('query',a.query_list),('reference',a.references)]:
        for s in records:qc.append({'source':origin,**sequence_stats(s)})
    _table(folder/'sequence_qc',qc)
    if rich:_table(folder/'sequence_qc',qc,'csv');_json(folder/'sequence_qc.json',qc)
    write_fasta([query],folder/'query_active.fasta')
    if a.references:write_fasta(a.references,folder/'references.fasta')
    all_hsp=[h.record() for h in a.hsps]
    selected_hsp=[h.record() for hits in a.selected.values() for h in hits]
    hsp_cols=list(all_hsp[0]) if all_hsp else ['query','subject','identity','align_length','qstart','qend','sstart','send','evalue','bitscore','qframe','sframe','qseq','sseq','matches','mismatches','gap_columns']
    _table(folder/'hsps_all',all_hsp,columns=hsp_cols)
    _table(folder/'hsps_selected',selected_hsp,columns=hsp_cols)
    u=units_for(a.mode,sequence_stats(query)['kind'])
    _table(folder/'coordinate_audit',coord_audit([h for group in a.selected.values() for h in group],a.mode),
           columns=['subject','qstart','qend','sstart','send','query_unit','reference_unit','alignment_unit',
                    'query_span','reference_span','query_aligned_chars','reference_aligned_chars',
                    'expected_query_span','expected_reference_span','query_consistent','reference_consistent',
                    'qframe','sframe'])
    _table(folder/'coordinate_units',[dict(mode=a.mode,query_native_unit=u.query,reference_native_unit=u.subject,
           hsp_alignment_columns=u.aligned,msa_alignment_columns=u.msa,query_translated=u.translated_query,
           reference_translated=u.translated_subject)],columns=['mode','query_native_unit','reference_native_unit',
           'hsp_alignment_columns','msa_alignment_columns','query_translated','reference_translated'])
    if rich:
        _table(folder/'hsps_selected',selected_hsp,'csv',hsp_cols)
        _json(folder/'hsps_selected.json',selected_hsp)
    if extended:
        _json(folder/'hsps_all.json',all_hsp)
        with (folder/'hsps_all.jsonl').open('w',encoding='utf-8') as f:
            for h in all_hsp:f.write(json.dumps(h,ensure_ascii=False,allow_nan=False)+'\n')
    refs={r.id:len(r.seq) for r in a.references}
    selected_raw=[h for hs in a.selected.values() for h in hs]
    coverage=coverage_data(selected_raw,len(query.seq),refs)
    cov_rows=[]
    for r in coverage:
        cov_rows.append({key:r.get(key) for key in ('subject','hsps','query_coverage','reference_length','reference_coverage','best_identity','best_bitscore','best_evalue')})
    _table(folder/'hits_by_reference',cov_rows,columns=['subject','hsps','query_coverage','reference_length','reference_coverage','best_identity','best_bitscore','best_evalue'])
    _table(folder/'coverage_regions',[{'subject':r['subject'],'side':side,'status':status,'start_1based':s,'end_1based':e,'length':e-s+1} for r in coverage for side in ('query','reference') for status,field in (('covered',f'{side}_segments'),('without_detected_hsp',f'{side}_uncovered')) for s,e in r[field]],columns=['subject','side','status','start_1based','end_1based','length'])
    if rich:_table(folder/'hits_by_reference',cov_rows,'csv');_json(folder/'coverage_regions.json',coverage)
    qbed=[];sbed=[]
    for h in selected_raw:
        q0,q1=h.qspan;s0,s1=h.sspan
        qbed.append((h.query,q0-1,q1,h.subject,min(1000,max(0,round(h.bitscore))),'-' if h.qstart>h.qend else '+'))
        sbed.append((h.subject,s0-1,s1,h.query,min(1000,max(0,round(h.bitscore))),'-' if h.sstart>h.send else '+'))
    _write_bed(folder/'hsp_query.bed',qbed)
    _write_bed(folder/'hsp_reference.bed',sbed)
    if extended:
        _write_bed(folder/'query_without_hsp.bed',[(query.id,s-1,e,'without_HSP',0,'.') for s,e in missing_intervals([h.qspan for h in selected_raw],len(query.seq))])
        _write_bed(folder/'references_without_hsp.bed',[(r['subject'],s-1,e,'without_HSP',0,'.') for r in coverage for s,e in r['reference_uncovered']])
    if rich:
        with (folder/'hsp_features.gff3').open('w',encoding='utf-8') as f:
            f.write('##gff-version 3\n')
            for i,h in enumerate(selected_raw,1):
                s,e=h.qspan
                strand='-' if h.qstart>h.qend else '+'
                attributes=f'ID=HSP_{i};Name={_gff(h.subject)};Target={_gff(h.subject)} {h.sspan[0]} {h.sspan[1]};Identity={h.identity:.5g}'
                f.write(f'{_gff(h.query)}\tAlignView\tmatch\t{s}\t{e}\t{h.bitscore:.4g}\t{strand}\t.\t{attributes}\n')
    if a.msa_records:
        msa=a.msa_records
        write_fasta(msa,folder/'msa_alignment.fasta')
        met=msa_metrics(msa)
        _table(folder/'msa_positions',met['positions'])
        write_fasta([Seq('CONSENSUS','CONSENSUS',met['consensus'])],folder/'msa_consensus.fasta')
        variable_positions=[p['position']-1 for p in met['positions'] if p['variable']]
        if variable_positions:write_fasta([Seq(x.id,x.description,''.join(x.seq[i] for i in variable_positions)) for x in msa],folder/'msa_variable_sites.fasta')
        if rich:
            _table(folder/'msa_positions',met['positions'],'csv')
            _json(folder/'msa_metrics.json',{k:v for k,v in met.items() if k!='positions'})
            _nexus(msa,folder/'msa_alignment.nex')
            _phylip(msa,folder/'msa_alignment.phy')
            _stockholm(msa,folder/'msa_alignment.sto')
            _clustal(msa,folder/'msa_alignment.aln')
        if len(msa)<=250:
            dm=distance_matrix(msa)
            _table(folder/'msa_pairwise_distances',[{'seq1':msa[i].id,'seq2':msa[j].id,'p_distance':dm[i][j]} for i in range(len(msa)) for j in range(i+1,len(msa))],columns=['seq1','seq2','p_distance'])
            if rich:_table(folder/'msa_pairwise_distances',[{'seq1':msa[i].id,'seq2':msa[j].id,'p_distance':dm[i][j]} for i in range(len(msa)) for j in range(i+1,len(msa))],'csv')
    for name,rows in [('motif_occurrences',a.motif_rows),('orfs',a.orf_rows)]:
        if rows:_table(folder/name,rows)
    if rich:
        src=a.work/'figures'
        if src.exists():shutil.copytree(src,folder/'figures',dirs_exist_ok=True)
        for name in ('report.html','figures_catalog.tsv','coverage.svg','profile.svg','dotplot.svg'):
            file=a.work/name
            if file.is_file():shutil.copy2(file,folder/name)
    tree=a.work/'tree_fasttree.nwk'
    if tree.exists() and rich:shutil.copy2(tree,folder/'tree_fasttree.nwk')
    params=a.settings()
    meta=dict(version=VERSION,created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              profile=profile,platform=platform.platform(),python=sys.version.split()[0],
              conda_environment=__import__('os').environ.get('CONDA_DEFAULT_ENV',''),
              query_id=query.id,source_query=str(a.query_path),source_references=str(a.refs_path) if a.refs_path else None,
              sha256_inputs={str(p):sha256(p) for p in (a.query_path,a.refs_path) if p},
              search_mode=a.mode,parameters=params,commands=a.tools.history,
              modules=a.tools.loaded,warning=_quality_report(qc),
              definitions={'identite':'Identité locale HSP (pas identité de séquence complète)',
                  'coverage':'Union des coordonnées QUERY couvertes / longueur QUERY',
                  'missing':'Sans HSP détecté, pas absence démontrée d’homologie',
                  'bed':'0-based half-open','tsv_gff3':'1-based inclusif',
                  'blastx':'Les coordonnées QUERY sont nucléotidiques et les coordonnées de la référence protéiques',
                  'msa_distance':'Distance p avec suppression des colonnes à gap ou ? par paire'})
    _json(folder/'provenance.json',meta)
    (folder/'README_exports.txt').write_text(
        'ALIGNVIEW 6 — EXPORTS\n\n'
        'Les fichiers BED utilisent des coordonnées 0-based / end-exclusive.\n'
        'coordinate_units.tsv décrit les unités nt/aa en fonction du mode BLAST.\n'
        'coordinate_audit.tsv compare les étendues nucléotidiques avec les résidus alignés.\n'
        'Les tableaux de positions et GFF3 utilisent des positions 1-based inclusives.\n'
        'Les identités proviennent des HSP locaux, et non de séquences globalement alignées.\n'
        'Les régions sans HSP ne prouvent pas une absence d’homologie.\n'
        'Les distances p du MSA ignorent les positions avec gaps pour la paire concernée.\n'
        'Les fichiers NEXUS/PHYLIP nécessitent un MSA et leur présence est conditionnelle.\n'
        'Les profils publication/complete incluent le rapport HTML autonome et le dossier figures/ si calculés.\n',encoding='utf-8')
    files=sorted((x for x in folder.rglob('*') if x.is_file()),key=lambda p:str(p.relative_to(folder)))
    with (folder/'MANIFEST.tsv').open('w',encoding='utf-8') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['file','bytes','sha256'])
        for file in files:writer.writerow([str(file.relative_to(folder)),file.stat().st_size,sha256(file)])
    return folder
