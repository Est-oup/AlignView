"""Explicit nucleotide workflows; VSEARCH scores are never BLAST E-values."""
from __future__ import annotations

import csv
import html
import json
import math
from pathlib import Path
import re
import shutil

from av_core import fasta_iter

ACTIONS = {
    'search': ('Global search against references', 'Compare the active query against nucleotide references'),
    'dereplicate': ('Dereplicate sequences', 'Merge identical sequences and retain their abundances'),
    'cluster': ('Cluster sequences', 'Create identity-based clusters and export their representatives'),
    'chimeras': ('Detect chimeras de novo', 'Dereplication followed by UCHIME; abundances matter'),
}
FIELDS = ('query','target','id','alnlen','qilo','qihi','tilo','tihi','qcov','tcov','qstrand','qrow','trow')


def publish_directory(staged,target):
    """Retain the previous valid directory if publishing the new one fails."""
    backup=target.with_name('.'+target.name+'_previous')
    if backup.exists():
        if not target.exists():backup.rename(target)
        else:shutil.rmtree(backup)
    if target.exists():target.rename(backup)
    try:staged.rename(target)
    except BaseException:
        if backup.exists():backup.rename(target)
        raise
    if backup.exists():shutil.rmtree(backup)


def commands(action, source, references, output, options):
    from av_tools import Command
    from av_parameters import validate
    validate(options,context='vsearch')
    output=Path(output)
    if action not in ACTIONS:raise ValueError('Analyse VSEARCH inconnue: '+action)
    limits=['--minseqlength',str(options.vsearch_min_length),'--maxseqlength',str(options.vsearch_max_length)]
    sizein=['--sizein'] if options.vsearch_sizein else []
    common=['--threads',str(options.threads),'--strand',options.vsearch_strand]
    if action=='search':
        if references is None:raise ValueError('VSEARCH search requires a reference FASTA.')
        args=['vsearch','--usearch_global',str(source),'--db',str(references),
            '--id',str(options.vsearch_id),'--iddef',str(options.vsearch_iddef),
            '--query_cov',str(options.vsearch_query_cov),'--target_cov',str(options.vsearch_target_cov),
            '--maxaccepts',str(options.vsearch_maxaccepts),'--maxrejects',str(options.vsearch_maxrejects),
            '--maxhits',str(options.hits),'--minwordmatches',str(options.vsearch_minwordmatches),
            '--qmask',options.vsearch_masking,'--dbmask',options.vsearch_masking,
            '--userout',str(output/'matches_raw.tsv'),'--userfields','+'.join(FIELDS),*common,*limits]
        if options.vsearch_n_mismatch:args+=['--n_mismatch']
        return [Command(args,'VSEARCH: global search','vsearch_search.log',allow_empty=True)]
    if action=='cluster':
        args=['vsearch','--cluster_fast',str(source),'--id',str(options.vsearch_id),
            '--iddef',str(options.vsearch_iddef),'--centroids',str(output/'centroids.fasta'),
            '--uc',str(output/'clusters.uc'),'--sizeout','--qmask',options.vsearch_masking,
            *sizein,*common,*limits]
        if options.vsearch_n_mismatch:args+=['--n_mismatch']
        return [Command(args,'VSEARCH: clustering','vsearch_cluster.log')]
    unique=output/'unique.fasta'
    derep=['vsearch','--derep_fulllength',str(source),'--output',str(unique),
           '--uc',str(output/'dereplication.uc'),'--sizeout',*sizein,'--strand',options.vsearch_strand,*limits]
    result=[Command(derep,'VSEARCH: dereplication','vsearch_dereplicate.log')]
    if action=='chimeras':
        result.append(Command(['vsearch','--uchime_denovo',str(unique),
            '--chimeras',str(output/'chimeras.fasta'),'--nonchimeras',str(output/'nonchimeras.fasta'),
            '--borderline',str(output/'borderline.fasta'),'--uchimeout',str(output/'uchime.tsv'),
            '--abskew',str(options.vsearch_abskew),'--minh',str(options.vsearch_minh),
            '--qmask',options.vsearch_masking,*limits],
            'VSEARCH: de novo chimera detection','vsearch_chimeras.log'))
    return result


def count_fasta(path):
    return sum(1 for _ in fasta_iter(path))


def summarize(action, output, input_count, eligible_count):
    output=Path(output)
    summary={'input_sequences':input_count,'eligible_sequences':eligible_count,
             'excluded_by_length':input_count-eligible_count}
    preview=[]
    if action=='search':
        raw=output/'matches_raw.tsv'
        if not raw.is_file():raise ValueError('VSEARCH: missing search table.')
        matches=0;targets=set()
        with raw.open(encoding='utf-8') as stream, (output/'matches.tsv').open('w',encoding='utf-8',newline='') as target:
            writer=csv.writer(target,delimiter='\t');writer.writerow(FIELDS)
            for values in csv.reader(stream,delimiter='\t'):
                if not values:continue
                if len(values)!=len(FIELDS):raise ValueError('VSEARCH: unexpected column count.')
                row=dict(zip(FIELDS,values))
                for key in ('id','qcov','tcov'):
                    value=float(row[key])
                    if not math.isfinite(value) or not 0<=value<=100:raise ValueError('VSEARCH: invalid percentage.')
                    row[key]=value
                if row['qstrand'] not in ('+','-'):raise ValueError('VSEARCH: invalid orientation.')
                matches+=1;targets.add(row['target']);writer.writerow(values)
                if len(preview)<200:
                    preview.append({k:row[k] for k in ('query','target','id','qcov','tcov','qstrand')})
        raw.unlink()
        summary.update(matches=matches,matched_references=len(targets))
    elif action in ('dereplicate','cluster'):
        filename='unique.fasta' if action=='dereplicate' else 'centroids.fasta'
        path=output/filename
        if not path.is_file():raise ValueError('VSEARCH: missing representatives.')
        representatives=0;total=0
        for record in fasta_iter(path):
            size=re.search(r';size=(\d+)(?:;|$)',record.id)
            abundance=int(size[1]) if size else 1
            representatives+=1;total+=abundance
            if len(preview)<200:preview.append({'representative':record.id,'length_nt':len(record.seq),'abundance':abundance})
        if not representatives:raise ValueError('VSEARCH: no representative sequences.')
        summary.update(representatives=representatives,total_abundance=total)
    else:
        for name in ('unique','chimeras','nonchimeras','borderline'):
            path=output/(name+'.fasta')
            if not path.is_file():raise ValueError('VSEARCH: sortie absente: '+path.name)
            summary[name]=count_fasta(path)
        if summary['unique']!=summary['chimeras']+summary['nonchimeras']+summary['borderline']:
            raise ValueError('VSEARCH: incomplete chimera classification.')
        for name in ('chimeras','borderline','nonchimeras'):
            for record in fasta_iter(output/(name+'.fasta')):
                if len(preview)<200:preview.append({'sequence':record.id,'classification':name})
    return summary,preview


def report_section(work):
    index=Path(work)/'vsearch/analyses.json'
    if not index.is_file():return ''
    analyses=json.loads(index.read_text(encoding='utf-8'))
    esc=lambda value:html.escape(str(value),quote=True)
    parts=['<section id="vsearch" class="section"><div class="panel"><h2>VSEARCH analyses</h2>',
        '<p>Nucleotide analyses. VSEARCH search uses global alignment; BLAST E-values and bitscores are not calculated.</p>',
        '<p class="hint">Previews are limited to 200 rows per analysis; full files are in the vsearch subdirectory of the project and exports.</p></div>']
    for action,result in analyses.items():
        parts.append('<article class="panel"><h2>'+esc(ACTIONS[action][0])+'</h2>')
        parts.append('<dl>'+''.join('<dt>'+esc(key)+'</dt><dd>'+esc(value)+'</dd>' for key,value in result['summary'].items())+'</dl>')
        if action=='search':parts.append('<p class="hint">Identity follows iddef; VSEARCH qcov/tcov are percentages and differ from a union of local HSPs.</p>')
        if action=='chimeras':parts.append('<p class="notice">Abundances and dataset composition affect UCHIME; no detected chimera is not an absolute validation.</p>')
        rows=result.get('preview',[])
        if rows:
            columns=list(rows[0]);parts.append('<div class="tablewrap"><table><thead><tr>'+''.join('<th>'+esc(c)+'</th>' for c in columns)+'</tr></thead><tbody>')
            parts.extend('<tr>'+''.join('<td>'+esc(row.get(c,''))+'</td>' for c in columns)+'</tr>' for row in rows)
            parts.append('</tbody></table></div>')
        else:parts.append('<p>No results retained for this preview.</p>')
        parts.append('<details><summary>Tool and parameters used</summary><pre>'+esc(json.dumps({'settings':result['settings'],'commands':result['commands']},ensure_ascii=False,indent=2))+'</pre></details></article>')
    parts.append('</section>')
    return ''.join(parts)
