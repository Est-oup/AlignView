"""Temporary workspace catalog: immutable analysis snapshots, optional persistence.

The application owns the catalog directory, independently of any loaded Session.
Opening an entry copies it into a working Session; calculations never write into
another entry. Explicit project saving includes the catalog through history/.
"""
from __future__ import annotations

import json
from pathlib import Path
import uuid


RESULT_VIEWS = {
    'search': ('Matches', 'Matches', ''),
    'search_batch': ('Reference comparisons', 'Batch', 'comparison_batch.json'),
    'msa': ('Alignment', 'Alignment', ''),
    'tree': ('Tree', 'Tree', ''),
    'translation': ('Six-frame translation', 'Artifact', 'six_frames.faa'),
    'orfs': ('ORFs', 'Artifact', 'orfs.tsv'),
    'motif': ('Motif positions', 'Artifact', 'motifs.tsv'),
    'extract': ('Extracted region', 'Artifact', ''),
    'dotplot': ('Dotplot', 'Artifact', 'dotplot.json'),
    'conservation': ('Conservation', 'Artifact', 'msa_positions.tsv'),
    'distances': ('Distances', 'Artifact', 'msa_pdistance.tsv'),
    'statistics': ('Sequence quality', 'Artifact', 'fasta_statistics.tsv'),
}


def entry_label(session, entry):
    """Keep user labels separate from immutable scientific snapshots."""
    return session.navigation_state.get(entry.get('id'),{}).get('label',entry.get('label','Result available'))


def rename_result(session, identifier, label):
    """Rename a completed result in this workspace without rewriting its files."""
    if not isinstance(label,str) or not label.strip() or not label.isprintable() or len(label.strip())>120:
        raise ValueError('Use a result name of 1–120 printable characters.')
    entry=next((row for row in entries(session,'result') if row['id']==identifier),None)
    if entry is None:raise ValueError('Choose a completed result to rename.')
    label=label.strip()
    if label!=entry_label(session,entry):
        session.navigation_state.setdefault(identifier,{})['label']=label
        session.dirty=True
    return label


def entries(session, kind=None):
    """Read display metadata only; opening validates the complete project."""
    rows=[]
    for path in session.history:
        path=Path(path)
        data=json.loads((path/'session.json').read_text(encoding='utf-8'))
        entry=dict(data.get('workspace_entry') or {})
        if not entry:
            entry=dict(id=path.name,kind='result',label='Previous analysis · '+data['active_query_id'],
                       view='Alignment' if data.get('msa') else 'Matches',artifact='')
        entry['label']=entry_label(session,entry)
        entry['path']=path
        if kind is None or entry['kind']==kind:rows.append(entry)
    return rows


def capture(session, kind, label=None, *, view=None, artifact='', reuse=''):
    """Retain a completed state without asking the user for a directory.

    A catalog root is attached by the interactive application. Batch callers do
    not pay for temporary snapshots and keep their explicit saving semantics.
    """
    root=getattr(session,'library_root',None)
    if root is None:return None
    from av_project import save_project
    if reuse:
        previous=next((row for row in entries(session,'data') if row.get('source_hash')==reuse),None)
        if previous:
            session.workspace_entry={key:value for key,value in previous.items() if key!='path'}
            return previous
    title,default_view,default_artifact=RESULT_VIEWS.get(kind,(kind.replace('_',' ').title(),'Artifact',''))
    previous=getattr(session,'workspace_entry',{})
    catalog=entries(session)
    inputs=[row['id'] for row in catalog if row['kind']=='data' and row.get('source_hash') in
            (session.input_hashes['query'],session.input_hashes.get('refs',''))]
    if kind in ('tree','conservation','distances') and session.result_sources.get('msa'):
        inputs=[session.result_sources['msa']]
    elif kind=='msa' and session.a.msa in ('hits','hsp','selection') and not session.single and session.result_sources.get('search'):
        inputs=[session.result_sources['search']]
    if kind=='data':inputs=[session.data_origin] if session.data_origin else []
    number=1+sum(row.get('analysis')==kind for row in catalog if row['kind']=='result')
    identifier=uuid.uuid4().hex
    entry=dict(id=identifier,kind='data' if kind=='data' else 'result',analysis=kind,
               label=(label or title+' · '+session.query.id)+(f' · #{number}' if kind!='data' else ''),
               view=view or default_view,artifact=artifact or default_artifact,
               source=session.source_query_path.name,source_hash=reuse,
               dataset_hash=session.input_hashes['query'],
               parent=inputs[0] if inputs else None,inputs=inputs,origin=previous.get('id'),query_id=session.query.id)
    destination=Path(root)/identifier
    session.workspace_entry=entry
    if kind=='search_batch':
        entry['members']=[run['result_id'] for run in session.comparison_batch['runs'] if run.get('result_id')]
    sources=dict(session.result_sources)
    if kind!='data':session.result_sources[kind]=identifier
    try:save_project(session,destination,include_history=False)
    except BaseException:
        session.workspace_entry=previous
        session.result_sources=sources
        raise
    session.history.append(destination)
    session.dirty=True
    return {**entry,'path':destination}


def capture_calculation(session, kind, details):
    if kind=='search':return  # Search also has to publish its filtered evidence.
    artifact=details.get('artifact','')
    if kind.startswith('vsearch_'):
        action=kind.removeprefix('vsearch_')
        return capture(session,kind,'VSEARCH '+action+' · '+session.source_query_path.name,
                       view='Vsearch',artifact='vsearch/'+action+'/result.json')
    label=None
    if kind=='msa':
        prep=details.get('preparation',{})
        engine=details.get('parameters',{}).get('aligner',session.a.aligner)
        label=f'Alignment · {len(session.msa_records)} sequences · {engine.upper()}'
        if prep.get('query_frame'):label+=f' · frame {prep["query_frame"]}'
    elif kind=='tree':
        alignment=next((row for row in entries(session) if row['id']==session.result_sources.get('msa')),None)
        label='Tree · '+(alignment['label'] if alignment else 'alignment')
    elif kind=='search_batch':
        label=f'Reference comparisons · {len(session.comparison_batch["query_ids"])} queries → {session.source_refs_path.name}'
    return capture(session,kind,label,artifact=artifact)
