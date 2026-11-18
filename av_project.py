"""Portable projects with original inputs, results and atomic directory publication."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import tempfile

from av_core import HSP, read_fasta, sha256
from av_version import VERSION

SCHEMA = 1
INTEGRITY_FILE = 'project_manifest.json'
IGNORE = shutil.ignore_patterns('database*', 'mmseqs_tmp', 'raw_hsps.tsv', 'history','pending_msa','pending_tree','pending_vsearch_*','.alignview-output-*')


def local_path(root, value):
    if not isinstance(value,str) or not value:raise ValueError('Invalid project file path.')
    path = (root/value).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError(f'Project file missing or outside the project: {value}')
    return path


def verify_integrity(root, data):
    """Check complete snapshots when present; schema-1 legacy projects still load.

    This detects accidental changes, not hostile edits by someone who can also
    replace the manifest. Original source paths are never needed for validation.
    """
    integrity=data.get('integrity')
    if integrity is None:return
    if not isinstance(integrity,dict) or integrity.get('schema')!=1 or integrity.get('manifest')!=INTEGRITY_FILE:
        raise ValueError('Invalid project integrity information.')
    try:
        manifest=json.loads(local_path(root,INTEGRITY_FILE).read_text(encoding='utf-8'))
        files=manifest['files']
        if manifest.get('schema')!=1 or not isinstance(files,dict) or 'session.json' not in files:
            raise ValueError('Invalid project integrity manifest.')
        for name,record in files.items():
            if not isinstance(record,dict) or not isinstance(record.get('sha256'),str):
                raise ValueError('Invalid project integrity entry: '+str(name))
            path=local_path(root,name)
            if path.stat().st_size!=record.get('bytes') or sha256(path)!=record['sha256']:
                raise ValueError('A saved project file has changed: '+name)
    except (OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError('Project integrity manifest is missing or invalid.') from exc


def write_integrity(root):
    """Include scientific results, metadata, logs and archived analyses."""
    files={}
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path==root/INTEGRITY_FILE:continue
        name=path.relative_to(root).as_posix()
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError('Project file points outside the project: '+name)
        files[name]={'bytes':path.stat().st_size,'sha256':sha256(path)}
    (root/INTEGRITY_FILE).write_text(json.dumps({'schema':1,'files':files},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')


def validate_reference_selection(selection, reference_ids):
    if not isinstance(selection,list) or any(not isinstance(identifier,str) or identifier not in reference_ids for identifier in selection) or len(set(selection))!=len(selection):
        raise ValueError('The project has an invalid reference selection.')
    return list(selection)


def read_project(folder):
    root = Path(folder).expanduser().resolve()
    if root.is_file():root=root.parent
    try:
        data=json.loads((root/'session.json').read_text(encoding='utf-8'))
        if not isinstance(data,dict) or data.get('schema')!=SCHEMA:
            raise ValueError('Unsupported project format version.')
        if not isinstance(data.get('inputs'),dict) or not data['inputs'].get('query') or not isinstance(data.get('args'),dict):
            raise ValueError('Invalid AlignView project: missing inputs or parameters.')
        if not isinstance(data.get('active_query_id'),str) or not isinstance(data.get('mode'),str):
            raise ValueError('Invalid AlignView project: missing query or mode.')
        # Input checks below retain their specific errors for older projects.
        if not isinstance(data.get('hsps',[]),list):raise ValueError('Invalid HSP list.')
        try:
            for row in data.get('hsps',[]):HSP(**row)
        except (TypeError,ValueError) as exc:raise ValueError('Invalid HSPs in the project.') from exc
        for name in ('query','references'):
            if data['inputs'].get(name):
                item=data['inputs'][name]
                path=local_path(root,item['path'])
                if sha256(path)!=item['sha256']:
                    raise ValueError(f'The input file has changed: {path.name}')
        if data.get('msa'):local_path(root,data['msa'])
        verify_integrity(root,data)
        return root,data
    except (KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError('Invalid AlignView project: incomplete session.json.') from exc


def save_project(session, destination, include_history=True):
    target=Path(destination).expanduser().resolve()
    if target.exists():raise ValueError(f'The directory already exists: {target}')
    selection=validate_reference_selection(getattr(session,'reference_selection',[]),{record.id for record in session.references})
    target.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.alignview-save-',dir=target.parent) as temporary:
        staged=Path(temporary)/'project'
        shutil.copytree(session.work,staged,ignore=IGNORE)
        inputs=staged/'inputs';inputs.mkdir(exist_ok=True)
        snapshots={}
        for name,source in (('query',session.query_path),('references',session.refs_path)):
            if source is None:continue
            filename=name+('.fasta.gz' if source.suffix.lower()=='.gz' else '.fasta')
            copied=inputs/filename
            shutil.copyfile(source,copied)
            original=session.source_query_path if name=='query' else session.source_refs_path
            snapshots[name]={'path':copied.relative_to(staged).as_posix(),'sha256':sha256(copied),'original_path':str(original)}
        msa=''
        if session.msa_file:
            copied=staged/'alignment_active.fasta'
            shutil.copyfile(session.msa_file,copied)
            msa=copied.name
        state=dict(schema=SCHEMA,version=VERSION,inputs=snapshots,
                   integrity={'schema':1,'manifest':INTEGRITY_FILE},
                   args={k:v for k,v in vars(session.a).items() if not k.startswith('_')},
                   active_query_id=session.query.id,mode=session.mode,search_completed=session.search_completed,
                   search_parameters=session.search_parameters,
                   calculation_runs=session.calculation_runs,
                   reference_selection=selection,workflow_goal=getattr(session,'workflow_goal',''),
                   sequence_selection=getattr(session,'sequence_selection',[]),
                   comparison_selection=getattr(session,'comparison_selection',[]),
                   comparison_batch=getattr(session,'comparison_batch',{}),
                   workspace_entry=getattr(session,'workspace_entry',{}),
                   navigation_state=getattr(session,'navigation_state',{}),
                   result_sources=getattr(session,'result_sources',{}),
                   vsearch_results=getattr(session,'vsearch_results',{}),
                   hsps=[asdict(h) for h in session.hsps],msa=msa,
                   dot=session.dot_data,motifs=session.motif_rows,orfs=session.orf_rows,
                   commands=session.tools.history,started_utc=session.start)
        (staged/'session.json').write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        if include_history:
            for i,archive in enumerate(getattr(session,'history',[]),1):
                shutil.copytree(archive,staged/'history'/f'{i:03d}',ignore=IGNORE)
        write_integrity(staged)
        # Directory rename is atomic on the destination filesystem.
        staged.rename(target)
    return target


def resume_project(folder, cli_args, ui):
    from alignview import Session
    root,data=read_project(folder)
    args=argparse.Namespace(**vars(cli_args))
    runtime={'query','references','resume','save','out','no_menu','simple_ui','no_color','no_open',
             'executable','tools_dir','module_policy','diagnose','diagnose_json','setup','install_tools','tools_prefix','export_dir','html_out','docs','lang','load_profile','write_profile'}
    overrides={name:getattr(cli_args,name) for name in getattr(cli_args,'_explicit',set()) if hasattr(cli_args,name)}
    for name,value in data.get('args',{}).items():
        if name not in runtime and hasattr(args,name):setattr(args,name,value)
    args.query=str(local_path(root,data['inputs']['query']['path']))
    args.references=str(local_path(root,data['inputs']['references']['path'])) if data['inputs'].get('references') else None
    args.query_id=data['active_query_id'];args.mode=data['mode']
    if args.preset=='sensible':args.preset='sensitive'  # Migrate older saved presets.
    original_filters=tuple(getattr(args,n) for n in ('hits','identity','coverage','bitscore'))
    if cli_args.load_profile:
        from av_parameters import load_profile
        load_profile(args,cli_args.load_profile)
        vars(args).update(getattr(cli_args,'_profile_overrides',{}))
    vars(args).update(overrides)
    from av_parameters import validate
    validate(args)
    session=Session(args,ui)
    try:
        shutil.copytree(root,session.work,dirs_exist_ok=True,ignore=IGNORE)
        # Detach from the original project so it can be moved or deleted afterwards.
        session.query_path=local_path(session.work,data['inputs']['query']['path'])
        session.refs_path=local_path(session.work,data['inputs']['references']['path']) if args.references else None
        session.source_query_path=Path(data['inputs']['query'].get('original_path',session.query_path))
        session.source_refs_path=Path(data['inputs']['references'].get('original_path',session.refs_path)) if args.references else None
        session.a.query=str(session.query_path);session.a.references=str(session.refs_path) if session.refs_path else None
        session.query_search_path=session.work/'query_active.fasta'
        from av_core import write_fasta
        write_fasta(session.query_list if session.single else [session.query],session.query_search_path)
        session.hsps=[HSP(**row) for row in data.get('hsps',[])]
        session.search_completed=data.get('search_completed',bool(session.hsps) or any('search' in row.get('stderr_log','') and row.get('exit_code')==0 for row in data.get('commands',[])))
        session.search_parameters=data.get('search_parameters',{})
        session.calculation_runs=data.get('calculation_runs',{})
        session.vsearch_results=data.get('vsearch_results',{})
        session.reference_selection=validate_reference_selection(data.get('reference_selection',[]),session.refmap)
        session.sequence_selection=validate_reference_selection(data.get('sequence_selection',[]),{r.id for r in session.query_list})
        session.comparison_selection=validate_reference_selection(data.get('comparison_selection',[]),{r.id for r in session.query_list})
        session.comparison_batch=data.get('comparison_batch') or {}
        session.workspace_entry=data.get('workspace_entry') or {}
        session.navigation_state=data.get('navigation_state') or {}
        session.result_sources=data.get('result_sources') or {}
        session.dot_data=data.get('dot');session.motif_rows=data.get('motifs',[]);session.orf_rows=data.get('orfs',[])
        session.tools.history=data.get('commands',[]);session.start=data.get('started_utc',session.start)
        if data.get('msa'):
            session.msa_file=local_path(session.work,data['msa']);session.msa_records=read_fasta(session.msa_file)
            if len({len(r.seq) for r in session.msa_records})!=1:raise ValueError('Project MSA is not aligned: sequences have different lengths.')
        # Result snapshots must be detached too: the user may move or remove
        # the original project while this reopened workspace remains active.
        session.history=[]
        history=root/'history'
        if history.is_dir():
            owned_history=Path(session.temporary.name)/'history_sources'
            for archive in sorted(history.iterdir()):
                if not archive.is_dir():continue
                if not archive.resolve().is_relative_to(root):
                    raise ValueError('Archived result points outside the project: '+archive.name)
                copied=owned_history/archive.name
                shutil.copytree(archive,copied,ignore=IGNORE)
                session.history.append(copied)
        session.apply_filters(preserve_msa=original_filters==tuple(getattr(args,n) for n in ('hits','identity','coverage','bitscore')))
        saved_goal=data.get('workflow_goal','')
        if saved_goal in ('compare','align','view','sequence'):session.workflow_goal=saved_goal
        elif not session.workflow_goal and session.msa_records:session.workflow_goal='view'
        session.last_saved=root;session.dirty=bool(cli_args.load_profile or (set(overrides)-runtime))
        ui.success(f'Project resumed: {session.query.id} · {len(session.hsps)} HSP · no calculations rerun')
        return session
    except BaseException:
        session.close()
        raise
