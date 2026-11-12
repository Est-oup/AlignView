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
IGNORE = shutil.ignore_patterns('database*', 'mmseqs_tmp', 'raw_hsps.tsv', 'history','pending_msa','pending_tree','pending_vsearch_*')


def local_path(root, value):
    path = (root/value).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError(f'Fichier de projet absent ou hors du projet : {value}')
    return path


def read_project(folder):
    root = Path(folder).expanduser().resolve()
    if root.is_file():root=root.parent
    try:
        data=json.loads((root/'session.json').read_text(encoding='utf-8'))
        if not isinstance(data,dict) or data.get('schema')!=SCHEMA:
            raise ValueError('Version du format de projet non prise en charge.')
        if not isinstance(data.get('inputs'),dict) or not data['inputs'].get('query') or not isinstance(data.get('args'),dict):
            raise ValueError('Projet AlignView invalide : entrées ou paramètres manquants.')
        if not isinstance(data.get('active_query_id'),str) or not isinstance(data.get('mode'),str):
            raise ValueError('Projet AlignView invalide : requête ou mode manquant.')
        if not isinstance(data.get('hsps',[]),list):raise ValueError('Liste de HSP invalide.')
        try:
            for row in data.get('hsps',[]):HSP(**row)
        except (TypeError,ValueError) as exc:raise ValueError('HSP invalides dans le projet.') from exc
        for name in ('query','references'):
            if data['inputs'].get(name):
                item=data['inputs'][name]
                path=local_path(root,item['path'])
                if sha256(path)!=item['sha256']:
                    raise ValueError(f'Le fichier d’entrée a changé : {path.name}')
        if data.get('msa'):local_path(root,data['msa'])
        return root,data
    except (KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError('Projet AlignView invalide : session.json est incomplet.') from exc


def save_project(session, destination, include_history=True):
    target=Path(destination).expanduser().resolve()
    if target.exists():raise ValueError(f'Le dossier existe déjà : {target}')
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
                   args={k:v for k,v in vars(session.a).items() if not k.startswith('_')},
                   active_query_id=session.query.id,mode=session.mode,search_completed=session.search_completed,
                   search_parameters=session.search_parameters,
                   calculation_runs=session.calculation_runs,
                   vsearch_results=getattr(session,'vsearch_results',{}),
                   hsps=[asdict(h) for h in session.hsps],msa=msa,
                   dot=session.dot_data,motifs=session.motif_rows,orfs=session.orf_rows,
                   commands=session.tools.history,started_utc=session.start)
        (staged/'session.json').write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        if include_history:
            for i,archive in enumerate(getattr(session,'history',[]),1):
                shutil.copytree(archive,staged/'history'/f'{i:03d}',ignore=IGNORE)
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
        session.dot_data=data.get('dot');session.motif_rows=data.get('motifs',[]);session.orf_rows=data.get('orfs',[])
        session.tools.history=data.get('commands',[]);session.start=data.get('started_utc',session.start)
        if data.get('msa'):
            session.msa_file=local_path(session.work,data['msa']);session.msa_records=read_fasta(session.msa_file)
            if len({len(r.seq) for r in session.msa_records})!=1:raise ValueError('MSA du projet non aligné : longueurs différentes.')
        session.history=sorted((root/'history').glob('*')) if (root/'history').is_dir() else []
        session.apply_filters(preserve_msa=original_filters==tuple(getattr(args,n) for n in ('hits','identity','coverage','bitscore')))
        session.last_saved=root;session.dirty=bool(cli_args.load_profile or (set(overrides)-runtime))
        ui.success(f'Projet repris : {session.query.id} · {len(session.hsps)} HSP · aucun calcul relancé')
        return session
    except BaseException:
        session.close()
        raise
