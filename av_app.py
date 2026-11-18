"""Application lifecycle: startup, history transitions, batch actions and recovery."""
from __future__ import annotations
import argparse
import json
import os
import shlex
import shutil
import sys
import tempfile
import uuid
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from av_core import read_fasta
from av_tools import Toolchain, CommandError
from av_ui import UI, HomeRequested
from av_docs import TOPICS, ALIASES as DOC_ALIASES, read_document

from av_start import startup, StartRequest, choose_file
from av_parameters import write_profile
from av_project import save_project, resume_project
from av_setup import setup_tools, diagnostic

from av_cli import parse_args
from av_session import Session

def show_tool_list():
    groups={
        'Homology and coverage':['MMseqs2 BLASTP for protein references','BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX','DIAMOND BLASTP / BLASTX',
                                  'Multiple HSPs, union coverage, gaps, substitutions and frames',
                                  'QUERY and reference maps, identity × coverage'],
        'Alignments':['MUSCLE 5, Clustal Omega and FastTree','Needle global, Water local','MAFFT hits / all / fragments / HSP',
                      'ALV and ANSI MSA viewer','Consensus, conservation, variable sites and p-distances'],
        'Sequences':['FASTA checks, GC, N and ambiguous symbols','Both-strand IUPAC motifs and region extraction',
                    'Six-frame translation, start-to-stop ORFs','Entropy, sliding profiles and k-mers'],
        'Visualization and reproducibility':['Conditional SVG atlas and FastTree tree with scale and support','Standalone HTML report, guided navigation and SVG exports','SVG coverage, dotplot, profiles and MSA',
                                           'CSV, TSV, JSON/JSONL, BED, GFF3, NEXUS, PHYLIP, Stockholm and SHA256','Slurm scripts, automatic modules and logs',
                                           'Offline manual: --docs; standalone HTML: --html-out report.html']}
    for group,items in groups.items():
        print('\n'+str(group))
        for item in items:print('  • '+str(item))


def _apply_start(args,request,demo_root):
    args._guided=True;args._intent=request.intent
    args._data_origin=request.origin
    args._comparison_queries=[]
    args._workflow_goal=request.intent if request.intent in ('compare','align','view','sequence') else ''
    args.query_id='';args.mode=None;args.msa='none';args.resume=None
    if request.intent=='resume':
        args.query=None;args.references=None;args.resume=str(request.query)
    elif request.intent=='demo':
        from av_demo import write_demo
        args.query=str(write_demo(demo_root));args.references=None;args.mode='msa';args.msa='view'
    else:
        args.query=str(request.query)
        args.references=str(request.references) if request.references else None
        if request.intent=='view':args.mode='msa';args.msa='view'
        elif request.intent=='align':args.mode='msa'


def _select_start(args,ui,demo_root):
    if args.demo:_apply_start(args,StartRequest('demo'),demo_root)
    if not args.query and not args.resume:
        if not ui.tty:raise ValueError('Use an interactive terminal or supply a FASTA, --demo or --resume.')
        command=getattr(args,'_command','')
        if command in ('view','align','compare'):
            first=choose_file(ui,'Alignment' if command=='view' else 'Query sequences' if command=='compare' else 'Sequences to align')
            if first is None:return False
            refs=choose_file(ui,'Reference sequences',first.parent) if command=='compare' else None
            if command=='compare' and refs is None:return False
            request=StartRequest(command,first,refs)
        else:request=startup(ui)
        if request is None:return False
        mode=args.mode;query_id=args.query_id
        _apply_start(args,request,demo_root)
        if 'mode' in args._explicit:args.mode=mode
        if 'query_id' in args._explicit:args.query_id=query_id
    if getattr(args,'_command','')=='compare' and not args.resume and not args.references:
        if not ui.tty or args.no_menu:raise ValueError('compare requires a reference FASTA.')
        reference=choose_file(ui,'Reference sequences',Path(args.query).expanduser().parent)
        if reference is None:return False
        args.references=str(reference)
    ui.workflow=()
    args._defer_query_choice=ui.tty and not args.no_menu and (getattr(args,'_command','')=='compare' or getattr(args,'_intent','')=='compare')
    return True


def _register_inputs(session):
    """Make imported files accessible as standalone datasets in this workspace."""
    from av_library import capture
    legacy=not session.workspace_entry
    if session.single:
        if not session.workspace_entry:
            capture(session,'data',session.source_query_path.name,view='Sequences',reuse=session.input_hashes['query'])
    for source,original in (() if session.single else ((session.query_path,session.source_query_path),(session.refs_path,session.source_refs_path))):
        fresh=argparse.Namespace(**vars(session.a))
        fresh.query=str(source);fresh.references=None;fresh.query_id='';fresh.mode='msa';fresh.msa='none'
        fresh._comparison_queries=[]
        data=Session(fresh,session.ui)
        try:
            data.data_origin=''
            data.source_query_path=original;data.library_root=session.library_root;data.history=session.history
            capture(data,'data',original.name,view='Sequences',reuse=data.input_hashes['query'])
        finally:data.close()
    if legacy:
        if session.search_completed:capture(session,'search','Imported matches · '+session.query.id)
        if session.msa_records:capture(session,'msa','Imported alignment · '+session.source_query_path.name)
        if (session.work/'tree_fasttree.nwk').is_file():capture(session,'tree','Imported tree · '+session.query.id)
        for action in session.vsearch_results:
            artifact='vsearch/'+action+'/result.json'
            if (session.work/artifact).is_file():
                capture(session,'vsearch_'+action,'Imported VSEARCH '+action,view='Vsearch',artifact=artifact)


def main(argv=None):
    if sys.version_info<(3,10):
        print('AlignView requires Python 3.10 or newer.',file=sys.stderr);return 2
    args=parse_args(argv)
    ui=UI(not args.no_color,args.simple_ui or args.no_menu)
    if args.write_profile:
        try:
            path=write_profile(args,args.write_profile);ui.success('Calculation profile saved: '+str(path));return 0
        except (OSError,ValueError) as exc:ui.error(str(exc));return 1
    if args.tools_prefix:os.environ['ALIGNVIEW_TOOLS_PREFIX']=str(Path(args.tools_prefix).expanduser().resolve())
    if args.docs:
        try:
            content=read_document(args.docs)
            if ui.tty and not args.no_menu:ui.page(content,TOPICS[DOC_ALIASES.get(args.docs,args.docs)][1])
            else:print(content,end='' if content.endswith('\n') else '\n')
            return 0
        except (ValueError,OSError) as exc:ui.error(str(exc));return 1
    if args.list_tools:show_tool_list();return 0
    if args.list_genetic_codes:
        from av_genetics import TABLES, context_dependent
        for identifier,(name,_,_) in TABLES.items():
            print(f'{identifier:>2}  {name}'+(' · context-dependent stops' if context_dependent(identifier) else ''))
        return 0
    if args.setup:
        try:setup_tools(ui,args.setup,args.tools_prefix,args.install_tools);return 0
        except (ValueError,CommandError,OSError) as exc:ui.error(str(exc));return 1
    if args.diagnose or args.diagnose_json:
        try:
            with tempfile.TemporaryDirectory(prefix='alignview_diagnostic_') as tmp:
                tools=Toolchain(Path(tmp),None,auto_modules=args.module_policy=='auto',executable_paths=args.executable,tool_dirs=args.tools_dir)
                rows=diagnostic(tools)
                if args.diagnose_json:print(json.dumps(rows,ensure_ascii=False,indent=2))
                else:
                    ui.banner('DIAGNOSTICS / SCIENTIFIC TOOLS')
                    for row in rows:
                        ui.info(f'{row["command"]}: '+(row['path']+' · '+row['version'] if row['available'] else f'missing → --setup {row["setup_profile"]} --install-tools'))
            return 0
        except (ValueError,OSError) as exc:ui.error(str(exc));return 1
    session=None;history=[];failed=False;trail=[];navigation={}
    with tempfile.TemporaryDirectory(prefix='alignview_workspace_') as workspace_name:
        workspace=Path(workspace_name)
        try:
            if not _select_start(args,ui,workspace/'demo'):return 0
            session=resume_project(args.resume,args,ui) if args.resume else Session(args,ui)
            for i,path in enumerate(session.history,1):
                archived=workspace/f'imported_{i:03d}'
                shutil.copytree(path,archived);history.append(archived)
            session.history=history
            interactive=ui.tty and not args.no_menu
            library=workspace/'catalog'
            if interactive:
                session.library_root=library
                navigation.update(session.navigation_state);session.navigation_state=navigation
                _register_inputs(session)
            if args.resume:session.resume_actions()
            if not args.resume:
                try:session.initial()
                except HomeRequested:pass
                except CommandError as exc:
                    if args.no_menu or not ui.tty:raise
                    ui.error(str(exc));ui.info('The project remains open. Use Commands → Tools and execution history to set up tools.')
                    ui.pause()
            while True:
                if session.a.export_dir:session.export_explicit_directory()
                if session.a.html_out:session.export_html(session.a.html_out)
                if not session.pending and not session.a.no_menu and ui.tty:
                    session.can_go_back=bool(trail);session.nav()
                if not session.pending:break
                request,value=session.pending;session.pending=None
                previous_saved=session.last_saved;previous_dirty=session.dirty
                queued=None
                reference_origin=None
                if request=='references' and isinstance(value,tuple):value,reference_origin=value
                if request=='analyze_entry':request='browse';queued='analysis'
                elif request=='query_analysis':
                    value,queued=value;request='query'
                returning=request=='back'
                if returning:
                    if not trail:continue
                    value=trail.pop();request='browse'
                # Transition copies are private Back navigation, not new results.
                archive=workspace/('navigation_'+uuid.uuid4().hex)
                save_project(session,archive,include_history=False)
                fresh=argparse.Namespace(**vars(session.a))
                fresh.resume=None;fresh.demo=False;fresh.save=False;fresh.out='';fresh.export_dir=''
                fresh.html_out='';fresh._explicit=set()
                fresh._command='';fresh._intent='load';fresh._workflow_goal=''
                fresh._comparison_queries=list(session.comparison_selection)
                fresh.msa='none';fresh.tree=False;fresh.pairwise='none';fresh.report=False;fresh.dotplot=False
                fresh.vsearch='none'
                if request in ('query','reconfigure','references'):
                    from av_project import read_project
                    root,state=read_project(archive)
                    fresh.query=str(root/state['inputs']['query']['path'])
                    fresh.references=str(root/state['inputs']['references']['path']) if state['inputs'].get('references') else None
                    fresh._guided=True;fresh._intent='load'
                    if request=='query':
                        fresh.query_id=value
                        from av_core import sequence_type
                        selected=next(record for record in session.query_list if record.id==value)
                        if sequence_type(selected.seq)!=sequence_type(session.query.seq):fresh.mode=None
                    elif request=='references':
                        fresh.references=str(value);fresh.mode=None;fresh.query_id=(session.comparison_selection or [session.query.id])[0];queued='compare'
                    else:fresh.mode=None
                elif request=='files':
                    _apply_start(fresh,value,workspace/'demo')
                    if fresh.resume:request='resume'
                elif request=='resume':fresh.resume=str(value)
                elif request=='browse':fresh.resume=str(value)
                replacement=None
                try:
                    replacement=resume_project(fresh.resume,fresh,ui) if fresh.resume else Session(fresh,ui)
                    restored_navigation=dict(replacement.navigation_state)
                    if interactive:
                        replacement.library_root=library
                        if request!='browse':
                            for i,path in enumerate(replacement.history,1):
                                imported=workspace/f'import_{len(history):03d}_{i:03d}'
                                shutil.copytree(path,imported);history.append(imported)
                            navigation.update(replacement.navigation_state)
                        replacement.history=history;replacement.navigation_state=navigation
                    if request in ('query','reconfigure','references'):
                        replacement.source_query_path=Path(state['inputs']['query'].get('original_path',replacement.source_query_path))
                        if replacement.refs_path and request!='references':
                            replacement.source_refs_path=Path(state['inputs']['references'].get('original_path',replacement.source_refs_path))
                        elif reference_origin:replacement.source_refs_path=reference_origin
                    if interactive and request!='browse':_register_inputs(replacement)
                    if request=='query' and session.single and session.msa_records:
                        target=replacement.work/'alignment_active.fasta'
                        shutil.copyfile(session.msa_file,target)
                        replacement.msa_file=target;replacement.msa_records=read_fasta(target)
                        replacement.result_sources=dict(session.result_sources)
                        replacement.calculation_runs=dict(session.calculation_runs)
                    if fresh.resume and Path(fresh.resume).resolve().is_relative_to(workspace.resolve()):
                        replacement.last_saved=None;replacement.dirty=True
                    if request=='browse':
                        replacement.last_saved=previous_saved;replacement.dirty=previous_dirty
                        saved_view=restored_navigation.get(replacement.workspace_entry.get('id'),{}).get('tab')
                        replacement._open_view=(saved_view if returning else None) or ('ResultHome' if replacement.workspace_entry.get('kind')=='result' else 'Home')
                    if queued=='analysis':replacement._open_analysis=True
                    elif queued:replacement._queued_action=queued
                    if not fresh.resume:
                        try:replacement.initial()
                        except HomeRequested:pass
                        except CommandError as exc:ui.error(str(exc));ui.pause()
                except (ValueError,CommandError,OSError,KeyboardInterrupt) as exc:
                    if replacement:replacement.close()
                    ui.error(str(exc) or 'Change canceled.');ui.pause();continue
                if not returning:trail.append(archive)
                if not interactive:history.append(archive)
                session.close();session=replacement;session.history=history
                if request!='browse':session.dirty=True
            session.save_and_exit()
            return 0
        except KeyboardInterrupt:
            ui.warn('Session interrupted.')
            failed=True
            return 130
        except (ValueError,CommandError,OSError) as exc:
            ui.error(str(exc));failed=True
            return 1
        finally:
            if session:
                if failed and (session.hsps or session.msa_records or session.tools.history):
                    try:
                        recovery=Path(tempfile.mkdtemp(prefix='alignview_recovery_'))/'project'
                        save_project(session,recovery)
                        ui.warn(f'Recovery copy retained: {recovery}\n  Resume: alignview --resume '+shlex.quote(str(recovery)))
                    except (ValueError,OSError) as exc:
                        # Keep the original workspace if even the recovery copy failed.
                        session.temporary._finalizer.detach()
                        ui.warn(f'Cannot create recovery copy ({exc}). Work retained in {session.work}')
                session.close() if not (failed and session.work.exists() and not session.temporary._finalizer.alive) else None
