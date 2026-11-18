"""Free workspace navigation. Scientific calculations stay in Session and av_core."""
from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import dataclass, replace
import io
from pathlib import Path
import re
import shutil

from av_core import sequence_stats, sequence_type
from av_parameters import snapshot
from av_start import choose_file, choose_sequence, choose_sequences, startup
from av_tools import CommandError
from av_ui import HomeRequested, menu_item as item
from av_advanced import AdvancedActions


@dataclass(frozen=True)
class WorkspaceState:
    sequence_count: int
    reference_count: int
    hit_count: int
    aligned: bool
    searched: bool
    can_align: bool
    align_reason: str
    search_stale: bool = False

    @classmethod
    def from_session(cls,s):
        records=[r for r in s.query_list if r.id in s.sequence_selection] if s.single and s.sequence_selection else s.query_list
        count=len(records) if s.single else 1+len(s.references)
        reason=''
        if count<2:reason='At least two sequences are required. Compare one query with references, or open a FASTA containing multiple sequences.'
        elif s.mode in ('tblastn','tblastx'):
            reason='Multiple translation frames are possible. Use detected segments in the advanced alignment options.'
        elif s.single and len({sequence_type(r.seq) for r in records})>1:
            reason='The file mixes sequence types. Separate nucleotide and protein inputs.'
        calculated=getattr(s,'search_parameters',{})
        current=snapshot('search',s.a,s.mode)
        stale=bool(calculated) and (any(name in calculated and calculated[name]!=value
            for name,value in current.items() if name!='candidate_limit') or
            current.get('candidate_limit',0)>calculated.get('candidate_limit',float('inf')))
        return cls(len(s.query_list),len(s.references),len(s.selected),bool(s.msa_records),
                   bool(s.hsps) or getattr(s,'search_completed',False),not reason,reason,stale)




class Workbench(AdvancedActions):
    def __init__(self,session):
        self.s=session
        self.ui=session.ui

    @property
    def state(self):return WorkspaceState.from_session(self.s)

    def refresh_context(self):
        s=self.s;state=self.state
        self.ui.workflow=()
        source=(f'{s.source_query_path.name} · {state.sequence_count} sequences · active: {s.query.id}' if s.single
                else f'Query: {s.query.id} · {state.reference_count} references')
        if s.comparison_selection and any(section in self.ui.route for section in ('References','Comparison')):
            names=', '.join(s.comparison_selection)
            source=f'{s.source_query_path.name} · {len(s.comparison_selection)} selected: {names}'
        results=[]
        if state.search_stale:results.append('New search settings prepared; existing matches are retained')
        if state.searched:results.append(f'{state.hit_count} matches')
        if state.aligned:results.append(f'alignment ready ({len(s.msa_records[0].seq):,} columns)')
        if getattr(s,'vsearch_results',{}):results.append(f'VSEARCH: {len(s.vsearch_results)} analyses')
        run=getattr(s,'calculation_runs',{}).get('msa',{}).get('parameters',{})
        if run and run!=snapshot('msa',s.a,s.mode,run.get('strategy','')):
            results.append('alignment uses the previous parameters')
        if not results:results.append('no calculation started')
        preparing=any(part in self.ui.route for part in ('Find matching sequences','Comparison','Target dataset'))
        comparing=not s.single or preparing
        if comparing:
            from av_workspace import comparison_context
            query_ids=(s.comparison_selection or ([] if s.single or s._comparison_needs_selection else [s.query.id])) if preparing else None
            self.ui.context=(*comparison_context(s,query_ids),' · '.join(results))
        else:self.ui.context=(source,' · '.join(results))
        self.ui.compact_context_lines=2 if comparing else 1

    def choose(self,title,items,subtitle='',**kwargs):
        self.refresh_context()
        # User-facing choices always have short, numbered keys, even when the
        # domain value is a method name such as "diamond-blastx".
        mapping={}
        if any(len(entry.key)>1 for entry in items):
            mapping={str(i):entry.key for i,entry in enumerate(items,1)}
            items=[replace(entry,key=str(i)) for i,entry in enumerate(items,1)]
        selected=self.ui.menu(title,items,subtitle,**kwargs)
        return mapping.get(selected,selected)

    def run_action(self,action,*,changes=False,artifact=None,title=None):
        try:
            result=self.ui.present(action)
            if changes:self.s.dirty=True
            if artifact and self.s.library_root is None and (self.s.work/artifact).is_file():
                self.ui.page((self.s.work/artifact).read_text(errors='replace'),title or 'Results')
            self.refresh_context()
            return result
        except KeyboardInterrupt:
            self.ui.notice='Action interrupted. The project remains open.'
        except (ValueError,CommandError,OSError) as exc:
            self.ui.error(str(exc))
            self.ui.page(str(exc),'Action not completed')
        return None

    def run(self):
        from av_workspace import Workspace
        self.ui.route=['Workspace']
        self.workspace=Workspace(self.s)
        queued=getattr(self.s,'_queued_action',None);self.s._queued_action=None
        if getattr(self.s,'_open_analysis',False):
            self.s._open_analysis=False
            try:queued=self.new_analysis()
            except (HomeRequested,KeyboardInterrupt):self.workspace.notice='Back to the workspace. No calculation started.'
        while True:
            before=self.s.workspace_entry.get('id')
            try:
                self.refresh_context()
                action=queued or self._workspace_event();queued=None
                self.workspace.remember()
                if action=='back':
                    if self.workspace.tab=='ResultHome':
                        self.workspace.switch('Home',record=False);continue
                    if not self.workspace.back():
                        if getattr(self.s,'can_go_back',False):self.s.pending=('back',None);return
                        self.workspace.switch('Home',record=False)
                    continue
                if action=='quit':
                    if self.confirm_exit():return
                    continue
                if action=='commands':action=self.command_palette()
                if action=='analysis':action=self.new_analysis()
                if self.s.pending:return
                if not action:continue
                handlers={'view':self.view_alignment,'align':self.align,'compare':self.compare,
                    'inspect':self.inspect_reference,'vsearch_results':self.vsearch_results,'tree':self.tree,'tree_view':self.show_tree,
                    'sequence':self.inspect_sequence,'query':self.choose_query,'use_query':lambda:self.choose_query(from_data=True),'files':self.open_files,
                    'history':self.open_history,'save':self.save,'checkpoint':self.s.save_checkpoint,
                    'report':self.s.open_report,'export':self.export,'advanced':self.advanced,
                    'html_export':self.s.copy_report,'analysis':self.new_analysis,
                    'home':lambda:self.workspace.switch('Home'),'manage':self.manage_files,
                    'inspect_data':self.inspect_data,'sequences_view':lambda:self.workspace.switch('Sequences'),
                    'open_result':self.open_result,
                    'rename_result':self.rename_result,
                    'datasets':lambda:self.workspace.switch('Data'),'results':lambda:self.workspace.switch('Results'),
                    'open_entry':self.open_entry,'artifact':self.read_artifact,'reuse':self.reuse_result,
                    'inputs':self.open_inputs,
                    'conservation':lambda:self.run_action(self.s.msa_analysis,changes=True),
                    'distances':lambda:self.run_action(self.s.msa_distances,changes=True),
                    'statistics':lambda:self.run_action(self.s.show_stats,changes=True),
                    'profile_scan':lambda:self.run_action(self.s.profiles,changes=True),
                    'composition':lambda:self.run_action(self.s.composition,changes=True),
                    'tools':self.tools_and_provenance,'help':self.documentation,
                    'align_selection':lambda:self.align('selection'),'trace':self.trace_alignment,
                    'vsearch':self.vsearch,'pairwise':self.pairwise,
                    'motif':lambda:self.run_action(self.s.motif,changes=True,artifact='motifs.tsv',title='Motif positions'),
                    'translate':lambda:self.run_action(self.s.translate_frames,changes=True),
                    'orfs':lambda:self.run_action(self.s.find_orfs,changes=True),
                    'extract':lambda:self.run_action(self.s.extract_region,changes=True),
                    'dotplot':lambda:self.run_action(self.s.make_dotplot,changes=True),
                    'translation_table':lambda:self.edit_parameters('translation'),
                    'pick_references':self.pick_references}
                self.run_action(handlers[action])
                if getattr(self.s,'_queued_action',None):queued=self.s._queued_action;self.s._queued_action=None
                if self.s.pending:return
                entry=self.s.workspace_entry
                if entry.get('id')!=before and entry.get('kind')=='result':self.workspace.switch('ResultHome')
                elif action=='pick_references':self.workspace.switch('Matches')
                self.workspace.remember()
                self.workspace.notice=self.ui.notice
            except HomeRequested:
                entry=self.s.workspace_entry
                if entry.get('id')!=before and entry.get('kind')=='result':self.workspace.switch('ResultHome')
                else:self.workspace.switch('Home',record=False)
                self.workspace.notice='Back to the workspace. Data and results are retained.'
            except KeyboardInterrupt:
                self.workspace.notice='Action canceled. The workspace remains open.'
            except (ValueError,CommandError,OSError) as exc:
                self.ui.error(str(exc));self.ui.page(str(exc),'Action not completed')

    def _workspace_event(self):
        from av_terminal import screen,read_key,draw
        from av_workspace import TABS
        w=self.workspace
        if w.tab in ('Home','ResultHome'):return self.landing_event()
        if w.tab in TABS:return self.catalog_event()
        if self.ui.simple:
            primary,label=w.primary()
            selected=self.choose('Workspace',[
                item('1',label,w._summary()),item('2','New analysis','Choose any operation'),
                item('3','Data','Browse loaded datasets'),item('4','Results','Browse completed analyses'),
                item('5','Add data','Load another FASTA'),item('6','Save / Export','Optional, for after closing'),
                item('7','Commands','All actions and settings'),item('8','Close workspace','End this session')])
            if selected is None:return 'back'
            if selected=='4':
                w.switch('Results')
                return None
            return {'1':primary,'2':'analysis','3':'datasets','5':'files','6':'save','7':'commands','8':'quit'}[selected]
        previous=None
        with screen(mouse=True):
            while True:
                signature=(self.ui.size(),w.tab,w.active,w.query,w.filtering,w.offset,w.notice,
                    tuple(getattr(self.s,'reference_selection',[])))
                if signature!=previous:
                    draw(w.render(*self.ui.size()));previous=signature
                key=read_key()
                if key=='tick':continue
                if key=='eof':return 'quit'
                click=re.fullmatch(r'mouse:0;(\d+);(\d+)M',key)
                if click:
                    x,y=map(int,click.groups())
                    if y==2:
                        if x<=10:w.switch('Data')
                        elif x<=25:w.switch('Results')
                        continue
                    if y==3:
                        return 'files' if x<15 else 'analysis' if x<32 else 'report' if x<49 else 'save'
                    if y in w.mouse_rows:
                        w.active=w.mouse_rows[y]
                        if x<=5 and w.tab in ('Sequences','Matches'):key=' '
                        else:return w.primary()[0]
                    else:continue
                if w.filtering:
                    if key=='esc':w.query='';w.filtering=False
                    elif key in ('\r','\n'):w.filtering=False
                    elif key in ('\x7f','\b'):w.query=w.query[:-1];w.active=0
                    elif len(key)==1 and key.isprintable():w.query+=key;w.active=0
                    continue
                rows=w.rows()
                if key in ('q','Q'):return 'quit' if w.tab in TABS else 'back'
                if key in (':','\x10'):return 'commands'
                if key in ('?','F1'):return 'help'
                if key=='\x13':return 'checkpoint'
                if key=='s':return 'save'
                if key=='o':return 'files'
                if key=='c':return 'compare'
                if key=='a':return 'align_selection' if self.s.reference_selection and w.tab=='Matches' else 'align'
                if key=='v':return 'view'
                if key=='n':return 'analysis'
                if key=='w':return 'report'
                if key=='h':w.switch('Results');continue
                if key=='F2':return 'home'
                if key=='d':w.switch('Data');continue
                if key=='i' and w.tab=='Sequences':return 'sequence'
                if key=='t':return 'tree_view' if (self.s.work/'tree_fasttree.nwk').is_file() else 'tree'
                if key=='r' and w.tab=='Alignment':return 'trace'
                if key=='u':return 'use_query' if w.tab=='Sequences' else 'reuse'
                if key in ('1','2'):w.switch(TABS[int(key)-1]);continue
                if key=='\t':w.switch('Results' if w.tab in ('Data','Sequences') else 'Data');continue
                if key in ('left','right'):return 'back' if key=='left' else w.primary()[0]
                if key in ('esc','F2'):
                    if w.query:w.query='';w.active=0
                    else:return 'back'
                    continue
                if key=='/':w.filtering=True;w.query='';w.active=0;continue
                if key==' ' and w.tab=='Matches' and rows:
                    name=rows[w.active]['subject']
                    if name in self.s.reference_selection:self.s.reference_selection.remove(name)
                    else:self.s.reference_selection.append(name)
                    self.s.save_settings();w.notice=f'{len(self.s.reference_selection)} reference(s) selected. Press a to align them.'
                    continue
                if key==' ' and w.tab=='Sequences' and rows:
                    row=rows[w.active];name=row['record'].id
                    if row['role']=='Reference':
                        w.notice='Choose references from the matches table, or open their dataset in Data.'
                    elif name in self.s.sequence_selection:self.s.sequence_selection.remove(name)
                    else:self.s.sequence_selection.append(name)
                    self.s.save_settings();continue
                delta=1 if key in ('down','j') else -1 if key in ('up','k') else w.page_size if key=='pagedown' else -w.page_size if key=='pageup' else 0
                if key.startswith('mouse:64;'):delta=-3
                elif key.startswith('mouse:65;'):delta=3
                if w.tab in ('Tree','Artifact','Vsearch'):w.offset=max(0,w.offset+delta)
                else:w.active=max(0,min(len(rows)-1,w.active+delta))
                if key=='home':w.active=w.offset=0
                elif key=='end':w.active=max(0,len(rows)-1);w.offset=10**9
                if key in ('\r','\n'):return w.primary()[0]

    def landing_event(self):
        title,options,summary,status=self.workspace.landing()
        from av_workspace import comparison_context
        self.ui.context=(summary,status) if self.s.single else comparison_context(self.s)
        self.ui.compact_context_lines=1 if self.s.single else 2
        self.ui.notice=''
        shortcuts={'a':'align','c':'compare','n':'analysis','o':'files','r':'results','d':'datasets',
                   's':'save','w':'report',':':'commands','\x10':'commands','\x13':'checkpoint','?':'help'}
        if self.workspace.tab=='Home':shortcuts.update(q='quit',Q='quit')
        with self.ui.section(self.s.source_query_path.name):
            selected=self.ui.menu(title,[item(c.key,c.label,c.detail,enabled=c.enabled,reason=c.reason) for c in options],
                subtitle=summary if not self.s.single else '',
                hint=(status if not self.s.single and self.workspace.tab=='ResultHome' else
                      '↑↓ choose · Enter open · Esc back · q '+('close' if self.workspace.tab=='Home' else 'file menu')),
                shortcuts=shortcuts)
        if selected is None:return 'back' if self.workspace.tab=='Home' else 'home'
        return next((c.action for c in options if c.key==selected),selected)

    def catalog_event(self):
        is_data=self.workspace.tab=='Data'
        rows=self.workspace.rows()
        if not is_data:
            members={identifier for row in rows for identifier in row.get('members',[])}
            rows=[row for row in rows if row.get('id') not in members]
        if not rows:
            self.ui.page('No completed analyses yet.\n\nReturn to the file menu and choose an analysis. Results will appear here after a successful run.','Previous results')
            return 'back'
        with self.ui.section('Loaded files' if is_data else 'Previous results'):
            selected=self.choose('Choose a file' if is_data else 'Choose a previous result',[
                item(str(i),row['label'],row.get('source','')+(' · open the file menu' if is_data else ' · open its result summary'))
                for i,row in enumerate(rows,1)],hint='Enter selects. Esc returns. No calculation is started.')
        if selected is None:return 'back'
        row=rows[int(selected)-1];self.workspace.active=int(selected)-1
        if row.get('path'):self.s.pending=('browse',row['path'])
        else:self.workspace.switch('Home')
        return None

    def manage_files(self):
        with self.ui.section('Files and results'):
            selected=self.choose('Files and previous results',[
                item('1','Choose a loaded file','Return to its file menu; keep every completed analysis.'),
                item('2','Load another FASTA','Add sequences without removing the files already loaded.'),
                item('3','Previous results','Reopen completed analyses without running them again.'),
                item('4','Save or export','Optional: retain a project after closing, or export a result.'),
                item('5','Settings and help','Calculation settings, available software and the manual.'),
            ])
        if selected=='1':self.workspace.switch('Data')
        elif selected=='2':self.open_files()
        elif selected=='3':self.workspace.switch('Results')
        elif selected=='4':self.save()
        elif selected=='5':
            with self.ui.section('Settings and help'):
                action=self.choose('Settings and help',[
                    item('1','Analysis settings','Engine options and calculation profiles.'),
                    item('2','Available tools and execution logs','Installed programs and commands actually run.'),
                    item('3','Read the manual','Offline explanations and examples.'),
                    item('4','All commands','Search all operations by name.'),
                ])
            if action=='1':self.advanced()
            elif action=='2':self.tools_and_provenance()
            elif action=='3':self.documentation()
            elif action=='4':
                chosen=self.command_palette()
                if chosen:self.s._queued_action=chosen

    def inspect_data(self):
        with self.ui.section('Inspect sequences'):
            selected=self.choose('Look at the loaded sequences',[
                item('1','Browse sequences','Read identifiers, lengths and sequence content.'),
                item('2','Read an existing alignment','Choose this only if this file is already aligned.',
                     enabled=self.workspace.can_read_alignment(),reason='Needs at least two already aligned sequences with equal-length rows.'),
            ])
        if selected=='1':self.workspace.switch('Sequences')
        elif selected=='2':self.view_alignment()

    def open_result(self):
        view=self.s.workspace_entry.get('view','')
        if view=='Alignment':self.view_alignment()
        elif view=='Matches':self.workspace.switch('Matches')
        elif view=='Tree':self.show_tree()
        elif view=='Vsearch':self.vsearch_results()
        elif view=='Batch':self.open_comparison_batch()
        else:self.read_artifact()

    def rename_result(self):
        from av_library import entry_label,rename_result
        entry=self.s.workspace_entry
        if entry.get('kind')!='result':
            self.ui.notice='Open a completed result first.';return
        current=entry_label(self.s,entry)
        with self.ui.section('Rename result'):
            label=self.ui.ask('Result name (q to cancel)',current,'1–120 characters; retained when you save the project.')
        if not label or label.strip().casefold()=='q':return
        label=rename_result(self.s,entry['id'],label)
        self.ui.notice='Result named: '+label

    def open_comparison_batch(self):
        from av_library import entries
        catalog={row['id']:row for row in entries(self.s)}
        runs=self.s.comparison_batch.get('runs',[])
        with self.ui.section('Comparison results'):
            selected=self.choose('Results by query sequence',[
                item(str(i),run['query_id']+' · '+run['status'].replace('_',' '),
                     (f'{run.get("matches",0)} matching references. Open this sequence’s comparison.' if run['status']=='completed' else
                      run.get('error','Not run because the comparison was interrupted.')))
                for i,run in enumerate(runs,1)],hint='Each query has its own identity, coverage and retained references.')
        if selected:
            run=runs[int(selected)-1];row=catalog.get(run.get('result_id'))
            if row:self.s.pending=('browse',row['path'])
            else:self.ui.page(run.get('error','This sequence was not compared.'),'Comparison '+run['status'].replace('_',' '))

    def command_palette(self):
        commands=self.workspace.commands()
        with self.ui.section('Commands'):
            selected=self.choose('Find a command',[
                item(str(i),c.label,c.detail,enabled=c.enabled,reason=c.reason)
                for i,c in enumerate(commands,1)],hint='Type to search commands. Enter runs the selected action.',searching=True)
        return commands[int(selected)-1].action if selected else None

    def new_analysis(self):
        """Offer a few explained groups; choose scope before a sequence operation."""
        if self.workspace.tab in ('Data','Results'):
            row=self.workspace.selected()
            if row and row.get('path'):
                self.s.pending=('analyze_entry',row['path']);return
        commands={command.action:command for command in self.workspace.commands()}
        nucleotide_records=[record for record in self.s.query_list if sequence_type(record.seq)=='DNA/RNA']
        with self.ui.section('Other analyses'):
            group=self.choose('Other analyses',[
                item('sequence','Sequence checks, motifs and regions','Check composition, find a motif or extract a sequence region.'),
                item('translation','Translate DNA or find ORFs','Produce protein sequences from a chosen nucleotide sequence.',
                     enabled=bool(nucleotide_records),reason='This file contains no nucleotide sequences to translate.'),
                item('alignment','Alignment statistics and tree','Read conservation, calculate distances or build a tree.',
                     enabled=bool(self.s.msa_records),reason='Create or open an alignment first, from the file menu.'),
                item('pair','Compare two sequences','A pairwise alignment or a dotplot.'),
                item('vsearch','Nucleotide collections','Search, cluster, dereplicate or check chimeras with VSEARCH.',
                     enabled=all(sequence_type(r.seq)=='DNA/RNA' for r in self.s.query_list),reason='Requires a nucleotide dataset.'),
            ],hint='Esc returns to the file menu. Choosing a group does not start a calculation.')
            groups={
                'sequence':('Sequence checks, motifs and regions',['statistics','motif','extract','profile_scan','composition']),
                'translation':('Translation and ORFs',['translate','orfs','translation_table']),
                'alignment':('Alignment analyses',['conservation','distances','tree']),
                'pair':('Compare two sequences',['pairwise','dotplot']),
            }
            if group in groups:
                title,actions=groups[group]
                extra={'profile_scan':item('profile_scan','Sliding-window profile','Composition or entropy along the selected sequence.'),
                       'composition':item('composition','Count sequence words','Count k-mers in the selected sequence.')}
                options=[]
                for action in actions:
                    if action in extra:options.append(extra[action]);continue
                    command=commands[action]
                    requires_nt=action in ('motif','translate','orfs','translation_table')
                    options.append(item(action,command.label,command.detail,
                        enabled=bool(nucleotide_records) if requires_nt else command.enabled,
                        reason='Choose a nucleotide dataset for this action.' if requires_nt else command.reason))
                selected=self.choose(title,options)
            else:selected=group
        if selected in ('motif','translate','orfs','extract','dotplot','profile_scan','composition'):
            row=self.workspace.selected() if self.workspace.tab=='Sequences' else None
            record=row['record'] if row and row.get('role')!='Reference' else None
            candidates=nucleotide_records if selected in ('motif','translate','orfs') else self.s.query_list
            if record and record not in candidates:record=None
            if record is None and len(self.s.query_list)>1:
                record=choose_sequence(self.ui,candidates,'Choose the sequence to analyse',self.s.query.id)
                if record is None:return
            if record and record.id!=self.s.query.id:
                self.s.pending=('query_analysis',(record.id,selected));return
        return selected

    def open_entry(self):
        row=self.workspace.selected()
        if row is None:return
        if row.get('path'):self.s.pending=('browse',row['path'])
        else:self.workspace.switch(row.get('view','Sequences'))

    def read_artifact(self):
        from av_project import local_path
        from av_library import entry_label
        artifact=self.s.workspace_entry.get('artifact','')
        if artifact:
            path=local_path(self.s.work,artifact)
            self.ui.page(path.read_text(encoding='utf-8',errors='replace'),entry_label(self.s,self.s.workspace_entry))
        else:self.vsearch_results()

    def open_inputs(self):
        from av_library import entries
        identifiers=self.s.workspace_entry.get('inputs',[])
        rows=[row for row in entries(self.s) if row['id'] in identifiers]
        if not rows:return
        if len(rows)==1:self.s.pending=('browse',rows[0]['path']);return
        selected=self.choose('Open input',[item(str(i),row['label']) for i,row in enumerate(rows,1)])
        if selected:self.s.pending=('browse',rows[int(selected)-1]['path'])

    def reuse_result(self):
        from av_project import local_path
        kind=self.s.workspace_entry.get('analysis','')
        artifact=self.s.workspace_entry.get('artifact','')
        if kind=='orfs':artifact='orfs.faa'
        if kind=='msa' and self.s.msa_file:path=self.s.msa_file
        elif artifact and Path(artifact).suffix.lower() in ('.fa','.fna','.faa','.fasta'):
            path=local_path(self.s.work,artifact)
        else:
            self.ui.page('This result does not contain a reusable sequence set. Open a translation, extracted region, ORFs or alignment.','Use as data');return
        from av_start import StartRequest
        self.s.pending=('files',StartRequest('load',path,origin=self.s.workspace_entry.get('id','')))

    def pick_references(self):
        from av_workspace import evidence_rows
        rows,_=evidence_rows(self.s)
        if not rows:
            self.ui.page('Run a comparison with retained matches first. Press c to review the search plan.','References required');return
        with self.ui.section('Reference selection'):
            choice=self.choose('Choose references to align',[
                item('1',f'Use all {len(rows)} matching references','Include the active query and every currently retained match'),
                item('2','Pick individual references in Matches','Space checks a reference; n prepares the checked alignment'),
            ],hint='This choice prepares the alignment scope; it does not start an aligner.')
        if choice=='1':
            self.s.reference_selection=[row['subject'] for row in rows];self.s.save_settings()
        if choice:
            self.workspace.switch('Matches')
            self.ui.notice=(f'{len(self.s.reference_selection)} references checked. Press a to prepare the alignment.'
                if choice=='1' else 'Use Space to check references, then press a to prepare the alignment.')

    def inspect_sequence(self):
        row=self.workspace.selected()
        if not row or 'record' not in row:return
        r=row['record'];stats=sequence_stats(r)
        if self.s.workflow_goal=='sequence':self.s._workflow_reviewed=True
        self.ui.page(r.id+'\n'+r.description+'\n\n'+f'{stats["kind"]} · {len(r.seq):,} characters\n\n'+
            '\n'.join(f'{i+1:>9}  {r.seq[i:i+80]}' for i in range(0,len(r.seq),80)),'Sequence details')

    def inspect_reference(self):
        row=self.workspace.selected()
        if not row or 'subject' not in row:return
        if row.get('global_match'):
            self.ui.page(f'Reference: {row["subject"]}\nVSEARCH global identity: {row["best_identity"]:.2f}%\n'+
                f'Query coverage: {row["query_coverage"]:.2f}%\nReference coverage: {row["reference_coverage"]:.2f}%\n\n'+
                'Native global-search evidence. No BLAST E-value or bitscore is assigned.','Reference evidence')
            return
        from av_units import units_for
        units=units_for(self.s.mode,sequence_type(self.s.query.seq))
        hsps=self.s.selected.get(row['subject'],[])
        target=self.s.refmap[row['subject']]
        output=io.StringIO()
        with redirect_stdout(output):
            print(f'Query: {self.s.query.id} · {len(self.s.query.seq):,} {units.query}\nQuery FASTA: {self.s.source_query_path}\n'+
                f'Target sequence: {target.id} · {len(target.seq):,} {units.subject}\nTarget FASTA: {self.s.source_refs_path}\n'+
                f'Method: {self.s.mode.upper()}\nStrongest-HSP identity: {row["best_identity"]:.2f}%\n'+
                f'Query coverage (interval union): {row["query_coverage"]:.2f}%\n')
            print(f'Aligned regions ({len(hsps)} local segments)\n'+
                'Coordinates are 1-based inclusive in the prepared ungapped sequences.\n'+
                f'Query positions: {units.query}; target positions: {units.subject}; aligned rows: {units.aligned}.')
            for index,hsp in enumerate(hsps,1):
                print(f'  {index}. Query {hsp.qstart} → {hsp.qend} {units.query}  /  Target {hsp.sstart} → {hsp.send} {units.subject}'+
                      f' · {hsp.identity:.1f}% identity · {hsp.align_length} aligned columns')
            print('\nRead each segment below: Q = query, S = target; | marks identical residues.\n'+
                  'Local segments can be separated by unaligned regions; this is not a full-length alignment.\n')
            for hsp in hsps:self.s.print_alignment(hsp)
        self.ui.page(output.getvalue(),'Aligned regions · '+target.id)

    def choose_query(self,from_data=False):
        selected=self.workspace.selected() if from_data and self.workspace.tab=='Sequences' else None
        record=selected.get('record') if selected and selected.get('role')!='Reference' else None
        if record is None:record=choose_sequence(self.ui,self.s.query_list,'Choose the active query',self.s.query.id)
        if record and record.id!=self.s.query.id:self.s.pending=('query',record.id)
        elif record:self.ui.notice='This query is already active.'

    def open_files(self):
        from av_start import StartRequest
        path=choose_file(self.ui,'Add sequences to this workspace')
        if path:self.s.pending=('files',StartRequest('load',path))

    def open_history(self):
        self.workspace.switch('Results')

    def trace_alignment(self):
        row=self.workspace.selected()
        if row and 'record' in row:self.s._trace_id=row['record'].id
        self.view_alignment()

    def compare(self):
        with self.ui.section('Find matching sequences'):
            self.refresh_context()
            return self.prepare_comparison()

    def prepare_comparison(self):
        if self.s.single or self.s._comparison_needs_selection:
            selected=choose_sequences(self.ui,self.s.query_list,self.s.comparison_selection or self.s.sequence_selection)
            if selected is None:return
            self.s.comparison_selection=[record.id for record in selected]
            self.s._comparison_needs_selection=False;self.s.save_settings()
            if not self.s.single and selected[0].id!=self.s.query.id:
                self.s.pending=('query_analysis',(selected[0].id,'compare'));return
        if self.s.single:
            self.choose_target()
            return
        return self.search()

    def choose_target(self):
        """A target change creates a fresh comparison, keeping earlier results."""
        with self.ui.section('Target dataset'):
            rows=self.workspace.catalog('data') if hasattr(self,'workspace') else []
            choice=self.choose('Choose the target dataset',[
                *[item(str(i),row['label'],
                       ('Search within the query file itself' if row.get('source_hash')==self.s.input_hashes['query'] else
                        'Search for matching sequences in this loaded FASTA')+' · '+row.get('source',''))
                  for i,row in enumerate(rows,1)],
                item('file','Load a target FASTA','Choose the file containing sequences to search against'),
            ],hint='Query sequences are searched AGAINST this target. The target is a local FASTA dataset.')
            if choice is None:return False
            if choice=='file':
                target=choose_file(self.ui,'Choose target FASTA')
                if target:self.s.pending=('references',target)
            else:
                from av_project import read_project,local_path
                root,state=read_project(rows[int(choice)-1]['path'])
                info=state['inputs']['query']
                self.s.pending=('references',(local_path(root,info['path']),Path(info.get('original_path',info['path']))))
            return bool(self.s.pending)

    def comparison_inputs(self,mode):
        ids=self.s.comparison_selection or [self.s.query.id]
        lines=['Method: '+mode.upper(),'','Query FASTA: '+str(self.s.source_query_path),
               f'Sequences to search ({len(ids)}):']
        selected=set(ids)
        lines.extend('  '+r.id+f' · {len(r.seq):,} characters' for r in self.s.query_list if r.id in selected)
        lines+=['','Target FASTA: '+str(self.s.source_refs_path),
                f'Search against all {len(self.s.references)} target sequences.',
                'Target sequence preview (up to 8):']
        lines.extend('  '+r.id+f' · {len(r.seq):,} characters · {sequence_type(r.seq)}' for r in self.s.references[:8])
        lines+=['','Each query is searched independently against this local target dataset.',
                'Choose Target dataset in the search menu to change it.']
        self.ui.page('\n'.join(lines),'Search inputs: queries → target')

    def explore(self):
        with self.ui.section('Results'):
            while True:
                state=self.state
                entries=[];actions={}
                if state.aligned:
                    entries.append(item('msa','Read the alignment',f'{len(self.s.msa_records)} sequences · navigation and position search'))
                    actions['msa']=self.view_alignment
                if state.searched:
                    entries.append(item('hits','Matches',f'{state.hit_count} references · identity and query coverage'))
                    actions['hits']=self.hits
                if state.aligned and len(self.s.msa_records)>=3:
                    existing=(self.s.work/'tree_fasttree.nwk').is_file()
                    entries.append(item('tree','View the tree' if existing else 'Create a tree',
                        'Exploratory FastTree tree · SVG figure available in the report'))
                    actions['tree']=self.show_tree if existing else self.tree
                if self.s.vsearch_results:
                    entries.append(item('vsearch','VSEARCH results','Global matches, clusters, unique sequences or chimeras'))
                    actions['vsearch']=self.vsearch_results
                if not entries:
                    self.ui.page('Press c to compare or a to align from the workspace.','No results available');return
                choice=self.choose('Results',entries,hint='Find Checks and statistics under Commands → Additional analyses and settings.')
                if choice is None:return
                self.run_action(actions[choice])

    def view_alignment(self):
        if not self.s.msa_records:
            lengths={len(r.seq) for r in self.s.query_list}
            if not self.s.single or len(self.s.query_list)<2 or len(lengths)!=1:
                self.ui.page('The loaded sequences are not an equal-length alignment. Press a to align them first.','No alignment available');return
            # Choosing View explicitly asserts alignment; opening a dataset does not.
            self.s.a.msa='view';self.s.load_msa()
        if hasattr(self,'workspace') and self.workspace.tab=='Alignment':
            selected=self.workspace.selected()
            if selected:self.s._viewer_focus_id=selected['record'].id
        self.s.show_msa_window()

    def sequences(self):
        with self.ui.section('Sequences'):
            while True:
                records=self.s.query_list if self.s.single else [self.s.query,*self.s.references]
                selected=choose_sequence(self.ui,records,'Browse a sequence')
                if selected is None:return
                stats=sequence_stats(selected)
                text=[selected.id,selected.description,'',f'Type: {stats["kind"]}',
                      f'Length: {len(selected.seq):,} characters','',
                      *[f'{i+1:>8}  {selected.seq[i:i+80]}' for i in range(0,len(selected.seq),80)]]
                self.ui.page('\n'.join(text),'Sequence: '+selected.id)

    def hits(self):
        with self.ui.section('Matches'):
            while True:
                rows=self.s.coverage_rows()
                if not rows:
                    self.ui.page('No matches retained.\n\nRelax the filters in Advanced → Calculation settings → Search,\nor press c in the workspace to rerun the comparison.','Search completed without matches')
                    return
                options=[item(str(i),r['subject'],f'{r["best_identity"]:.1f}% local identity · {r["query_coverage"]:.1f}% of the query · {r["hsps"]} segments') for i,r in enumerate(rows,1)]
                choice=self.choose('Matches',options,f'{len(rows)} retained references · limit: {self.s.a.hits}',
                    hint='Choose a reference to inspect its aligned segments. / searches its identifier.')
                if choice is None:return
                selected=rows[int(choice)-1]
                output=io.StringIO()
                with redirect_stdout(output):
                    print(f'Reference: {selected["subject"]}\nIdentity: {selected["best_identity"]:.1f}% on the best local segment')
                    print(f'Query coverage: {selected["query_coverage"]:.1f}% (union of segments)')
                    for hsp in self.s.selected[selected['subject']]:self.s.print_alignment(hsp)
                self.ui.page(output.getvalue(),'Aligned segments: '+selected['subject'])


    def analyze(self):
        """Compatibility entry point for additional analyses."""
        return self.additional_analyses()

    def align(self,strategy=None):
        if not self.state.can_align and not (self.s.selected and self.s.mode in ('tblastn','tblastx')):
            self.ui.page(self.state.align_reason,'Cannot align these inputs');return False
        with self.ui.section('Alignment'):
            strategy=strategy or ('selection' if self.s.single and self.s.sequence_selection else 'all' if self.s.single or self.s.mode.startswith('global-') else 'hits' if self.s.selected else 'all')
            if self.s.single and strategy=='selection' and len(self.s.sequence_selection)<2:
                self.ui.page('Check at least two sequences with Space, or clear the selection to align the whole dataset.','Alignment selection');return False
            if self.s.mode in ('tblastn','tblastx') and strategy!='selection':strategy='hsp'
            labels={'all':'All sequences','hits':'Query and retained references','fragments':'Add a partial sequence','hsp':'Detected segments only','selection':f'Query and {len(self.s.reference_selection)} checked references'}
            if self.s.single:labels['selection']=f'{len(self.s.sequence_selection)} checked sequences'
            while True:
                engine='mafft' if strategy=='fragments' else self.s.a.aligner
                choice=self.choose('Create an alignment',[
                    item('1','Run the alignment',f'{labels[strategy]} · {engine.upper()} · {self.s.a.threads} CPU'),
                    item('2','Advanced options','Engine, sequences, settings, commands and profiles'),
                ],subtitle=self.s.tools.engine_info('muscle5' if engine=='muscle' and not self.s.tools.has('muscle') else engine),
                    hint='The result opens directly in the viewer.')
                if choice is None:return False
                if choice=='2':strategy=self.alignment_options(strategy);continue
                aliases=('muscle','muscle5') if engine=='muscle' else (engine,)
                if not self.ensure_engine(aliases,engine):continue
                return bool(self.run_action(lambda:self.s.msa_run(strategy),changes=True))


    def recommended_mode(self):
        from av_start import suggest_modes
        from collections import Counter
        qtype=sequence_type(self.s.query.seq)
        rtype=Counter(sequence_type(r.seq) for r in self.s.references).most_common(1)[0][0]
        mode,_,_=suggest_modes(qtype,rtype)
        if mode in ('blastp','blastx') and not self.s.tools.has(mode) and self.s.tools.has('diamond'):
            mode='diamond-'+mode
        return mode


    def search(self):
        with self.ui.section('Comparison'):
            mode=self.s.mode if self.s.mode and self.s.mode not in ('msa','global-dna','global-aa') else self.recommended_mode()
            if not mode:mode=self.choose_search_method()
            if not mode:return False
            while True:
                query_ids=self.s.comparison_selection or [self.s.query.id]
                query_label=query_ids[0] if len(query_ids)==1 else f'{len(query_ids)} selected sequences'
                choice=self.choose('Find matching sequences',[
                    item('1','Run '+mode.upper()+' search',f'{query_label} → {self.s.source_refs_path.name} · {len(self.s.references)} target sequences'),
                    item('2','Sequences to search',self.s.source_query_path.name+' · '+', '.join(query_ids),enabled=len(self.s.query_list)>1,
                         reason='The file contains a single query.'),
                    item('3','VSEARCH parameters' if mode=='vsearch' else 'Search threshold',
                         f'Identity ≥ {self.s.a.vsearch_id:.0%}' if mode=='vsearch' else f'{self.s.a.preset} · E-value ≤ {self.s.a.evalue:g}'),
                    item('4','Maximum matches per query',str(self.s.a.hits)),
                    item('5','Advanced options','Method, filters, translation, CPUs, settings and commands'),
                    item('6','Target dataset: '+self.s.source_refs_path.name,'Change the FASTA searched for matching sequences'),
                    item('7','Review search inputs','Full file paths, selected sequence IDs and search direction'),
                ],subtitle=self.s.tools.engine_info('diamond' if mode.startswith('diamond-') else 'mmseqs' if mode=='mmseqs-blastp' else mode),
                    hint='Search direction: query sequences → target FASTA. Esc returns without running.')
                if choice is None:return False
                if choice=='6':
                    if self.choose_target():return False
                    continue
                if choice=='7':self.comparison_inputs(mode);continue
                if choice=='2':
                    selected=choose_sequences(self.ui,self.s.query_list,query_ids)
                    if selected is None:continue
                    self.s.comparison_selection=[record.id for record in selected];self.s.save_settings()
                    if selected[0].id!=self.s.query.id:self.s.pending=('query_analysis',(selected[0].id,'compare'));return False
                    continue
                if choice=='3':
                    if mode=='vsearch':self.edit_parameters('vsearch',strategy='search')
                    else:self.choose_preset()
                    continue
                if choice=='4':self.s._change_numeric('hits','Maximum number of references',int,1);continue
                if choice=='5':
                    if mode=='vsearch':self.edit_parameters('vsearch',strategy='search')
                    else:mode=self.search_options(mode)
                    continue
                program='diamond' if mode.startswith('diamond-') else 'mmseqs' if mode=='mmseqs-blastp' else mode
                if not self.ensure_engine((program,),program if program in ('diamond','mmseqs','vsearch') else 'blast'):continue
                if program not in ('diamond','mmseqs','vsearch') and not self.ensure_engine(('makeblastdb',),'blast'):continue
                if len(query_ids)>1:
                    self.s.pending=('browse',self.s.compare_many(query_ids,mode));return True
                if mode=='vsearch':
                    self.s.vsearch_run('search');self.s.dirty=True;return True
                previous=(self.s.mode,self.s.a.mode);previous_hsps=self.s.hsps
                try:
                    self.s.mode=self.s.a.mode=mode
                    self.s.check_type();self.s.search();self.s.dirty=True
                except BaseException:
                    if self.s.hsps is previous_hsps:self.s.mode,self.s.a.mode=previous
                    else:self.s.dirty=True
                    raise
                self.ui.notice='Search retained in Results. Inspect matches or choose another analysis.'
                return True

    def ensure_engine(self,aliases,profile):
        if any(self.s.tools.available(command) for command in aliases):return True
        from av_setup import setup_tools
        with self.ui.section('Missing tool'):
            while True:
                choice=self.choose('Install '+profile.upper(),[
                    item('1','Install in the AlignView environment','Use Conda, Mamba or Micromamba'),
                    item('2','See other options','HPC module, explicit path or manual installation'),
                ],hint='Required command: '+(' or '.join(aliases))+'. No calculation has started.')
                if choice is None:return False
                if choice=='2':
                    self.ui.page('HPC module: load the module provided by the cluster.\n\nLocal path: --tools-dir /path/bin\nExplicit executable: --executable '+aliases[0]+'=/path/'+aliases[0]+
                        '\n\nInstallation: alignview --setup '+profile+' --install-tools','Install the required engine')
                    continue
                setup_tools(self.ui,profile)
                if self.ui.yesno('Run this installation?',False):
                    setup_tools(self.ui,profile,install=True);self.s.tools.refresh_managed()
                    return any(self.s.tools.available(command) for command in aliases)

    def vsearch(self,action=None):
        from av_vsearch import ACTIONS
        with self.ui.section('VSEARCH'):
            if action is None:
                action=self.choose('Which VSEARCH analysis?',[
                    item(key,title,description,enabled=key!='search' or not self.s.single,
                         reason='Load references for a global search.')
                    for key,(title,description) in ACTIONS.items()],
                    hint='DNA/RNA only · global search distinct from BLAST · no E-values or bitscores')
            if action is None:return False
            while True:
                choice=self.choose('Prepare VSEARCH: '+action,[
                    item('1','Run the analysis',ACTIONS[action][0]),
                    item('2','VSEARCH parameters','Identity, strands, coverage, abundances and limits'),
                    item('3','Inspect commands and input preparation','Engine, paths, normalization and actual parameters'),
                    item('4','Tools and installation','Inspect VSEARCH and its version'),
                    item('5','Maximum number of references',str(self.s.a.hits),enabled=action=='search',
                         reason='The reported-reference limit applies to global searches.'),
                ],subtitle=self.s.tools.engine_info('vsearch'),
                    hint=('The search uses the active query: '+self.s.query.id if action=='search' else
                          f'All {len(self.s.query_list)} sequences from the first FASTA; gaps removed and U converted to T.'))
                if choice is None:return False
                if choice=='2':self.edit_parameters('vsearch',strategy=action)
                elif choice=='3':self.preview(lambda:self.s.vsearch_run(action,preview=True))
                elif choice=='4':self.software()
                elif choice=='5':self.s._change_numeric('hits','Maximum number of references',int,1)
                elif self.ensure_engine(('vsearch',),'vsearch'):
                    self.run_action(lambda:self.s.vsearch_run(action),changes=True)
                    return True

    def vsearch_results(self):
        from av_vsearch import ACTIONS
        if not self.s.vsearch_results:return
        labels={'input_sequences':'Input sequences','eligible_sequences':'Within length limits',
                'excluded_by_length':'Excluded by length','matches':'Matches',
                'matched_references':'Matching references','representatives':'Representatives',
                'total_abundance':'Total abundance','unique':'Unique sequences','chimeras':'Chimeras',
                'nonchimeras':'Non-chimeras','borderline':'Borderline'}
        while True:
            action=self.choose('VSEARCH results',[item(key,ACTIONS[key][0],
                ' · '.join(labels.get(name,name)+': '+str(value) for name,value in result['summary'].items()))
                for key,result in self.s.vsearch_results.items()])
            if action is None:return
            result=self.s.vsearch_results[action]
            lines=[ACTIONS[action][0],'Tool: VSEARCH','']
            lines.extend(labels.get(key,key)+': '+str(value) for key,value in result['summary'].items())
            lines+=['','Preview limited to 200 rows; complete files: vsearch/'+action+'/','']
            rows=result.get('preview',[])
            if rows:
                columns=list(rows[0]);lines.append(' | '.join(columns))
                lines.extend(' | '.join(str(row.get(key,'')) for key in columns) for row in rows)
            lines+=['','Settings and commands: Commands → Tools and execution history.']
            self.s._workflow_reviewed=True
            self.ui.page('\n'.join(lines),ACTIONS[action][0])


    def tree(self):
        if not self.s.msa_records or len(self.s.msa_records)<3:
            self.ui.page('Create or open an alignment with at least three sequences first.','Alignment required');return
        with self.ui.section('Exploratory tree'):
            while True:
                choice=self.choose('Create a tree',[
                    item('1','Infer the tree','FastTree · using the available alignment'),
                    item('2','Advanced options','Model, support, seed and command'),
                ],subtitle=self.s.tools.engine_info('FastTree' if self.s.tools.has('FastTree') else 'fasttree'),
                    hint='Approximate exploratory tree. Local support values are not bootstraps.')
                if choice is None:return
                if choice=='2':self.tree_options();continue
                if not self.ensure_engine(('FastTree','fasttree'),'fasttree'):continue
                self.s.make_tree();return

    def show_tree(self):
        from av_tree import tree_text
        self.ui.page(tree_text((self.s.work/'tree_fasttree.nwk').read_text()),'Exploratory tree')


    def data(self):
        with self.ui.section('Files and query'):
            while True:
                can_change=len(self.s.query_list)>1
                choice=self.choose('Files and query',[
                    item('1','Choose the active sequence' if self.s.single else 'Change the active query',
                         f'{len(self.s.query_list)} sequences ; targeted analyses use the active sequence',enabled=can_change,
                         reason='The loaded file contains only one sequence.'),
                    item('2','Open other files','Choose an alignment, sequences or references'),
                    item('3','Return to a previous analysis',f'{len(self.s.history)} analyses in history',enabled=bool(self.s.history),
                         reason='History records analyses when files or the active query change.'),
                    item('4','Open a saved project','Restore files, parameters and calculations'),
                    item('5','Browse sequences','Identifiers, length and FASTA content'),
                ],hint='Changing inputs keeps the current analysis in this session\'s history.')
                if choice is None:return
                if choice=='1':
                    record=choose_sequence(self.ui,self.s.query_list,'Choose the active sequence' if self.s.single else 'Choose the active query',self.s.query.id)
                    if record and record.id!=self.s.query.id:self.s.pending=('query',record.id);return
                elif choice=='2':
                    request=startup(self.ui)
                    if request:self.s.pending=('files',request);return
                elif choice=='3':
                    from av_project import read_project
                    choices=[]
                    for i,path in enumerate(self.s.history,1):
                        _,data=read_project(path)
                        choices.append(item(str(i),data['active_query_id'],data['mode']+' · '+data['started_utc']))
                    selected=self.choose('Previous analyses',choices)
                    if selected:self.s.pending=('resume',self.s.history[int(selected)-1]);return
                elif choice=='4':
                    folder=self.ui.ask('Project directory (q to cancel)')
                    if folder and folder.casefold()!='q':
                        from av_project import read_project
                        root,_=read_project(folder);self.s.pending=('resume',root);return
                elif choice=='5':self.sequences()

    def save(self):
        with self.ui.section('Save'):
            while True:
                choice=self.choose('Save and share',[
                    item('1','Save the project for later','Optional: retain every dataset and result after closing'),
                    item('2','Preview HTML report','Open a browser preview without saving or exporting'),
                    item('3','Export files','FASTA alignment, TSV matches, SVG tree or complete folder'),
                    item('4','Export HTML report','Keep a standalone report after closing'),
                ],hint='Saving is never required to navigate or run another analysis. Export destinations are relative to '+str(Path.cwd())+'.')
                if choice is None:return
                self.run_action({'1':self.s.save_checkpoint,'2':self.s.open_report,'3':self.export,'4':self.s.copy_report}[choice])

    def export(self):
        with self.ui.section('Export'):
            native_search=not self.state.searched and 'search' in self.s.vsearch_results
            choice=self.choose('What would you like to export?',[
                item('msa','FASTA alignment','Use the alignment in another program',enabled=bool(self.s.msa_records),reason='Create or open an alignment.'),
                item('hits','TSV global matches' if native_search else 'TSV matches',
                     'Native VSEARCH global matches, coordinates and aligned rows' if native_search else 'One row per reference, identity and coverage',
                     enabled=self.state.searched or native_search,reason='Run a comparison.'),
                item('tree','SVG tree','Vector figure to open or publish',enabled=(self.s.work/'tree_fasttree.nwk').is_file(),reason='Create a tree in Results.'),
                item('all','Results folder','Tables, alignments, figures, HTML and VSEARCH results'),
            ])
            if choice is None:return
            stem='alignview_'+re.sub('[^A-Za-z0-9_.-]','_',self.s.query.id)[:35]
            suffix={'msa':'_alignment.fasta','hits':'_matches.tsv','tree':'_tree.svg','all':'_exports'}[choice]
            destination=self.ui.ask('Destination (q to cancel)',stem+suffix)
            if not destination or destination.casefold()=='q':return
            target=Path(destination).expanduser().resolve()
            if target.exists():raise ValueError('The destination already exists: '+str(target))
            if choice!='all':self.s.export_artifact(choice,target)
            else:
                previous=self.s.a.export_profile
                try:
                    self.s.a.export_profile='publication';self.s.export_summary()
                    target.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copytree(self.s.work/'exports',target)
                finally:self.s.a.export_profile=previous
            self.ui.notice='Results exported to '+str(target)

    def settings(self):
        return self.advanced()


    def confirm_exit(self):
        if not getattr(self.s,'dirty',True):
            self.s.exit_decision='saved' if self.s.last_saved and Path(self.s.last_saved).is_dir() else 'clean'
            return True
        from av_library import entries
        if not self.s.last_saved and not entries(self.s,'result') and not self.s.calculation_runs:
            self.s.exit_decision='clean';return True
        choice=self.choose('Close workspace?',[
            item('1','Continue working','Return to the workspace'),
            item('2','Close without saving','Remove temporary data and results; saved copies remain intact'),
            item('3','Save a copy and close','Optional: retain every dataset and result for another session'),
        ])
        if choice=='3':
            if self.s.save_checkpoint():self.s.exit_decision='saved';return True
        elif choice=='2':self.s.exit_decision='discard';return True
        return False
