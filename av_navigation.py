"""Task-based navigation. Scientific calculations stay in Session and av_core."""
from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import dataclass, replace
import io
from pathlib import Path
import re
import shutil

from av_core import sequence_stats, sequence_type
from av_parameters import snapshot
from av_start import choose_file, choose_sequence, startup
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
        count=len(s.query_list) if s.single else 1+len(s.references)
        reason=''
        if count<2:reason='At least two sequences are required. Open another file in Files and query.'
        elif s.mode in ('tblastn','tblastx'):
            reason='Multiple translation frames are possible. Use detected segments in the advanced alignment options.'
        elif s.single and len({sequence_type(r.seq) for r in s.query_list})>1:
            reason='The file mixes sequence types. Separate nucleotide and protein inputs.'
        calculated=getattr(s,'search_parameters',{})
        current=snapshot('search',s.a,s.mode)
        stale=bool(calculated) and any(name in calculated and calculated[name]!=value
            for name,value in current.items() if name!='candidate_limit')
        return cls(len(s.query_list),len(s.references),len(s.selected),bool(s.msa_records),
                   bool(s.hsps) or getattr(s,'search_completed',False),not reason,reason,stale)


HOME_ITEMS = (
    ('1','Results','Read the alignment, matches and tree'),
    ('2','Compare with references','Choose the query and find its matches'),
    ('3','Align sequences','Prepare an alignment and open its viewer'),
    ('4','Save and share','Project, HTML report and result files'),
    ('5','Files and query','Switch sequences, open files or resume a project'),
    ('6','Advanced','Additional analyses, statistics, settings and tools'),
)


class Workbench(AdvancedActions):
    def __init__(self,session):
        self.s=session
        self.ui=session.ui

    @property
    def state(self):return WorkspaceState.from_session(self.s)

    def refresh_context(self):
        s=self.s;state=self.state
        source=(f'{s.source_query_path.name} · {state.sequence_count} sequences · active: {s.query.id}' if s.single
                else f'Query: {s.query.id} · {state.reference_count} references')
        results=[]
        if state.search_stale:results.append('Search needs rerunning: parameters changed')
        if state.searched:results.append(f'{state.hit_count} correspondances')
        if state.aligned:results.append(f'alignment ready ({len(s.msa_records[0].seq):,} columns)')
        if getattr(s,'vsearch_results',{}):results.append(f'VSEARCH: {len(s.vsearch_results)} analyses')
        run=getattr(s,'calculation_runs',{}).get('msa',{}).get('parameters',{})
        if run and run!=snapshot('msa',s.a,s.mode,run.get('strategy','')):
            results.append('alignment uses the previous parameters')
        if not results:results.append('no calculation started')
        saved='unsaved changes' if getattr(s,'dirty',True) else 'no changes'
        self.ui.context=(source,' · '.join([*results,saved]))

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
            if artifact and (self.s.work/artifact).is_file():
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
        self.ui.route=['Home']
        self.ui.positions.pop(('Home','Project'),None)
        while True:
            state=self.state
            ready=state.aligned or state.searched or bool(self.s.vsearch_results)
            hint=('Results ready: press Enter to open them.' if ready else
                  'Choose Compare with references or Align sequences.')
            try:
                entries=[item(*entry) for entry in HOME_ITEMS]
                entries[0]=replace(entries[0],enabled=ready,reason='No results: run a comparison or an alignment.')
                entries[2]=replace(entries[2],enabled=state.can_align or bool(self.s.selected),reason=state.align_reason)
                choice=self.choose('Project',entries,hint=hint,default='1' if ready else '2',
                    shortcuts={'v':'view','s':'save','c':'data','t':'software','h':'help','x':'exit'})
                if choice in (None,'exit'):
                    if self.confirm_exit():return
                    continue
                handlers={'1':self.explore,'2':self.compare,'3':self.align,'4':self.save,
                          '5':self.data,'6':self.settings,'view':self.view_alignment,'save':self.s.save_checkpoint,
                          'data':self.data,'software':self.software,'help':self.help}
                self.run_action(handlers[choice])
                if self.s.pending:return
            except HomeRequested:
                self.ui.positions[('Home','Project')]='1' if self.state.aligned or self.state.searched or self.s.vsearch_results else '2'
                self.ui.notice='Back to Home. Data and results are retained.'
            except KeyboardInterrupt:
                self.ui.notice='Action cancelled. The project remains open.'
            except (ValueError,CommandError,OSError) as exc:
                self.ui.error(str(exc));self.ui.page(str(exc),'Action not completed')

    def compare(self):
        if self.s.single:
            with self.ui.section('References'):
                reference=choose_file(self.ui,'Choose references')
                if reference:self.s.pending=('references',reference)
            return
        return self.search()

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
                    self.ui.page('Run a comparison or alignment from Home.','No results available');return
                choice=self.choose('Results',entries,hint='Detailed checks and statistics are in Advanced.')
                if choice is None:return
                self.run_action(actions[choice])

    def view_alignment(self):
        state=self.state
        if not state.aligned and not (self.s.single and len(self.s.query_list)>=2 and len({len(r.seq) for r in self.s.query_list})==1):
            self.ui.page('No alignment is available.\n\nChoose Align sequences.','Open an alignment')
            return
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
                    self.ui.page('No matches retained.\n\nRelax the filters in Advanced → Calculation settings → Search,\nor run a more sensitive search from Home.','Search completed without matches')
                    return
                options=[item(str(i),r['subject'],f'{r["best_identity"]:.1f}% local identity · {r["query_coverage"]:.1f}% of the query · {r["hsps"]} segments') for i,r in enumerate(rows,1)]
                choice=self.choose('Matches',options,f'{len(rows)} retained references · limite: {self.s.a.hits}',
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

    def align(self):
        with self.ui.section('Alignment'):
            strategy=('all' if self.s.single or self.s.mode.startswith('global-') else 'hits' if self.s.selected else 'all')
            if self.s.mode in ('tblastn','tblastx'):strategy='hsp'
            labels={'all':'All sequences','hits':'Query and retained references','fragments':'Add a partial sequence','hsp':'Detected segments only'}
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
                if mode=='vsearch':return self.vsearch(action='search')
                choice=self.choose('Compare with references',[
                    item('1','Run comparison',f'{mode.upper()} · query: {self.s.query.id} · {len(self.s.references)} references'),
                    item('2','Choose query',self.s.query.id,enabled=len(self.s.query_list)>1,
                         reason='The file contains a single query.'),
                    item('3','Search threshold',f'{self.s.a.preset} · E-value ≤ {self.s.a.evalue:g}'),
                    item('4','Maximum number of references',str(self.s.a.hits)),
                    item('5','Advanced options','Method, filters, translation, CPUs, settings and commands'),
                ],subtitle=self.s.tools.engine_info('diamond' if mode.startswith('diamond-') else 'mmseqs' if mode=='mmseqs-blastp' else mode),
                    hint='The search compares only the displayed query with the references.')
                if choice is None:return False
                if choice=='2':
                    record=choose_sequence(self.ui,self.s.query_list,'Choose the active query',self.s.query.id)
                    if record and record.id!=self.s.query.id:self.s.pending=('query',record.id);return False
                    continue
                if choice=='3':self.choose_preset();continue
                if choice=='4':self.s._change_numeric('hits','Maximum number of references',int,1);continue
                if choice=='5':mode=self.search_options(mode);continue
                program='diamond' if mode.startswith('diamond-') else 'mmseqs' if mode=='mmseqs-blastp' else mode
                if not self.ensure_engine((program,),program if program in ('diamond','mmseqs') else 'blast'):continue
                if program not in ('diamond','mmseqs') and not self.ensure_engine(('makeblastdb',),'blast'):continue
                previous=(self.s.mode,self.s.a.mode);previous_hsps=self.s.hsps
                try:
                    self.s.mode=self.s.a.mode=mode
                    self.s.check_type();self.s.search();self.s.dirty=True
                except BaseException:
                    if self.s.hsps is previous_hsps:self.s.mode,self.s.a.mode=previous
                    else:self.s.dirty=True
                    raise
                self.hits()
                return True

    def ensure_engine(self,aliases,profile):
        if any(self.s.tools.available(command) for command in aliases):return True
        from av_setup import setup_tools
        with self.ui.section('Missing tool'):
            while True:
                choice=self.choose('Installer '+profile.upper(),[
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
                    self.vsearch_results();return True

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
            lines+=['','Settings and commands: Advanced → Tools and provenance.']
            self.ui.page('\n'.join(lines),ACTIONS[action][0])


    def tree(self):
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
                self.s.make_tree();self.show_tree();return

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
                    item('1','Save the project','Resume data, settings and all calculations'),
                    item('2','Create an HTML report','A standalone file to open in a browser'),
                    item('3','Export files','FASTA alignment, TSV matches, SVG tree or complete folder'),
                ],hint='Destinations are relative to '+str(Path.cwd())+'. Existing files are never overwritten.')
                if choice is None:return
                self.run_action({'1':self.s.save_checkpoint,'2':self.s.open_report,'3':self.export}[choice])

    def export(self):
        with self.ui.section('Export'):
            choice=self.choose('What would you like to export?',[
                item('msa','FASTA alignment','Use the alignment in another program',enabled=bool(self.s.msa_records),reason='Create or open an alignment.'),
                item('hits','TSV matches','One row per reference, identity and coverage',enabled=self.state.searched,reason='Run a comparison.'),
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
        choice=self.choose('Quitter AlignView ?',[
            item('1','Save and quit','Choose a destination and retain the complete project'),
            item('2','Quit without saving','Discard changes; previously saved copies remain intact'),
            item('3','Continue working','Return to the project'),
        ])
        if choice=='1':
            if self.s.save_checkpoint():self.s.exit_decision='saved';return True
        elif choice=='2':self.s.exit_decision='discard';return True
        return False
