"""Session data and scientific operations, independent of entry-point orchestration."""
from __future__ import annotations
import argparse
import collections
import datetime
import json
import math
import os
import re
import shlex
import shutil
import sys
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from av_core import Seq, HSP, read_fasta, write_fasta, sequence_type, sequence_stats, merge_intervals, interval_length, select_hits, coverage_data, translate_nt, reverse_complement, orfs, find_motifs, sliding_scores, msa_metrics, distance_matrix, dotplot, sha256, export_table, normalize_search_record, validate_search_results, validate_alignment_results
from av_tools import Toolchain, CommandError
from av_ui import UI, BOLD, GREEN, CYAN, YELLOW, RED, GREY, HomeRequested
from av_visual import dotplot_svg, plot_svg
from av_report import build_report
from av_tree import TREE_FIGURE, tree_svg

from av_units import units_for
from av_exports import generate_exports
from av_start import suggest_modes
from av_version import VERSION
from av_parameters import PARAMETERS, validate as validate_parameters, snapshot as calculation_settings, extra_cli
from av_project import save_project
from av_setup import diagnostic
from av_msa_tui import MSAViewer, msa_static

from av_cli import HELP_TEXT

def show_file(file:Path,ui:UI,max_lines=160):
    if not file.exists():ui.warn('File not found');return
    ui.page(file.read_text(errors='replace'),file.name)


def text_segments(spans):return ', '.join(f'{a:,}–{b:,} ({b-a+1:,})' for a,b in spans) if spans else 'none'


def ascii_bar(spans,length,width=62,color=True):
    if length<=0:return ''
    out=[]
    for i in range(width):
        a=i*length/width+1;b=(i+1)*length/width
        covered=any(s<=b and e>=a for s,e in spans)
        out.append('█' if covered else '·')
    return ''.join((GREEN+x if x=='█' else GREY+x)+'\033[0m' for x in out) if color else ''.join(out)


class Session:
    def __init__(self,args,ui):
        self.a=args;self.ui=ui;self.exit_decision=None;self.pending=None;self.history=[];self._viewer=None
        self.reference_selection=[];self._alignment_aliases={};self._expected_alignment_records=None
        self.sequence_selection=[];self.workspace_entry={};self.navigation_state={};self.library_root=None
        self.comparison_selection=list(getattr(args,'_comparison_queries',[]));self.comparison_batch={}
        self._comparison_needs_selection=False
        self.result_sources={}
        self.data_origin=getattr(args,'_data_origin','')
        self.dirty=True;self.last_saved=None;self.search_completed=False;self.search_parameters={};self.calculation_runs={}
        self.query_path=Path(args.query).expanduser().resolve()
        self.refs_path=Path(args.references).expanduser().resolve() if args.references else None
        for path in (self.query_path,self.refs_path):
            if path and not path.is_file():raise ValueError(f'File not found: {path}')
        self.source_query_path=self.query_path;self.source_refs_path=self.refs_path
        self.temporary=tempfile.TemporaryDirectory(prefix='alignview_')
        self.work=Path(self.temporary.name)/'results';self.work.mkdir()
        # Immutable input copies keep a session independent from source-file edits.
        for name in ('query','refs'):
            source=getattr(self,name+'_path')
            if source:
                folder=Path(self.temporary.name)/'input_sources'/name;folder.mkdir(parents=True)
                copied=folder/source.name;shutil.copyfile(source,copied)
                setattr(self,name+'_path',copied)
        self.input_hashes={'query':sha256(self.query_path),'refs':sha256(self.refs_path) if self.refs_path else ''}
        self.query_list=read_fasta(self.query_path)
        self.references=read_fasta(self.refs_path) if self.refs_path else []
        self.single=not self.refs_path
        self.query=self.query_list[0]
        self.refmap={r.id:r for r in self.references}
        self.query_search_path=self.query_path
        if args.query_id or (not self.single and len(self.query_list)>1):
            candidate=next((r for r in self.query_list if r.id==args.query_id),None) if args.query_id else None
            if args.query_id and candidate is None:
                raise ValueError(f'Unknown query identifier: {args.query_id}')
            if candidate is None and (args.no_menu or not ui.tty):
                raise ValueError('Multiple query sequences: provide --query-id in non-interactive mode.')
            if candidate is None:
                if getattr(args,'_defer_query_choice',False):
                    candidate=self.query_list[0];self._comparison_needs_selection=True
                else:
                    labels=[(str(i),r.id,f'{len(r.seq):,} characters · {sequence_type(r.seq)}','SEQUENCES') for i,r in enumerate(self.query_list,1)]
                    choice=ui.menu('Choose the active query',labels,
                                   f'{len(self.query_list)} sequences · analyses focus on one query')
                    if choice is None:raise KeyboardInterrupt
                    candidate=self.query_list[int(choice)-1]
            self.query=candidate
            self.query_search_path=self.work/'query_active.fasta'
            write_fasta([self.query],self.query_search_path)
            if not self._comparison_needs_selection:ui.info(f'Selected query: {self.query.id} · {len(self.query.seq):,} characters')
        self.input_preparation={'query_gap_characters_removed':0,'reference_gap_characters_removed':0,
            'coordinate_system':'1-based positions in prepared ungapped sequences'}
        if not self.single:
            raw_query=self.query;raw_references=self.references
            self.query=normalize_search_record(raw_query)
            self.references=[normalize_search_record(r) for r in raw_references]
            self.refmap={r.id:r for r in self.references}
            self.input_preparation.update(
                query_gap_characters_removed=len(raw_query.seq)-len(self.query.seq),
                reference_gap_characters_removed=sum(len(a.seq)-len(b.seq) for a,b in zip(raw_references,self.references)),
                normalization='Remove -/.; replace ? by N/X; nucleotide U becomes T')
            if self.input_preparation['query_gap_characters_removed'] or self.input_preparation['reference_gap_characters_removed']:
                ui.info('Comparison inputs prepared without gaps. Coordinates refer to these ungapped sequences.')
        # Plain prepared FASTA supports compressed source inputs for every engine.
        self.query_search_path=self.work/'query_active.fasta'
        write_fasta([self.query] if not self.single else self.query_list,self.query_search_path)
        self.refs_search_path=self.work/'references_active.fasta' if self.refs_path else None
        if self.refs_search_path:write_fasta(self.references,self.refs_search_path)
        self.tools=Toolchain(self.work,ui,auto_modules=args.module_policy=='auto',executable_paths=args.executable,tool_dirs=args.tools_dir)
        self.hsps=[];self.selected={};self.msa_records=None;self.msa_file=None;self.vsearch_results={}
        self.dot_data=None;self.motif_rows=[];self.orf_rows=[]
        self.created=[];self.mode=args.mode or ('msa' if self.single else '')
        command=getattr(args,'_command','')
        if not self.single and not args.mode and command in ('compare','align') and not self._comparison_needs_selection:
            qkind=sequence_type(self.query.seq);rkind=collections.Counter(sequence_type(r.seq) for r in self.references).most_common(1)[0][0]
            suggested,global_mode,_=suggest_modes(qkind,rkind)
            chosen=global_mode if command=='align' else suggested
            if not chosen:raise ValueError('Cannot infer a compatible method. Choose --mode explicitly, or prepare matching sequence alphabets.')
            self.mode=args.mode=chosen
        intent=getattr(args,'_workflow_goal','') or getattr(args,'_intent','') or command
        self.workflow_goal=intent if intent in ('compare','align','view','sequence') else ('compare' if not self.single else 'view' if args.msa=='view' else '')
        self._workflow_reviewed=False
        self.start=datetime.datetime.now(datetime.timezone.utc).isoformat()
        if not self._comparison_needs_selection:self.check_type()
        self.save_settings()
        self.dirty=False

    def close(self):
        self.temporary.cleanup()

    def check_type(self):
        qtype=sequence_type(self.query.seq)
        if self.single:return
        ref_type=collections.Counter(sequence_type(r.seq) for r in self.references).most_common(1)[0][0]
        types={
            'vsearch':('DNA/RNA','DNA/RNA'),
            'blastn':('DNA/RNA','DNA/RNA'),'blastp':('protein','protein'),
            'blastx':('DNA/RNA','protein'),'mmseqs-blastp':('protein','protein'),'tblastn':('protein','DNA/RNA'),
            'tblastx':('DNA/RNA','DNA/RNA'),'diamond-blastp':('protein','protein'),
            'diamond-blastx':('DNA/RNA','protein'),'global-dna':('DNA/RNA','DNA/RNA'),
            'global-aa':('protein','protein')}
        if self.mode in types:
            want=types[self.mode]
            if (qtype,ref_type)!=want:
                self.ui.warn(f'Detected FASTA types: QUERY={qtype}, REFS≈{ref_type} ; expected for {self.mode}: {want[0]}/{want[1]}.')
                if self.mode.startswith('global-') and qtype!='unknown' and ref_type!='unknown':
                    raise ValueError('Global MAFFT: incompatible DNA/protein types; use homogeneous FASTA input or translate the sequences.')

    def settings(self):
        keys=('mode','task','evalue','identity','coverage','bitscore','hits','threads','sensitivity',
              'preset','msa','mafft','aligner','tree','frame','gapopen','gapextend','matrix','free_ends','kmer')
        d={k:getattr(self.a,k) for k in keys};d['mode']=self.mode
        d.update({p.name:getattr(self.a,p.name,p.default) for p in PARAMETERS})
        d['native_defaults']='Options not supplied use defaults from the engine version recorded in commands.json.'
        d['calculation_runs']=self.calculation_runs
        d['vsearch_analyses']={action:{k:value for k,value in result.items() if k!='preview'} for action,result in self.vsearch_results.items()}
        u=units_for(self.mode,sequence_type(self.query.seq))
        d['units']={'query_coordinates':u.query,'reference_coordinates':u.subject,'alignment_columns':u.aligned,'msa_columns':u.msa}
        d['timestamp_utc']=self.start;d['query']=str(self.query_path);d['active_query_id']=self.query.id
        d['references']=str(self.refs_path) if self.refs_path else None
        d['search_completed']=self.search_completed
        d['search_run']=dict(self.search_parameters)
        d['reference_selection']=list(self.reference_selection)
        d['input_preparation']=dict(self.input_preparation)
        d['workflow_goal']=self.workflow_goal
        d['original_query']=str(self.source_query_path);d['original_references']=str(self.source_refs_path) if self.source_refs_path else None
        return d

    def save_settings(self):
        self.dirty=True
        (self.work/'parameters.json').write_text(json.dumps(self.settings(),ensure_ascii=False,indent=2)+'\n')
        (self.work/'inputs_sha256.tsv').write_text('file\tsha256\n'+
             '\n'.join(f'{p}\t{self.input_hashes[name]}' for name,p in [('query',self.query_path),('refs',self.refs_path)] if p)+'\n')

    def export_hsps(self):
        columns=('query','subject','identity','align_length','qstart','qend','sstart','send',
                 'evalue','bitscore','qframe','sframe','qseq','sseq','matches','mismatches','gap_columns')
        all_records=[h.record() for h in self.hsps]
        selected=[h.record() for group in self.selected.values() for h in group]
        export_table(self.work/'hsps_all.tsv',all_records,columns)
        export_table(self.work/'hsps_selected.tsv',selected,columns)
        export_table(self.work/'coverage.tsv',[
            {**r,'query_segments':str(r['query_segments']),'query_uncovered':str(r['query_uncovered']),
             'reference_segments':str(r['reference_segments']),'reference_uncovered':str(r['reference_uncovered'])}
            for r in self.coverage_rows()],sep='\t')
        (self.work/'hsps.json').write_text(json.dumps(selected,ensure_ascii=False,indent=2))
        with (self.work/'hsp_query.bed').open('w') as f:
            f.write('# Query BED 0-based half-open; BlastX query coordinates = nucleotides\n')
            for h in [h for g in self.selected.values() for h in g]:
                a,b=h.qspan; f.write(f'{h.query}\t{a-1}\t{b}\t{h.subject}\t{min(1000,round(h.bitscore))}\t{"-" if h.qstart>h.qend else "+"}\n')
        with (self.work/'hsp_reference.bed').open('w') as f:
            f.write('# Reference BED 0-based half-open; TBLASTN/TBLASTX reference coordinates = nucleotides\n')
            for h in [h for g in self.selected.values() for h in g]:
                a,b=h.sspan;f.write(f'{h.subject}\t{a-1}\t{b}\t{h.query}\t{min(1000,round(h.bitscore))}\t{"-" if h.sstart>h.send else "+"}\n')

    def coverage_rows(self):
        if self.single:return []
        return coverage_data([h for group in self.selected.values() for h in group],
                             len(self.query.seq),{r.id:len(r.seq) for r in self.references})

    def ref_select(self):
        if self.single:return None,0
        labels=list(self.selected) if self.selected else [r.id for r in self.references]
        if self.a.no_menu or not self.ui.tty:
            n=self.a.ref
            if not 1<=n<=len(labels):raise ValueError('Reference index out of range: '+str(n))
            return self.refmap[labels[n-1]],n
        options=[]
        for i,name in enumerate(labels,1):
            group=self.selected.get(name,[])
            detail=f'{group[0].identity:.1f}% · E={group[0].evalue:.2g}' if group else f'{len(self.refmap[name].seq):,} characters'
            options.append((str(i),name,detail,'REFERENCES'))
        choice=self.ui.menu('CHOOSE A REFERENCE',options,f'{len(labels):,} references · / search')
        if choice is None:return None,0
        n=int(choice);self.a.ref=n
        return self.refmap[labels[n-1]],n

    def show_dashboard(self):
        self.ui.banner('SESSION / settings and data in memory')
        self.ui.heading('SESSION STATUS')
        for k,v in [
            ('Query',f'{self.query_path.name} · {len(self.query_list)} seq. · {len(self.query.seq):,} positions'),
            ('References',f'{self.refs_path.name} · {len(self.references):,} seq.' if self.refs_path else 'none — single FASTA'),
            ('Program',self.mode.upper()),('Thresholds',f'E ≤ {self.a.evalue:g} · id ≥ {self.a.identity}% · coverage QUERY ≥ {self.a.coverage}%'),
            ('Search',f'{len(self.hsps)} HSP · {len(self.selected)} retained references'),
            ('Alignment',self.msa_file.name if self.msa_file else 'not generated'),
            ('Dotplot','generated' if self.dot_data else 'not generated'),
            ('Launch directory',str(Path.cwd())),
            ('Saving','temporary session · not saved' if not self.last_saved else 'unsaved changes' if self.dirty else 'project saved')]:
            print('  '+self.ui.c(GREY,f'{str(k):<23}')+str(v))
        self.ui.info('Loaded modules: '+(', '.join(m for _,m in self.tools.loaded) or str('none; using PATH')))
        if self.last_saved:self.ui.info('Last saved copy: '+str(self.last_saved))


    def _change_numeric(self,key,label,cast,minimum=None,maximum=None):
        previous=getattr(self.a,key)
        answer=self.ui.ask(label,str(previous),'Enter = keep · q = cancel')
        if answer.lower() in ('q','cancel'):return False
        try:
            value=cast(answer)
            if not math.isfinite(float(value)):raise ValueError
            if minimum is not None and value<minimum:raise ValueError
            if maximum is not None and value>maximum:raise ValueError
        except ValueError:
            self.ui.warn('Invalid value: parameter unchanged.');return False
        setattr(self.a,key,value);self.save_settings()
        self.ui.success(f'{label}: {previous} → {value}')
        return value!=previous

    def choose_preset(self):
        from av_navigation import Workbench
        return Workbench(self).choose_preset()

    def settings_hub(self,initial=False):
        from av_navigation import Workbench
        return Workbench(self).search_settings()


    def search(self):
        if self.mode.startswith('global-') or self.single:return
        self.ui.heading('SEARCH / EXECUTION')
        validate_parameters(self.a,context='search',mode=self.mode)
        max_target=self.a.search_candidates or max(100,self.a.hits*10)
        self.ui.info(f'{self.mode.upper()} · {self.a.threads} threads')
        self.ui.info(f'Query: {self.query.id} · {self.source_query_path}')
        self.ui.info(f'Target FASTA: {self.source_refs_path} · {len(self.references)} sequences')
        self._calculation_start=len(self.tools.history)
        candidates=self.tools.search(self.mode,self.query_search_path,self.refs_search_path,self.a.threads,self.a.evalue,
                                    max_target,self.a.sensitivity,self.a.task,options=self.a)
        prepared=validate_search_results(candidates,self.query,self.references,self.mode)
        self.hsps=prepared
        self.comparison_batch={}
        self.search_completed=True
        self.search_parameters=calculation_settings('search',self.a,self.mode)
        self.record_calculation('search',self.search_parameters,preparation=dict(self.input_preparation))
        self.apply_filters()
        from av_library import capture
        capture(self,'search',f'Matches · {self.query.id} → {self.source_refs_path.name}')
        self.ui.success(f'{len(self.hsps)} HSPs found · {len(self.selected)} retained references')

    def search_preview(self,mode=None):
        mode=mode or self.mode
        commands=self.tools.search_commands(mode,self.query_search_path,self.refs_search_path,self.a.threads,
            self.a.evalue,self.a.search_candidates or max(100,self.a.hits*10),self.a.sensitivity,self.a.task,self.a)
        return self.tools.describe(commands)

    def compare_many(self,query_ids,mode):
        """Run independent searches and retain a grouped, resumable summary."""
        query_ids=list(query_ids)
        known={record.id:record for record in self.query_list}
        if (len(query_ids)<2 or any(not isinstance(identifier,str) or identifier not in known for identifier in query_ids)
                or len(set(query_ids))!=len(query_ids)):
            raise ValueError('Choose at least two distinct query identifiers from this FASTA.')
        if len({sequence_type(known[identifier].seq) for identifier in query_ids})!=1:
            raise ValueError('Compare nucleotide and protein queries in separate selections.')
        from av_cli import MODES
        if mode not in MODES or mode.startswith('global-'):raise ValueError('Choose a reference-search method for grouped comparisons.')
        if self.single:raise ValueError('Choose a target FASTA before running comparisons.')
        if self.library_root is None:raise ValueError('Grouped comparisons require an interactive workspace.')
        def prepared(identifier):
            args=argparse.Namespace(**vars(self.a))
            args.query=str(self.query_path);args.references=str(self.refs_path)
            args.query_id=identifier;args.mode=mode;args.msa='none';args.tree=False
            args._comparison_queries=[identifier];args.resume=None
            result=Session(args,self.ui)
            result.source_query_path=self.source_query_path;result.source_refs_path=self.source_refs_path
            result.library_root=self.library_root;result.history=self.history
            result.tools.env.update(self.tools.env)
            return result
        group=prepared(query_ids[0]);runs=[];history_count=len(self.history)
        try:
            group.comparison_selection=query_ids
            for index,identifier in enumerate(query_ids,1):
                self.ui.stage(f'Comparing {index}/{len(query_ids)} · {identifier}')
                leaf=None
                attempt=group.work/'attempts'/f'{index:04d}';attempt.mkdir(parents=True)
                try:
                    leaf=prepared(identifier)
                    if mode=='vsearch':
                        leaf.vsearch_run('search')
                        matches=leaf.vsearch_results['search']['summary'].get('matches',0)
                    else:
                        leaf.search();matches=len(leaf.selected)
                    runs.append(dict(query_id=identifier,status='completed',result_id=leaf.workspace_entry['id'],matches=matches))
                except KeyboardInterrupt:
                    runs.append(dict(query_id=identifier,status='canceled',error='Comparison interrupted.'))
                    break
                except (CommandError,ValueError,OSError) as exc:
                    runs.append(dict(query_id=identifier,status='failed',error=str(exc)))
                    self.ui.warn(f'{identifier}: {exc}')
                finally:
                    if leaf:
                        (attempt/'commands.json').write_text(json.dumps(leaf.tools.history,ensure_ascii=False,indent=2),encoding='utf-8')
                        if (leaf.work/'logs').is_dir():shutil.copytree(leaf.work/'logs',attempt/'logs')
                        leaf.close()
            for identifier in query_ids[len(runs):]:runs.append(dict(query_id=identifier,status='not_run'))
            group.comparison_batch=dict(query_ids=query_ids,reference_count=len(self.references),mode=mode,runs=runs)
            (group.work/'comparison_batch.json').write_text(json.dumps(group.comparison_batch,ensure_ascii=False,indent=2),encoding='utf-8')
            group._calculation_start=0
            settings=calculation_settings('vsearch' if mode=='vsearch' else 'search',self.a,mode)
            group.record_calculation('search_batch',settings,query_ids=query_ids,
                completed=sum(run['status']=='completed' for run in runs),failed=sum(run['status']=='failed' for run in runs))
            return Path(self.history[-1])
        finally:
            if len(self.history)!=history_count:self.dirty=True
            group.close()

    def record_calculation(self,kind,parameters,**details):
        self.calculation_runs[kind]=dict(parameters=parameters,**details,
            commands=[dict(row) for row in self.tools.history[self._calculation_start:]])
        self.save_settings()
        from av_library import capture_calculation
        capture_calculation(self,kind,self.calculation_runs[kind])

    def invalidate_derived(self,clear_msa=False):
        """Remove only generated, now stale outputs in this owned temporary workspace."""
        names=['report.html','coverage.svg','profile.svg','figures_catalog.tsv','figures','exports','tree_fasttree.nwk']
        self.calculation_runs.pop('tree',None)
        if clear_msa:
            names+=['dotplot.svg','dotplot.json','msa_positions.tsv','msa_pdistance.tsv','msa_consensus.fasta','msa_conservation.svg','alignment_active.fasta']
            self.msa_records=None;self.msa_file=None;self.dot_data=None;self._viewer=None
            self.calculation_runs.pop('msa',None)
        for name in names:
            path=self.work/name
            if path.is_dir():shutil.rmtree(path)
            elif path.is_file():path.unlink()

    def apply_filters(self,preserve_msa=False):
        selected=select_hits(self.hsps,len(self.query.seq),self.a.identity,self.a.coverage,
                                  self.a.hits,self.a.bitscore)
        if selected!=self.selected and not preserve_msa:self.invalidate_derived(clear_msa=True)
        self.selected=selected
        self.export_hsps();self.save_settings()

    def show_hits(self):
        if self.single or self.mode.startswith('global-'):
            self.ui.warn('No local search: choose BLAST / DIAMOND to obtain matches.');return
        rows=self.coverage_rows();self.ui.heading(f'RESULTS  /  {len(rows)} REFERENCES')
        print(f'  {"#":>2}  {"Reference":<30} {"Id.%":>7} {"Cov.%":>7} {"HSP":>5} {"E-value":>11} {"Bits":>8}')
        for i,r in enumerate(rows,1):
            print(f'  {i:>2}  {r["subject"][:30]:<30} {r["best_identity"]:7.1f} {r["query_coverage"]:7.1f} {r["hsps"]:5d} {r["best_evalue"]:11.2g} {r["best_bitscore"]:8.1f}')
        if not rows:self.ui.warn('No hits after filtering. Relax thresholds or change the method.')
        units=units_for(self.mode,sequence_type(self.query.seq))
        self.ui.info(f'Identity = % per HSP; coverage = union on QUERY ({units.query}) ; HSP alignment ({units.aligned}).')

    def show_coverage(self,by_ref=False):
        if not self.selected:self.ui.warn('No retained HSPs; coverage unavailable.');return
        rows=self.coverage_rows();q_len=len(self.query.seq)
        self.ui.heading('COVERAGE MAPS / REFERENCES' if by_ref else 'COVERAGE MAP / QUERY')
        if not by_ref:
            all_spans=merge_intervals([h.qspan for grp in self.selected.values() for h in grp],q_len)
            pct=100*interval_length(all_spans)/q_len
            print(f'  {self.ui.c(GREEN,"UNION TOTAL".ljust(24))}  {ascii_bar(all_spans,q_len,57,self.ui.color)}  {pct:5.1f}%')
        for r in rows:
            segments=r['reference_segments'] if by_ref else r['query_segments']
            length=r['reference_length'] if by_ref else q_len
            pct=r['reference_coverage'] if by_ref else r['query_coverage']
            print(f'  {r["subject"][:24]:24}  {ascii_bar(segments,length,57,self.ui.color)}  {pct or 0:5.1f}%')
            miss=r['reference_uncovered'] if by_ref else r['query_uncovered']
            print('     '+self.ui.c(GREY,'No HSP: '+text_segments(miss)))
        units=units_for(self.mode,sequence_type(self.query.seq))
        self.ui.info(f'█ HSP · gray = no HSP · QUERY in {units.query} · REFERENCE in {units.subject} · aligned columns in {units.aligned}')
        if 'blastx' in self.mode or self.mode=='tblastx':
            self.ui.warn('BLASTX/TBLASTX: QUERY uses nucleotides; internal alignments use amino acids.')
        if by_ref and self.mode in ('tblastn','tblastx'):
            self.ui.warn('TBLASTN/TBLASTX: reference coordinates use nucleotides.')

    def inspect_hit(self):
        if not self.selected:self.ui.warn('No HSPs to explore.');return
        self.ui.heading('INSPECTOR / REFERENCE AND HSP')
        ref,_=self.ref_select()
        if ref is None:return
        rows=self.selected.get(ref.id,[])
        print(); self.ui.info(f'{ref.id} · target length {len(ref.seq):,} · {len(rows)} HSPs')
        for i,h in enumerate(rows,1):
            print(f'  {self.ui.c(CYAN,str(i))}  Q {h.qstart}–{h.qend}  /  S {h.sstart}–{h.send} '
                  f'· {h.identity:.1f}% identity · E={h.evalue:.2g} · bits={h.bitscore:.1f} · frame Q {h.qframe}/S {h.sframe}')
        if not rows:return
        val=self.ui.ask('HSP number (0 = all)','0')
        try:ix=int(val)
        except ValueError:ix=0
        if ix>len(rows) or ix<0:self.ui.warn('Invalid segment number.');return
        import contextlib,io
        output=io.StringIO()
        with contextlib.redirect_stdout(output):
            for h in rows if ix==0 else [rows[ix-1]]:self.print_alignment(h)
        self.ui.page(output.getvalue(),f'Alignments HSP / {ref.id}')

    def print_alignment(self,h:HSP,width=64):
        self.ui.heading(f'HSP  /  {h.subject}  ·  {h.identity:.1f}% identity')
        units=units_for(self.mode,sequence_type(self.query.seq))
        self.ui.info(f'QUERY {h.query}: {h.qstart} → {h.qend} {units.query} · TARGET {h.subject}: {h.sstart} → {h.send} {units.subject}')
        self.ui.info(f'{h.align_length} alignment columns ({units.aligned}) · E={h.evalue:g}')
        self.ui.info(f'{h.matches} identities · {h.mismatches} substitutions · {h.gap_columns} columns with gaps · frame Q={h.qframe}, S={h.sframe}')
        qa=h.qseq;sa=h.sseq
        if not qa or not sa:self.ui.warn('The HSP alignment has no qseq/sseq.');return
        nt_q=self.mode in ('blastx','diamond-blastx','tblastx')
        nt_s=self.mode in ('tblastn','tblastx')
        qstep=3 if nt_q else 1;sstep=3 if nt_s else 1
        qsign=1 if h.qend>=h.qstart else -1;ssign=1 if h.send>=h.sstart else -1
        qp=h.qstart;sp=h.sstart
        for i in range(0,min(len(qa),len(sa)),width):
            q=qa[i:i+width];s=sa[i:i+width]
            qn=sum(c not in '-.' for c in q);sn=sum(c not in '-.' for c in s)
            endq=qp+qsign*(qn*qstep-1) if qn else qp
            ends=sp+ssign*(sn*sstep-1) if sn else sp
            marks=''.join(' ' if x in '-.' or y in '-.' else '|' if x.upper()==y.upper() else '.' for x,y in zip(q,s))
            if self.ui.color:
                colors=[RED if c in '-.' or t in '-.' else GREEN if c.upper()==t.upper() else YELLOW for c,t in zip(q,s)]
                cq=''.join(self.ui.c(color,c) for color,c in zip(colors,q))
                cs=''.join(self.ui.c(color,c) for color,c in zip(colors,s))
            else:cq=q;cs=s
            print(f'  Q {str(qp) if qn else "—":>8}  {cq}  {endq if qn else "—"}')
            print(f'  {"":>11}  {marks}')
            print(f'  S {str(sp) if sn else "—":>8}  {cs}  {ends if sn else "—"}\n')
            qp+=qsign*qn*qstep;sp+=ssign*sn*sstep
        self.ui.info('| identity · . substitution · blank gap. Coordinates are 1-based inclusive in the prepared ungapped query and target.')
        if self.ui.color:self.ui.info('Green = identity, yellow = substitution, red = gap.')

    def hsp_diagnostics(self):
        self.ui.heading('DIAGNOSTICS / SUBSTITUTIONS, GAPS AND FRAMES')
        for ref,items in self.selected.items():
            groups=collections.Counter(h.qframe for h in items)
            spans=merge_intervals([h.qspan for h in items],len(self.query.seq))
            total_cols=sum(len(h.qseq) for h in items)
            gaps=sum(h.gap_columns for h in items);subs=sum(h.mismatches for h in items)
            self.ui.info(f'{ref} · {len(items)} HSP · {interval_length(spans)} unique QUERY positions')
            print(f'      {total_cols} cumulative columns (HSPs can overlap), {subs} substitutions, {gaps} gap columns ; frames {dict(groups)}')
        self.ui.warn('Summing substitutions across overlapping HSPs can count positions more than once.')

    def identity_coverage(self):
        self.ui.heading('PLOT / LOCAL IDENTITY VS QUERY COVERAGE')
        rows=self.coverage_rows()
        if not rows:self.ui.warn('No hits');return
        W=56;H=16;board=[[' ' for _ in range(W+1)] for _ in range(H+1)]
        for i in range(H+1):board[i][0]='│'
        for j in range(W+1):board[H][j]='─'
        board[H][0]='└'
        for i,r in enumerate(rows,1):
            x=min(W,max(1,round(r['query_coverage']/100*W)))
            y=min(H-1,max(0,H-1-round(r['best_identity']/100*(H-1))))
            board[y][x]=str(i%10)
        for row in board:print('     '+''.join(row))
        print('     0% coverage QUERY'+' '*(W-23)+'100%')
        for i,r in enumerate(rows,1):print(f'    {i:>2}. {r["subject"]}: {r["best_identity"]:.1f}% identity HSP / {r["query_coverage"]:.1f}% coverage')
        self.ui.warn('These metrics do not represent whole-sequence percent identity.')

    def run_pairwise(self,algorithm,reference=None,preview=False):
        if self.single:self.ui.warn('Provide a second FASTA for Needle/Water.');return
        self.ui.heading(f'{algorithm.upper()} / PAIRWISE ALIGNMENT')
        ref=reference
        if ref is None:ref,_=self.ref_select()
        if ref is None:return
        rows=self.selected.get(ref.id,[])
        h=rows[0] if rows else None
        qframe=self.a.frame or (h.qframe if h else 1)
        sframe=self.a.sframe or (h.sframe if h and h.sframe else 1)
        if self.mode in ('blastx','diamond-blastx','tblastx','tblastn'):
            self.ui.warn(f'Full-sequence translation: Q frame {qframe}, REF frame {sframe} where applicable.')
        prefix=self.work/f'{algorithm}_{re.sub("[^A-Za-z0-9_.-]","_",ref.id)}'
        self._calculation_start=len(self.tools.history)
        result=self.tools.emboss(algorithm,self.query,ref,prefix,mode=self.mode,
                        gapopen=self.a.gapopen,gapextend=self.a.gapextend,
                        free_ends=self.a.free_ends,matrix=self.a.matrix,frame=qframe,sframe=sframe,options=self.a,preview=preview)
        if preview:
            translated=self.mode in ('blastx','diamond-blastx','tblastx','tblastn')
            return (f'Sequences: {self.query.id} / {ref.id}.\n'+
                (f'Full translation: QUERY frame {qframe}, table {self.a.genetic_code} ; reference frame {sframe}, table {self.a.reference_genetic_code or self.a.genetic_code} ; stop codons retained.\n' if translated else '')+result)
        file,_=result
        self.a.pairwise='global' if algorithm=='needle' else 'local'
        self.record_calculation(algorithm,calculation_settings(algorithm,self.a,self.mode),
            reference=ref.id,query_frame=qframe,reference_frame=sframe,genetic_code=self.a.genetic_code,
            reference_genetic_code=self.a.reference_genetic_code or self.a.genetic_code,artifact=file.name)
        self.dirty=True
        if self.library_root is None:show_file(file,self.ui,max_lines=180)
        self.created.append(str(file.name))

    def msa_run(self,strategy,preview=False):
        engine='mafft' if strategy=='fragments' else self.a.aligner
        effective=argparse.Namespace(**vars(self.a));effective.aligner=engine
        validate_parameters(effective,context='msa')
        self.ui.heading(f'{engine.upper()} / MULTIPLE ALIGNMENT')
        self._calculation_start=len(self.tools.history)
        self._msa_preparation={}
        self._alignment_aliases={};self._expected_alignment_records=None
        pending=self.work/'pending_msa'
        if pending.exists():shutil.rmtree(pending)
        pending.mkdir()
        if self.single:
            selected=[r for r in self.query_list if r.id in self.sequence_selection] if strategy=='selection' else self.query_list
            if len(selected)<2:self.ui.warn('An MSA requires at least two sequences.');return
            selected=[normalize_search_record(r) for r in selected]
            input_file=self.work/'msa_prepared.fasta';write_fasta(selected,input_file)
            self._expected_alignment_records=selected
            self._msa_preparation['normalization']='Remove source gaps; ? becomes N/X; nucleotide U becomes T'
            self._msa_preparation['record_ids']=[r.id for r in selected]
        else:
            if strategy=='hsp' or (strategy=='selection' and self.mode in ('tblastn','tblastx')):
                return self._msa_hsps(preview=preview,subjects=self.reference_selection if strategy=='selection' else None)
            if self.mode in ('tblastn','tblastx'):
                self.ui.warn('Full MAFFT is unavailable for TBLASTN/TBLASTX: multiple frames/ORFs are possible. Extract protein ORFs first.');return
            if strategy=='selection':
                if not self.reference_selection:raise ValueError('Select references in Matches before aligning the selection.')
                selected=[self.refmap[name] for name in self.reference_selection]
            else:
                selected=self.references if strategy=='all' or self.mode.startswith('global-') or (strategy=='fragments' and not self.selected) else [self.refmap[k] for k in self.selected]
            if not selected:self.ui.warn('No retained references; use all or relax thresholds.');return
            if self.mode in ('blastx','diamond-blastx'):
                frames={h.qframe for v in self.selected.values() for h in v if h.qframe}
                frame=self.a.frame or (next(iter(self.selected.values()))[0].qframe if self.selected else None)
                if frame is None:
                    frame=1
                    self.ui.warn('No frame detected or forced: explicitly using QUERY frame +1 (change with --frame).')
                if frame not in (-3,-2,-1,1,2,3):self.ui.warn('Invalid frame');return
                if len(frames)>1:self.ui.warn(f'Different frames {sorted(frames)}: only QUERY frame {frame} is used.')
                tr=translate_nt(self.query.seq,frame,self.a.genetic_code).replace('*','X')
                self._msa_preparation.update(query_frame=frame,genetic_code=self.a.genetic_code,stop_handling='stop codons replaced with X')
                query_record=Seq('QUERY',f'QUERY_frame_{frame}',tr)
                self.ui.info('QUERY DNA translated to protein (stop codons replaced with X for MAFFT).')
            else:query_record=Seq('QUERY','QUERY',self.query.seq)
            ref_records=[Seq('REF_'+str(i),f'REF_{i} {r.id}',r.seq) for i,r in enumerate(selected,1)]
            targets=[self.query,*selected];counts=collections.Counter(r.id for r in targets)
            reserved={r.id for r in targets if counts[r.id]==1};used=set()
            def display_id(record,role):
                name=record.id
                if counts[name]>1:
                    name=role+':'+name
                    while name in reserved or name in used:name=role+':'+name
                used.add(name);return name
            self._alignment_aliases={'QUERY':(display_id(self.query,'query'),self.query.description)}
            for alias,record in zip(ref_records,selected):
                self._alignment_aliases[alias.id]=(display_id(record,'ref'),record.description)
            self._msa_preparation['reference_alignment_ids']={r.id:self._alignment_aliases[a.id][0] for a,r in zip(ref_records,selected)}
            self._expected_alignment_records=[query_record,*ref_records]
            self._msa_preparation['reference_ids']=[r.id for r in selected]
            self._msa_preparation['sequence_aliases']={key:value[0] for key,value in self._alignment_aliases.items()}
            input_file=self.work/'msa_input.fasta'
            if strategy=='fragments':
                commands=[]
                if len(ref_records)==1:
                    ref_aln=pending/'msa_reference_alignment.fasta';write_fasta(ref_records,ref_aln)
                else:
                    refs_file=self.work/'msa_refs.fasta';write_fasta(ref_records,refs_file)
                    ref_aln=pending/'msa_reference_alignment.fasta'
                    commands.append(self.tools.mafft_command(refs_file,ref_aln,self.a.mafft,self.a.threads,options=self.a))
                frag=self.work/'msa_fragment.fasta';write_fasta([query_record],frag)
                out=pending/'alignment_fragments.fasta'
                commands.append(self.tools.mafft_command(ref_aln,out,'auto',self.a.threads,fragments=frag,options=self.a))
                if preview:return self.msa_preparation_text()+'\nAlign references, then add QUERY; addfragments requires MAFFT.\n'+self.tools.describe(commands)
                self.tools.execute(commands)
                return self._set_msa(out,'fragments')
            write_fasta([query_record,*ref_records],input_file)
        if self.a.aligner=='mafft' and self.a.mafft in ('global','local') and len(selected)>200:
            self.ui.warn('G-INS-i/L-INS-i is expensive above 200 sequences: use --mafft auto.');return
        out=pending/('alignment_'+strategy+'.fasta')
        if preview:return self.msa_preparation_text()+'\n'+self.tools.describe([self.tools.msa_command(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)])
        self.tools.msa(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)
        return self._set_msa(out,strategy)

    def msa_preparation_text(self):
        if self._msa_preparation.get('query_frame'):
            return f'Preparation: QUERY translated in frame {self._msa_preparation["query_frame"]}, NCBI table {self.a.genetic_code}, stop codons replaced with X.'
        return 'Preparation: full-length sequences; no additional translation.'

    def _msa_hsps(self,preview=False,subjects=None):
        groups=self.selected
        if subjects is not None:
            if not subjects or any(name not in groups for name in subjects):
                raise ValueError('Each checked reference requires retained local HSPs for translated alignment.')
            groups={name:groups[name] for name in subjects}
        if not groups:self.ui.warn('No HSPs for HSP alignment.');return
        anchor=next(iter(groups.values()))[0]
        self._msa_preparation=dict(anchor_subject=anchor.subject,anchor_query_span=anchor.qspan,
            hsp_gap_handling='Gaps removed before realignment',stop_handling='stop codons replaced with X',
            inherited_search_parameters=dict(self.search_parameters))
        anchor_intervals=anchor.qspan
        query_id='QUERY_HSP'
        while query_id in self.selected:query_id='query:'+query_id
        kept=[Seq(query_id,f'{query_id} {self.query.id}:{anchor.qstart}-{anchor.qend}',anchor.qseq.replace('-','').replace('*','X'))]
        manifest=[]
        for name,group in groups.items():
            choices=[]
            for h in group:
                a=max(anchor_intervals[0],h.qspan[0]);b=min(anchor_intervals[1],h.qspan[1])
                ov=max(0,b-a+1)/max(1,anchor_intervals[1]-anchor_intervals[0]+1)
                choices.append((ov,h.bitscore,h))
            overlap,_,best=max(choices,key=lambda it:(it[0],it[1]))
            accepted=overlap>=self.a.hsp_overlap
            manifest.append(dict(subject=name,anchor_overlap_pct=round(100*overlap,2),included=accepted,
                                 qstart=best.qstart,qend=best.qend))
            if accepted:kept.append(Seq(name,name,best.sseq.replace('-','').replace('*','X')))
        if not preview:export_table(self.work/'msa_hsp_selection.tsv',manifest)
        if len(kept)<2:self.ui.warn(f'No HSP covers at least {self.a.hsp_overlap:.0%} of the anchor region.');return
        self.ui.info(f'{len(kept)-1} compatible HSPs, {len(manifest)-(len(kept)-1)} excluded.')
        input_file=self.work/'msa_hsp_input.fasta';out=self.work/'pending_msa'/'alignment_hsp.fasta'
        self._alignment_aliases={f'HSP_{i}':(r.id,r.description) for i,r in enumerate(kept,1)}
        prepared=[Seq(alias,alias,r.seq) for alias,r in zip(self._alignment_aliases,kept)]
        self._expected_alignment_records=prepared
        self._msa_preparation['sequence_aliases']={'QUERY':query_id,**{key:value[0] for key,value in self._alignment_aliases.items()}}
        self._msa_preparation['reference_alignment_ids']={r.id:r.id for r in kept[1:]}
        self._msa_preparation['reference_ids']=[r.id for r in kept[1:]]
        self._msa_preparation['checked_reference_scope']=list(subjects) if subjects is not None else None
        write_fasta(prepared,input_file)
        if preview:return f'Preparation: anchor {anchor.subject} Q {anchor.qstart}–{anchor.qend}; overlap threshold {self.a.hsp_overlap:.0%} ; {len(kept)} sequences ; gaps removed and stop codons converted to X.\n'+self.tools.describe([self.tools.msa_command(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)])
        self.tools.msa(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)
        self.ui.warn('Exploratory MSA: HSPs can represent non-equivalent domains.')
        return self._set_msa(out,'hsp')

    def _set_msa(self,out,strategy):
        records=read_fasta(out)
        expected=self._expected_alignment_records
        if expected:
            try:
                records=validate_alignment_results(records,expected,
                    nucleotide=all(sequence_type(r.seq)=='DNA/RNA' for r in expected))
            except ValueError as exc:
                raise CommandError('Invalid alignment output: '+str(exc)) from exc
        elif len({len(record.seq) for record in records})!=1:
            raise CommandError('The engine produced sequences with different lengths.')
        aliases=self._alignment_aliases
        if aliases:
            records=[Seq(aliases[r.id][0],aliases[r.id][1],r.seq) if r.id in aliases else r for r in records]
        # Validate completely before replacing the previous successful alignment.
        write_fasta(records,out)
        target=self.work/out.name
        if out!=target:out.replace(target)
        out=target
        self.dirty=True
        self.invalidate_derived()
        self.msa_file=out;self.msa_records=records
        self.a.msa=strategy
        self._viewer=None
        self.created.append(out.name)
        self.record_calculation('msa',calculation_settings('msa',self.a,self.mode,strategy),
            preparation=getattr(self,'_msa_preparation',{}),output_sha256=sha256(out))
        engine='mafft' if strategy=='fragments' else self.a.aligner
        self.ui.success(f'MSA {engine}/{strategy} · {len(self.msa_records)} sequences · {len(self.msa_records[0].seq)} columns')
        if not self.ui.tty or self.a.no_menu:self.show_msa_window()
        return out

    def load_msa(self):
        if self.msa_records:return self.msa_records
        if self.single and (len(self.query_list)>=2 or self.a.msa=='view') and len({len(r.seq) for r in self.query_list})==1:
            # Equal length is not proof of alignment. --msa view is an explicit user assertion.
            if self.a.msa!='view' and self.ui.tty and not self.a.no_menu:
                self.ui.warn('All sequences have the same length; this does not prove they are aligned.')
                if not self.ui.yesno('Is this FASTA already aligned?',False):
                    self.ui.info('Choose Align sequences to align this file.')
                    return None
            self.msa_records=self.query_list;self.msa_file=self.query_path
            from av_library import capture
            capture(self,'msa','Imported alignment · '+self.source_query_path.name)
            return self.msa_records
        self.ui.warn('No MSA available; align sequences first.')
        return None

    def show_msa_window(self):
        """Full-screen viewer: two-way scrolling, variable sites and search."""
        seqs=self.load_msa()
        if not seqs:return
        ncol=len(seqs[0].seq)
        if any(len(r.seq)!=ncol for r in seqs):
            self.ui.warn('Different lengths: FASTA is not aligned. Run an aligner first.');return
        if self.a.no_menu or not self.ui.tty:
            print(str(msa_static(seqs,width=75,limit=40,color=False)))
            return
        self._workflow_reviewed=True
        if self._viewer is None or self._viewer.seqs is not seqs:
            self._viewer=MSAViewer(seqs,color=self.ui.color)
            self._viewer.row=next((i for i,record in enumerate(seqs) if record.id==getattr(self,'_viewer_focus_id',self.query.id)),0)
            state=self.navigation_state.get(self.workspace_entry.get('id'),{}).get('viewer',{})
            for name in ('row','col','xoff','yoff','only_variable','row_filter','trace_id','difference','region'):
                if name in state:setattr(self._viewer,name,state[name])
            self._viewer.checked=set(state.get('checked',[]))
        if getattr(self,'_trace_id',None):
            self._viewer.trace_id=self._trace_id;self._viewer.difference=True
            self._trace_id=None
        self._viewer.workflow_lines=[]
        self._viewer.result=None
        result=self._viewer.run()
        if result=='home':raise HomeRequested
        if result=='web':
            self.make_report(open_preview=True)
        elif result=='export':
            self.copy_report()
        elif result=='alignment_export':self.export_alignment_selection(self._viewer)
        elif result=='save':self.save_checkpoint()

    def msa_analysis(self):
        seqs=self.load_msa()
        if not seqs:return
        metric=msa_metrics(seqs)
        self.ui.heading('CONSERVATION / VARIABLE SITES')
        self.ui.info(f'{metric["n"]} sequences, {metric["length"]:,} columns, {metric["variable"]:,} variables')
        self.ui.info(f'{metric["parsimony_informative"]:,} parsimony-informative, {metric["gapped"]:,} with gaps, {metric["invariant"]:,} strictly invariant')
        export_table(self.work/'msa_positions.tsv',metric['positions'])
        write_fasta([Seq('CONSENSUS','CONSENSUS',metric['consensus'])],self.work/'msa_consensus.fasta')
        scores=[dict(center=r['position'],value=r['conservation']) for r in metric['positions']]
        (self.work/'msa_conservation.svg').write_text(plot_svg(scores,'MSA conservation (%)'))
        for r in metric['positions'][:20]:
            print(f'  {r["position"]:>5}  {r["consensus"]}  conserv. {r["conservation"]:5.1f}% · gap {r["gap_pct"]:5.1f}% · entropy {r["entropy"]:.2f}')
        if len(metric['positions'])>20:self.ui.info('All columns are in msa_positions.tsv.')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('conservation',{},engine='AlignView Python',version=VERSION)

    def msa_distances(self):
        seqs=self.load_msa()
        if not seqs:return
        matrix=distance_matrix(seqs);labels=[s.id for s in seqs]
        with (self.work/'msa_pdistance.tsv').open('w') as f:
            f.write('sequence\t'+'\t'.join(labels)+'\n')
            for name,row in zip(labels,matrix):f.write(name+'\t'+'\t'.join(f'{v:.6g}' if v is not None else 'NA' for v in row)+'\n')
        self.ui.heading('P-DISTANCES / GAPS EXCLUDED')
        self.ui.info('Row and column numbers refer to these sequences :')
        for i,name in enumerate(labels,1):print(f'  {i:>3}  {name}')
        print('\n       '+' '.join(f'{i:>7}' for i in range(1,len(labels)+1)))
        for i,row in enumerate(matrix,1):
            print(f'  {i:>3}  '+' '.join('     NA' if value is None else f'{value:7.3f}' for value in row))
        self.ui.info('Uncorrected p-distance; gapped columns are excluded for each pair. Maximum 250 sequences.')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('distances',{},engine='AlignView Python',version=VERSION)

    def vsearch_run(self,action,preview=False):
        from av_vsearch import commands, summarize, publish_directory
        source=[self.query] if action=='search' else self.query_list
        references=self.references if action=='search' else []
        if action=='search' and not references:raise ValueError('VSEARCH search: load a reference FASTA.')
        if any(sequence_type(record.seq)!='DNA/RNA' for record in [*source,*references]):
            raise ValueError('VSEARCH requires nucleotide sequences; proteins are not supported.')
        def normalized(records):
            return [Seq(r.id,r.id,r.seq.upper().replace('-','').replace('.','').replace('U','T')) for r in records]
        normalized_source=normalized(source)
        eligible=sum(self.a.vsearch_min_length<=len(r.seq)<=self.a.vsearch_max_length for r in normalized_source)
        if not eligible:raise ValueError('VSEARCH: no sequences meet the selected length limits.')
        pending=self.work/('pending_vsearch_'+action)
        if pending.exists():shutil.rmtree(pending)
        pending.mkdir()
        input_file=pending/'input.fasta';write_fasta(normalized_source,input_file)
        refs=pending/'references.fasta' if references else None
        if refs:write_fasta(normalized(references),refs)
        prepared=commands(action,input_file,refs,pending,self.a)
        preparation=f'{len(source)} sequences · gaps removed and U converted to T · {eligible} within the length limits.'
        if preview:return preparation+'\n\n'+self.tools.describe(prepared)
        self.ui.info(preparation)
        if action=='search':
            self.ui.info(f'Query: {self.query.id} · {self.source_query_path}')
            self.ui.info(f'Target FASTA: {self.source_refs_path} · {len(self.references)} sequences')
        self._calculation_start=len(self.tools.history)
        self.tools.execute(prepared)
        summary,rows=summarize(action,pending,len(source),eligible)
        settings=calculation_settings('vsearch',self.a,'vsearch' if action=='search' else self.mode,strategy=action)
        settings['action']=action
        result={'action':action,'engine':'VSEARCH','settings':settings,'summary':summary,'preview':rows,
                'preparation':{'gaps_removed':True,'rna_to_dna':True,'query_id':self.query.id if action=='search' else None},
                'commands':[dict(row) for row in self.tools.history[self._calculation_start:]]}
        (pending/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        parent=self.work/'vsearch';parent.mkdir(exist_ok=True)
        target=parent/action
        publish_directory(pending,target)
        if action=='search':
            # Publish the new comparison only after execution and validation
            # succeed. Earlier local attempts remain in the workspace catalog.
            self.mode=self.a.mode='vsearch'
            self.hsps=[];self.selected={};self.search_completed=False
            self.search_parameters={};self.reference_selection=[];self.comparison_batch={}
            self.calculation_runs.pop('search',None)
            self.invalidate_derived(clear_msa=True)
            for kind in ('search','msa','tree'):self.result_sources.pop(kind,None)
            self.export_hsps()
        self.vsearch_results[action]=result
        (parent/'analyses.json').write_text(json.dumps(self.vsearch_results,ensure_ascii=False,indent=2),encoding='utf-8')
        self.record_calculation('vsearch_'+action,settings,summary=summary)
        self.dirty=True
        self.ui.success('VSEARCH analysis completed: '+action)
        for key,value in summary.items():self.ui.info(f'{key}: {value}')
        self.ui.info('Files: vsearch/'+action+'/ · results included in the next HTML report.')
        return result

    def make_tree(self,preview=False):
        """FastTree on an existing MSA; exploratory inference without molecular dating."""
        if not self.msa_records or not self.msa_file:
            self.ui.warn('Create an MSA before running FastTree.');return
        if len(self.msa_records)<3:
            self.ui.warn('FastTree requires at least 3 sequences.');return
        nucleotide=all(sequence_type(r.seq)=='DNA/RNA' for r in self.msa_records)
        pending=self.work/'pending_tree';pending.mkdir(exist_ok=True)
        output=pending/'tree_fasttree.nwk'
        tree_input=pending/'alignment_safe_ids.fasta'
        aliases={f'TREE_{i}':record.id for i,record in enumerate(self.msa_records,1)}
        write_fasta([Seq(alias,alias,record.seq) for alias,record in zip(aliases,self.msa_records)],tree_input)
        if preview:return self.tools.describe([self.tools.fasttree_command(tree_input,output,nucleotide,self.a)])
        self._calculation_start=len(self.tools.history)
        self.tools.fasttree(tree_input,output,nucleotide,options=self.a)
        # Native Newick labels are restricted; restore arbitrary user IDs safely.
        from av_tree import relabel_newick
        restored=relabel_newick(output.read_text(encoding='utf-8'),aliases)
        figure=tree_svg(restored)
        output.write_text(restored,encoding='utf-8')
        target=self.work/output.name;output.replace(target);output=target
        figure_dir=self.work/'figures';figure_dir.mkdir(exist_ok=True)
        (figure_dir/(TREE_FIGURE+'.svg')).write_text(figure,encoding='utf-8')
        self.a.tree=True
        self.record_calculation('tree',calculation_settings('tree',self.a,self.mode),
            model=(('JC' if nucleotide else 'JTT') if self.a.fasttree_model=='auto' else self.a.fasttree_model),
            alignment_sha256=sha256(self.msa_file),output_sha256=sha256(output),sequence_aliases=aliases)
        self.dirty=True
        self.created.append(output.name)
        self.ui.success('Exploratory tree created: tree_fasttree.nwk (Newick)')
        self.ui.info('Figure SVG: figures/'+TREE_FIGURE+'.svg · included in the next HTML report.')
        self.ui.warn('FastTree provides an approximate phylogeny; check the model, support and MSA quality before publication.')

    def select_aligner(self):
        from av_navigation import Workbench
        return Workbench(self).choose_aligner()

    def show_stats(self):
        self.ui.heading('FASTA QUALITY CHECK')
        seqs=self.query_list if self.single else [self.query,*self.references]
        rows=[sequence_stats(r) for r in seqs]
        export_table(self.work/'fasta_statistics.tsv',rows)
        print(f'  {"Sequence":<27} {"Type":<11} {"Length":>10} {"GC %":>7} {"N %":>6} {"H(bits)":>8}')
        for r in rows:
            print(f'  {r["id"][:27]:<27} {r["kind"]:<11} {r["length"]:>10,} {r["gc_pct"] if r["gc_pct"] is not None else "—":>7} {r["n_pct"]:>6} {r["entropy"]:>8}')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('statistics',{},engine='AlignView Python',version=VERSION)

    def motif(self):
        self.ui.heading('IUPAC MOTIF SEARCH / BOTH STRANDS')
        seqs=[self.query]
        motif=self.ui.ask('IUPAC motif (e.g. ATG, NNNGTT)',self.a.motif or 'ATG')
        if not motif:return
        data=[]
        for r in seqs:data.extend(find_motifs(r,motif,True))
        self.a.motif=motif
        self.motif_rows=data;export_table(self.work/'motifs.tsv',data,['seqid','start','end','strand','matched'])
        self._calculation_start=len(self.tools.history)
        self.record_calculation('motif',{'motif':motif,'query_id':self.query.id},engine='AlignView Python',version=VERSION)
        self.ui.info(f'{len(data)} occurrence(s) · motif {motif}')
        for r in data[:50]:print(f'  {r["seqid"][:23]:23} {r["start"]:>8}–{r["end"]:<8} {r["strand"]}  {r["matched"]}')
        if len(data)>50:self.ui.warn('Display limited to 50 occurrences; full export is in motifs.tsv.')

    def find_orfs(self):
        self.ui.heading('ORF FINDER / SIX FRAMES')
        if sequence_type(self.query.seq)!='DNA/RNA':self.ui.warn('Nucleotide QUERY required.');return
        n=self.ui.ask('Minimum length in amino acids','60',f'ORF: NCBI table {self.a.genetic_code} · starts {self.a.orf_starts}')
        try:min_aa=int(n)
        except ValueError:self.ui.warn('Invalid length');return
        if min_aa<1:self.ui.warn('A positive length is required');return
        rows=orfs(self.query,min_aa,include_partial=False,table=self.a.genetic_code,starts=self.a.orf_starts);self.orf_rows=rows
        export_table(self.work/'orfs.tsv',rows,['seqid','start','end','frame','strand','length_aa','complete'])
        write_fasta([Seq(f'ORF_{i}',f'ORF_{i} {r["seqid"]}:{r["start"]}-{r["end"]} frame={r["frame"]}',r['protein'])
                     for i,r in enumerate(rows,1)],self.work/'orfs.faa')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('orfs',calculation_settings('translation',self.a,self.mode),
            engine='AlignView Python',version=VERSION,minimum_aa=min_aa,output_sha256=sha256(self.work/'orfs.faa'))
        self.ui.info(f'{len(rows)} complete ORFs (≥{min_aa} aa)')
        for r in rows[:35]:print(f'  {r["start"]:>9}–{r["end"]:<9}  frame {r["frame"]:>2}  {r["length_aa"]:>5} aa  {r["strand"]}')
        self.ui.warn('This simple detector does not replace Prodigal/GeneMark or structural annotation.')

    def translate_frames(self):
        self.ui.heading(f'SIX-FRAME TRANSLATION / NCBI TABLE {self.a.genetic_code}')
        if sequence_type(self.query.seq)!='DNA/RNA':self.ui.warn('Nucleotide QUERY required.');return
        seqs=[]
        for frame in (1,2,3,-1,-2,-3):
            aa=translate_nt(self.query.seq,frame,self.a.genetic_code)
            seqs.append(Seq(f'frame_{frame}',f'QUERY frame {frame} table={self.a.genetic_code}',aa))
            print(f'  {frame:>2}  {len(aa)} aa  {aa.count("*")} STOP  {aa[:95]}'+('…' if len(aa)>95 else ''))
        write_fasta(seqs,self.work/'six_frames.faa')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('translation',calculation_settings('translation',self.a,self.mode),
            engine='AlignView Python',version=VERSION,output_sha256=sha256(self.work/'six_frames.faa'))
        self.ui.info('Translation available in Results. Use as data to continue with these proteins.')

    def extract_region(self):
        self.ui.heading('EXTRACT A REGION / QUERY')
        try:
            start=int(self.ui.ask('Start (1-based)','1'))
            end=int(self.ui.ask('End (inclusive)',str(min(300,len(self.query.seq)))))
        except ValueError:self.ui.warn('Invalid coordinates');return
        if not 1<=start<=end<=len(self.query.seq):self.ui.warn(f'Coordinates outside range 1–{len(self.query.seq)}');return
        s=self.query.seq[start-1:end];rc=self.ui.yesno('Reverse complement?')
        if rc:
            if sequence_type(self.query.seq)!='DNA/RNA':self.ui.warn('Reverse complement requires nucleic acids');return
            s=reverse_complement(s)
        out=self.work/f'query_{start}_{end}.fasta'
        write_fasta([Seq(f'{self.query.id}_{start}_{end}',f'{self.query.id}_{start}_{end}'+('_RC' if rc else ''),s)],out)
        self._calculation_start=len(self.tools.history)
        self.record_calculation('extract',{'start':start,'end':end,'reverse_complement':rc,'query_id':self.query.id},
            engine='AlignView Python',version=VERSION,artifact=out.name)
        self.ui.success(f'{len(s)} characters extracted: {out.name}')
        for i in range(0,min(200,len(s)),70):print('  '+s[i:i+70])

    def profiles(self):
        self.ui.heading('SLIDING-WINDOW PROFILES')
        metric=self.ui.ask('Metric: gc, entropy or gap','gc' if sequence_type(self.query.seq)=='DNA/RNA' else 'entropy')
        try:
            w=int(self.ui.ask('Window',str(min(101,len(self.query.seq)))))
            st=int(self.ui.ask('Step',str(max(1,w//4))))
            if metric=='gc' and sequence_type(self.query.seq)!='DNA/RNA':self.ui.warn('GC does not apply to proteins.');return
            data=sliding_scores(self.query.seq,w,st,metric)
        except ValueError as e:self.ui.warn(str(e));return
        export_table(self.work/f'profile_{metric}.tsv',data)
        (self.work/f'profile_{metric}.svg').write_text(plot_svg(data,f'Profile {metric} (window {w}, step {st})'))
        self._calculation_start=len(self.tools.history)
        self.record_calculation('profile',{'metric':metric,'window':w,'step':st},artifact=f'profile_{metric}.tsv')
        numeric=[r['value'] for r in data if r['value'] is not None]
        self.ui.info(f'{len(data)} windows · mean {sum(numeric)/len(numeric):.2f}' if numeric else 'No data')
        for r in data[:30]:print(f'  {r["start"]:>7}–{r["end"]:<7}  {r["value"]}')
        self.ui.info(f'Chart: profile_{metric}.svg (retained after saving)')

    def composition(self):
        self.ui.heading('COMPOSITION / QUERY K-MERS')
        k=self.ui.ask('k (1–8)','3' if sequence_type(self.query.seq)=='DNA/RNA' else '2')
        try:k=int(k)
        except ValueError:self.ui.warn('Invalid k');return
        if not 1<=k<=8:self.ui.warn('k must be between 1 and 8');return
        c=collections.Counter(self.query.seq[i:i+k] for i in range(len(self.query.seq)-k+1)
                              if '-' not in self.query.seq[i:i+k] and '?' not in self.query.seq[i:i+k])
        rows=[dict(kmer=key,count=count,frequency=count/max(1,sum(c.values()))) for key,count in c.most_common()]
        export_table(self.work/f'kmers_{k}.tsv',rows)
        self._calculation_start=len(self.tools.history)
        self.record_calculation('composition',{'k':k},artifact=f'kmers_{k}.tsv')
        n=sum(c.values());self.ui.info(f'{len(c)} distinct words of length {k} across {n:,} windows')
        for row in rows[:30]:
            bar='▰'*min(30,round(30*row['count']/max(1,rows[0]['count'])))
            print(f'  {row["kmer"]:<12} {row["count"]:>8,}  {100*row["frequency"]:>6.2f}%  {self.ui.c(GREEN,bar)}')

    def _compute_dotplot(self,ref,k):
        """Dotplots compare sequences using the same alphabet only.

        For BLASTX/TBLASTN/TBLASTX, project nucleotide FASTA input into amino acids
        using the highest-scoring HSP. Axes then use translated amino acids,
        rather than original DNA positions.
        """
        u=units_for(self.mode,sequence_type(self.query.seq))
        q=self.query.seq;r=ref.seq;frame_q=None;frame_s=None
        hs=self.selected.get(ref.id,[])
        best=hs[0] if hs else None
        if u.translated_query:
            frame_q=self.a.frame or (best.qframe if best else None)
            if not frame_q:raise ValueError('Translated dotplot requires an HSP QUERY frame or explicit --frame.')
            q=translate_nt(q,frame_q,self.a.genetic_code).replace('*','X')
        if u.translated_subject:
            frame_s=best.sframe if best else None
            if not frame_s:raise ValueError('Translated dotplot requires a known reference frame from an HSP.')
            r=translate_nt(r,frame_s,self.a.reference_genetic_code or self.a.genetic_code).replace('*','X')
        dtypeq='aa' if u.translated_query else u.query
        dtypes='aa' if u.translated_subject else u.subject
        if dtypeq!=dtypes:raise ValueError('Dotplot: incompatible alphabets; translate sequences before comparison.')
        data=dotplot(q,r,k)
        data.update(xunit=dtypeq,yunit=dtypes,qframe=frame_q,sframe=frame_s,
                    note='Dotplot of translated amino-acid sequences; Q/S frames '+str(frame_q)+'/'+str(frame_s) if frame_q or frame_s else 'Dotplot in '+dtypeq)
        return data

    def make_dotplot(self):
        if self.single and len(self.query_list)<2:self.ui.warn('Dotplot requires at least two sequences.');return
        self.ui.heading('DOTPLOT / EXACT WORD MATCHES')
        if self.single:
            from av_start import choose_sequence
            candidates=[r for r in self.query_list if r.id!=self.query.id]
            ref=choose_sequence(self.ui,candidates,'Choose a sequence for the dotplot')
        else:ref,_=self.ref_select()
        if ref is None:return
        default_k='3' if (self.mode in ('blastx','tblastn','tblastx','blastp','diamond-blastp','diamond-blastx','mmseqs-blastp') or sequence_type(self.query.seq)=='protein') else str(self.a.kmer)
        k_str=self.ui.ask('Exact-match word length k',default_k,'DNA: 8–12; proteins: 3–5. Axis units are shown in the figure.')
        try:k=int(k_str)
        except ValueError:self.ui.warn('Invalid k');return
        if not 2<=k<=60:self.ui.warn('k must be between 2 and 60');return
        self.dot_data=self._compute_dotplot(ref,k)
        (self.work/'dotplot.svg').write_text(dotplot_svg(self.dot_data,f'{self.query.id} ({self.dot_data["xunit"]})',f'{ref.id} ({self.dot_data["yunit"]})'),encoding='utf-8')
        (self.work/'dotplot.json').write_text(json.dumps(self.dot_data))
        self._calculation_start=len(self.tools.history)
        self.record_calculation('dotplot',{'k':k,'query_id':self.query.id,'reference_id':ref.id},
            engine='AlignView Python',version=VERSION)
        self.ui.success(f'{len(self.dot_data["points"]):,} matches (k={k}) ; dotplot.svg')
        self.ui.info(self.dot_data['note'])
        if self.dot_data['truncated']:self.ui.warn('Point limit reached; increase k.')

    def make_report(self, open_preview=True):
        """Build a fully standalone HTML report without a temporary server."""
        self.ui.heading('REPORT / BROWSER EXPLORATION')
        refs=[r for r in self.query_list if r.id!=self.query.id] if self.single else self.references
        if self.workspace_entry.get('analysis')=='search_batch':
            from av_report import build_comparison_summary
            report=build_comparison_summary(self.work,self.comparison_batch,self.settings())
        else:report=build_report(self.work,self.query,refs,
                [h for group in self.selected.values() for h in group],self.mode,self.settings(),
                self.tools.history,self.msa_records,self.dot_data,
                sha_query=sha256(self.query_path),sha_refs=sha256(self.refs_path) if self.refs_path else '',
                motifs=self.motif_rows,orf_rows=self.orf_rows)
        self.ui.success('Report created · standalone HTML file · no network connection required.')
        if open_preview and self.ui.tty and not self.a.no_menu:
            self._open_local_html(report,temporary=True)
            self.ui.info('Session preview: no export required. Use Export HTML to retain a copy after closing.')
        else:
            self.ui.info('Temporary HTML prepared; retain it with --save, --html-out or saving on exit.')
        return report

    def _html_target(self):
        safe=re.sub(r'[^A-Za-z0-9_.-]+','_',self.query.id).strip('._')[:45] or 'query'
        return (Path.cwd()/f'alignview_{safe}_report.html').resolve()

    def export_html(self,destination):
        """Persist one self-contained HTML with no implicit project save."""
        target=Path(destination).expanduser().resolve()
        if target.suffix.lower()!='.html':raise ValueError('The .html extension is required.')
        if target.exists():raise ValueError('The HTML already exists: '+str(target)+'. Choose another name.')
        report=self.make_report(open_preview=False)
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as output:output.write(report.read_bytes())
        self.ui.success('Standalone HTML created: '+str(target))
        return target

    def resume_actions(self):
        """Only explicit actions run after resuming; saved action flags stay inert."""
        explicit=getattr(self.a,'_explicit',set())
        if 'vsearch' in explicit and self.a.vsearch!='none':self.vsearch_run(self.a.vsearch)
        if 'msa' in explicit and self.a.msa!='none':
            if self.single and self.a.msa in ('hits','hsp','fragments'):
                raise ValueError('This MSA strategy requires a reference FASTA.')
            if self.a.msa=='view':self.load_msa()
            else:
                strategy='all' if self.single else 'hits' if self.a.msa in ('align','full','global') else self.a.msa
                self.msa_run(strategy)
        if 'pairwise' in explicit and self.a.pairwise!='none':
            self.run_pairwise('needle' if self.a.pairwise=='global' else 'water')
        if 'dotplot' in explicit and self.a.dotplot:self._dot_batch()
        if 'motif' in explicit and self.a.motif:
            rows=[]
            for record in (self.query_list if self.single else [self.query]):
                rows.extend(find_motifs(record,self.a.motif,True))
            self.motif_rows=rows;export_table(self.work/'motifs.tsv',rows,['seqid','start','end','strand','matched'])
            self.dirty=True
        if 'tree' in explicit and self.a.tree:
            if not self.msa_records:raise ValueError('Create an MSA before --tree: add --msa align, or use the menus.')
            self.make_tree()
        if 'report' in explicit and self.a.report:self.make_report(open_preview=False)

    def _publish_html(self,prepared=None):
        """Only on explicit user request; does not export the complete session."""
        # Rebuild the report on each request to include new HSP/MSA results.
        report=prepared or self.make_report(open_preview=False)
        default=self._html_target()
        self.ui.info('The HTML will be retained independently of the session (explicit export).')
        dest=self.ui.ask('Destination HTML file',str(default),
                         'Permanent local path: independent of /tmp or an HTTP server')
        if not dest:return None
        target=Path(os.path.expanduser(dest)).resolve()
        if target.suffix.lower()!='.html':self.ui.warn('The .html extension is required.');return None
        if target.exists() and not self.ui.yesno('The report exists. Replace the file?',False):
            self.ui.info('Copy cancelled.');return None
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(report,target)
        self.ui.success(f'Standalone HTML created: {target}')
        self.ui.info(f'Open locally: {target.as_uri()}')
        return target

    def _open_local_html(self,target,temporary=False):
        remote=bool(os.environ.get('SSH_CONNECTION') or os.environ.get('SSH_TTY') or os.environ.get('SLURM_JOB_ID'))
        local_gui=bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform=='darwin' or os.name=='nt')
        if not remote and local_gui and not self.a.no_open:
            import webbrowser
            try:
                if webbrowser.open(target.as_uri(),new=2):
                    self.ui.success('Opened in the local browser.');return
            except (OSError,webbrowser.Error):pass
            self.ui.warn('Cannot open the browser; use the path above.')
        if remote:
            import socket,getpass
            user=getpass.getuser();host=socket.gethostname()
            self.ui.info('VS Code Remote SSH: open the file through the remote explorer.')
            self.ui.info('Or retrieve the file from a terminal on your computer if the host is reachable :')
            print('    scp '+shlex.quote(f'{user}@{host}:{target}')+' .')
            self.ui.info('On a private Slurm node, use the login node (scp -J) or VS Code Remote.')
        self.ui.info('Preview available until the session closes.' if temporary else 'The local HTML remains available after the session closes.')

    def open_report(self):
        self.make_report(open_preview=True)

    def copy_report(self):
        self.ui.heading('EXPORT / STANDALONE HTML')
        target=self._publish_html()
        if target:self._open_local_html(target)

    def export_summary(self):
        stats=[sequence_stats(r) for r in (self.query_list if self.single else [self.query,*self.references])]
        export_table(self.work/'fasta_statistics.tsv',stats)
        if not self.single:self.export_hsps()
        if self.msa_records:
            m=msa_metrics(self.msa_records)
            export_table(self.work/'msa_positions.tsv',m['positions'])
            write_fasta([Seq('CONSENSUS','CONSENSUS',m['consensus'])],self.work/'msa_consensus.fasta')
        was_dirty=self.dirty
        self.save_settings()
        self.dirty=was_dirty
        # Publication/complete profiles include the report and exportable figures.
        if self.a.export_profile in ('publication','complete'):
            self.make_report(open_preview=False)
        (self.work/'commands.json').write_text(json.dumps(self.tools.history,ensure_ascii=False,indent=2))
        summary=dict(mode=self.mode,query=self.query.id,reference_count=len(self.references),
                     hsps=len(self.hsps),hits=len(self.selected),
                     query_length=len(self.query.seq),msa=str(self.msa_file.name if self.msa_file else ''),
                     generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        (self.work/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
        export_dir=generate_exports(self,self.a.export_profile)
        self.ui.success(f'Scientific exports ready: {len(list(export_dir.iterdir()))} files · profile {self.a.export_profile}')

    def export_artifact(self,kind,destination):
        """Write one requested scientific result without changing session state."""
        import csv
        import io
        if kind=='msa':
            if not self.msa_records:raise ValueError('No alignment available.')
            content=''.join('>'+r.id+'\n'+r.seq+'\n' for r in self.msa_records)
        elif kind=='hits':
            if not self.search_completed and not self.hsps:
                if 'search' not in self.vsearch_results:raise ValueError('No comparison available.')
                content=(self.work/'vsearch/search/matches.tsv').read_text(encoding='utf-8')
            else:
                rows=self.coverage_rows()
                columns=list(rows[0]) if rows else ['subject','best_identity','query_coverage','hsps']
                output=io.StringIO();writer=csv.DictWriter(output,fieldnames=columns,delimiter='\t',lineterminator='\n')
                writer.writeheader();writer.writerows(rows);content=output.getvalue()
        elif kind=='tree':
            from av_tree import tree_svg

            content=str(tree_svg((self.work/'tree_fasttree.nwk').read_text()))
        else:raise ValueError('Unknown export type: '+kind)
        target=Path(destination).expanduser().resolve();target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('x',encoding='utf-8') as stream:stream.write(content)
        return target

    def export_alignment_selection(self,viewer):
        records,columns=viewer.export_selection()
        if not records or not columns:raise ValueError('The current selection is empty.')
        suffix=f'{columns[0]}_{columns[-1]}'
        selected=self.ui.ask('Export selection folder (q to cancel)',f'alignview_selection_{suffix}')
        if not selected or selected.casefold()=='q':return
        target=Path(selected).expanduser().resolve()
        target.parent.mkdir(parents=True,exist_ok=True);target.mkdir(exist_ok=False)
        try:
            write_fasta(records,target/'alignment.fasta')
            (target/'columns.tsv').write_text('output_column\tmsa_column\n'+
                ''.join(f'{i}\t{j}\n' for i,j in enumerate(columns,1)),encoding='utf-8')
            (target/'selection.json').write_text(json.dumps(dict(
                source_alignment_sha256=sha256(self.msa_file),sequence_ids=[r.id for r in records],
                original_msa_columns=columns,trace_sequence=viewer.trace_id,
                variable_sites_only=viewer.only_variable),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        except BaseException:
            shutil.rmtree(target);raise
        self.ui.success('Alignment selection exported: '+str(target))
        return target

    def export_center(self):
        from av_navigation import Workbench
        return Workbench(self).export()

    def export_explicit_directory(self):
        """Batch export requested independently of --save."""
        if not self.a.export_dir:return
        self.export_summary()
        target=Path(self.a.export_dir).expanduser().resolve()
        if target.exists():raise ValueError(f'Export directory already exists: {target}')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copytree(self.work/'exports',target)
        self.ui.success(f'Requested exports: {target}')

    def slurm_script(self):
        self.ui.heading('SLURM / GENERATE AN SBATCH SCRIPT')
        threads=self.a.threads
        jobname=re.sub('[^A-Za-z0-9_-]','_',self.query.id)[:30]
        location=f'alignview_{jobname}_slurm'
        launcher=Path(sys.argv[0]).resolve()
        command=([sys.executable,str(launcher)] if launcher.suffix=='.pyz' else ['alignview'] if launcher.name=='alignview' else [str(ROOT/'alignview.sh')])+['inputs/query'+('.fasta.gz' if self.query_path.suffix=='.gz' else '.fasta')]
        if self.refs_path:command.append('inputs/references'+('.fasta.gz' if self.refs_path.suffix=='.gz' else '.fasta'))
        if not self.single:
            command+=['--mode',self.mode]
            if len(self.query_list)>1:command+=['--query-id',self.query.id]
        if self.a.msa!='none':command+=['--msa',self.a.msa,'--aligner',self.a.aligner]
        if self.a.module_policy!='auto':command+=['--module-policy',self.a.module_policy]
        for item in self.a.executable:command+=['--executable',item]
        for item in self.a.tools_dir:command+=['--tools-dir',item]
        command+=['--export-profile',self.a.export_profile,'--preset',self.a.preset]
        command+=['--bitscore',str(self.a.bitscore),'--sensitivity',self.a.sensitivity,
                  '--mafft',self.a.mafft,'--gapopen',str(self.a.gapopen),
                  '--gapextend',str(self.a.gapextend)]
        if self.a.task:command+=['--task',self.a.task]
        if self.a.frame:command+=['--frame',str(self.a.frame)]
        if self.a.matrix:command+=['--matrix',self.a.matrix]
        if self.a.free_ends:command+=['--free-ends']
        if self.a.pairwise!='none':command+=['--pairwise',self.a.pairwise,'--ref',str(self.a.ref)]
        if self.a.tree:command+=['--tree']
        if self.a.dotplot:command+=['--dotplot','--kmer',str(self.a.kmer)]
        if self.a.motif:command+=['--motif',self.a.motif]
        if self.a.vsearch!='none':command+=['--vsearch',self.a.vsearch]
        if self.a.tools_prefix:command+=['--tools-prefix',self.a.tools_prefix]
        command+=extra_cli(self.a)
        command+=['--threads',str(threads),'--hits',str(self.a.hits),'--evalue',str(self.a.evalue),
                  '--identity',str(self.a.identity),'--coverage',str(self.a.coverage),'--no-menu','--report','--save','--out',location]
        contents=(f'#!/bin/bash\n#SBATCH --job-name=alignview_{jobname}\n#SBATCH --cpus-per-task={threads}\n'
                  '#SBATCH --mem=8G\n#SBATCH --time=04:00:00\n#SBATCH --output=alignview_%j.log\n'
                  'set -euo pipefail\n# Run from the saved project directory.\n'
                  'cd -- "$(dirname -- "${BASH_SOURCE[0]}")"\n# Adjust the partition and resources for the cluster.\n'
                  +shlex.join(command)+'\n')
        file=self.work/'alignview_job.sbatch';file.write_text(contents)
        self.dirty=True
        self.ui.success(f'File prepared: {file.name}')
        show_file(file,self.ui)
        self.ui.info('No automatic submission: inspect the script, save the project, then run sbatch.')

    def diagnostics(self):
        rows=diagnostic(self.tools)
        for row in rows:
            self.ui.info(f'{row["command"]}: '+(f'{row["path"]} · {row["version"]}' if row['available'] else f'missing → alignview --setup {row["setup_profile"]} --install-tools'))

    def guide(self):
        self.ui.heading('METHOD GUIDE / SCIENTIFIC LIMITS')
        for key,text in HELP_TEXT.items():
            print();print('  '+self.ui.c(CYAN+BOLD,key.upper()))
            self.ui.wrapped(text,'    ')

    def edit_filters(self):
        self.ui.heading('FILTERS / RERANK HSPS WITHOUT RERUNNING BLAST')
        self.ui.info('Reapply identity, coverage or bitscore filters to existing HSPs.')
        changes={}
        for k,label in [('identity','Minimum identity (%)'),('coverage','Coverage QUERY union (%)'),
                        ('bitscore','Bitscore minimal'),('hits','Maximum number of references')]:
            current=getattr(self.a,k)
            value=self.ui.ask(label,str(current))
            try:changes[k]=int(value) if k=='hits' else float(value)
            except ValueError:self.ui.warn('Invalid value; threshold unchanged')
        if not 0<=changes.get('identity',self.a.identity)<=100 or not 0<=changes.get('coverage',self.a.coverage)<=100 or changes.get('hits',self.a.hits)<1 or changes.get('bitscore',self.a.bitscore)<0:
            self.ui.warn('Thresholds out of range: use valid values.');return
        for key,value in changes.items():setattr(self.a,key,value)
        self.apply_filters()
        self.ui.success(f'{len(self.selected)} retained references without rerunning the search.')
        self.show_hits()

    def rerun(self):
        if self.single or self.mode.startswith('global-'):self.ui.warn('This mode has no local search.');return
        if self.ui.yesno('Rerun the search and replace this session\'s BLAST results?'):
            self.msa_records=None;self.msa_file=None;self.dot_data=None;self.search()

    def single_intro(self):
        self.show_dashboard()
        if self.a.msa=='view':self.show_msa_window()
        elif self.a.msa in ('align','full','all','global'):self.msa_run('all')

    def initial(self):
        command=getattr(self.a,'_command','')
        if self.single and (self.a.no_menu or not self.ui.tty):
            if command in ('view','report') and len({len(r.seq) for r in self.query_list})!=1:
                raise ValueError('This task requires an aligned FASTA with equal sequence lengths. Use align first.')
            if command=='align' and len(self.query_list)<2:
                raise ValueError('align requires at least two sequences.')
        if getattr(self.a,'_command','') in ('view','align','compare') and self.ui.tty and not self.a.no_menu:
            self.a._guided=True;self.a._intent=self.a._command
        explicit_actions=self.a.msa!='none' or self.a.tree or self.a.vsearch!='none' or self.a.pairwise!='none' or self.a.dotplot or self.a.motif or self.a.report
        if self.ui.tty and not self.a.no_menu and (getattr(self.a,'_guided',False) or not self.mode or not explicit_actions):
            from av_navigation import Workbench
            workbench=Workbench(self)
            intent=getattr(self.a,'_intent','')
            if self.single and self.a.msa=='view':
                if len({len(r.seq) for r in self.query_list})!=1:
                    self.a.msa='none'
                    self.ui.notice='Different lengths: choose Align sequences.'
                    return
                self.load_msa()
            elif not self.single and not self.mode:
                self.mode=workbench.recommended_mode() or ''
                self.a.mode=self.mode
            if intent=='view':self.show_msa_window()
            elif intent=='demo':self._open_view='Home'
            elif intent=='align':workbench.run_action(workbench.align)
            elif intent=='compare':workbench.run_action(workbench.compare)
            return
        self.ui.banner('ALIGNMENTS / GUIDED EXPLORATION')
        if self.a.vsearch!='none':
            self.vsearch_run(self.a.vsearch)
            if self.a.msa=='view':self.load_msa()
            elif self.a.msa!='none':self.msa_run('all')
            if self.a.tree:self.make_tree()
            if self.a.report:self.make_report()
            return
        if self.single:
            self.single_intro()
            if self.a.motif:
                data=[]
                for s in self.query_list:data.extend(find_motifs(s,self.a.motif,True))
                self.motif_rows=data
                export_table(self.work/'motifs.tsv',data,['seqid','start','end','strand','matched'])
            if self.a.dotplot:self._dot_batch()
            if self.a.tree:
                if not self.msa_records:self.msa_run('all')
                self.make_tree()
            if self.a.report:self.make_report()
            return
        if not self.mode:
            if self.a.no_menu:raise ValueError('Choose a search method with --mode, or use alignview compare for automatic type-based selection.')
            from av_navigation import Workbench
            self.mode=Workbench(self).recommended_mode() or 'blastn'
            self.a.mode=self.mode
        if self.mode.startswith('global-'):
            self.ui.info('GLOBAL mode: no BLAST search; alignment runs on request.')
            if self.a.msa in ('align','full','all','global','hits'):self.msa_run('all')
        else:
            self.search()
            if self.a.msa in ('align','full','hits','global'):self.msa_run('hits')
            elif self.a.msa in ('all','fragments','hsp'):self.msa_run(self.a.msa)
        if self.a.pairwise!='none':self.run_pairwise('needle' if self.a.pairwise=='global' else 'water')
        if self.a.motif:
            data=find_motifs(self.query,self.a.motif,True)
            self.motif_rows=data;export_table(self.work/'motifs.tsv',data,['seqid','start','end','strand','matched'])
            self.ui.success(f'Motif {self.a.motif}: {len(data)} occurrence(s)')
        if self.a.dotplot:self._dot_batch()
        if self.a.tree:
            if not self.msa_records:self.msa_run('all')
            self.make_tree()
        if self.a.report:self.make_report()

    def _dot_batch(self):
        refs=self.references if self.references else [r for r in self.query_list if r.id!=self.query.id]
        if not refs:return
        ref=self.refmap[next(iter(self.selected))] if self.selected else refs[0]
        self.dot_data=self._compute_dotplot(ref,self.a.kmer)
        (self.work/'dotplot.svg').write_text(dotplot_svg(self.dot_data,f'{self.query.id} ({self.dot_data["xunit"]})',f'{ref.id} ({self.dot_data["yunit"]})'))
        (self.work/'dotplot.json').write_text(json.dumps(self.dot_data))

    def exit_dialog(self):
        from av_navigation import Workbench
        return Workbench(self).confirm_exit()


    def nav(self):
        from av_navigation import Workbench
        return Workbench(self).run()


    def show_alv(self):
        if not self.msa_file:self.load_msa()
        if not self.msa_file:self.ui.warn('No multiple alignment available');return
        print(self.tools.alv(self.msa_file)[:150000])

    def list_outputs(self):
        self.ui.heading('FILES / TEMPORARY WORKSPACE')
        for file in sorted(self.work.rglob('*')):
            if file.is_file() and 'database' not in str(file):
                print(f'  {str(file.relative_to(self.work)):<48} {file.stat().st_size/1024:9.1f} KiB')
        self.ui.info('Temporary work: save or export explicitly to retain results.')

    def save_checkpoint(self):
        if self.last_saved and Path(self.last_saved).is_dir():
            if not self.dirty:
                self.ui.notice='Project is already saved: '+str(self.last_saved);return Path(self.last_saved)
            previous=Path(self.last_saved)
            stem=re.sub(r'_snapshot_\d+$','',previous.name)
            number=2
            while (previous.parent/f'{stem}_snapshot_{number:03d}').exists():number+=1
            selected=previous.parent/f'{stem}_snapshot_{number:03d}'
        else:
            safe=re.sub(r'[^A-Za-z0-9_.-]','_',self.source_query_path.stem)[:50] or 'project'
            default=self.a.out or 'alignview_'+safe
            selected=self.ui.ask('Project directory (q to cancel)',default)
            if not selected or str(selected).casefold()=='q':return None
        self.export_summary()
        target=save_project(self,selected)
        self.last_saved=target;self.dirty=False
        self.ui.success(f'Project saved: {target}')
        self.ui.info('Resume: alignview resume '+shlex.quote(str(target)))
        return target

    def data_hub(self):
        from av_navigation import Workbench
        return Workbench(self).data()

    def save_and_exit(self):
        if self.exit_decision=='clean':
            self.ui.info('Session closed.');return
        if self.exit_decision=='saved':
            self.ui.info(f'Project retained: {self.last_saved}');return
        if self.exit_decision=='discard':
            self.ui.info('Session closed without saving temporary results.');return
        if not self.a.save:
            if self.a.no_menu or not self.ui.tty:
                if self.a.resume and not self.dirty:
                    self.ui.info('Project retained: '+str(self.last_saved or self.a.resume));return
                if self.a.out:self.ui.info('--out names a directory; add --save to retain the project there.')
                if self.a.report and not (self.a.html_out or self.a.export_dir):
                    self.ui.warn('Report not retained: add --html-out report.html, --save or --export-dir DIRECTORY.')
                self.ui.info('Finished: temporary results removed (use --save to retain them).')
                return
            self.ui.heading('END SESSION / SAVE')
            self.ui.info('Analysis files are temporary. Unsaved files will be removed.')
            if not self.ui.yesno('Keep this session\'s files?'):
                self.ui.info('Session ended; temporary files removed.');return
        self.export_summary()
        target=self.a.out or f'alignview_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
        target_path=Path(target).expanduser().resolve()
        if target_path.exists():
            self.ui.warn(f'Directory already exists: {target_path}')
            if self.a.no_menu or not self.ui.tty:raise ValueError(f'Output directory already exists: {target_path}. Choose --out with an unused path.')
            n=2
            while target_path.with_name(f'{target_path.name}_{n}').exists():n+=1
            proposed=target_path.with_name(f'{target_path.name}_{n}')
            selected=self.ui.ask('New directory (q = cancel saving)',str(proposed))
            if selected.lower()=='q':
                if self.ui.yesno('Cancel saving and remove temporary files?',False):return
                selected=str(proposed)
            target_path=Path(selected).expanduser().resolve()
            if target_path.exists():
                self.ui.warn('Directory already exists: automatically choosing an unused alternative directory.')
                target_path=proposed
        target_path.parent.mkdir(parents=True,exist_ok=True)
        save_project(self,target_path)
        self.ui.success(f'Results saved: {target_path}')
        if (target_path/'report.html').exists():self.ui.info('Open report.html in a browser; no network connection is required.')
