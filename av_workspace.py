"""Two persistent destinations, object readers and contextual analyses.

The workspace reads session state. Calculations remain explicit controller actions;
redrawing, filtering and switching views never invoke a scientific engine.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re

from av_core import sequence_type
from av_terminal import fit
from av_version import VERSION

TABS = ('Data', 'Results')
VIEWS = ('Home', 'ResultHome', *TABS, 'Sequences', 'Alignment', 'Matches', 'Tree', 'Artifact', 'Vsearch')


def comparison_context(session, query_ids=None):
    """The direction of a search, using original filenames rather than scratch paths."""
    if query_ids is None:
        if session.workspace_entry.get('analysis')=='search_batch':
            query_ids=session.comparison_batch.get('query_ids',[])
        else:query_ids=[] if session.single or session._comparison_needs_selection else [session.query.id]
    identifiers=list(query_ids)
    source=session.source_query_path.name
    if len(identifiers)==1:query=f'Query: {identifiers[0]} · {source}'
    elif identifiers:query=f'Queries: {source} · {len(identifiers)} selected sequences'
    else:query=f'Query file: {source} · choose from {len(session.query_list)} sequences'
    target=(f'Target: {session.source_refs_path.name} · {len(session.references)} sequences'
            if session.source_refs_path else 'Target: not chosen · select a FASTA to search against')
    return query,target


@dataclass(frozen=True)
class Command:
    key: str
    label: str
    detail: str
    action: str
    enabled: bool = True
    reason: str = ''


def evidence_rows(session):
    """Keep local HSP evidence and native global-search metrics distinct."""
    if session.search_completed or session.hsps:
        return session.coverage_rows(), 'Local search'
    result=session.vsearch_results.get('search')
    if not result:return [], 'Comparison'
    path=session.work/'vsearch/search/matches.tsv'
    if path.is_file():
        with path.open(encoding='utf-8') as stream:raw=list(csv.DictReader(stream,delimiter='\t'))
    else:raw=result.get('preview',[])
    rows=[]
    for row in raw:
        target=row.get('target','')
        if not target:continue
        rows.append(dict(subject=target,best_identity=float(row.get('id',0)),
            query_coverage=float(row.get('qcov',0)),reference_coverage=float(row.get('tcov',0)),
            hsps=None,best_evalue=None,best_bitscore=None,global_match=True))
    return rows,'VSEARCH global search'


class Workspace:
    """Small view state independent from calculation/project state."""
    def __init__(self,session):
        self.session=session
        saved=session.navigation_state.get(self.object_id,{})
        requested=getattr(session,'_open_view',None)
        self.tab=requested if requested in VIEWS else ('ResultHome' if session.workspace_entry.get('kind')=='result' else 'Home')
        if getattr(session.a,'_intent','')=='demo' and not requested:self.tab='Home'
        self.active=0;self.query='';self.filtering=False;self.offset=0
        self.notice='';self.page_size=8;self._rows=[]
        self.views={};self.trail=[];self._catalog_cache=None
        self.mouse_rows={};self.body_hits={}
        self.views=dict(saved.get('views',{}))
        if 'sequence_selection' in saved:
            session.sequence_selection=[name for name in saved['sequence_selection'] if any(r.id==name for r in session.query_list)]
        if 'reference_selection' in saved:
            session.reference_selection=[name for name in saved['reference_selection'] if name in session.refmap]
        self._restore_view()

    @property
    def object_id(self):
        return self.session.workspace_entry.get('id','input:'+self.session.input_hashes['query'])

    def _restore_view(self):
        state=self.views.get(self.tab,{})
        self.active=state.get('active',0);self.offset=state.get('offset',0)
        self.query=state.get('query','');self.filtering=False

    def remember(self):
        self.views[self.tab]=dict(active=self.active,offset=self.offset,query=self.query)
        state=self.session.navigation_state.setdefault(self.object_id,{})
        state.update(tab=self.tab,views=self.views,sequence_selection=list(self.session.sequence_selection),
                     reference_selection=list(self.session.reference_selection))
        if self.session._viewer is not None:
            viewer=self.session._viewer
            state['viewer']={name:getattr(viewer,name) for name in
                ('row','col','xoff','yoff','only_variable','row_filter','trace_id','difference','region') if hasattr(viewer,name)}
            state['viewer']['checked']=sorted(viewer.checked)

    def back(self):
        if self.trail:
            self.switch(self.trail.pop(),record=False)
            return True
        return False

    def catalog(self,kind):
        from av_library import entries,entry_label
        signature=tuple(self.session.history)
        if self._catalog_cache is None or self._catalog_cache[0]!=signature:
            self._catalog_cache=(signature,entries(self.session))
        return [dict(row,label=entry_label(self.session,row)) for row in self._catalog_cache[1] if row['kind']==kind]

    def can_read_alignment(self):
        s=self.session
        return bool(s.msa_records) or (s.single and len(s.query_list)>=2 and len({len(r.seq) for r in s.query_list})==1)

    def landing(self):
        """Explain the loaded file before offering a small set of explicit actions."""
        s=self.session
        kinds={sequence_type(record.seq) for record in s.query_list}
        kind=next(iter(kinds)) if len(kinds)==1 else 'mixed-type'
        lengths=[len(record.seq.replace('-','').replace('.','')) for record in s.query_list]
        span=f'{min(lengths):,}' if min(lengths)==max(lengths) else f'{min(lengths):,}–{max(lengths):,}'
        unit='nt' if kind=='DNA/RNA' else 'aa' if kind=='protein' else 'characters'
        summary=f'{len(lengths)} {kind} sequence'+('s' if len(lengths)!=1 else '')+f' · {span} {unit}'
        if self.tab=='ResultHome':
            entry=s.workspace_entry
            names={'msa':'Alignment','search':'Sequence matches','tree':'Tree','motif':'Motif search',
                   'translation':'Translation','orfs':'ORF search','extract':'Region extraction',
                   'statistics':'Sequence checks','conservation':'Alignment conservation','distances':'Sequence distances',
                   'profile':'Sequence profile','composition':'Sequence composition','dotplot':'Dotplot',
                   'needle':'Global pairwise alignment','water':'Local pairwise alignment'}
            names['search_batch']='Sequence comparisons'
            title=names.get(entry.get('analysis'), 'Analysis')+' ready'
            if entry.get('analysis')=='msa' and s.msa_records:
                summary=f'{len(s.msa_records)} aligned sequences · {len(s.msa_records[0].seq):,} columns'
            elif entry.get('analysis')=='search':
                summary=f'{len(s.selected)} matching targets · {len(s.hsps)} local aligned segments'
            elif entry.get('analysis')=='search_batch':
                runs=s.comparison_batch.get('runs',[])
                completed=sum(run['status']=='completed' for run in runs)
                summary=f'{completed}/{len(runs)} queries completed · {len(s.references)} reference sequences'
            elif entry.get('analysis') in ('tree','conservation','distances') and s.msa_records:
                summary=f'Based on an alignment of {len(s.msa_records)} sequences'
            elif entry.get('analysis')=='statistics':summary=f'{len(s.query_list) if s.single else 1+len(s.references)} sequences inspected'
            elif entry.get('analysis') in ('motif','translation','orfs','extract','profile','composition'):
                summary='Input sequence: '+s.query.id
            repeat={'msa':'align','search':'compare','search_batch':'compare','tree':'tree','translation':'translate',
                    'motif':'motif','orfs':'orfs','extract':'extract','dotplot':'dotplot'}.get(entry.get('analysis'))
            options=[
                Command('1','Open this result','Read the completed result; no calculation is repeated.','open_result'),
                Command('2','Try different settings' if repeat else 'Choose another analysis',
                        'Keep this result and prepare a separate attempt.',repeat or 'analysis'),
                Command('3','Back to file actions','Return to the menu for '+s.source_query_path.name+'.','home'),
                Command('4','Previous results','Reopen any earlier completed analysis.','results'),
                Command('5','Save or export','Optional: keep work after closing or create reusable files.','save'),
                Command('6','Rename this result','Give this attempt a name you can find in Previous results.','rename_result'),
            ]
            from av_library import entry_label
            return title,options,summary,entry_label(s,entry)
        from av_navigation import WorkspaceState
        state=WorkspaceState.from_session(s)
        align_detail=(f'Compare positions across {len(s.sequence_selection)} checked sequences.' if s.sequence_selection else
                      f'Compare positions across the {len(s.query_list)} loaded sequences.' if s.single else
                      'Choose which loaded sequences to align together.')
        options=[
            Command('1','Look at the sequences','Browse sequence content, or explicitly read an existing alignment.','inspect_data'),
            Command('2','Align sequences',align_detail,'align',state.can_align,state.align_reason),
            Command('3','Find matching sequences','Choose the sequences to search and the target FASTA to search against.','compare'),
            Command('4','Other analyses','Motifs, translation, regions, alignment statistics and trees.','analysis'),
            Command('5','Files and previous results','Change the loaded file, reopen results, or save/export.','manage'),
        ]
        completed=[row for row in self.catalog('result') if row.get('dataset_hash')==s.input_hashes['query'] or
                   (not row.get('dataset_hash') and row.get('source')==s.source_query_path.name)]
        members={identifier for row in completed for identifier in row.get('members',[])}
        completed=[row for row in completed if row.get('id') not in members]
        status=f'{len(completed)} completed result(s) available in Files and previous results.' if completed else 'File loaded. No analysis has been run.'
        return 'File loaded — choose an action',options,summary,status

    def switch(self,name,record=True):
        if name not in VIEWS:raise ValueError('Unknown workspace view: '+name)
        if name==self.tab:return
        self.remember()
        if record:self.trail.append(self.tab)
        self.tab=name;self._restore_view()

    def rows(self):
        s=self.session
        if self.tab in TABS:
            rows=self.catalog('data' if self.tab=='Data' else 'result')
            if self.tab=='Data' and not rows:
                rows=[dict(label=s.source_query_path.name,source=s.source_query_path.name,view='Sequences',path=None)]
        elif self.tab=='Matches':rows=evidence_rows(s)[0]
        elif self.tab=='Sequences':
            roles=not s.single
            rows=[dict(record=r,role=('Query' if r.id==s.query.id else 'Query candidate') if roles else 'Input sequence') for r in s.query_list]
            if not s.single:rows.extend(dict(record=r,role='Reference') for r in s.references)
        elif self.tab=='Alignment':rows=[dict(record=r) for r in (s.msa_records or [])]
        else:rows=[]
        if self.query:
            term=self.query.casefold()
            rows=[r for r in rows if term in (r.get('label','')+' '+r.get('source','') if 'label' in r else r['subject'] if 'subject' in r else r['record'].id+' '+r['record'].description).casefold()]
        self._rows=rows;self.active=max(0,min(self.active,len(rows)-1))
        return rows

    def selected(self):
        rows=self.rows()
        return rows[self.active] if rows else None

    def primary(self):
        s=self.session
        if self.tab in ('Home','ResultHome'):return 'open_result' if self.tab=='ResultHome' else 'inspect_data','Open result' if self.tab=='ResultHome' else 'Look at the sequences'
        if self.tab in TABS:return 'open_entry','Open selected '+('dataset' if self.tab=='Data' else 'result')
        if self.tab in ('Artifact','Vsearch'):return 'artifact','Read result'
        if self.tab=='Alignment':return ('view','Open alignment' if s.msa_records else 'Read file as alignment') if self.can_read_alignment() else ('align','Create alignment')
        if self.tab=='Matches':
            rows=self.rows()
            if rows:return 'inspect','Inspect global match' if rows[self.active].get('global_match') else 'View aligned regions'
            if s.vsearch_results and not s.search_completed:return 'vsearch_results','Open VSEARCH results'
            return ('compare','Run comparison') if not s.single else ('compare','Add references and compare')
        if self.tab=='Tree':return ('tree_view','Open tree') if (s.work/'tree_fasttree.nwk').is_file() else ('tree','Create tree')
        return 'sequence','Inspect selected sequence'

    def commands(self):
        s=self.session;aligned=bool(s.msa_records);has_refs=bool(s.references)
        nucleic=sequence_type(s.query.seq)=='DNA/RNA'
        all_nucleic=all(sequence_type(r.seq)=='DNA/RNA' for r in s.query_list)
        from av_navigation import WorkspaceState
        state=WorkspaceState.from_session(s)
        translated=s.mode in ('tblastn','tblastx')
        can_align=state.can_align or (translated and bool(s.selected))
        align_reason=state.align_reason if not translated else 'Run a local comparison first; only detected segments can be aligned across reference frames.'
        chosen=len(getattr(s,'reference_selection',[]))
        can_select=chosen>0 and (state.can_align or (translated and all(name in s.selected for name in s.reference_selection)))
        translation_relevant=nucleic or s.mode in ('tblastn','tblastx')
        return [
            Command('analysis','New analysis','Choose an operation using the loaded data','analysis'),
            Command('data','Browse loaded datasets','All files stay available in this session','datasets'),
            Command('results','Browse results',f'{len(self.catalog("result"))} completed analyses; no saving required','results'),
            Command('rename_result','Rename this result','Give the current attempt a recognizable name','rename_result',
                    s.workspace_entry.get('kind')=='result','Open a completed result first.'),
            Command('open_alignment','Open alignment','Full-screen navigation, trace sequence and region export','view',self.can_read_alignment(),'Load equal-length aligned sequences, or create an alignment first.'),
            Command('compare','Find matching sequences','Choose query sequences and a target dataset; review the search direction','compare'),
            Command('align','Create an alignment',f'{chosen} selected references' if chosen else 'Choose all sequences or the retained references','align',can_align,align_reason),
            Command('selection','Align selected references',f'{chosen} references checked in Matches','align_selection',can_select,'Use Space in Matches to select references with compatible sequence or HSP data.'),
            Command('tree','Create an exploratory tree','FastTree on the current alignment','tree',aligned and len(s.msa_records)>=3,'At least three aligned sequences are required.'),
            Command('query','Choose active query',f'{len(s.query_list)} query candidates','query',len(s.query_list)>1,'Only one query is loaded.'),
            Command('open','Add sequences','Load another FASTA; keep all earlier data and results','files'),
            Command('reuse','Use result as sequence data','Continue with a translation, ORFs, extracted region or alignment','reuse'),
            Command('inputs','Open result inputs','Return to the exact dataset or alignment that produced this result','inputs',bool(s.workspace_entry.get('inputs')),'Open a calculated result first.'),
            Command('save','Save project for later','Optional: retain the whole workspace after closing','checkpoint'),
            Command('report','Preview HTML report','Open in a browser without choosing an export destination','report'),
            Command('html','Export HTML report','Keep a standalone file after closing','html_export'),
            Command('export','Export data or figures','Alignment FASTA, evidence TSV, tree SVG or result folder','export'),
            Command('vsearch','VSEARCH nucleotide analyses','Global search, dereplication, clustering and chimera detection','vsearch',all_nucleic,'VSEARCH requires nucleotide sequences.'),
            Command('pairwise','Compare two sequences','Needle global alignment or Water local alignment','pairwise',has_refs,'Load a reference file.'),
            Command('motif','Find a nucleotide motif','IUPAC positions on both strands of the active query','motif',nucleic,'Select a nucleotide query.'),
            Command('translate','Translate all six frames',f'NCBI genetic code {s.a.genetic_code}','translate',nucleic,'Select a nucleotide query.'),
            Command('orfs','Find simple ORFs','Exploratory start-to-stop coding regions','orfs',nucleic,'Select a nucleotide query.'),
            Command('extract','Extract a query region','Inclusive native coordinates, optional reverse complement','extract'),
            Command('dotplot','Create a dotplot','Exact word matches between two sequences','dotplot',has_refs or len(s.query_list)>1,'Load at least two sequences.'),
            Command('conservation','Alignment conservation','Consensus, variable sites and conservation','conservation',aligned,'Open or create an alignment.'),
            Command('distances','Sequence distances','Pairwise p-distances on the alignment','distances',aligned,'Open or create an alignment.'),
            Command('statistics','Sequence quality','Lengths, composition and ambiguous symbols','statistics'),
            Command('translation_table','Choose translation table','NCBI genetic codes and ORF start codons','translation_table',translation_relevant,'Load a nucleotide query or translated nucleotide references.'),
            Command('advanced','Additional analyses and settings','VSEARCH, translation, diagnostics, profiles and HPC','advanced'),
            Command('tools','Tools and execution history','Availability, versions, commands and logs','tools'),
            Command('help','Read the manual','Task guides and scientific interpretation','help'),
        ]

    def render(self,width=100,height=30):
        if self.tab in ('Home','ResultHome'):
            title,options,summary,status=self.landing()
            old_context,old_route=self.session.ui.context,self.session.ui.route
            old_compact=self.session.ui.compact_context_lines
            self.session.ui.context=(summary,status) if self.session.single else comparison_context(self.session)
            self.session.ui.compact_context_lines=1 if self.session.single else 2
            self.session.ui.route=['Workspace',self.session.source_query_path.name]
            from av_ui import menu_item
            try:
                items=[menu_item(c.key,c.label,c.detail,enabled=c.enabled,reason=c.reason) for c in options]
                frame=self.session.ui._menu_frame(title,items,
                    subtitle=summary if not self.session.single else '',
                    hint=(status if not self.session.single and self.tab=='ResultHome' else
                          '↑↓ choose · Enter open · Esc back · q '+('close' if self.tab=='Home' else 'file menu')),
                    dimensions=(width,height)).lines
                return [fit(line,width-1) for line in frame][:max(1,height-1)]
            finally:
                self.session.ui.context,self.session.ui.route=old_context,old_route
                self.session.ui.compact_context_lines=old_compact
        s=self.session;rows=self.rows();compact=height<20
        title=s.source_query_path.stem
        saved='Saved copy available' if s.last_saved else 'Temporary session · saving optional'
        current='Data' if self.tab in ('Data','Sequences') else 'Results'
        tabs=' '+'    '.join(f'{i} [{name}]' if name==current else f'{i} {name}' for i,name in enumerate(TABS,1))
        detail='' if self.tab in TABS else ' / '+self.tab
        inputs=([f' {line}' for line in comparison_context(s,[s.query.id])] if not s.single else
                [f' Input FASTA: {len(s.query_list)} sequences  ·  active sequence: {s.query.id}'])
        head=[f' AlignView {VERSION}  /  {title}  ·  {saved}',
              tabs+detail,
              ' o Add data   n New analysis   w Preview HTML   s Save / Export',
              *(inputs if self.tab not in TABS else [' All datasets and results remain accessible until you close this session.']),
              ' ─'+'─'*max(1,width-5)]
        action,label=self.primary()
        foot=[f' Enter {label}  ·  Esc Back  ·  : Commands  ·  q '+('Close' if self.tab in TABS else 'Back'),
              (' Space Check ref  ·  a Align checked  ·  / Filter' if self.tab=='Matches' else
               ' Space Select sequence  ·  a Align selection  ·  u Use sequence  ·  / Filter' if self.tab=='Sequences' else
               ' / Filter results  ·  n New analysis  ·  1 Data  ·  2 Results' if self.tab in TABS else
               ' n New analysis  ·  u Use as data  ·  1 Data  ·  2 Results')]
        if self.filtering:foot.append(' Filter: '+self.query+'▏  ·  Enter apply / Esc clear')
        else:foot.append(' '+(self.notice or s.ui.notice or 'Explore any analysis. Earlier results are kept automatically during this session.'))
        if not compact:head.append(' '+self._summary())
        available=max(1,height-1-len(head)-len(foot))
        content=self._body(rows,width,available)
        self.mouse_rows={len(head)+relative+1:index for relative,index in self.body_hits.items() if len(head)+relative<height-1-len(foot)}
        self.page_size=max(1,available-3)
        content=content[:available]+['']*max(0,available-len(content))
        return [fit(line,width-1) for line in [*head,*content,*foot]][:max(1,height-1)]

    def _summary(self):
        s=self.session
        if self.tab=='Data':return f'{len(self.rows())} dataset(s). Enter opens sequences; o adds another FASTA.'
        if self.tab=='Results':return f'{len(self.rows())} result(s). Open an earlier attempt without removing newer work.'
        if self.tab in ('Artifact','Vsearch'):
            from av_library import entry_label
            return entry_label(s,s.workspace_entry)
        if self.tab=='Alignment':
            return (f'{len(s.msa_records)} sequences × {len(s.msa_records[0].seq):,} columns. Enter opens the alignment workspace.'
                    if s.msa_records else 'No alignment yet. Open an aligned FASTA or create an alignment.')
        if self.tab=='Matches':
            if not self.rows() and s.vsearch_results and not s.search_completed:return f'{len(s.vsearch_results)} completed VSEARCH analyses · Enter opens their native results'
            kind=evidence_rows(s)[1]
            return f'{kind} · {len(self.rows())} references · {len(getattr(s,"reference_selection",[]))} selected for alignment'
        if self.tab=='Tree':return 'Exploratory topology. Drawing root is arbitrary; local support is not bootstrap support.'
        if self.session.single:return f'{len(self.session.query_list)} sequences · {len(s.sequence_selection)} checked for alignment. n opens analyses.'
        return 'Query = one active sequence. References = the database to compare against. u changes the active query.'

    def _body(self,rows,width,available):
        s=self.session
        self.body_hits={}
        if self.tab in ('Artifact','Vsearch'):
            from av_project import local_path
            artifact=s.workspace_entry.get('artifact','')
            if artifact:
                path=local_path(s.work,artifact)
                with path.open(encoding='utf-8',errors='replace') as stream:
                    from itertools import islice
                    preview=list(islice(stream,max(0,self.offset),self.offset+max(1,available-1)))
                return [' '+line.rstrip() for line in preview]+[' Enter reads the complete result; u reuses sequence outputs.']
            return self._empty('Analysis result','Enter reads this result. n starts another analysis.')
        if self.tab in TABS:
            if not rows:return self._empty('No results yet.','Press n to choose an analysis. No project save is needed.')
            count=max(1,available-2);start=min(max(0,self.active-count//2),max(0,len(rows)-count))
            body=[]
            for i,row in enumerate(rows[start:start+count],start):
                marker='›' if i==self.active else ' '
                self.body_hits[len(body)]=i
                body.append(f' {marker} {row["label"]}')
                if width>=90 and len(rows)<count//2:
                    body.append('     '+row.get('source','')+' · '+row.get('view','Sequences'))
            body.append(f' {len(rows)} '+('datasets' if self.tab=='Data' else 'results')+' · ↑↓ select · Enter open')
            return body
        if self.tab=='Tree':
            path=s.work/'tree_fasttree.nwk'
            if not path.is_file():return self._empty('Build a tree from the current alignment.',
                'Enter: prepare FastTree' if s.msa_records and len(s.msa_records)>=3 else 'Create an alignment with at least three sequences first.')
            from av_tree import tree_text
            lines=tree_text(path.read_text()).splitlines()[4:]
            self.offset=max(0,min(self.offset,max(0,len(lines)-available+1)))
            return [' Topology preview · Enter opens the complete tree',*[' '+line for line in lines[self.offset:self.offset+max(1,available-1)]]]
        if not rows:
            if self.query:return self._empty('No matching identifiers.','Esc clears the filter.')
            if self.tab=='Matches' and s.search_completed:
                if s.hsps:
                    return self._empty('No matches pass the current filters.',
                        f'{len(s.hsps)} local segment'+(' was' if len(s.hsps)==1 else 's were')+' found. Press c to review thresholds and filters.')
                return self._empty('Search completed: no local aligned segments found.',
                    'Press c to review the method, target and search threshold. No hit does not rule out homology.')
            if self.tab=='Matches' and s.vsearch_results:
                return ['', ' Native VSEARCH results', '', *[' '+name+': '+', '.join(f'{key}={value}' for key,value in result.get('summary',{}).items()) for name,result in s.vsearch_results.items()], '', ' Enter opens the result tables and files.']
            if self.tab=='Matches':return self._empty('Find the references that match your query.',
                'Press c to prepare a comparison. Native VSEARCH and local BLAST results stay distinct.')
            if self.tab=='Alignment' and self.can_read_alignment():return self._empty('The loaded file contains equal-length sequences.',
                'Enter reads the file as an alignment; a creates a new alignment. Equal length alone is not proof of alignment.')
            if self.tab=='Alignment':return self._empty('Inspect homologous sequences position by position.',
                'Press a to prepare an alignment. Existing alignments can be opened without an engine.')
            return self._empty('No records are loaded.','Press o to open a dataset.')
        count=max(1,available-3)
        start=min(max(0,self.active-count//2),max(0,len(rows)-count))
        chosen=set(getattr(s,'reference_selection',[]))
        if self.tab=='Matches':
            idw=max(12,width-54)
            body=[f'     {"Reference":<{idw}} {"Query cov":>9} {"Identity*":>9} {"HSPs":>5} {"E-value":>10}']
            for i,row in enumerate(rows[start:start+count],start):
                marker='›' if i==self.active else ' '
                pick='✓' if row['subject'] in chosen else ' '
                ev='—' if row.get('best_evalue') is None else f'{row["best_evalue"]:.1g}'
                hsps='—' if row.get('hsps') is None else str(row['hsps'])
                self.body_hits[len(body)]=i
                body.append(f' {marker} [{pick}] {fit(row["subject"],idw):<{idw}} {row["query_coverage"]:8.1f}% {row["best_identity"]:8.1f}% {hsps:>5} {ev:>10}')
            body.append(' Enter inspects global metrics. Identity*: native VSEARCH global identity.' if rows[0].get('global_match') else
                        ' Enter reads aligned regions. Identity*: strongest local segment, not whole-sequence identity.')
        elif self.tab=='Alignment':
            idw=max(12,min(28,width//3));columns=max(1,width-idw-7)
            body=[f' {"Sequence":<{idw}}  Alignment preview · first {columns} columns']
            for i,row in enumerate(rows[start:start+count],start):
                self.body_hits[len(body)]=i
                r=row['record'];body.append(f' {"›" if i==self.active else " "} {fit(r.id,idw):<{idw}}  {r.seq[:columns]}')
            body.append(' Enter: full viewer with trace selection, position search and scoped export')
        else:
            idw=max(14,width-43);body=[f'   {"Identifier":<{idw}} {"Length":>10} {"Type":<9} Role']
            for i,row in enumerate(rows[start:start+count],start):
                self.body_hits[len(body)]=i
                r=row['record'];pick='✓' if r.id in s.sequence_selection and row['role']!='Reference' else ' '
                body.append(f' {"›" if i==self.active else " "} {pick} {fit(r.id,idw):<{idw}} {len(r.seq):10,} {sequence_type(r.seq):<9} {row["role"]}')
            body.append(' Enter details · Space select · a align · n all analyses')
        body.append(f' Showing {start+1}–{min(len(rows),start+count)} of {len(rows)} · ↑↓ / PgUp PgDn navigate')
        return body

    @staticmethod
    def _empty(title,detail):return ['', ' '+title,'',' '+detail]
