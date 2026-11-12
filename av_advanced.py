"""Advanced analyses and settings, kept out of the everyday workbench."""
from __future__ import annotations

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import re
import shutil

from av_core import sequence_type
from av_parameters import relevant, snapshot, load_profile, write_profile, BY_NAME
from av_defaults import SEARCH_PRESETS
from av_start import choose_file, choose_sequence, startup
from av_ui import menu_item as item


class AdvancedActions:
    def alignment_scope(self):
        state=self.state
        return self.choose('Sequences to align',[
            item('hits','Retained references','The query and matching reference sequences',enabled=bool(self.s.selected) and state.can_align,
                 reason='Requires search results with retained references and a compatible frame.'),
            item('all','All references','Also include references without a match',enabled=state.can_align,reason=state.align_reason),
            item('fragments','Add a partial sequence','Add the query to a reference alignment using MAFFT',enabled=state.can_align,reason=state.align_reason),
            item('hsp','Detected segments only','Exploratory comparison of local regions, excluding flanks',enabled=bool(self.s.selected),
                 reason='Run a search before comparing its segments.'),
        ])

    def choose_aligner(self):
        choice=self.choose('Alignment engine',[
            item('mafft','MAFFT','General-purpose choice; automatic strategy by default'),
            item('muscle','MUSCLE 5','Alternative for DNA and proteins'),
            item('clustalo','Clustal Omega','Alternative for multiple alignments'),
        ],f'Current engine: {self.s.a.aligner.upper()}')
        if choice:
            self.s.a.aligner=choice;self.s.save_settings();self.ui.notice='Selected engine: '+choice.upper()

    def choose_search_method(self):
        from alignview import MODES
        from av_start import suggest_modes
        from collections import Counter
        qtype=sequence_type(self.s.query.seq)
        rtype=Counter(sequence_type(r.seq) for r in self.s.references).most_common(1)[0][0]
        default,_,_=suggest_modes(qtype,rtype)
        allowed={'blastn':('blastn','tblastx','vsearch'),'blastp':('blastp','diamond-blastp','mmseqs-blastp'),
                 'blastx':('blastx','diamond-blastx'),'tblastn':('tblastn',)}.get(default,())
        return self.choose('Search method',[item(mode,MODES[mode][0],MODES[mode][1],
            enabled=not allowed or mode in allowed,reason=f'Incompatible with the detected types: {qtype} → {rtype}.')
            for mode in MODES if not mode.startswith('global-')],f'{qtype} → {rtype} · recommandation: {default or "check required"}')

    def pairwise(self):
        with self.ui.section('Two sequences'):
            choice=self.choose('How should the two sequences be compared?',[
                item('needle','End to end','Global alignment using Needle, including sequence ends'),
                item('water','Best local segment','Optimal local alignment using Water'),
            ])
            if choice and self.ensure_engine((choice,),'emboss'):
                ref,_=self.s.ref_select()
                if ref is None:return
                while True:
                    action=self.choose('Prepare '+choice.upper(),[
                        item('1','Run comparison',f'{self.s.query.id} / {ref.id}'),
                        item('2','Calculation parameters',self.parameter_summary(choice)),
                        item('3','Inspect commands','Matrix, penalties, sequence ends and translation'),
                        item('4','Calculation profiles','Load or save parameters'),
                    ],subtitle=self.s.tools.engine_info(choice),hint='Needle: AlignView penalizes end gaps by default; enable free ends to relax them.' if choice=='needle' else 'Water optimizes a local segment; end-gap penalties do not apply.')
                    if action is None:return
                    if action=='2':self.edit_parameters(choice)
                    elif action=='3':self.preview(lambda:self.s.run_pairwise(choice,reference=ref,preview=True))
                    elif action=='4':self.profiles()
                    else:self.run_action(lambda:self.s.run_pairwise(choice,reference=ref),changes=True);return

    def other_analyses(self):
        with self.ui.section('Other analyses'):
            while True:
                nucleic=sequence_type(self.s.query.seq)=='DNA/RNA'
                choice=self.choose('Sequence analyses',[
                    item('1','Find a motif','IUPAC motif and positions on both strands',enabled=nucleic,reason='IUPAC motif analysis requires a DNA/RNA sequence.'),
                    item('2','Translation and coding regions','Translate, find ORFs or extract a region',enabled=nucleic,reason='Translation requires a DNA/RNA sequence.'),
                    item('3','Composition and profiles','GC, entropy and k-mer frequencies'),
                    item('4','Dotplot','Visualize identical words shared by two sequences',enabled=len(self.s.query_list)>1 or not self.s.single,
                         reason='Load at least two sequences to compare them.'),
                ])
                if choice is None:return
                actions={'1':lambda:self.run_action(self.s.motif,changes=True,artifact='motifs.tsv',title='Motif positions'),
                         '2':self.coding_tools,'3':self.composition_tools,
                         '4':lambda:self.run_action(self.s.make_dotplot,changes=True)}
                self.run_action(actions[choice])

    def coding_tools(self):
        choice=self.choose('Translation and regions',[
            item('1','Translate all six frames',f'Table NCBI {self.s.a.genetic_code} ; export protein sequences'),
            item('2','Find simple ORFs','Start-to-stop segments; exploratory analysis'),
            item('3','Extract a region','Choose coordinates and orientation'),
            item('4','Choose the translation table',f'QUERY: table NCBI {self.s.a.genetic_code} · starts ORF: {self.s.a.orf_starts}'),
        ])
        if choice=='4':self.edit_parameters('translation')
        elif choice:self.run_action({'1':self.s.translate_frames,'2':self.s.find_orfs,'3':self.s.extract_region}[choice],changes=True)

    def composition_tools(self):
        choice=self.choose('Composition',[
            item('1','Overall quality and composition','Length, GC and ambiguous symbols'),
            item('2','Sequence profile','Sliding-window GC, complexity or gaps'),
            item('3','K-mer frequencies','Count exact k-mers'),
        ])
        if choice:self.run_action({'1':self.s.show_stats,'2':self.s.profiles,'3':self.s.composition}[choice],changes=choice!='1')


    def documentation(self):
        from av_docs import documentation_menu
        documentation_menu(self.ui)

    def search_settings(self,mode=None):
        with self.ui.section('Search'):
            while True:
                choice=self.choose('Search settings',[
                    item('1','E-value threshold preset',self.s.a.preset+' ; DIAMOND sensitivity where applicable'),
                    item('2','Maximum number of references',str(self.s.a.hits)+' maximum retained references'),
                    item('3','Maximum E-value',f'{self.s.a.evalue:g} ; plus petit = plus strict'),
                    item('4','CPU count',str(self.s.a.threads)),
                    item('5','Result filters','Identity, coverage and bitscore without rerunning the search'),
                    item('6','Engine parameters','Seeds, masking, scores and sensitivity, depending on the method'),
                ],hint='Filters reuse existing results. Search parameter changes apply to the next calculation.')
                if choice is None:return
                if choice=='1':self.choose_preset()
                elif choice in ('2','3','4'):
                    name,label,cast,minimum={'2':('hits','Maximum number of references',int,1),
                        '3':('evalue','Maximum E-value',float,1e-300),'4':('threads','CPU count',int,1)}[choice]
                    changed=self.s._change_numeric(name,label,cast,minimum)
                    if changed and name=='hits' and self.state.searched:self.s.apply_filters()
                elif choice=='5':self.filters()
                elif choice=='6':
                    method=mode or self.s.mode
                    if method in ('msa','global-dna','global-aa',''):
                        if not self.s.single:method=self.recommended_mode()
                        else:
                            self.ui.page('Choose a method in Compare with references after adding references.','Search method required');continue
                    self.edit_parameters('search',method)

    def alignment_settings(self):
        while True:
            choice=self.choose('Alignment settings',[
                item('1','Choose the engine',self.s.a.aligner.upper()),
                item('2','Calculation parameters',self.parameter_summary('msa')),
                item('3','Calculation profiles','Load or save settings'),
            ])
            if choice is None:return
            if choice=='1':self.choose_aligner()
            elif choice=='2':self.edit_parameters('msa')
            else:self.profiles()

    def edit_parameters(self,context,mode=None,strategy=''):
        mode=mode or self.s.mode
        engine=('MAFFT' if strategy=='fragments' else self.s.a.aligner.upper()) if context=='msa' else mode.upper() if context=='search' else 'FASTTREE' if context=='tree' else context.upper()
        with self.ui.section('Parameters '+engine):
            while True:
                fields=relevant(context,self.s.a,mode,strategy=strategy)
                choice=self.choose('Calculation parameters · '+engine,
                    [item(str(i),p.label,p.display(self.s.a)+' · '+p.effect) for i,p in enumerate(fields,1)]+
                    [item('r','Restore initial values','Reset only this calculation\'s parameters')],
                    hint=(f'Traduction: table NCBI {self.s.a.genetic_code}. ' if mode in ('blastx','diamond-blastx','tblastn','tblastx') else '')+
                        'F1: parameter effect. /: search. Engine default means the option is not supplied and depends on the engine version.')
                if choice is None:return
                if choice=='r':
                    for p in fields:setattr(self.s.a,p.name,p.default)
                    self.s.save_settings();self.ui.notice='Initial values restored for '+engine;continue
                p=fields[int(choice)-1]
                old=getattr(self.s.a,p.name)
                if p.name in ('genetic_code','reference_genetic_code'):
                    from av_genetics import TABLES
                    options=[item(str(identifier),f'Table {identifier} · {data[0]}','NCBI translation code') for identifier,data in TABLES.items()]
                    if p.name=='reference_genetic_code':options.insert(0,item('auto','Use the query table','Inherit the query translation table'))
                    picked=self.choose(p.label,options)
                    if picked is None:continue
                    value=None if picked=='auto' else int(picked)
                elif p.choices or p.kind is bool:
                    values=p.choices if p.choices else (True,False)
                    picked=self.choose(p.label,[item(str(i),('Yes' if value else 'No') if p.kind is bool else str(value),
                        'Initial value' if value==p.default else p.effect) for i,value in enumerate(values,1)]+
                        [item('r','Initial value','Restore '+str(p.default if p.default is not None else 'automatic'))],p.effect)
                    if picked is None:continue
                    value=p.default if picked=='r' else values[int(picked)-1]
                else:
                    answer=self.ui.ask(p.label,str(old) if old is not None else 'default',
                        p.effect+'\nEnter: keep · default: reset · q: cancel')
                    if answer.casefold() in ('q','cancel','cancel'):continue
                    try:
                        value=p.default if answer.casefold() in ('default','defaut','default','native','auto') else p.kind(answer)
                        p.check(value)
                    except ValueError as exc:
                        self.ui.page(str(exc)+'\n\nThe parameter is unchanged.','Invalid value');continue
                setattr(self.s.a,p.name,value);self.s.save_settings()
                self.ui.notice=p.label+': '+p.display(self.s.a)+' ; applies to the next calculation.'

    def profiles(self):
        while True:
            choice=self.choose('Calculation profiles',[
                item('1','Load a JSON profile','Apply validated settings; previous calculations retain their provenance'),
                item('2','Save current settings','Create a JSON profile reusable with --load-profile on Linux or HPC'),
                item('3','View current profile','View settings for all engines and filters'),
            ],hint='A profile contains settings without sequences or executable paths. Loading one does not run an engine.')
            if choice is None:return
            if choice=='3':
                from av_parameters import profile
                self.ui.page(json.dumps(profile(self.s.a),ensure_ascii=False,indent=2),'Profile actuel');continue
            filename=self.ui.ask('File JSON (q to cancel)', 'alignview-profile.json' if choice=='2' else '')
            if not filename or filename.casefold()=='q':continue
            try:
                if choice=='2':write_profile(self.s.a,filename);self.ui.notice='Profile saved: '+filename
                else:
                    old=tuple(getattr(self.s.a,n) for n in ('hits','identity','coverage','bitscore'))
                    load_profile(self.s.a,filename)
                    if self.state.searched and old!=tuple(getattr(self.s.a,n) for n in ('hits','identity','coverage','bitscore')):self.s.apply_filters()
                    else:self.s.save_settings()
                    self.ui.notice='Profile loaded. The new settings apply to the next calculation.'
            except (OSError,ValueError) as exc:self.ui.page(str(exc),'Profile not applied')

    def provenance(self):
        runs=getattr(self.s,'calculation_runs',{})
        text=['Current settings can differ from those that produced the results.',
              'Each successful calculation records its settings and the versions and commands actually executed.',
              'null / native: option not supplied; consult the engine version and its documentation.',
              f'Current tables: QUERY {self.s.a.genetic_code} ; references {self.s.a.reference_genetic_code}.',
              '',json.dumps(runs,ensure_ascii=False,indent=2) if runs else 'No detailed calculation provenance (older project or imported alignment).']
        if self.s.vsearch_results:text+=['',json.dumps({key:result['settings'] for key,result in self.s.vsearch_results.items()},ensure_ascii=False,indent=2)]
        if self.s.tools.history:text+=['',json.dumps(self.s.tools.history,ensure_ascii=False,indent=2)]
        self.ui.page('\n'.join(text),'Completed calculations and provenance')

    def software(self):
        from av_setup import software_hub
        with self.ui.section('Logiciels'):software_hub(self.ui,self.s.tools)

    def help(self):
        from av_docs import read_document
        self.ui.page(read_document('glossary'),'Understanding AlignView')

    def search_results(self):
        with self.ui.section('Search results'):
            while True:
                count=len(self.s.selected)
                choice=self.choose('Matches and coverage',[
                    item('1','Choose a match',f'{count} references ; search and read their segments'),
                    item('2','Query coverage','Detected and unmatched regions',enabled=bool(count),reason='No references retained after filtering.'),
                    item('3','Reference coverage','Detected regions on each reference sequence',enabled=bool(count),reason='No references retained after filtering.'),
                    item('4','Identity and coverage','Compare local similarity and the fraction of the query covered',enabled=bool(count),reason='No references retained after filtering.'),
                    item('5','Segment quality','Mismatches, gaps and reading frames',enabled=bool(count),reason='No references retained after filtering.'),
                ])
                if choice is None:return
                if choice=='1':self.hits()
                elif choice=='3':
                    output=io.StringIO()
                    with redirect_stdout(output):self.s.show_coverage(by_ref=True)
                    self.ui.page(output.getvalue(),'Reference coverage')
                else:self.run_action({'2':self.s.show_coverage,'4':self.s.identity_coverage,'5':self.s.hsp_diagnostics}[choice])

    def alignment_statistics(self):
        with self.ui.section('Alignment statistics'):
            while True:
                count=len(self.s.msa_records or [])
                choice=self.choose('Analyze the alignment',[
                    item('1','Conservation and consensus','Identify conserved, variable and ambiguous columns'),
                    item('2','Sequence distances','Pairwise p-distances, excluding gaps in each pair'),

                ])
                if choice is None:return
                self.run_action({'1':self.s.msa_analysis,'2':self.s.msa_distances}[choice],changes=True)

    def filters(self):
        while True:
            choice=self.choose('Filter matches',[
                item('1','Minimum local identity',f'{self.s.a.identity:g}% per aligned segment'),
                item('2','Minimum query coverage',f'{self.s.a.coverage:g}% using the union of segments'),
                item('3','Bitscore minimal',str(self.s.a.bitscore)),
                item('4','Reset filters','Reset identity, coverage and bitscore to 0'),
            ],hint='Thresholds filter existing matches. A changed selection invalidates the previous derived alignment.')
            if choice is None:return
            if choice=='4':
                self.s.a.identity=self.s.a.coverage=self.s.a.bitscore=0
                changed=True
            else:
                name,label,maximum={'1':('identity','Minimum identity (%)',100),
                    '2':('coverage','Coverage minimale (%)',100),'3':('bitscore','Bitscore minimal',None)}[choice]
                changed=self.s._change_numeric(name,label,float,0,maximum)
            if changed:
                if self.state.searched:self.s.apply_filters()
                else:self.s.save_settings()
                self.ui.notice=f'Filters applied: {len(self.s.selected)} retained references.'

    def choose_preset(self):
        choice=self.choose('E-value threshold preset',[
            item('standard','Standard','E-value 1e-5 · initial setting'),
            item('sensitive','More permissive','E-value 1e-2 ; DIAMOND very-sensitive when used'),
            item('strict','Stricter','E-value 1e-10 ; DIAMOND sensitive when used'),
        ],f'The limit remains {self.s.a.hits} references, regardless of the preset.')
        if not choice:return False
        old=(self.s.a.evalue,self.s.a.sensitivity)
        self.s.a.preset=choice
        for name,value in SEARCH_PRESETS[choice].items():setattr(self.s.a,name,value)
        self.s.save_settings();self.ui.notice=f'Profile {choice} applied. Limit preserved: {self.s.a.hits} references.'
        return old!=(self.s.a.evalue,self.s.a.sensitivity)

    def parameter_summary(self,context,mode=None,strategy=''):
        fields=relevant(context,self.s.a,mode or self.s.mode,strategy=strategy)
        custom=[p.label+': '+p.display(self.s.a) for p in fields if getattr(self.s.a,p.name)!=p.default]
        if custom:return ' · '.join(custom[:2])+(' · other modified settings' if len(custom)>2 else '')
        if context=='msa' and (strategy=='fragments' or self.s.a.aligner=='mafft'):
            return 'MAFFT '+self.s.a.mafft+' · native gap costs · 1000 refinement iterations in global/local'
        return 'AlignView initial settings and engine defaults; open to inspect'

    def preview(self,prepare):
        result=self.run_action(prepare)
        if result:self.ui.page(result,'Commands and calculation preparation')

    def advanced(self):
        with self.ui.section('Advanced'):
            while True:
                choice=self.choose('Advanced',[
                    item('1','Additional analyses','VSEARCH, pairwise comparison, motifs and translation'),
                    item('2','Checks and statistics','Sequence quality, aligned segments, conservation and distances'),
                    item('3','Calculation settings','Search, alignment, translation tables and profiles'),
                    item('4','Tools and provenance','Available software, versions, settings and executed commands'),
                    item('5','HPC calculation','Prepare a Slurm script without submitting a job'),
                    item('6','Help','Documentation and result interpretation'),
                ])
                if choice is None:return
                self.run_action({'1':self.additional_analyses,'2':self.diagnostics,'3':self.calculation_options,
                    '4':self.tools_and_provenance,'5':self.s.slurm_script,'6':self.help_options}[choice])

    def additional_analyses(self):
        while True:
            choice=self.choose('Additional analyses',[
                item('1','VSEARCH analyses','Global search, dereplication, clustering and chimeras',
                     enabled=all(sequence_type(r.seq)=='DNA/RNA' for r in self.s.query_list),reason='VSEARCH requires nucleotide sequences.'),
                item('2','Compare two sequences','Needle: global; Water: local',enabled=not self.s.single,reason='Load a reference file.'),
                item('3','Sequences, motifs and translation','ORFs, genetic codes, extraction, composition and dotplot'),
            ])
            if choice is None:return
            self.run_action({'1':self.vsearch,'2':self.pairwise,'3':self.other_analyses}[choice])
            if self.s.pending:return

    def diagnostics(self):
        while True:
            choice=self.choose('Checks and statistics',[
                item('1','Sequence quality','Lengths, composition and ambiguous characters'),
                item('2','Search details','Coverage, identity and aligned segment quality',
                     enabled=self.state.searched,reason='Run a comparison.'),
                item('3','Alignment statistics','Conservation, consensus and distances',
                     enabled=self.state.aligned,reason='Create or open an alignment.'),
            ])
            if choice is None:return
            self.run_action({'1':self.s.show_stats,'2':self.search_results,'3':self.alignment_statistics}[choice])

    def calculation_options(self):
        while True:
            choice=self.choose('Calculation settings',[
                item('1','Search','Thresholds, filters, reference limit and CPUs'),
                item('2','Alignment','Multiple alignment engine and settings'),
                item('3','Translation and ORFs',f'Table NCBI {self.s.a.genetic_code} · starts {self.s.a.orf_starts}'),
                item('4','Calculation profiles','Load or save reusable settings'),
            ])
            if choice is None:return
            self.run_action({'1':self.search_settings,'2':self.alignment_settings,
                '3':lambda:self.edit_parameters('translation'),'4':self.profiles}[choice])

    def tools_and_provenance(self):
        while True:
            choice=self.choose('Tools and provenance',[
                item('1','Available software','Paths, versions and installation'),
                item('2','Completed calculations and provenance','Actual settings and commands'),
                item('3','Project data','Source files, active query and result status'),
            ])
            if choice is None:return
            self.run_action({'1':self.software,'2':self.provenance,'3':self.s.show_dashboard}[choice])

    def help_options(self):
        while True:
            choice=self.choose('Help',[
                item('1','Understand results','Search, alignment, coverage and identity'),
                item('2','Complete manual','Guides, examples, formats and troubleshooting'),

            ])
            if choice is None:return
            self.run_action({'1':self.help,'2':self.documentation}[choice])

    def alignment_options(self,strategy):
        while True:
            choice=self.choose('Advanced alignment options',[
                item('1','Choose the engine',self.s.a.aligner.upper(),enabled=strategy!='fragments',reason='Adding fragments requires MAFFT.'),
                item('2','Choose sequences',strategy,enabled=not self.s.single,reason='The single file is aligned in full.'),
                item('3','Calculation parameters',self.parameter_summary('msa',strategy=strategy)),
                item('4','Inspect commands','Sequence preparation and engine arguments'),
                item('5','Calculation profiles','Load or save your settings'),
            ])
            if choice is None:return strategy
            if choice=='1':self.choose_aligner()
            elif choice=='2':strategy=self.alignment_scope() or strategy
            elif choice=='3':self.edit_parameters('msa',strategy=strategy)
            elif choice=='4':self.preview(lambda:self.s.msa_run(strategy,preview=True))
            else:self.profiles()

    def search_options(self,mode):
        while True:
            choice=self.choose('Advanced comparison options',[
                item('1','Choose the method',mode.upper()),
                item('2','Search settings','E-value, filters, references and CPUs'),
                item('3','Engine parameters',self.parameter_summary('search',mode)),
                item('4','Inspect commands','Indexing, options and versions'),
                item('5','Calculation profiles','Load or save your settings'),
            ])
            if choice is None:return mode
            if choice=='1':
                mode=self.choose_search_method() or mode
                if mode=='vsearch':return mode
            elif choice=='2':self.search_settings(mode)
            elif choice=='3':self.edit_parameters('search',mode)
            elif choice=='4':self.preview(lambda:self.s.search_preview(mode))
            else:self.profiles()

    def tree_options(self):
        while True:
            choice=self.choose('Advanced tree options',[
                item('1','Calculation parameters',self.parameter_summary('tree')),
                item('2','Inspect the command','FastTree options before calculation'),
                item('3','Calculation profiles','Load or save your settings'),
            ])
            if choice is None:return
            if choice=='1':self.edit_parameters('tree')
            elif choice=='2':self.preview(lambda:self.s.make_tree(preview=True))
            else:self.profiles()
