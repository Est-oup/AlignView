#!/usr/bin/env python3
"""AlignView — explorateur d'alignements ADN/protéines pour Linux et clusters HPC.

Bibliothèque standard Python uniquement. Les logiciels scientifiques sont optionnels et
chargés à la demande depuis PATH / modules HPC. Tout travail est temporaire jusqu'à
confirmation explicite de sauvegarde en fin de session.
"""
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
import textwrap
import traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from av_core import (Seq,HSP,read_fasta,write_fasta,sequence_type,sequence_stats,merge_intervals,
    missing_intervals,interval_length,select_hits,coverage_data,group_hsps,translate_nt,
    reverse_complement,orfs,find_motifs,sliding_scores,msa_metrics,distance_matrix,dotplot,
    sha256,export_table)
from av_tools import Toolchain,CommandError,GROUPS
from av_ui import UI,BLUE,BOLD,GREEN,CYAN,YELLOW,RED,GREY,HomeRequested
from av_visual import coverage_svg,dotplot_svg,plot_svg
from av_report_v6 import build_report_v6
from av_tree import TREE_FIGURE, tree_svg
from av_docs import TOPICS, ALIASES as DOC_ALIASES, read_document
from av_vsearch import ACTIONS as VSEARCH_ACTIONS
from av_i18n import set_language, language, translate as tr
from av_units import units_for, coord_audit
from av_exports import generate_exports,PRESETS
from av_start import startup, StartRequest, suggest_modes, choose_file, choose_sequence
from av_version import VERSION
from av_defaults import DEFAULT_HITS, SEARCH_PRESETS
from av_parameters import (PARAMETERS, add_arguments, validate as validate_parameters,
    snapshot as calculation_settings, load_profile, write_profile, extra_cli)
from av_project import save_project, resume_project
from av_setup import PACKAGES, setup_tools, software_hub, diagnostic
from av_msa_tui import MSAViewer, msa_static

MODES={
    'vsearch':('VSEARCH','Recherche globale de similarité entre séquences nucléotidiques'),
    'blastn':('BLASTN','ADN → ADN : recherche locale nucléotidique'),
    'blastp':('BLASTP','Protéine → protéine : recherche locale'),
    'blastx':('BLASTX','ADN → protéine : 6 cadres de traduction de QUERY'),
    'tblastn':('TBLASTN','Protéine → ADN : 6 cadres côté références'),
    'tblastx':('TBLASTX','ADN → ADN : traduction des 2 côtés'),
    'diamond-blastp':('DIAMOND BLASTP','Recherche protéine → protéine, optimisée pour gros jeux'),
    'diamond-blastx':('DIAMOND BLASTX','ADN → protéine, DIAMOND, QUERY traduite'),
    'mmseqs-blastp':('MMseqs2','Protéine → protéine : recherche locale rapide, HSP inspectables'),
    'global-dna':('GLOBAL ADN','Needle et MAFFT sans recherche de similarité'),
    'global-aa':('GLOBAL PROTÉINES','Needle et MAFFT sans recherche de similarité'),
}
HELP_TEXT={
 'blast':'BLAST/DIAMOND : homologues locaux (HSP). Ils peuvent ignorer les extrémités et des régions divergentes. Une absence de HSP n’est pas une preuve d’absence d’homologie.',
 'coverage':'Couverture = longueur de l’union des intervalles QUERY touchés par les HSP / longueur totale QUERY. Plusieurs HSP qui se chevauchent ne sont pas comptés deux fois. Pour BLASTX les positions QUERY sont nucléotidiques.',
 'identity':'Identité = correspondances exactes / longueur de l’alignement local (HSP). Ce n’est PAS une identité globale de la séquence. Les gaps et substitutions sont distingués.',
 'needle':'EMBOSS Needle implémente un alignement global de type Needleman–Wunsch. Il impose un chemin de bout en bout et permet de voir les extrémités, parfois au prix de gaps artificiels si les séquences sont non colinéaires.',
 'water':'EMBOSS Water implémente Smith–Waterman, alignement local optimal ; ses résultats peuvent différer de BLAST qui est heuristique.',
 'mafft':'MAFFT global (G-INS-i) / local (L-INS-i) sont des stratégies d’alignement multiple, et non des recherches d’homologie. MAFFT --auto choisit suivant les données. ADD-FRAGMENTS ajoute une requête partielle à un alignement de références.',
 'hsp':'MSA des HSP : segments alignables localement, sans les régions flanquantes. AlignView écarte les HSP ne recouvrant pas suffisamment la région ancre. Ne pas interpréter ce MSA comme des séquences intégrales.',
 'dot':'Dotplot = correspondances exactes de k-mers à leurs positions. Une diagonale suggère une colinéarité, plusieurs diagonales peuvent indiquer des répétitions ou duplications ; cela ne mesure pas une identité globale.',
 'orf':'ORF : segment entre un codon de départ et un STOP, dans la table NCBI choisie. Départs ATG par défaut ou initiateurs de la table. Les arrêts contextuels nécessitent une annotation CDS ; ce détecteur simple n’est pas un prédicteur de gènes.',
 'msa':'Sites variables = au moins deux résidus distincts non-gap dans une colonne ; informatifs en parcimonie = au moins deux états représentés chacun ≥ 2 fois. Distances p ignorent les colonnes avec gap dans une paire.',
 'save':'Tous les fichiers calculés sont temporaires. L’enregistrement est proposé une seule fois en fin de session, sauf si --save est indiqué explicitement. L’arrêt brutal de la session peut faire perdre les résultats.',
}


class CLIParser(argparse.ArgumentParser):
    def error(self,message):super().error(tr(message))


def build_parser():
    p=CLIParser(prog='alignview',formatter_class=argparse.RawDescriptionHelpFormatter,
        description='AlignView — exploration, alignements, visualisation et diagnostic en terminal.',
        epilog='Exemples :\n  ./alignview.sh query.fna refs.faa --mode blastx\n  ./alignview.sh sequences.fasta\n  ./alignview.sh query.fna refs.faa --mode diamond-blastx --report --save\n  ./alignview.sh query.fna refs.fna --mode blastn --msa full --save --no-menu')
    section_0=p.add_argument_group('Entrées et reprise')
    section_1=p.add_argument_group('Recherche et filtres')
    section_2=p.add_argument_group('Alignements et analyses')
    section_3=p.add_argument_group('Rapports, exports et sauvegardes')
    section_4=p.add_argument_group('Moteurs scientifiques et installation')
    section_5=p.add_argument_group('Interface et documentation')
    section_0.add_argument('query',nargs='?',help='FASTA de la requête, ou FASTA unique à explorer')
    section_0.add_argument('--query-id',default='',help='Identifiant FASTA à utiliser comme requête quand le premier fichier en contient plusieurs')
    section_0.add_argument('references',nargs='?',help='FASTA des références')
    section_0.add_argument('--mode',choices=['msa',*MODES],help='Programme de recherche ou alignement')
    section_1.add_argument('--task',default='',help='Tâche NCBI compatible avec le programme')
    section_1.add_argument('--evalue',type=float,default=1e-5,help='E-value maximale (1e-5 par défaut)')
    section_1.add_argument('--identity',type=float,default=0,help='Identité locale minimale %% (0)')
    section_1.add_argument('--coverage',type=float,default=0,help='Couverture QUERY par union des HSP %% (0)')
    section_1.add_argument('--bitscore',type=float,default=0,help='Bitscore local minimal (0)')
    section_1.add_argument('--hits',type=int,default=DEFAULT_HITS,help='Nombre maximal de références retenues (100 par défaut)')
    section_1.add_argument('--preset',choices=['standard','sensible','strict'],default='standard',help='Profil de recherche : standard, sensible ou strict ; les options explicites sont prioritaires')
    section_1.add_argument('--threads',type=int,default=4,help='Nombre de threads (4)')
    section_1.add_argument('--sensitivity',choices=['fast','sensitive','more-sensitive','very-sensitive','ultra-sensitive'],default='sensitive')
    section_2.add_argument('--msa',choices=['none','view','align','full','hits','all','global','fragments','hsp'],default='none')
    section_2.add_argument('--mafft',choices=['auto','global','local','fast'],default='auto')
    section_2.add_argument('--aligner',choices=['mafft','muscle','clustalo'],default='mafft',help='Moteur MSA : MAFFT, MUSCLE 5, Clustal Omega')
    section_2.add_argument('--tree',action='store_true',help='Inférer un arbre exploratoire FastTree à partir du MSA calculé')
    section_2.add_argument('--pairwise',choices=['none','global','local'],default='none')
    section_2.add_argument('--ref',type=int,default=1,help='Indice d’une référence (1)')
    section_2.add_argument('--frame',type=int,choices=[-3,-2,-1,1,2,3],help='Forcer un cadre de traduction')
    section_2.add_argument('--gapopen',type=float,default=10)
    section_2.add_argument('--gapextend',type=float,default=.5)
    section_2.add_argument('--matrix',default='',help='Matrice de substitution EMBOSS')
    section_2.add_argument('--free-ends',action='store_true',help='Ne pas pénaliser les gaps terminaux Needle')
    section_3.add_argument('--report',action='store_true',help='Préparer le HTML et ses SVG ; conserver avec --save, --export-dir ou --html-out')
    section_3.add_argument('--html-out',metavar='FICHIER.html',default='',help='Enregistrer un HTML autonome à cet endroit, même sans --save ; refuse d’écraser')
    section_3.add_argument('--no-open',action='store_true',help='Ne pas ouvrir automatiquement le HTML dans le navigateur')
    section_3.add_argument('--preview-port',type=int,default=8765,help=argparse.SUPPRESS)
    section_2.add_argument('--dotplot',action='store_true',help='Générer dotplot QUERY / meilleure référence')
    section_2.add_argument('--kmer',type=int,default=9,help='Longueur des mots du dotplot (9)')
    section_2.add_argument('--motif',default='',help='Motif IUPAC sur la requête')
    section_2.add_argument('--vsearch',choices=['none',*VSEARCH_ACTIONS],default='none',help='Analyse nucléotidique VSEARCH : recherche, déduplication, regroupement ou chimères')
    section_5.add_argument('--list-genetic-codes',action='store_true',help='Afficher les tables de traduction NCBI disponibles, puis quitter')
    section_3.add_argument('--out',default='',help='Dossier du projet avec --save ; chemin relatif au dossier de lancement, créé à la fin')
    section_3.add_argument('--save',action='store_true',help='Sauvegarde explicite à la fin (batch compris)')
    section_5.add_argument('--no-menu',action='store_true',help='Sans menu interactif')
    section_5.add_argument('--simple-ui',action='store_true',help='Menus numérotés classiques')
    section_5.add_argument('--no-color',action='store_true',help='Sans couleurs ANSI')
    section_4.add_argument('--executable',action='append',default=[],metavar='COMMANDE=/chemin',
                   help='Chemin exact vers un outil externe, répétable, ex. blastx=/opt/ncbi/bin/blastx')
    section_4.add_argument('--tools-dir',action='append',default=[],metavar='DOSSIER',
                   help='Ajouter des exécutables locaux au PATH d’AlignView (répétable)')
    section_4.add_argument('--module-policy',choices=['auto','never'],default='auto',
                   help='auto : modules HPC si commande absente ; never : PATH et chemins seulement')
    section_3.add_argument('--export-profile',choices=list(PRESETS),default='standard',
                   help='standard / publication / complete : fichiers préparés en cas de sauvegarde')
    section_3.add_argument('--export-dir',default='',metavar='DOSSIER',
                   help='Exporter explicitement les fichiers tabulaires dans ce dossier, même sans --save')
    section_4.add_argument('--diagnose',action='store_true',help='Diagnostic PATH/modules puis sortie')
    section_5.add_argument('--list-tools',action='store_true',help='Lister les fonctionnalités')
    section_5.add_argument('--docs',nargs='?',const='index',choices=[*TOPICS,*DOC_ALIASES],help='Consulter un guide hors ligne, sans ouvrir de données (défaut : index)')
    section_0.add_argument('--resume',metavar='PROJET',help='Reprendre un projet sauvegardé sans relancer ses calculs')
    section_0.add_argument('--demo',action='store_true',help='Explorer un exemple embarqué sans outil externe')
    section_4.add_argument('--setup',choices=sorted(PACKAGES),help='Préparer l’installation de moteurs via Conda/Mamba')
    section_4.add_argument('--install-tools',action='store_true',help='Exécuter le plan --setup dans un environnement dédié')
    section_4.add_argument('--tools-prefix',help='Préfixe Conda dédié (ou ALIGNVIEW_TOOLS_PREFIX)')
    section_4.add_argument('--diagnose-json',action='store_true',help='Diagnostic machine-readable JSON, sans chargement de modules')
    section_5.add_argument('--version',action='version',version='AlignView '+VERSION)
    section_5.add_argument('--lang',choices=('fr','en'),default='en' if os.environ.get('ALIGNVIEW_LANG')=='en' else 'fr',help='Langue de l’interface, du HTML et des guides (ou ALIGNVIEW_LANG)')
    add_arguments(p)
    for group in p._action_groups:
        group.title=tr(group.title)
        for action in group._group_actions:
            if action.help and action.help!=argparse.SUPPRESS:action.help=tr(action.help)
    p.description=tr(p.description);p.epilog=tr(p.epilog)
    return p


def parse_args(argv=None):
    raw=list(argv if argv is not None else sys.argv[1:])
    selected='en' if os.environ.get('ALIGNVIEW_LANG')=='en' else 'fr'
    for i,token in enumerate(raw):
        if token.startswith('--lang='):selected=token.split('=',1)[1]
        elif token=='--lang' and i+1<len(raw):selected=raw[i+1]
    set_language(selected if selected in ('fr','en') else 'fr')
    p=build_parser()
    a=p.parse_args(argv)
    source=list(argv if argv is not None else sys.argv[1:])
    explicit={p._option_string_actions[token.split('=',1)[0]].dest for token in source
        if token.startswith('--') and token.split('=',1)[0] in p._option_string_actions}
    a._explicit=explicit
    if a.load_profile:
        try:load_profile(a,a.load_profile,explicit)
        except ValueError as exc:p.error(str(exc))
    presets=SEARCH_PRESETS
    flags={'evalue':'--evalue','sensitivity':'--sensitivity'}
    for key,value in presets[a.preset].items():
        if (not a.load_profile or 'preset' in explicit) and not any(x==flags[key] or x.startswith(flags[key]+'=') for x in source):setattr(a,key,value)
    a._profile_overrides={name:getattr(a,name) for name in explicit if name in {p.name for p in PARAMETERS}|{'hits','identity','coverage','bitscore','preset'}}
    if a.load_profile and 'preset' in explicit:a._profile_overrides.update(evalue=a.evalue,sensitivity=a.sensitivity)
    try:validate_parameters(a)
    except ValueError as exc:p.error(str(exc))
    # Sans fichiers, l'interface démarre sur l'écran d'accueil (sauf en batch).
    if not a.query and a.no_menu and not (a.diagnose or a.diagnose_json or a.list_tools or a.list_genetic_codes or a.docs or a.setup or a.resume or a.demo or a.write_profile):p.error('Préciser un FASTA en mode --no-menu')
    if a.mode=='vsearch' and a.vsearch=='none':a.vsearch='search'
    if a.vsearch=='search' and a.mode is None and not a.resume:a.mode='vsearch'
    if a.vsearch!='none' and a.mode not in (None,'msa','vsearch'):p.error('Utiliser --vsearch sans mode BLAST/DIAMOND ; les analyses ont des métriques différentes.')
    if not all(math.isfinite(x) for x in (a.evalue,a.identity,a.coverage,a.bitscore,a.gapopen,a.gapextend)) or a.ref<1 or a.gapopen<0 or a.gapextend<0 or a.evalue<=0 or a.hits<1 or a.threads<1 or a.kmer<2 or a.kmer>60 or not 0<=a.identity<=100 or not 0<=a.coverage<=100 or a.bitscore<0 or not 0<=a.preview_port<=65535:
        p.error('Valeurs invalides : evalue > 0, hits/threads > 0, 2 ≤ kmer ≤ 60, identité/couverture entre 0 et 100')
    if a.install_tools and not a.setup:p.error('--install-tools nécessite --setup PROFIL')
    if a.resume and (a.query or a.references or a.demo):p.error('--resume ne se combine pas avec des fichiers ou --demo')
    if a.resume and explicit & {'mode','query_id'}:p.error('Pour changer de méthode ou de séquence après --resume, utiliser les menus Analyser / Données.')
    if a.html_out and Path(a.html_out).suffix.lower()!='.html':p.error('--html-out attend un nom de fichier avec l’extension .html')
    if a.references and a.mode=='msa':p.error('--mode msa uniquement avec un seul FASTA')
    if a.query and not a.references and a.mode not in (None,'msa'):p.error('Ce mode nécessite un FASTA de références.')
    if a.query and not a.references and a.msa in ('hits','hsp','fragments'):p.error('Cette stratégie MSA nécessite un FASTA de références.')
    if a.references and a.msa=='view':p.error('--msa view nécessite un seul FASTA déjà aligné.')
    return a


def show_file(file:Path,ui:UI,max_lines=160):
    if not file.exists():ui.warn('Fichier absent');return
    ui.page(file.read_text(errors='replace'),file.name)


def text_segments(spans):return ', '.join(f'{a:,}–{b:,} ({b-a+1:,})' for a,b in spans) if spans else 'aucune'


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
        self.dirty=True;self.last_saved=None;self.search_completed=False;self.search_parameters={};self.calculation_runs={}
        self.query_path=Path(args.query).expanduser().resolve()
        self.refs_path=Path(args.references).expanduser().resolve() if args.references else None
        for path in (self.query_path,self.refs_path):
            if path and not path.is_file():raise ValueError(f'Fichier introuvable : {path}')
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
                raise ValueError(f'Identifiant QUERY introuvable : {args.query_id}')
            if candidate is None and (args.no_menu or not ui.tty):
                raise ValueError('FASTA QUERY multiple : préciser --query-id en mode non interactif.')
            if candidate is None:
                labels=[(str(i),r.id,f'{len(r.seq):,} caractères · {sequence_type(r.seq)}','SÉQUENCES') for i,r in enumerate(self.query_list,1)]
                choice=ui.menu('Choisir la requête active',labels,
                               f'{len(self.query_list)} séquences · analyses ciblées sur une requête')
                if choice is None:raise KeyboardInterrupt
                candidate=self.query_list[int(choice)-1]
            self.query=candidate
            self.query_search_path=self.work/'query_active.fasta'
            write_fasta([self.query],self.query_search_path)
            ui.info(f'Requête sélectionnée : {self.query.id} · {len(self.query.seq):,} caractères')
        # Normalized plain FASTA inputs also support gzip sources for every engine.
        self.query_search_path=self.work/'query_active.fasta'
        write_fasta([self.query] if not self.single else self.query_list,self.query_search_path)
        self.refs_search_path=self.work/'references_active.fasta' if self.refs_path else None
        if self.refs_search_path:write_fasta(self.references,self.refs_search_path)
        self.tools=Toolchain(self.work,ui,auto_modules=args.module_policy=='auto',executable_paths=args.executable,tool_dirs=args.tools_dir)
        self.hsps=[];self.selected={};self.msa_records=None;self.msa_file=None;self.vsearch_results={}
        self.dot_data=None;self.motif_rows=[];self.orf_rows=[]
        self.created=[];self.mode=args.mode or ('msa' if self.single else '')
        self.start=datetime.datetime.now(datetime.timezone.utc).isoformat()
        self.check_type()
        self.save_settings()

    def close(self):
        self.temporary.cleanup()

    def check_type(self):
        qtype=sequence_type(self.query.seq)
        if self.single:return
        ref_type=collections.Counter(sequence_type(r.seq) for r in self.references).most_common(1)[0][0]
        types={
            'vsearch':('ADN/ARN','ADN/ARN'),
            'blastn':('ADN/ARN','ADN/ARN'),'blastp':('protéine','protéine'),
            'blastx':('ADN/ARN','protéine'),'mmseqs-blastp':('protéine','protéine'),'tblastn':('protéine','ADN/ARN'),
            'tblastx':('ADN/ARN','ADN/ARN'),'diamond-blastp':('protéine','protéine'),
            'diamond-blastx':('ADN/ARN','protéine'),'global-dna':('ADN/ARN','ADN/ARN'),
            'global-aa':('protéine','protéine')}
        if self.mode in types:
            want=types[self.mode]
            if (qtype,ref_type)!=want:
                self.ui.warn(f'Types FASTA détectés : QUERY={qtype}, REFS≈{ref_type} ; attendus pour {self.mode} : {want[0]}/{want[1]}.')
                if self.mode.startswith('global-') and qtype!='inconnu' and ref_type!='inconnu':
                    raise ValueError('MAFFT global : types ADN/protéine incompatibles ; utiliser un FASTA homogène ou traduire les séquences.')

    def settings(self):
        keys=('mode','task','evalue','identity','coverage','bitscore','hits','threads','sensitivity',
              'preset','msa','mafft','aligner','tree','frame','gapopen','gapextend','matrix','free_ends','kmer')
        d={k:getattr(self.a,k) for k in keys};d['mode']=self.mode
        d.update({p.name:getattr(self.a,p.name,p.default) for p in PARAMETERS})
        d['native_defaults']='Options non transmises : défauts de la version du moteur enregistrée dans commands.json.'
        d['calculation_runs']=self.calculation_runs
        d['language']=language()
        d['vsearch_analyses']={action:{k:value for k,value in result.items() if k!='preview'} for action,result in self.vsearch_results.items()}
        u=units_for(self.mode,sequence_type(self.query.seq))
        d['units']={'query_coordinates':u.query,'reference_coordinates':u.subject,'alignment_columns':u.aligned,'msa_columns':u.msa}
        d['timestamp_utc']=self.start;d['query']=str(self.query_path);d['active_query_id']=self.query.id
        d['references']=str(self.refs_path) if self.refs_path else None
        d['search_completed']=self.search_completed
        d['search_run']=dict(self.search_parameters)
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
            f.write('# Query BED 0-based half-open; BlastX query coordinates = nucléotides\n')
            for h in [h for g in self.selected.values() for h in g]:
                a,b=h.qspan; f.write(f'{h.query}\t{a-1}\t{b}\t{h.subject}\t{min(1000,round(h.bitscore))}\t{"-" if h.qstart>h.qend else "+"}\n')
        with (self.work/'hsp_reference.bed').open('w') as f:
            f.write('# Reference BED 0-based half-open; TBLASTN/TBLASTX reference coordinates = nucléotides\n')
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
            if not 1<=n<=len(labels):raise ValueError('Indice de référence hors plage : '+str(n))
            return self.refmap[labels[n-1]],n
        options=[]
        for i,name in enumerate(labels,1):
            group=self.selected.get(name,[])
            detail=f'{group[0].identity:.1f}% · E={group[0].evalue:.2g}' if group else f'{len(self.refmap[name].seq):,} caractères'
            options.append((str(i),name,detail,'RÉFÉRENCES'))
        choice=self.ui.menu('CHOISIR UNE RÉFÉRENCE',options,f'{len(labels):,} références · / rechercher')
        if choice is None:return None,0
        n=int(choice);self.a.ref=n
        return self.refmap[labels[n-1]],n

    def show_dashboard(self):
        self.ui.banner('SESSION  /  paramètres et données en mémoire')
        self.ui.heading('ÉTAT DE LA SESSION')
        for k,v in [
            ('Requête',f'{self.query_path.name} · {len(self.query_list)} séq. · {len(self.query.seq):,} positions'),
            ('Références',f'{self.refs_path.name} · {len(self.references):,} séq.' if self.refs_path else 'aucune — un FASTA'),
            ('Programme',self.mode.upper()),('Seuils',f'E ≤ {self.a.evalue:g} · id ≥ {self.a.identity}% · couverture QUERY ≥ {self.a.coverage}%'),
            ('Recherche',f'{len(self.hsps)} HSP · {len(self.selected)} références retenues'),
            ('Alignement',self.msa_file.name if self.msa_file else 'non généré'),
            ('Dotplot','généré' if self.dot_data else 'non généré'),
            ('Dossier de lancement',str(Path.cwd())),
            ('Sauvegarde','modifications à enregistrer' if self.dirty else 'projet enregistré')]:
            print('  '+self.ui.c(GREY,f'{tr(k):<23}')+tr(v))
        self.ui.info('Modules chargés : '+(', '.join(m for _,m in self.tools.loaded) or tr('aucun, PATH utilisé')))
        if self.last_saved:self.ui.info('Dernière copie enregistrée : '+str(self.last_saved))


    def _change_numeric(self,key,label,cast,minimum=None,maximum=None):
        previous=getattr(self.a,key)
        answer=self.ui.ask(label,str(previous),'Entrée = conserver · q = annuler')
        if answer.lower() in ('q','annuler'):return False
        try:
            value=cast(answer)
            if not math.isfinite(float(value)):raise ValueError
            if minimum is not None and value<minimum:raise ValueError
            if maximum is not None and value>maximum:raise ValueError
        except ValueError:
            self.ui.warn('Valeur invalide : paramètre inchangé.');return False
        setattr(self.a,key,value);self.save_settings()
        self.ui.success(f'{label} : {previous} → {value}')
        return value!=previous

    def choose_preset(self):
        from av_navigation import Workbench
        return Workbench(self).choose_preset()

    def settings_hub(self,initial=False):
        from av_navigation import Workbench
        return Workbench(self).search_settings()


    def search(self):
        if self.mode.startswith('global-') or self.single:return
        self.ui.heading('RECHERCHE / EXÉCUTION')
        validate_parameters(self.a,context='search',mode=self.mode)
        max_target=self.a.search_candidates or max(100,self.a.hits*10)
        self.ui.info(f'{self.mode.upper()} · {len(self.references)} références · {self.a.threads} threads')
        self._calculation_start=len(self.tools.history)
        self.hsps=self.tools.search(self.mode,self.query_search_path,self.refs_search_path,self.a.threads,self.a.evalue,
                                    max_target,self.a.sensitivity,self.a.task,options=self.a)
        self.search_completed=True
        self.search_parameters=calculation_settings('search',self.a,self.mode)
        self.record_calculation('search',self.search_parameters)
        self.apply_filters()
        self.ui.success(f'{len(self.hsps)} HSP trouvés · {len(self.selected)} références retenues')

    def search_preview(self,mode=None):
        mode=mode or self.mode
        commands=self.tools.search_commands(mode,self.query_search_path,self.refs_search_path,self.a.threads,
            self.a.evalue,self.a.search_candidates or max(100,self.a.hits*10),self.a.sensitivity,self.a.task,self.a)
        return self.tools.describe(commands)

    def record_calculation(self,kind,parameters,**details):
        self.calculation_runs[kind]=dict(parameters=parameters,**details,
            commands=[dict(row) for row in self.tools.history[self._calculation_start:]])
        self.save_settings()

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
            self.ui.warn('Aucune recherche locale : choisir BLAST / DIAMOND pour obtenir des hits.');return
        rows=self.coverage_rows();self.ui.heading(f'RÉSULTATS  /  {len(rows)} RÉFÉRENCES')
        print(f'  {"#":>2}  {"Référence":<30} {"Id.%":>7} {"Cov.%":>7} {"HSP":>5} {"E-value":>11} {"Bits":>8}')
        for i,r in enumerate(rows,1):
            print(f'  {i:>2}  {r["subject"][:30]:<30} {r["best_identity"]:7.1f} {r["query_coverage"]:7.1f} {r["hsps"]:5d} {r["best_evalue"]:11.2g} {r["best_bitscore"]:8.1f}')
        if not rows:self.ui.warn('Aucun hit après filtrage. Réduire les seuils ou changer de méthode.')
        units=units_for(self.mode,sequence_type(self.query.seq))
        self.ui.info(f'Identité = % sur HSP ; couverture = union sur QUERY ({units.query}) ; alignement HSP ({units.aligned}).')

    def show_coverage(self,by_ref=False):
        if not self.selected:self.ui.warn('Aucun HSP retenu ; couverture indisponible.');return
        rows=self.coverage_rows();q_len=len(self.query.seq)
        self.ui.heading('CARTES DE COUVERTURE / RÉFÉRENCES' if by_ref else 'CARTE DE COUVERTURE / QUERY')
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
            print('     '+self.ui.c(GREY,'Sans HSP : '+text_segments(miss)))
        units=units_for(self.mode,sequence_type(self.query.seq))
        self.ui.info(f'█ HSP · gris = sans HSP · QUERY en {units.query} · RÉFÉRENCE en {units.subject} · colonnes alignées en {units.aligned}')
        if 'blastx' in self.mode or self.mode=='tblastx':
            self.ui.warn('BLASTX/TBLASTX : QUERY en nucléotides, alignements internes en acides aminés.')
        if by_ref and self.mode in ('tblastn','tblastx'):
            self.ui.warn('TBLASTN/TBLASTX : coordonnées référence en nucléotides.')

    def inspect_hit(self):
        if not self.selected:self.ui.warn('Aucun HSP à explorer.');return
        self.ui.heading('INSPECTEUR / RÉFÉRENCE ET HSP')
        ref,_=self.ref_select()
        if ref is None:return
        rows=self.selected.get(ref.id,[])
        print(); self.ui.info(f'{ref.id} · longueur REF {len(ref.seq):,} · {len(rows)} HSP')
        for i,h in enumerate(rows,1):
            print(f'  {self.ui.c(CYAN,str(i))}  Q {h.qstart}–{h.qend}  /  S {h.sstart}–{h.send} '
                  f'· {h.identity:.1f}% identiques · E={h.evalue:.2g} · bits={h.bitscore:.1f} · frame Q {h.qframe}/S {h.sframe}')
        if not rows:return
        val=self.ui.ask('Numéro de HSP (0 = tous)','0')
        try:ix=int(val)
        except ValueError:ix=0
        if ix>len(rows) or ix<0:self.ui.warn('Index incorrect');return
        import contextlib,io
        output=io.StringIO()
        with contextlib.redirect_stdout(output):
            for h in rows if ix==0 else [rows[ix-1]]:self.print_alignment(h)
        self.ui.page(output.getvalue(),f'Alignements HSP / {ref.id}')

    def print_alignment(self,h:HSP,width=64):
        self.ui.heading(f'HSP  /  {h.subject}  ·  {h.identity:.1f}% identité')
        self.ui.info(f'QUERY {h.qstart} → {h.qend} · REF {h.sstart} → {h.send} · {h.align_length} colonnes · E={h.evalue:g}')
        self.ui.info(f'{h.matches} identiques · {h.mismatches} substitutions · {h.gap_columns} colonnes avec gap · cadre Q={h.qframe}, S={h.sframe}')
        qa=h.qseq;sa=h.sseq
        if not qa or not sa:self.ui.warn('L’alignement HSP ne contient pas qseq/sseq.');return
        nt_q=self.mode in ('blastx','diamond-blastx','tblastx')
        nt_s=self.mode in ('tblastn','tblastx')
        qstep=3 if nt_q else 1;sstep=3 if nt_s else 1
        qsign=1 if h.qend>=h.qstart else -1;ssign=1 if h.send>=h.sstart else -1
        qp=h.qstart;sp=h.sstart
        for i in range(0,min(len(qa),len(sa)),width):
            q=qa[i:i+width];s=sa[i:i+width]
            qn=sum(c!='-' for c in q);sn=sum(c!='-' for c in s)
            endq=qp+qsign*(qn*qstep-1) if qn else qp
            ends=sp+ssign*(sn*sstep-1) if sn else sp
            marks=''.join('|' if x==y and x!='-' else ' ' if x=='-' or y=='-' else '.' for x,y in zip(q,s))
            if self.ui.color:
                cq=''.join(self.ui.c(GREEN,c) if c==t and c!='-' else self.ui.c(RED,c) if c=='-' or t=='-' else self.ui.c(YELLOW,c) for c,t in zip(q,s))
                cs=''.join(self.ui.c(GREEN,t) if c==t and c!='-' else self.ui.c(RED,t) if c=='-' or t=='-' else self.ui.c(YELLOW,t) for c,t in zip(q,s))
            else:cq=q;cs=s
            print(f'  Q {qp:>8}  {cq}  {endq}')
            print(f'  {"":>11}  {marks}')
            print(f'  S {sp:>8}  {cs}  {ends}\n')
            qp+=qsign*qn*qstep;sp+=ssign*sn*sstep
        self.ui.info('Vert = identité, jaune = substitution, rouge = gap. Coordonnées locales aux HSP.')

    def hsp_diagnostics(self):
        self.ui.heading('DIAGNOSTIC / SUBSTITUTIONS, GAPS ET CADRES')
        for ref,items in self.selected.items():
            groups=collections.Counter(h.qframe for h in items)
            spans=merge_intervals([h.qspan for h in items],len(self.query.seq))
            total_cols=sum(len(h.qseq) for h in items)
            gaps=sum(h.gap_columns for h in items);subs=sum(h.mismatches for h in items)
            self.ui.info(f'{ref} · {len(items)} HSP · {interval_length(spans)} positions QUERY uniques')
            print(f'      {total_cols} colonnes cumulées (HSP, chevauchements possibles), {subs} substitutions, {gaps} colonnes gap ; cadres {dict(groups)}')
        self.ui.warn('La somme des substitutions sur HSP se chevauchant peut compter certaines positions plusieurs fois.')

    def identity_coverage(self):
        self.ui.heading('DIAGRAMME / IDENTITÉ LOCALE VS COUVERTURE QUERY')
        rows=self.coverage_rows()
        if not rows:self.ui.warn('Aucun hit');return
        W=56;H=16;board=[[' ' for _ in range(W+1)] for _ in range(H+1)]
        for i in range(H+1):board[i][0]='│'
        for j in range(W+1):board[H][j]='─'
        board[H][0]='└'
        for i,r in enumerate(rows,1):
            x=min(W,max(1,round(r['query_coverage']/100*W)))
            y=min(H-1,max(0,H-1-round(r['best_identity']/100*(H-1))))
            board[y][x]=str(i%10)
        for row in board:print('     '+''.join(row))
        print('     0% couverture QUERY'+' '*(W-23)+'100%')
        for i,r in enumerate(rows,1):print(f'    {i:>2}. {r["subject"]}: {r["best_identity"]:.1f}% identité HSP / {r["query_coverage"]:.1f}% couverture')
        self.ui.warn('Ces deux métriques ne sont pas un pourcentage global d’identité de séquence.')

    def run_pairwise(self,algorithm,reference=None,preview=False):
        if self.single:self.ui.warn('Fournir un second FASTA pour Needle/Water.');return
        self.ui.heading(f'{algorithm.upper()} / ALIGNEMENT PAIRWISE')
        ref=reference
        if ref is None:ref,_=self.ref_select()
        if ref is None:return
        rows=self.selected.get(ref.id,[])
        h=rows[0] if rows else None
        qframe=self.a.frame or (h.qframe if h else 1)
        sframe=self.a.sframe or (h.sframe if h and h.sframe else 1)
        if self.mode in ('blastx','diamond-blastx','tblastx','tblastn'):
            self.ui.warn(f'Traduction de la séquence entière : Q frame {qframe}, REF frame {sframe} si applicable.')
        prefix=self.work/f'{algorithm}_{re.sub("[^A-Za-z0-9_.-]","_",ref.id)}'
        self._calculation_start=len(self.tools.history)
        result=self.tools.emboss(algorithm,self.query,ref,prefix,mode=self.mode,
                        gapopen=self.a.gapopen,gapextend=self.a.gapextend,
                        free_ends=self.a.free_ends,matrix=self.a.matrix,frame=qframe,sframe=sframe,options=self.a,preview=preview)
        if preview:
            translated=self.mode in ('blastx','diamond-blastx','tblastx','tblastn')
            return (f'Séquences : {self.query.id} / {ref.id}.\n'+
                (f'Traduction complète : QUERY cadre {qframe}, table {self.a.genetic_code} ; référence cadre {sframe}, table {self.a.reference_genetic_code or self.a.genetic_code} ; STOP conservés.\n' if translated else '')+result)
        file,_=result
        self.a.pairwise='global' if algorithm=='needle' else 'local'
        self.record_calculation(algorithm,calculation_settings(algorithm,self.a,self.mode),
            reference=ref.id,query_frame=qframe,reference_frame=sframe,genetic_code=self.a.genetic_code,
            reference_genetic_code=self.a.reference_genetic_code or self.a.genetic_code)
        self.dirty=True
        show_file(file,self.ui,max_lines=180)
        self.created.append(str(file.name))

    def msa_run(self,strategy,preview=False):
        engine='mafft' if strategy=='fragments' else self.a.aligner
        effective=argparse.Namespace(**vars(self.a));effective.aligner=engine
        validate_parameters(effective,context='msa')
        self.ui.heading(f'{engine.upper()} / ALIGNEMENT MULTIPLE')
        self._calculation_start=len(self.tools.history)
        self._msa_preparation={}
        pending=self.work/'pending_msa';pending.mkdir(exist_ok=True)
        if self.single:
            selected=self.query_list
            if len(selected)<2:self.ui.warn('Un MSA requiert au moins deux séquences.');return
            input_file=self.query_search_path
        else:
            if strategy=='hsp':return self._msa_hsps(preview=preview)
            if self.mode in ('tblastn','tblastx'):
                self.ui.warn('MAFFT full non proposé pour TBLASTN/TBLASTX : plusieurs cadres/ORF possibles. Extraire les ORF protéiques au préalable.');return
            selected=self.references if strategy=='all' or self.mode.startswith('global-') or (strategy=='fragments' and not self.selected) else [self.refmap[k] for k in self.selected]
            if not selected:self.ui.warn('Aucune référence retenue ; utiliser all ou réduire les seuils.');return
            if self.mode in ('blastx','diamond-blastx'):
                frames={h.qframe for v in self.selected.values() for h in v if h.qframe}
                frame=self.a.frame or (next(iter(self.selected.values()))[0].qframe if self.selected else None)
                if frame is None:
                    frame=1
                    self.ui.warn('Aucun cadre détecté ou forcé : utilisation explicite du cadre QUERY +1 (--frame pour changer).')
                if frame not in (-3,-2,-1,1,2,3):self.ui.warn('Cadre invalide');return
                if len(frames)>1:self.ui.warn(f'Cadres différents {sorted(frames)} : seule QUERY frame {frame} est utilisée.')
                tr=translate_nt(self.query.seq,frame,self.a.genetic_code).replace('*','X')
                self._msa_preparation.update(query_frame=frame,genetic_code=self.a.genetic_code,stop_handling='STOP remplacés par X')
                query_record=Seq('QUERY',f'QUERY_frame_{frame}',tr)
                self.ui.info('QUERY ADN traduite en protéine (STOP remplacés par X pour MAFFT).')
            else:query_record=Seq('QUERY','QUERY',self.query.seq)
            ref_records=[Seq('REF_'+str(i),f'REF_{i} {r.id}',r.seq) for i,r in enumerate(selected,1)]
            input_file=self.work/'msa_input.fasta'
            if strategy=='fragments':
                commands=[]
                if len(ref_records)==1:
                    ref_aln=self.work/'msa_reference_alignment.fasta';write_fasta(ref_records,ref_aln)
                else:
                    refs_file=self.work/'msa_refs.fasta';write_fasta(ref_records,refs_file)
                    ref_aln=self.work/'msa_reference_alignment.fasta'
                    commands.append(self.tools.mafft_command(refs_file,ref_aln,self.a.mafft,self.a.threads,options=self.a))
                frag=self.work/'msa_fragment.fasta';write_fasta([query_record],frag)
                out=pending/'alignment_fragments.fasta'
                commands.append(self.tools.mafft_command(ref_aln,out,'auto',self.a.threads,fragments=frag,options=self.a))
                if preview:return self.msa_preparation_text()+'\nAligner les références, puis ajouter QUERY ; moteur MAFFT imposé par addfragments.\n'+self.tools.describe(commands)
                self.tools.execute(commands)
                return self._set_msa(out,'fragments')
            write_fasta([query_record,*ref_records],input_file)
        if self.a.aligner=='mafft' and self.a.mafft in ('global','local') and len(selected)>200:
            self.ui.warn('G-INS-i/L-INS-i coûteux > 200 séquences : utiliser --mafft auto.');return
        out=pending/('alignment_'+strategy+'.fasta')
        if preview:return self.msa_preparation_text()+'\n'+self.tools.describe([self.tools.msa_command(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)])
        self.tools.msa(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)
        return self._set_msa(out,strategy)

    def msa_preparation_text(self):
        if self._msa_preparation.get('query_frame'):
            return f'Préparation : QUERY traduite en cadre {self._msa_preparation["query_frame"]}, table NCBI {self.a.genetic_code}, STOP remplacés par X.'
        return 'Préparation : séquences intégrales ; aucune traduction supplémentaire.'

    def _msa_hsps(self,preview=False):
        if not self.selected:self.ui.warn('Aucun HSP pour MAFFT HSP.');return
        anchor=next(iter(self.selected.values()))[0]
        self._msa_preparation=dict(anchor_subject=anchor.subject,anchor_query_span=anchor.qspan,
            hsp_gap_handling='Gaps retirés avant réalignement',stop_handling='STOP remplacés par X',
            inherited_search_parameters=dict(self.search_parameters))
        anchor_intervals=anchor.qspan
        kept=[Seq('QUERY_HSP',f'QUERY_HSP_{anchor.qstart}_{anchor.qend}',anchor.qseq.replace('-','').replace('*','X'))]
        manifest=[]
        for name,group in self.selected.items():
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
        if len(kept)<2:self.ui.warn(f'Aucun HSP ne recouvre au moins {self.a.hsp_overlap:.0%} de la région ancre.');return
        self.ui.info(f'{len(kept)-1} HSP compatibles, {len(manifest)-(len(kept)-1)} exclus.')
        input_file=self.work/'msa_hsp_input.fasta';out=self.work/'pending_msa'/'alignment_hsp.fasta'
        write_fasta(kept,input_file)
        if preview:return f'Préparation : ancre {anchor.subject} Q {anchor.qstart}–{anchor.qend} ; seuil {self.a.hsp_overlap:.0%} ; {len(kept)} séquences ; gaps retirés et STOP → X.\n'+self.tools.describe([self.tools.msa_command(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)])
        self.tools.msa(input_file,out,engine=self.a.aligner,strategy=self.a.mafft,threads=self.a.threads,options=self.a)
        self.ui.warn('MSA exploratoire : les HSP peuvent représenter des domaines non équivalents.')
        return self._set_msa(out,'hsp')

    def _set_msa(self,out,strategy):
        records=read_fasta(out)
        if len({len(record.seq) for record in records})!=1:
            raise CommandError('Le moteur a produit des séquences de longueurs différentes.')
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
        self.ui.success(f'MSA {engine}/{strategy} · {len(self.msa_records)} séquences · {len(self.msa_records[0].seq)} colonnes')
        self.show_msa_window();return out

    def load_msa(self):
        if self.msa_records:return self.msa_records
        if self.single and (len(self.query_list)>=2 or self.a.msa=='view') and len({len(r.seq) for r in self.query_list})==1:
            # Même longueur != preuve d'alignement. Le mode --msa view est une affirmation explicite.
            if self.a.msa!='view' and self.ui.tty and not self.a.no_menu:
                self.ui.warn('Toutes les séquences ont la même longueur ; cela ne prouve pas qu’elles sont déjà alignées.')
                if not self.ui.yesno('Ce FASTA est-il déjà aligné ?',False):
                    self.ui.info('Choisir Analyser → Créer un alignement pour aligner ce fichier.')
                    return None
            self.msa_records=self.query_list;self.msa_file=self.query_path
            self.dirty=True
            return self.msa_records
        self.ui.warn('Aucun MSA disponible ; aligner les séquences avec MAFFT.')
        return None

    def show_msa_window(self):
        """Visionneuse plein écran : défilement bi-directionnel, sites variables, recherche."""
        seqs=self.load_msa()
        if not seqs:return
        ncol=len(seqs[0].seq)
        if any(len(r.seq)!=ncol for r in seqs):
            self.ui.warn('Longueurs différentes : FASTA non aligné. Utiliser MAFFT au préalable.');return
        if not self.ui.tty:
            print(tr(msa_static(seqs,width=75,limit=40,color=False)))
            return
        if self._viewer is None or self._viewer.seqs is not seqs:
            self._viewer=MSAViewer(seqs,color=self.ui.color)
            self._viewer.row=next((i for i,record in enumerate(seqs) if record.id==self.query.id),0)
        result=self._viewer.run()
        if result=='home':raise HomeRequested
        if result=='web':
            self.make_report(open_preview=True)
        elif result=='export':
            self.copy_report()

    def msa_analysis(self):
        seqs=self.load_msa()
        if not seqs:return
        metric=msa_metrics(seqs)
        self.ui.heading('CONSERVATION / SITES VARIABLES')
        self.ui.info(f'{metric["n"]} séquences, {metric["length"]:,} colonnes, {metric["variable"]:,} variables')
        self.ui.info(f'{metric["parsimony_informative"]:,} informatives en parcimonie, {metric["gapped"]:,} avec gap, {metric["invariant"]:,} strictement invariantes')
        export_table(self.work/'msa_positions.tsv',metric['positions'])
        write_fasta([Seq('CONSENSUS','CONSENSUS',metric['consensus'])],self.work/'msa_consensus.fasta')
        scores=[dict(center=r['position'],value=r['conservation']) for r in metric['positions']]
        (self.work/'msa_conservation.svg').write_text(plot_svg(scores,'Conservation MSA (%)'))
        for r in metric['positions'][:20]:
            print(f'  {r["position"]:>5}  {r["consensus"]}  conserv. {r["conservation"]:5.1f}% · gap {r["gap_pct"]:5.1f}% · entropie {r["entropy"]:.2f}')
        if len(metric['positions'])>20:self.ui.info('Toutes les colonnes dans msa_positions.tsv.')

    def msa_distances(self):
        seqs=self.load_msa()
        if not seqs:return
        matrix=distance_matrix(seqs);labels=[s.id for s in seqs]
        with (self.work/'msa_pdistance.tsv').open('w') as f:
            f.write('sequence\t'+'\t'.join(labels)+'\n')
            for name,row in zip(labels,matrix):f.write(name+'\t'+'\t'.join(f'{v:.6g}' if v is not None else 'NA' for v in row)+'\n')
        self.ui.heading('DISTANCES P / GAPS EXCLUS')
        self.ui.info('Les numéros des lignes et colonnes renvoient à ces séquences :')
        for i,name in enumerate(labels,1):print(f'  {i:>3}  {name}')
        print('\n       '+' '.join(f'{i:>7}' for i in range(1,len(labels)+1)))
        for i,row in enumerate(matrix,1):
            print(f'  {i:>3}  '+' '.join('     NA' if value is None else f'{value:7.3f}' for value in row))
        self.ui.info('Distance p non corrigée ; les colonnes à gap sont exclues par paire. Maximum 250 séquences.')

    def vsearch_run(self,action,preview=False):
        from av_vsearch import commands, summarize, publish_directory
        source=[self.query] if action=='search' else self.query_list
        references=self.references if action=='search' else []
        if action=='search' and not references:raise ValueError('VSEARCH recherche : charger un FASTA de références.')
        if any(sequence_type(record.seq)!='ADN/ARN' for record in [*source,*references]):
            raise ValueError('VSEARCH nécessite des séquences nucléotidiques ; les protéines ne sont pas prises en charge.')
        def normalized(records):
            return [Seq(r.id,r.id,r.seq.upper().replace('-','').replace('.','').replace('U','T')) for r in records]
        normalized_source=normalized(source)
        eligible=sum(self.a.vsearch_min_length<=len(r.seq)<=self.a.vsearch_max_length for r in normalized_source)
        if not eligible:raise ValueError('VSEARCH : aucune séquence dans les limites de longueur choisies.')
        pending=self.work/('pending_vsearch_'+action)
        if pending.exists():shutil.rmtree(pending)
        pending.mkdir()
        input_file=pending/'input.fasta';write_fasta(normalized_source,input_file)
        refs=pending/'references.fasta' if references else None
        if refs:write_fasta(normalized(references),refs)
        prepared=commands(action,input_file,refs,pending,self.a)
        preparation=f'{len(source)} séquences · gaps retirés et U→T · {eligible} dans les limites de longueur.'
        if preview:return preparation+'\n\n'+self.tools.describe(prepared)
        self.ui.info(preparation)
        self._calculation_start=len(self.tools.history)
        self.tools.execute(prepared)
        summary,rows=summarize(action,pending,len(source),eligible)
        settings=calculation_settings('vsearch',self.a,self.mode,strategy=action)
        settings['action']=action
        result={'action':action,'engine':'VSEARCH','settings':settings,'summary':summary,'preview':rows,
                'preparation':{'gaps_removed':True,'rna_to_dna':True,'query_id':self.query.id if action=='search' else None},
                'commands':[dict(row) for row in self.tools.history[self._calculation_start:]]}
        (pending/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
        parent=self.work/'vsearch';parent.mkdir(exist_ok=True)
        target=parent/action
        publish_directory(pending,target)
        self.vsearch_results[action]=result
        (parent/'analyses.json').write_text(json.dumps(self.vsearch_results,ensure_ascii=False,indent=2),encoding='utf-8')
        self.record_calculation('vsearch_'+action,settings,summary=summary)
        self.dirty=True
        self.ui.success('Analyse VSEARCH terminée : '+action)
        for key,value in summary.items():self.ui.info(f'{key} : {value}')
        self.ui.info('Fichiers : vsearch/'+action+'/ · résultats intégrés au prochain HTML.')
        return result

    def make_tree(self,preview=False):
        """FastTree sur MSA préexistant ; méthode exploratoire, non datée."""
        if not self.msa_records or not self.msa_file:
            self.ui.warn('Créer un MSA avant FastTree.');return
        if len(self.msa_records)<3:
            self.ui.warn('FastTree nécessite au moins 3 séquences.');return
        nucleotide=all(sequence_type(r.seq)=='ADN/ARN' for r in self.msa_records)
        pending=self.work/'pending_tree';pending.mkdir(exist_ok=True)
        output=pending/'tree_fasttree.nwk'
        if preview:return self.tools.describe([self.tools.fasttree_command(self.msa_file,output,nucleotide,self.a)])
        self._calculation_start=len(self.tools.history)
        self.tools.fasttree(self.msa_file,output,nucleotide,options=self.a)
        # Validate and draw before replacing a previous successful tree.
        figure=tree_svg(output.read_text(encoding='utf-8'))
        target=self.work/output.name;output.replace(target);output=target
        figure_dir=self.work/'figures';figure_dir.mkdir(exist_ok=True)
        (figure_dir/(TREE_FIGURE+'.svg')).write_text(figure,encoding='utf-8')
        self.a.tree=True
        self.record_calculation('tree',calculation_settings('tree',self.a,self.mode),
            model=(('JC' if nucleotide else 'JTT') if self.a.fasttree_model=='auto' else self.a.fasttree_model),
            alignment_sha256=sha256(self.msa_file),output_sha256=sha256(output))
        self.dirty=True
        self.created.append(output.name)
        self.ui.success('Arbre exploratoire créé : tree_fasttree.nwk (Newick)')
        self.ui.info('Figure SVG : figures/'+TREE_FIGURE+'.svg · intégrée au prochain rapport HTML.')
        self.ui.warn('FastTree fournit une phylogénie approximative ; vérifier modèle, support et qualité du MSA pour publication.')

    def select_aligner(self):
        from av_navigation import Workbench
        return Workbench(self).choose_aligner()

    def show_stats(self):
        self.ui.heading('CONTRÔLE QUALITÉ FASTA')
        seqs=self.query_list if self.single else [self.query,*self.references]
        rows=[sequence_stats(r) for r in seqs]
        export_table(self.work/'fasta_statistics.tsv',rows)
        print(f'  {"Séquence":<27} {"Type":<11} {"Taille":>10} {"GC %":>7} {"N %":>6} {"H(bits)":>8}')
        for r in rows:
            print(f'  {r["id"][:27]:<27} {r["kind"]:<11} {r["length"]:>10,} {r["gc_pct"] if r["gc_pct"] is not None else "—":>7} {r["n_pct"]:>6} {r["entropy"]:>8}')

    def motif(self):
        self.ui.heading('RECHERCHE DE MOTIF IUPAC / DEUX BRINS')
        seqs=[self.query]
        motif=self.ui.ask('Motif IUPAC (ex. ATG, NNNGTT)','ATG')
        if not motif:return
        data=[]
        for r in seqs:data.extend(find_motifs(r,motif,True))
        self.motif_rows=data;export_table(self.work/'motifs.tsv',data,['seqid','start','end','strand','matched'])
        self.ui.info(f'{len(data)} occurrence(s) · motif {motif}')
        for r in data[:50]:print(f'  {r["seqid"][:23]:23} {r["start"]:>8}–{r["end"]:<8} {r["strand"]}  {r["matched"]}')
        if len(data)>50:self.ui.warn('Affichage limité à 50 occurrences, export complet dans motifs.tsv.')

    def find_orfs(self):
        self.ui.heading('ORF FINDER / SIX CADRES')
        if sequence_type(self.query.seq)!='ADN/ARN':self.ui.warn('QUERY nucléotidique requise.');return
        n=self.ui.ask('Longueur minimale en acides aminés','60',f'ORF : table NCBI {self.a.genetic_code} · départs {self.a.orf_starts}')
        try:min_aa=int(n)
        except ValueError:self.ui.warn('Longueur invalide');return
        if min_aa<1:self.ui.warn('Longueur positive requise');return
        rows=orfs(self.query,min_aa,include_partial=False,table=self.a.genetic_code,starts=self.a.orf_starts);self.orf_rows=rows
        export_table(self.work/'orfs.tsv',rows,['seqid','start','end','frame','strand','length_aa','complete'])
        write_fasta([Seq(f'ORF_{i}',f'ORF_{i} {r["seqid"]}:{r["start"]}-{r["end"]} frame={r["frame"]}',r['protein'])
                     for i,r in enumerate(rows,1)],self.work/'orfs.faa')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('orfs',calculation_settings('translation',self.a,self.mode),
            engine='AlignView Python',version=VERSION,minimum_aa=min_aa,output_sha256=sha256(self.work/'orfs.faa'))
        self.ui.info(f'{len(rows)} ORF complètes (≥{min_aa} aa)')
        for r in rows[:35]:print(f'  {r["start"]:>9}–{r["end"]:<9}  cadre {r["frame"]:>2}  {r["length_aa"]:>5} aa  {r["strand"]}')
        self.ui.warn('Ce détecteur simple ne remplace pas Prodigal/GeneMark ni une annotation structurale.')

    def translate_frames(self):
        self.ui.heading(f'TRADUCTION SIX CADRES / TABLE NCBI {self.a.genetic_code}')
        if sequence_type(self.query.seq)!='ADN/ARN':self.ui.warn('QUERY nucléotidique requise.');return
        seqs=[]
        for frame in (1,2,3,-1,-2,-3):
            aa=translate_nt(self.query.seq,frame,self.a.genetic_code)
            seqs.append(Seq(f'frame_{frame}',f'QUERY frame {frame} table={self.a.genetic_code}',aa))
            print(f'  {frame:>2}  {len(aa)} aa  {aa.count("*")} STOP  {aa[:95]}'+('…' if len(aa)>95 else ''))
        write_fasta(seqs,self.work/'six_frames.faa')
        self._calculation_start=len(self.tools.history)
        self.record_calculation('translation',calculation_settings('translation',self.a,self.mode),
            engine='AlignView Python',version=VERSION,output_sha256=sha256(self.work/'six_frames.faa'))
        self.ui.info('Traductions complètes dans six_frames.faa (sur demande de sauvegarde).')

    def extract_region(self):
        self.ui.heading('EXTRAIRE UNE RÉGION / QUERY')
        try:
            start=int(self.ui.ask('Début (1-indexée)','1'))
            end=int(self.ui.ask('Fin (incluse)',str(min(300,len(self.query.seq)))))
        except ValueError:self.ui.warn('Bornes invalides');return
        if not 1<=start<=end<=len(self.query.seq):self.ui.warn(f'Bornes hors plage 1–{len(self.query.seq)}');return
        s=self.query.seq[start-1:end];rc=self.ui.yesno('Complément inverse ?')
        if rc:
            if sequence_type(self.query.seq)!='ADN/ARN':self.ui.warn('Complément inverse réservé aux acides nucléiques');return
            s=reverse_complement(s)
        out=self.work/f'query_{start}_{end}.fasta'
        write_fasta([Seq(f'{self.query.id}_{start}_{end}',f'{self.query.id}_{start}_{end}'+('_RC' if rc else ''),s)],out)
        self.ui.success(f'{len(s)} caractères extraits : {out.name}')
        for i in range(0,min(200,len(s)),70):print('  '+s[i:i+70])

    def profiles(self):
        self.ui.heading('PROFILS EN FENÊTRES GLISSANTES')
        metric=self.ui.ask('Métrique : gc, entropy ou gap','gc' if sequence_type(self.query.seq)=='ADN/ARN' else 'entropy')
        try:
            w=int(self.ui.ask('Fenêtre',str(min(101,len(self.query.seq)))))
            st=int(self.ui.ask('Pas',str(max(1,w//4))))
            if metric=='gc' and sequence_type(self.query.seq)!='ADN/ARN':self.ui.warn('Le GC ne s’applique pas aux protéines.');return
            data=sliding_scores(self.query.seq,w,st,metric)
        except ValueError as e:self.ui.warn(str(e));return
        export_table(self.work/f'profile_{metric}.tsv',data)
        (self.work/f'profile_{metric}.svg').write_text(plot_svg(data,f'Profil {metric} (fenêtre {w}, pas {st})'))
        numeric=[r['value'] for r in data if r['value'] is not None]
        self.ui.info(f'{len(data)} fenêtres · moyenne {sum(numeric)/len(numeric):.2f}' if numeric else 'Aucune donnée')
        for r in data[:30]:print(f'  {r["start"]:>7}–{r["end"]:<7}  {r["value"]}')
        self.ui.info(f'Graphique : profile_{metric}.svg (conservé après sauvegarde)')

    def composition(self):
        self.ui.heading('COMPOSITION / K-MERS DE QUERY')
        k=self.ui.ask('k (1–8)','3' if sequence_type(self.query.seq)=='ADN/ARN' else '2')
        try:k=int(k)
        except ValueError:self.ui.warn('k invalide');return
        if not 1<=k<=8:self.ui.warn('k entre 1 et 8');return
        c=collections.Counter(self.query.seq[i:i+k] for i in range(len(self.query.seq)-k+1)
                              if '-' not in self.query.seq[i:i+k] and '?' not in self.query.seq[i:i+k])
        rows=[dict(kmer=key,count=count,frequency=count/max(1,sum(c.values()))) for key,count in c.most_common()]
        export_table(self.work/f'kmers_{k}.tsv',rows)
        n=sum(c.values());self.ui.info(f'{len(c)} mots distincts de taille {k} sur {n:,} fenêtres')
        for row in rows[:30]:
            bar='▰'*min(30,round(30*row['count']/max(1,rows[0]['count'])))
            print(f'  {row["kmer"]:<12} {row["count"]:>8,}  {100*row["frequency"]:>6.2f}%  {self.ui.c(GREEN,bar)}')

    def _compute_dotplot(self,ref,k):
        """Seules des séquences dans le même alphabet sont comparées dans le dotplot.

        En BLASTX/TBLASTN/TBLASTX, projection en AA du/des FASTA nucléotidiques
        selon le HSP de meilleur score. Les axes deviennent des aa traduits,
        et ne sont PAS des positions ADN originales.
        """
        u=units_for(self.mode,sequence_type(self.query.seq))
        q=self.query.seq;r=ref.seq;frame_q=None;frame_s=None
        hs=self.selected.get(ref.id,[])
        best=hs[0] if hs else None
        if u.translated_query:
            frame_q=self.a.frame or (best.qframe if best else None)
            if not frame_q:raise ValueError('Dotplot traduit : un HSP avec cadre QUERY ou --frame explicite est nécessaire.')
            q=translate_nt(q,frame_q,self.a.genetic_code).replace('*','X')
        if u.translated_subject:
            frame_s=best.sframe if best else None
            if not frame_s:raise ValueError('Dotplot traduit : un HSP avec cadre de référence connu est nécessaire.')
            r=translate_nt(r,frame_s,self.a.reference_genetic_code or self.a.genetic_code).replace('*','X')
        dtypeq='aa' if u.translated_query else u.query
        dtypes='aa' if u.translated_subject else u.subject
        if dtypeq!=dtypes:raise ValueError('Dotplot : alphabets incompatibles ; traduire les séquences avant comparaison.')
        data=dotplot(q,r,k)
        data.update(xunit=dtypeq,yunit=dtypes,qframe=frame_q,sframe=frame_s,
                    note='Dotplot de séquences traduites en aa ; cadres Q/S '+str(frame_q)+'/'+str(frame_s) if frame_q or frame_s else 'Dotplot en '+dtypeq)
        return data

    def make_dotplot(self):
        if self.single and len(self.query_list)<2:self.ui.warn('Dotplot : au moins deux séquences nécessaires.');return
        self.ui.heading('DOTPLOT / COMPARAISON DE MOTIFS EXACTS')
        ref,_=self.ref_select() if not self.single else (self.query_list[1],2)
        if ref is None:return
        default_k='3' if (self.mode in ('blastx','tblastn','tblastx','blastp','diamond-blastp','diamond-blastx','mmseqs-blastp') or sequence_type(self.query.seq)=='protéine') else str(self.a.kmer)
        k_str=self.ui.ask('Longueur k du motif exact',default_k,'ADN : 8–12 ; protéines : 3–5. Les axes sont précisés dans la figure.')
        try:k=int(k_str)
        except ValueError:self.ui.warn('k invalide');return
        if not 2<=k<=60:self.ui.warn('k entre 2 et 60');return
        self.dot_data=self._compute_dotplot(ref,k)
        (self.work/'dotplot.svg').write_text(dotplot_svg(self.dot_data,f'{self.query.id} ({self.dot_data["xunit"]})',f'{ref.id} ({self.dot_data["yunit"]})'),encoding='utf-8')
        (self.work/'dotplot.json').write_text(json.dumps(self.dot_data))
        self.ui.success(f'{len(self.dot_data["points"]):,} correspondances (k={k}) ; dotplot.svg')
        self.ui.info(self.dot_data['note'])
        if self.dot_data['truncated']:self.ui.warn('Limite de points atteinte ; augmenter k.')

    def make_report(self, open_preview=True):
        """Construit un rapport HTML entièrement autonome ; aucun serveur éphémère."""
        self.ui.heading('RAPPORT / EXPLORATION DANS LE NAVIGATEUR')
        refs=self.query_list[1:] if self.single else self.references
        report=build_report_v6(self.work,self.query,refs,
                [h for group in self.selected.values() for h in group],self.mode,self.settings(),
                self.tools.history,self.msa_records,self.dot_data,
                sha_query=sha256(self.query_path),sha_refs=sha256(self.refs_path) if self.refs_path else '',
                motifs=self.motif_rows,orf_rows=self.orf_rows)
        self.ui.success('Rapport construit · fichier HTML autonome · aucune connexion réseau nécessaire.')
        if open_preview and self.ui.tty and not self.a.no_menu:
            self.open_report()
        else:
            self.ui.info('HTML temporaire préparé ; sauvegarde avec --save / demande à la fermeture.')
        return report

    def _html_target(self):
        safe=re.sub(r'[^A-Za-z0-9_.-]+','_',self.query.id).strip('._')[:45] or 'query'
        return (Path.cwd()/f'alignview_{safe}_report.html').resolve()

    def export_html(self,destination):
        """Persist one self-contained HTML with no implicit project save."""
        target=Path(destination).expanduser().resolve()
        if target.suffix.lower()!='.html':raise ValueError('Extension .html requise.')
        if target.exists():raise ValueError('Le HTML existe déjà : '+str(target)+'. Choisir un autre nom.')
        report=self.make_report(open_preview=False)
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as output:output.write(report.read_bytes())
        self.ui.success('HTML indépendant créé : '+str(target))
        return target

    def resume_actions(self):
        """Only explicit actions run after resuming; saved action flags stay inert."""
        explicit=getattr(self.a,'_explicit',set())
        if 'vsearch' in explicit and self.a.vsearch!='none':self.vsearch_run(self.a.vsearch)
        if 'msa' in explicit and self.a.msa!='none':
            if self.single and self.a.msa in ('hits','hsp','fragments'):
                raise ValueError('Cette stratégie MSA nécessite un FASTA de références.')
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
            if not self.msa_records:raise ValueError('Créer un MSA avant --tree : ajouter --msa align, ou utiliser les menus.')
            self.make_tree()
        if 'report' in explicit and self.a.report:self.make_report(open_preview=False)

    def _publish_html(self):
        """Uniquement à l'initiative de l'utilisateur ; n'exporte pas la session complète."""
        # Reconstruire le rapport à chaque demande pour inclure les nouveaux HSP/MSA.
        report=self.make_report(open_preview=False)
        default=self._html_target()
        self.ui.info('Le HTML sera conservé indépendamment de la session (action explicite).')
        dest=self.ui.ask('Fichier HTML de destination',str(default),
                         'Chemin local permanent : ne dépend pas de /tmp ni d’un serveur HTTP')
        if not dest:return None
        target=Path(os.path.expanduser(dest)).resolve()
        if target.suffix.lower()!='.html':self.ui.warn('Extension .html requise.');return None
        if target.exists() and not self.ui.yesno('Le rapport existe. Remplacer le fichier ?',False):
            self.ui.info('Copie annulée.');return None
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(report,target)
        self.ui.success(f'HTML indépendant créé : {target}')
        self.ui.info(f'Ouvrir en local : {target.as_uri()}')
        return target

    def _open_local_html(self,target):
        remote=bool(os.environ.get('SSH_CONNECTION') or os.environ.get('SSH_TTY') or os.environ.get('SLURM_JOB_ID'))
        local_gui=bool(os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY') or sys.platform=='darwin' or os.name=='nt')
        if not remote and local_gui and not self.a.no_open:
            import webbrowser
            try:
                if webbrowser.open(target.as_uri(),new=2):
                    self.ui.success('Ouvert dans le navigateur local.');return
            except (OSError,webbrowser.Error):pass
            self.ui.warn('Ouverture du navigateur impossible ; utiliser le chemin ci-dessus.')
        if remote:
            import socket,getpass
            user=getpass.getuser();host=socket.gethostname()
            self.ui.info('Accès via VS Code Remote SSH : ouvrir le fichier avec l’explorateur distant.')
            self.ui.info('Ou récupérer le fichier depuis un terminal sur TON ordinateur (si hôte joignable) :')
            print('    scp '+shlex.quote(f'{user}@{host}:{target}')+' .')
            self.ui.info('Sur un nœud SLURM privé, utiliser le nœud login (scp -J) ou VS Code Remote.')
        self.ui.info('Aucune URL HTTP avec jeton : le rapport ne renverra plus de 404 à la fermeture.')

    def open_report(self):
        self.ui.heading('OUVRIR / RAPPORT AUTONOME')
        target=self._publish_html()
        if target:self._open_local_html(target)

    def copy_report(self):
        self.ui.heading('EXPORTER / UN SEUL FICHIER HTML')
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
        self.save_settings()
        # Les profils publication/complete comprennent le rapport ET les figures exportables.
        if self.a.export_profile in ('publication','complete'):
            self.make_report(open_preview=False)
        (self.work/'commands.json').write_text(json.dumps(self.tools.history,ensure_ascii=False,indent=2))
        summary=dict(mode=self.mode,query=self.query.id,reference_count=len(self.references),
                     hsps=len(self.hsps),hits=len(self.selected),
                     query_length=len(self.query.seq),msa=str(self.msa_file.name if self.msa_file else ''),
                     generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        (self.work/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
        export_dir=generate_exports(self,self.a.export_profile)
        self.ui.success(f'Exports scientifiques prêts : {len(list(export_dir.iterdir()))} fichiers · profil {self.a.export_profile}')

    def export_center(self):
        from av_navigation import Workbench
        return Workbench(self).export()

    def export_explicit_directory(self):
        """Pour le batch, export demandé indépendamment de --save."""
        if not self.a.export_dir:return
        self.export_summary()
        target=Path(self.a.export_dir).expanduser().resolve()
        if target.exists():raise ValueError(f'Répertoire d’export déjà présent : {target}')
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copytree(self.work/'exports',target)
        self.ui.success(f'Exports demandés : {target}')

    def slurm_script(self):
        self.ui.heading('SLURM / GÉNÉRER UN SCRIPT SBATCH')
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
        command+=['--lang',self.a.lang]
        if self.a.tools_prefix:command+=['--tools-prefix',self.a.tools_prefix]
        command+=extra_cli(self.a)
        command+=['--threads',str(threads),'--hits',str(self.a.hits),'--evalue',str(self.a.evalue),
                  '--identity',str(self.a.identity),'--coverage',str(self.a.coverage),'--no-menu','--report','--save','--out',location]
        contents=(f'#!/bin/bash\n#SBATCH --job-name=alignview_{jobname}\n#SBATCH --cpus-per-task={threads}\n'
                  '#SBATCH --mem=8G\n#SBATCH --time=04:00:00\n#SBATCH --output=alignview_%j.log\n'
                  'set -euo pipefail\n# Exécuter depuis le dossier du projet sauvegardé.\n'
                  'cd -- "$(dirname -- "${BASH_SOURCE[0]}")"\n# Ajuster la partition et les ressources selon le cluster.\n'
                  +shlex.join(command)+'\n')
        file=self.work/'alignview_job.sbatch';file.write_text(contents)
        self.dirty=True
        self.ui.success(f'Fichier préparé : {file.name}')
        show_file(file,self.ui)
        self.ui.info('Aucune soumission automatique : vérifier puis sbatch après sauvegarde.')

    def diagnostics(self):
        rows=diagnostic(self.tools)
        for row in rows:
            self.ui.info(f'{row["command"]}: '+(f'{row["path"]} · {row["version"]}' if row['available'] else f'absent → alignview --setup {row["setup_profile"]} --install-tools'))

    def guide(self):
        self.ui.heading('GUIDE DES MÉTHODES / LIMITES SCIENTIFIQUES')
        for key,text in HELP_TEXT.items():
            print();print('  '+self.ui.c(CYAN+BOLD,key.upper()))
            self.ui.wrapped(text,'    ')

    def edit_filters(self):
        self.ui.heading('FILTRES / RECLASSER LES HSP SANS RELANCER BLAST')
        self.ui.info('Réapplique identité, couverture ou bitscore aux HSP déjà récupérés.')
        changes={}
        for k,label in [('identity','Identité minimale (%)'),('coverage','Couverture QUERY union (%)'),
                        ('bitscore','Bitscore minimal'),('hits','Nombre maximal de références')]:
            current=getattr(self.a,k)
            value=self.ui.ask(label,str(current))
            try:changes[k]=int(value) if k=='hits' else float(value)
            except ValueError:self.ui.warn('Valeur non valide ; seuil inchangé')
        if not 0<=changes.get('identity',self.a.identity)<=100 or not 0<=changes.get('coverage',self.a.coverage)<=100 or changes.get('hits',self.a.hits)<1 or changes.get('bitscore',self.a.bitscore)<0:
            self.ui.warn('Seuils hors plage : utiliser des valeurs valides.');return
        for key,value in changes.items():setattr(self.a,key,value)
        self.apply_filters()
        self.ui.success(f'{len(self.selected)} références retenues sans nouvelle recherche.')
        self.show_hits()

    def rerun(self):
        if self.single or self.mode.startswith('global-'):self.ui.warn('Pas de recherche locale pour ce mode.');return
        if self.ui.yesno('Relancer la recherche et remplacer les résultats BLAST de cette session ?'):
            self.msa_records=None;self.msa_file=None;self.dot_data=None;self.search()

    def single_intro(self):
        self.show_dashboard()
        if self.a.msa=='view':self.show_msa_window()
        elif self.a.msa in ('align','full','all','global'):self.msa_run('all')

    def initial(self):
        if self.ui.tty and not self.a.no_menu and (getattr(self.a,'_guided',False) or not self.mode):
            if self.single and self.a.msa=='view':
                if len({len(r.seq) for r in self.query_list})!=1:
                    self.a.msa='none'
                    self.ui.notice='Longueurs différentes : créer d’abord un alignement dans Analyser.'
                    self.save_settings();return
                self.load_msa()
            elif not self.single and not self.mode:
                from av_navigation import Workbench
                self.mode=Workbench(self).recommended_mode() or ''
                self.a.mode=self.mode
            self.save_settings()
            return
        self.ui.banner('ALIGNEMENTS / EXPLORATION GUIDÉE')
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
            if self.a.no_menu:raise ValueError('Mode requis en non-interactif : --mode blastx, blastn, etc.')
            from av_navigation import Workbench
            self.mode=Workbench(self).recommended_mode() or 'blastn'
            self.a.mode=self.mode
        if self.mode.startswith('global-'):
            self.ui.info('Mode GLOBAL : pas de recherche BLAST, alignement sur demande.')
            if self.a.msa in ('align','full','all','global','hits'):self.msa_run('all')
        else:
            self.search()
            if self.a.msa in ('align','full','hits','global'):self.msa_run('hits')
            elif self.a.msa in ('all','fragments','hsp'):self.msa_run(self.a.msa)
        if self.a.pairwise!='none':self.run_pairwise('needle' if self.a.pairwise=='global' else 'water')
        if self.a.motif:
            data=find_motifs(self.query,self.a.motif,True)
            self.motif_rows=data;export_table(self.work/'motifs.tsv',data,['seqid','start','end','strand','matched'])
            self.ui.success(f'Motif {self.a.motif} : {len(data)} occurrence(s)')
        if self.a.dotplot:self._dot_batch()
        if self.a.tree:
            if not self.msa_records:self.msa_run('all')
            self.make_tree()
        if self.a.report:self.make_report()

    def _dot_batch(self):
        refs=self.references if self.references else self.query_list[1:]
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
        if not self.msa_file:self.ui.warn('Alignement multiple absent');return
        print(self.tools.alv(self.msa_file)[:150000])

    def list_outputs(self):
        self.ui.heading('FICHIERS / ESPACE TEMPORAIRE')
        for file in sorted(self.work.rglob('*')):
            if file.is_file() and 'database' not in str(file):
                print(f'  {str(file.relative_to(self.work)):<48} {file.stat().st_size/1024:9.1f} Kio')
        self.ui.info('Travail temporaire : aucune conservation sans export explicite ou confirmation à la sortie.')

    def save_checkpoint(self):
        default=self.a.out or f'alignview_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
        selected=self.ui.ask('Dossier du projet (q = annuler)',default)
        if selected.lower()=='q':return None
        self.export_summary()
        target=save_project(self,selected)
        self.last_saved=target;self.dirty=False
        self.ui.success(f'Projet enregistré : {target}')
        self.ui.info('Reprendre : alignview --resume '+shlex.quote(str(target)))
        return target

    def data_hub(self):
        from av_navigation import Workbench
        return Workbench(self).data()

    def save_and_exit(self):
        if self.exit_decision=='saved':
            self.ui.info(f'Projet conservé : {self.last_saved}');return
        if self.exit_decision=='discard':
            self.ui.info('Session quittée sans sauvegarder les résultats temporaires.');return
        if not self.a.save:
            if not self.ui.tty:
                if self.a.out:self.ui.info('--out désigne un dossier ; ajouter --save pour y conserver le projet.')
                if self.a.report and not (self.a.html_out or self.a.export_dir):
                    self.ui.warn('Rapport non conservé : ajouter --html-out rapport.html, --save ou --export-dir DOSSIER.')
                self.ui.info('Fin : résultats temporaires effacés (--save pour les conserver).')
                return
            self.ui.heading('FIN DE SESSION / ENREGISTREMENT')
            self.ui.info('Les analyses sont temporaires. Sans sauvegarde, elles seront supprimées.')
            if not self.ui.yesno('Conserver les fichiers de cette session ?'):
                self.ui.info('Session terminée ; fichiers temporaires supprimés.');return
        self.export_summary()
        target=self.a.out or f'alignview_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}'
        target_path=Path(target).expanduser().resolve()
        if target_path.exists():
            self.ui.warn(f'Dossier déjà présent : {target_path}')
            if not self.ui.tty:raise ValueError(f'Répertoire de sortie déjà présent : {target_path}. Choisir --out avec un chemin libre.')
            n=2
            while target_path.with_name(f'{target_path.name}_{n}').exists():n+=1
            proposed=target_path.with_name(f'{target_path.name}_{n}')
            selected=self.ui.ask('Nouveau dossier (q = annuler la sauvegarde)',str(proposed))
            if selected.lower()=='q':
                if self.ui.yesno('Abandonner la sauvegarde et supprimer les temporaires ?',False):return
                selected=str(proposed)
            target_path=Path(selected).expanduser().resolve()
            if target_path.exists():
                self.ui.warn('Dossier déjà présent : choix automatique du dossier alternatif libre.')
                target_path=proposed
        target_path.parent.mkdir(parents=True,exist_ok=True)
        save_project(self,target_path)
        self.ui.success(f'Résultats enregistrés : {target_path}')
        if (target_path/'report.html').exists():self.ui.info('Ouvrir report.html dans le navigateur, sans connexion nécessaire.')


def show_tool_list():
    groups={
        'Homologie et couverture':['MMseqs2 BLASTP pour références protéiques','BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX','DIAMOND BLASTP / BLASTX',
                                  'HSP multiples, union de couverture, gaps, substitutions et cadres',
                                  'Cartes QUERY et références, identité × couverture'],
        'Alignements':['MUSCLE 5, Clustal Omega et FastTree','Needle global, Water local','MAFFT hits / all / fragments / HSP',
                      'ALV et visionneuse MSA ANSI','Consensus, conservation, sites variables, distances p'],
        'Séquences':['Contrôle FASTA, GC, N, ambiguïtés','Motifs IUPAC deux brins, extraction de région',
                    'Traduction six cadres, ORF ATG→STOP','Entropie, profils glissants, k-mers'],
        'Visualisation & reproductibilité':['Atlas SVG conditionnel et arbre FastTree avec échelle et supports','Rapport HTML autonome, navigation guidée et exports SVG','SVG couverture, dotplot, profils et MSA',
                                           'CSV, TSV, JSON/JSONL, BED, GFF3, NEXUS, PHYLIP, Stockholm et SHA256','Script Slurm, modules automatiques et logs',
                                           'Documentation hors ligne : --docs ; HTML seul : --html-out rapport.html']}
    for group,items in groups.items():
        print('\n'+tr(group))
        for item in items:print('  • '+tr(item))


def _apply_start(args,request,demo_root):
    args._guided=True
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
        if not ui.tty:raise ValueError('Lancer dans un terminal interactif ou fournir un FASTA, --demo ou --resume.')
        request=startup(ui)
        if request is None:return False
        _apply_start(args,request,demo_root)
    return True


def main(argv=None):
    if sys.version_info<(3,10):
        print('AlignView requiert Python ≥ 3.10.',file=sys.stderr);return 2
    args=parse_args(argv)
    ui=UI(not args.no_color,args.simple_ui)
    if args.write_profile:
        try:
            path=write_profile(args,args.write_profile);ui.success('Profil de calcul enregistré : '+str(path));return 0
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
            print(f'{identifier:>2}  {name}'+(' · arrêts contextuels' if context_dependent(identifier) else ''))
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
                    ui.banner('DIAGNOSTIC / OUTILS SCIENTIFIQUES')
                    for row in rows:
                        ui.info(f'{row["command"]}: '+(row['path']+' · '+row['version'] if row['available'] else f'absent → --setup {row["setup_profile"]} --install-tools'))
            return 0
        except (ValueError,OSError) as exc:ui.error(str(exc));return 1
    session=None;history=[];failed=False
    with tempfile.TemporaryDirectory(prefix='alignview_workspace_') as workspace_name:
        workspace=Path(workspace_name)
        try:
            if not _select_start(args,ui,workspace/'demo'):return 0
            session=resume_project(args.resume,args,ui) if args.resume else Session(args,ui)
            for i,path in enumerate(session.history,1):
                archived=workspace/f'imported_{i:03d}'
                shutil.copytree(path,archived);history.append(archived)
            session.history=history
            if args.resume:session.resume_actions()
            if not args.resume:
                try:session.initial()
                except HomeRequested:pass
                except CommandError as exc:
                    if args.no_menu or not ui.tty:raise
                    ui.error(str(exc));ui.info('Le projet reste ouvert. Réglages permet de préparer les logiciels ; Analyser relance le calcul.')
                    ui.pause()
            while True:
                if session.a.export_dir:session.export_explicit_directory()
                if session.a.html_out:session.export_html(session.a.html_out)
                if not session.a.no_menu and ui.tty:session.nav()
                if not session.pending:break
                request,value=session.pending;session.pending=None
                archive=workspace/f'analysis_{len(history)+1:03d}'
                save_project(session,archive,include_history=False)
                fresh=argparse.Namespace(**vars(session.a))
                fresh.resume=None;fresh.demo=False;fresh.save=False;fresh.out='';fresh.export_dir=''
                fresh.html_out='';fresh._explicit=set()
                fresh.msa='none';fresh.tree=False;fresh.pairwise='none';fresh.report=False;fresh.dotplot=False
                fresh.vsearch='none'
                if request in ('query','reconfigure'):
                    from av_project import read_project
                    root,state=read_project(archive)
                    fresh.query=str(root/state['inputs']['query']['path'])
                    fresh.references=str(root/state['inputs']['references']['path']) if state['inputs'].get('references') else None
                    if request=='query':fresh.query_id=value
                    else:fresh.mode=None
                elif request=='files':
                    _apply_start(fresh,value,workspace/'demo')
                elif request=='resume':fresh.resume=str(value)
                replacement=None
                try:
                    replacement=resume_project(fresh.resume,fresh,ui) if fresh.resume else Session(fresh,ui)
                    if request in ('query','reconfigure'):
                        replacement.source_query_path=Path(state['inputs']['query'].get('original_path',replacement.source_query_path))
                        if replacement.refs_path:
                            replacement.source_refs_path=Path(state['inputs']['references'].get('original_path',replacement.source_refs_path))
                    if request=='query' and session.single and session.msa_records:
                        target=replacement.work/'alignment_active.fasta'
                        shutil.copyfile(session.msa_file,target)
                        replacement.msa_file=target;replacement.msa_records=read_fasta(target)
                    if fresh.resume and Path(fresh.resume).resolve().is_relative_to(workspace.resolve()):
                        replacement.last_saved=None;replacement.dirty=True
                    if not fresh.resume:
                        try:replacement.initial()
                        except HomeRequested:pass
                        except CommandError as exc:ui.error(str(exc));ui.pause()
                    # Imported history must outlive the project from which it came.
                    for i,path in enumerate(replacement.history,1):
                        imported=workspace/f'import_{len(history):03d}_{i:03d}'
                        shutil.copytree(path,imported);history.append(imported)
                except (ValueError,CommandError,OSError,KeyboardInterrupt) as exc:
                    if replacement:replacement.close()
                    ui.error(str(exc) or 'Changement annulé.');ui.pause();continue
                history.append(archive)
                session.close();session=replacement;session.history=history
            session.save_and_exit()
            return 0
        except KeyboardInterrupt:
            ui.warn('Session interrompue.')
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
                        ui.warn(f'Récupération conservée : {recovery}\n  Reprendre : alignview --resume '+shlex.quote(str(recovery)))
                    except (ValueError,OSError) as exc:
                        # Keep the original workspace if even the recovery copy failed.
                        session.temporary._finalizer.detach()
                        ui.warn(f'Copie de récupération impossible ({exc}). Travail conservé dans {session.work}')
                session.close() if not (failed and session.work.exists() and not session.temporary._finalizer.alive) else None

if __name__=='__main__':raise SystemExit(main())
