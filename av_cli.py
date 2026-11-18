"""Task commands, argument validation and concise/full CLI help."""
from __future__ import annotations
import argparse
import math
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from av_docs import TOPICS, ALIASES as DOC_ALIASES
from av_vsearch import ACTIONS as VSEARCH_ACTIONS

from av_exports import PRESETS
from av_version import VERSION
from av_defaults import DEFAULT_HITS, SEARCH_PRESETS
from av_parameters import PARAMETERS, add_arguments, validate as validate_parameters, load_profile
from av_setup import PACKAGES


MODES={
    'vsearch':('VSEARCH','Global similarity search between nucleotide sequences'),
    'blastn':('BLASTN','DNA → DNA: local nucleotide search'),
    'blastp':('BLASTP','Protein → protein: local search'),
    'blastx':('BLASTX','DNA → protein: six-frame QUERY translation'),
    'tblastn':('TBLASTN','Protein → DNA: six-frame reference translation'),
    'tblastx':('TBLASTX','DNA → DNA: translation on both sides'),
    'diamond-blastp':('DIAMOND BLASTP','Protein → protein search optimized for large datasets'),
    'diamond-blastx':('DIAMOND BLASTX','DNA → protein, DIAMOND, translated QUERY'),
    'mmseqs-blastp':('MMseqs2','Protein → protein: fast local search with inspectable HSPs'),
    'global-dna':('GLOBAL DNA','Needle and MAFFT without a similarity search'),
    'global-aa':('GLOBAL PROTEINS','Needle and MAFFT without a similarity search'),
}
HELP_TEXT={
 'blast':'BLAST/DIAMOND find local high-scoring segments (HSPs). They can omit ends and divergent regions. No detected HSP is not proof of non-homology.',
 'coverage':'Coverage is the length of the union of QUERY intervals covered by HSPs divided by QUERY length. Overlaps are counted once. BLASTX QUERY positions use nucleotides.',
 'identity':'Identity is exact matches divided by local alignment length (HSP). It is not whole-sequence identity. Gaps and substitutions are reported separately.',
 'needle':'EMBOSS Needle performs global Needleman–Wunsch alignment. Its end-to-end path includes sequence ends and can introduce artificial gaps when sequences are not collinear.',
 'water':'EMBOSS Water performs optimal local Smith–Waterman alignment. Results can differ from heuristic BLAST searches.',
 'mafft':'MAFFT global (G-INS-i) and local (L-INS-i) are multiple-alignment strategies, not homology searches. MAFFT --auto chooses a strategy from the input. ADD-FRAGMENTS adds a partial query to an aligned reference set.',
 'hsp':'HSP MSA aligns local segments without flanking regions. AlignView excludes HSPs with insufficient overlap with the anchor region. Do not interpret this MSA as full-length sequences.',
 'dot':'Dotplots show exact k-mer matches at their positions. A diagonal suggests collinearity; multiple diagonals can indicate repeats or duplications. They do not measure whole-sequence identity.',
 'orf':'An ORF is a segment between a start and a stop codon under the selected NCBI code. Starts use ATG by default, or the code-specific initiators. Context-dependent stops require CDS annotation; this simple detector is not a gene predictor.',
 'msa':'Variable sites contain at least two distinct non-gap residues. Parsimony-informative sites contain at least two states, each appearing at least twice. Pairwise p-distances exclude gapped columns.',
 'save':'Calculated files are temporary until explicitly saved or exported. Abrupt termination can lose unsaved results.',
}


class CLIParser(argparse.ArgumentParser):
    """Keep everyday help short without hiding supported native controls."""
    def format_help(self):
        if getattr(self,'full',False):return super().format_help()
        common={'help','help_all','query','references','docs','version','no_menu','save','out','html_out','export_dir'}
        command=getattr(self,'command','')
        common|=({'query_id','mode','preset','hits','evalue','threads'} if command=='compare' else
                 {'aligner','mafft','threads','tree'} if command=='align' else
                 {'resume'} if command=='resume' else
                 {'setup','install_tools','tools_prefix','diagnose_json','tools_dir','executable'} if command=='tools' else
                 {'query_id','aligner','threads','tree','resume','demo'} if not command else set())
        previous=[(action,action.help) for action in self._actions]
        previous_usage=self.usage
        self.usage='%(prog)s [INPUT [REFERENCES]] [options]'
        try:
            for action,_ in previous:
                if action.dest not in common:action.help=argparse.SUPPRESS
            return super().format_help()
        finally:
            self.usage=previous_usage
            for action,help_text in previous:action.help=help_text


class CompleteHelp(argparse.Action):
    def __call__(self,parser,namespace,values,option_string=None):
        parser.full=True;parser.print_help();parser.exit()


def build_parser(*,full=False,command=""):
    p=CLIParser(prog='alignview',formatter_class=argparse.RawDescriptionHelpFormatter,
        description=('AlignView: inspect alignments, compare sequence evidence and share reproducible results.\n\n'
            'Commands:\n'
            '  view FILE                  Open an existing alignment\n'
            '  align FILE [REFERENCES]     Create an alignment\n'
            '  compare QUERY REFERENCES   Find matching sequences\n'
            '  report FILE_OR_PROJECT     Write standalone HTML and exit\n'
            '  resume PROJECT             Reopen a saved workspace\n'
            '  tools                      Check scientific engines'),
        epilog='Launch without arguments for the workspace launcher. Use --help-all for all engine parameters.\n'
            'Examples: alignview view aligned.fasta\n'
            '          alignview compare query.fna refs.faa --no-menu --save --out comparison\n'
            '          alignview report comparison --html-out report.html')
    p.full=full;p.command=command
    if command:p.prog='alignview '+command
    p.add_argument('--help-all',action=CompleteHelp,nargs=0,help='Show every engine-specific parameter')
    section_0=p.add_argument_group('Input and resuming')
    section_1=p.add_argument_group('Search and filtering')
    section_2=p.add_argument_group('Alignments and analyses')
    section_3=p.add_argument_group('Reports, exports and saving')
    section_4=p.add_argument_group('Scientific engines and installation')
    section_5=p.add_argument_group('Interface and documentation')
    section_0.add_argument('query',nargs='?',help='Query FASTA, or a single FASTA to explore')
    section_0.add_argument('--query-id',default='',help='FASTA identifier to use as the query when the first file contains several sequences')
    section_0.add_argument('references',nargs='?',help='Reference FASTA')
    section_0.add_argument('--mode',choices=['msa',*MODES],help='Search or alignment method')
    section_1.add_argument('--task',default='',help='NCBI task compatible with the selected program')
    section_1.add_argument('--evalue',type=float,default=1e-5,help='Maximum E-value (default: 1e-5)')
    section_1.add_argument('--identity',type=float,default=0,help='Minimum local identity %% (0)')
    section_1.add_argument('--coverage',type=float,default=0,help='Query coverage from the union of HSPs %% (0)')
    section_1.add_argument('--bitscore',type=float,default=0,help='Minimum local bitscore (0)')
    section_1.add_argument('--hits',type=int,default=DEFAULT_HITS,help='Maximum retained references (default: 100)')
    section_1.add_argument('--preset',choices=['standard','sensitive','strict'],default='standard',help='Search preset: standard, sensitive or strict; explicit options take priority')
    section_1.add_argument('--threads',type=int,default=4,help='Number of threads (4)')
    section_1.add_argument('--sensitivity',choices=['fast','sensitive','more-sensitive','very-sensitive','ultra-sensitive'],default='sensitive')
    section_2.add_argument('--msa',choices=['none','view','align','full','hits','all','global','fragments','hsp'],default='none')
    section_2.add_argument('--mafft',choices=['auto','global','local','fast'],default='auto')
    section_2.add_argument('--aligner',choices=['mafft','muscle','clustalo'],default='mafft',help='MSA engine: MAFFT, MUSCLE 5, Clustal Omega')
    section_2.add_argument('--tree',action='store_true',help='Infer an exploratory FastTree tree from the MSA')
    section_2.add_argument('--pairwise',choices=['none','global','local'],default='none')
    section_2.add_argument('--ref',type=int,default=1,help='Reference index (1)')
    section_2.add_argument('--frame',type=int,choices=[-3,-2,-1,1,2,3],help='Select a translation frame')
    section_2.add_argument('--gapopen',type=float,default=10)
    section_2.add_argument('--gapextend',type=float,default=.5)
    section_2.add_argument('--matrix',default='',help='EMBOSS substitution matrix')
    section_2.add_argument('--free-ends',action='store_true',help='Do not penalize terminal gaps in Needle')
    section_3.add_argument('--report',action='store_true',help='Prepare HTML and SVG figures; retain them using --save, --export-dir or --html-out')
    section_3.add_argument('--html-out',metavar='FILE.html',default='',help='Save a standalone HTML here, even without --save; refuses to overwrite')
    section_3.add_argument('--no-open',action='store_true',help='Do not automatically open HTML in the browser')
    section_3.add_argument('--preview-port',type=int,default=8765,help=argparse.SUPPRESS)
    section_2.add_argument('--dotplot',action='store_true',help='Generate a query / best-reference dotplot')
    section_2.add_argument('--kmer',type=int,default=9,help='Dotplot word length (9)')
    section_2.add_argument('--motif',default='',help='IUPAC motif in the query')
    section_2.add_argument('--vsearch',choices=['none',*VSEARCH_ACTIONS],default='none',help='VSEARCH nucleotide analysis: search, dereplication, clustering or chimeras')
    section_5.add_argument('--list-genetic-codes',action='store_true',help='List available NCBI translation tables and exit')
    section_3.add_argument('--out',default='',help='Project directory with --save; relative to the launch directory, created at the end')
    section_3.add_argument('--save',action='store_true',help='Explicitly save at the end, including batch runs')
    section_5.add_argument('--no-menu',action='store_true',help='Without interactive menus')
    section_5.add_argument('--simple-ui',action='store_true',help='Classic numbered menus')
    section_5.add_argument('--no-color',action='store_true',help='Without ANSI colors')
    section_4.add_argument('--executable',action='append',default=[],metavar='COMMAND=/path',
                   help='Exact external-tool path, repeatable, e.g. blastx=/opt/ncbi/bin/blastx')
    section_4.add_argument('--tools-dir',action='append',default=[],metavar='DIRECTORY',
                   help='Add local executable directories to AlignView\'s PATH (repeatable)')
    section_4.add_argument('--module-policy',choices=['auto','never'],default='auto',
                   help='auto: try HPC modules for missing commands; never: use existing paths and the managed environment')
    section_3.add_argument('--export-profile',choices=list(PRESETS),default='standard',
                   help='standard / publication / complete: files prepared when saving')
    section_3.add_argument('--export-dir',default='',metavar='DIRECTORY',
                   help='Explicitly export files to this directory, even without --save')
    section_4.add_argument('--diagnose',action='store_true',help='Report available executables and versions, then exit')
    section_5.add_argument('--list-tools',action='store_true',help='List features')
    section_5.add_argument('--docs',nargs='?',const='index',metavar='TOPIC',choices=[*TOPICS,*DOC_ALIASES],help='Read an offline guide without opening data (default: index)')
    section_0.add_argument('--resume',metavar='PROJECT',help='Resume a saved project without rerunning its calculations')
    section_0.add_argument('--demo',action='store_true',help='Explore a bundled example without external tools')
    section_4.add_argument('--setup',choices=sorted(PACKAGES),help='Prepare engine installation using Conda/Mamba')
    section_4.add_argument('--install-tools',action='store_true',help='Run the --setup plan in a dedicated environment')
    section_4.add_argument('--tools-prefix',help='Dedicated Conda prefix (or ALIGNVIEW_TOOLS_PREFIX)')
    section_4.add_argument('--diagnose-json',action='store_true',help='Machine-readable JSON diagnostics without loading modules')
    section_5.add_argument('--version',action='version',version='AlignView '+VERSION)
    add_arguments(p)
    return p


def parse_args(argv=None):
    source=list(argv if argv is not None else sys.argv[1:])
    verbs={'view','align','compare','report','resume','tools'}
    command=source.pop(0) if source and source[0] in verbs else ''
    if command=='resume':source=['--resume',*source]
    if command in ('view','align') and not any(token=='--msa' or token.startswith('--msa=') for token in source):
        source+=['--msa','view' if command=='view' else 'all']
    p=build_parser(command=command)
    a=p.parse_args(source)
    a._command=command;a._guided=command in ('view','align','compare');a._intent=command
    if command=='tools' and (a.query or a.references or a.resume or a.demo):
        p.error('tools does not accept datasets; use view, align, compare or resume.')
    if command=='compare' and a.no_menu and not a.references:
        p.error('compare requires a query FASTA and a reference FASTA in --no-menu mode.')
    if command=='report':
        if a.references:p.error('report accepts one alignment file or project directory.')
        if not (a.query or a.resume):p.error('report requires an alignment file or project directory.')
        a.no_menu=True
        if a.query and Path(a.query).expanduser().is_dir():
            if a.references:p.error('report accepts one alignment file or project directory.')
            a.resume=a.query;a.query=None
        elif a.query:a.msa='view'
        if not a.html_out:
            item=a.resume or a.query
            if not item:p.error('report requires an alignment file or project directory.')
            a.html_out=Path(item).stem+'_alignview.html'
    elif command=='tools' and not a.setup:a.diagnose=True
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
    # Without input files, open the launcher unless running in batch mode.
    if not a.query and a.no_menu and not (a.diagnose or a.diagnose_json or a.list_tools or a.list_genetic_codes or a.docs or a.setup or a.resume or a.demo or a.write_profile):p.error('Provide a FASTA in --no-menu mode')
    if a.mode=='vsearch' and a.vsearch=='none':a.vsearch='search'
    if a.vsearch=='search' and a.mode is None and not a.resume:a.mode='vsearch'
    if a.vsearch!='none' and a.mode not in (None,'msa','vsearch'):p.error('Use --vsearch without a BLAST/DIAMOND mode; these analyses use different metrics.')
    if not all(math.isfinite(x) for x in (a.evalue,a.identity,a.coverage,a.bitscore,a.gapopen,a.gapextend)) or a.ref<1 or a.gapopen<0 or a.gapextend<0 or a.evalue<=0 or a.hits<1 or a.threads<1 or a.kmer<2 or a.kmer>60 or not 0<=a.identity<=100 or not 0<=a.coverage<=100 or a.bitscore<0 or not 0<=a.preview_port<=65535:
        p.error('Invalid values: evalue > 0, hits/threads > 0, 2 ≤ kmer ≤ 60, identity/coverage between 0 and 100')
    if a.install_tools and not a.setup:p.error('--install-tools requires --setup PROFILE')
    if a.resume and (a.query or a.references or a.demo):p.error('--resume cannot be combined with input files or --demo')
    if a.resume and explicit & {'mode','query_id'}:p.error('To change the method or query after --resume, use Find matching sequences / Files and query.')
    if a.html_out and Path(a.html_out).suffix.lower()!='.html':p.error('--html-out requires a filename ending in .html')
    if a.references and a.mode=='msa':p.error('--mode msa requires a single FASTA')
    if a.query and not a.references and a.mode not in (None,'msa'):p.error('This mode requires a reference FASTA.')
    if a.query and not a.references and a.msa in ('hits','hsp','fragments'):p.error('This MSA strategy requires a reference FASTA.')
    if a.references and a.msa=='view':p.error('--msa view requires a single aligned FASTA.')
    return a
