"""Typed calculation settings, shared by CLI, menus, profiles and provenance.

None/native means that no flag is sent: the installed engine owns the default.
No free-form command fragments are accepted here.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re

from av_version import VERSION
from av_genetics import TABLES


@dataclass(frozen=True)
class Parameter:
    name: str
    label: str
    effect: str
    default: object = None
    kind: type = float
    minimum: float | None = None
    maximum: float | None = None
    choices: tuple = ()

    @property
    def flag(self):return '--'+self.name.replace('_','-')

    def check(self,value):
        if value is None and self.default is None:return value
        if self.kind is bool:
            if type(value) is not bool:raise ValueError('expected true or false')
        elif self.kind in (int,float):
            if isinstance(value,bool) or not isinstance(value,(int,float)):raise ValueError('expected a number')
            if not math.isfinite(value) or (self.kind is int and int(value)!=value):raise ValueError('invalid number')
            if self.minimum is not None and value<self.minimum:raise ValueError(f'minimum: {self.minimum:g}')
            if self.maximum is not None and value>self.maximum:raise ValueError(f'maximum: {self.maximum:g}')
            value=self.kind(value)
        elif not isinstance(value,str):raise ValueError('expected text')
        if self.choices and value not in self.choices:raise ValueError('choix: '+', '.join(map(str,self.choices)))
        if self.kind is str and not self.choices and value and not re.fullmatch(r'[A-Za-z0-9_.+-]+',value):
            raise ValueError('use a matrix or task name without spaces')
        return value

    def display(self,args):
        value=getattr(args,self.name,self.default)
        if value is None:
            if self.name=='search_candidates':return f'Auto: {max(100,args.hits*10)} engine candidates'
            if self.name in ('frame','sframe'):return 'Meilleur HSP, sinon +1 (AlignView)'
            if self.name=='needle_endopen':return f'Comme gapopen: {args.gapopen:g} (AlignView)'
            if self.name=='needle_endextend':return f'Comme gapextend: {args.gapextend:g} (AlignView)'
            if self.name=='mafft_maxiterate':return 'Auto: choix MAFFT ; global/local: 1000 ; fast: 0'
        if value is None or value in ('native',''):return 'Installed engine default'
        if self.kind is bool:return 'Yes' if value else 'No'
        return str(value)


def P(name,label,effect,default=None,kind=float,lo=None,hi=None,choices=()):
    return Parameter(name,label,effect,default,kind,lo,hi,tuple(choices))


PARAMETERS = [
    P('genetic_code','Query translation table (NCBI)','Code used for internal query translation and translated query searches.',1,int,choices=tuple(TABLES)),
    P('reference_genetic_code','Reference translation table','Unset: use the query table; passed to TBLASTN/TBLASTX and internal reference translation.',None,int,choices=tuple(TABLES)),
    P('orf_starts','ORF start codons','atg: ATG only; table: all initiators from the selected NCBI table.','atg',str,choices=('atg','table')),
    P('vsearch_id','Minimum VSEARCH identity','Fraction from 0 to 1; distinct from the BLAST local identity filter.',.97,float,0,1),
    P('vsearch_iddef','VSEARCH identity definition','0 to 4; 2 excludes terminal gaps. See the VSEARCH documentation.',2,int,choices=(0,1,2,3,4)),
    P('vsearch_strand','VSEARCH strands','both: search and cluster both strands; de novo UCHIME uses the plus strand.','both',str,choices=('plus','both')),
    P('vsearch_query_cov','Minimum VSEARCH query coverage','Fraction from 0 to 1, calculated by VSEARCH.',0.,float,0,1),
    P('vsearch_target_cov','Minimum VSEARCH reference coverage','Fraction from 0 to 1, calculated by VSEARCH.',0.,float,0,1),
    P('vsearch_maxaccepts','Accepted candidates before stopping','0: no accepted-candidate limit; distinct from the reported reference limit.',0,int,0),
    P('vsearch_maxrejects','Rejected candidates before stopping','0: do not stop after rejections; also set maxaccepts and minwordmatches to 0 for exhaustive search.',32,int,0),
    P('vsearch_minwordmatches','Minimum shared words for VSEARCH','0 disables this prefilter; useful for short sequences, at greater computational cost.',10,int,0),
    P('vsearch_masking','VSEARCH masking','dust: low complexity; none: no masking. Internal input sequences are normalized to uppercase.','dust',str,choices=('none','dust')),
    P('vsearch_n_mismatch','Count ambiguous symbols as mismatches','Prevent N from counting as a match to every compatible nucleotide.',False,bool),
    P('vsearch_sizein','Read header abundances','Use existing ;size=N; annotations. Otherwise each input record counts as one sequence.',False,bool),
    P('vsearch_min_length','Minimum VSEARCH sequence length','Shorter sequences are excluded and counted in the summary.',1,int,1),
    P('vsearch_max_length','Maximum VSEARCH sequence length','Longer sequences are excluded and counted in the summary.',50000,int,1),
    P('vsearch_abskew','UCHIME abundance skew','Minimum parent abundance relative to the candidate chimera.',2.,float,1),
    P('vsearch_minh','Minimum UCHIME chimera score','Higher values make detection stricter, with fewer candidates detected.',.28,float,1e-12),
    P('threads','CPU','Resources for engines that support threads; standard FastTree remains single-threaded.',4,int,1),
    P('evalue','Maximum E-value','Threshold sent to the engine. Changing it requires a new search.',1e-5,float,1e-300),
    P('search_candidates','Engine candidate limit','Auto = max(100, 10 × retained references). This candidate limit can affect results and differs from the 100 displayed hits.',None,int,1),
    P('task','Task BLAST','Changes seeds and scores. Unset: BLASTN uses megablast; other programs use their native task.','',str),
    P('sensitivity','DIAMOND sensitivity','More sensitive search requires more work; this setting is specific to DIAMOND.','sensitive',str,choices=('fast','sensitive','more-sensitive','very-sensitive','ultra-sensitive')),
    P('blast_word_size','BLAST seed size','Smaller seeds increase sensitivity and cost; allowed values depend on the task.',None,int,2),
    P('blast_reward','DNA match score','Reward for an identical nucleotide; depends on mismatch and gap costs.',None,int,1),
    P('blast_penalty','DNA mismatch score','Negative mismatch score; affects scoring and statistics.',None,int,hi=-1),
    P('blast_matrix','BLAST protein matrix','Examples: BLOSUM62, BLOSUM45, PAM30. Gap costs must be compatible.',None,str),
    P('blast_gapopen','BLAST gap opening','BLAST-specific cost, independent of Needle/Water. TBLASTX is ungapped.',None,int,0),
    P('blast_gapextend','BLAST gap extension','BLAST-specific cost. Only supported matrix/cost combinations are accepted.',None,int,0),
    P('blast_masking','Masquage DUST / SEG','Mask low complexity. Disabling masking may reveal nonspecific matches.','native',str,choices=('native','yes','no')),
    P('blast_soft_masking','Masquage souple BLAST','true: mask seeds; false: also mask extension. The default depends on the program.','native',str,choices=('native','true','false')),
    P('blast_strand','DNA query strand','both: both strands; plus/minus: restrict the search to one strand.','both',str,choices=('both','plus','minus')),
    P('blast_max_hsps','Maximum HSPs per reference','Limit reported segments and potentially the coverage calculated from them.',None,int,1),
    P('blast_comp_stats','BLAST composition statistics','0: disabled; 1/2/3: NCBI corrections depending on the engine. Affects scores and E-values.','native',str,choices=('native','0','1','2','3')),
    P('diamond_matrix','DIAMOND matrix','Protein matrix; gap costs must be compatible.',None,str),
    P('diamond_gapopen','DIAMOND gap opening','DIAMOND-specific cost, independent of other engines.',None,int,0),
    P('diamond_gapextend','DIAMOND gap extension','Penalize gap extension.',None,int,0),
    P('diamond_masking','DIAMOND masking','0: disabled; 1: tantan; seg: SEG (DIAMOND ≥ 2.0.12).','native',str,choices=('native','0','1','seg')),
    P('diamond_comp_stats','DIAMOND composition correction','0: disabled; 1: composition correction; 3/4: matrix adjustments, depending on the version.','native',str,choices=('native','0','1','3','4')),
    P('diamond_max_hsps','Maximum DIAMOND HSPs','0: all; multiple HSPs may increase coverage and computational cost.',None,int,0),
    P('diamond_min_orf','Minimum DIAMOND BLASTX ORF','AlignView uses 1 rather than DIAMOND\'s automatic choice. 0: native automatic behavior.',1,int,0),
    P('mmseqs_sensitivity','MMseqs2 sensitivity (-s)','From 1 (fast) to 7.5 (sensitive). Independent of DIAMOND sensitivity.',None,float,1,7.5),
    P('mmseqs_mask','MMseqs2 masking','1: low complexity masked; 0: disabled.','native',str,choices=('native','0','1')),
    P('mmseqs_comp_bias','MMseqs2 composition correction','0/1: disable/enable composition-bias correction.','native',str,choices=('native','0','1')),
    P('aligner','MSA engine','Algorithms and their parameters are not interchangeable.','mafft',str,choices=('mafft','muscle','clustalo')),
    P('mafft','MAFFT strategy','auto: native choice; global: G-INS-i; local: L-INS-i; fast: FFT-NS-2.','auto',str,choices=('auto','global','local','fast')),
    P('mafft_op','MAFFT gap opening (--op)','Group-to-group alignment gap penalty; independent of EMBOSS.',None,float,0),
    P('mafft_ep','MAFFT offset (--ep)','MAFFT offset parameter; different from EMBOSS gapextend.',None,float,0),
    P('mafft_maxiterate','MAFFT refinement iterations','Global/local: AlignView defaults to 1000 iterations. Choose global, local or fast to set this number explicitly.',None,int,0,1000),
    P('muscle_method','MUSCLE 5 algorithm','align: PPP; super5: intended for larger datasets. MUSCLE 3 is unsupported.','align',str,choices=('align','super5')),
    P('muscle_perturb','MUSCLE perturbation seed','0: standard HMM; another seed explores alignment variability.',0,int,0),
    P('muscle_perm','Guide-tree permutation','Change the guide-tree traversal. One permutation produces one MSA.','none',str,choices=('none','abc','acb','bca')),
    P('muscle_refineiters','MUSCLE 5 refinement iterations','Refinement iteration count; cost and results may change.',None,int,0),
    P('clustalo_iterations','Clustal Omega iterations','Combined guide-tree/HMM cycles. 0: no iterative refinement.',0,int,0),
    P('clustalo_full','Full Clustal Omega distances','Full distance matrix instead of mBed; quadratic cost.',False,bool),
    P('clustalo_full_iter','Full distances during iterations','Use a full matrix during refinement; requires iterations.',False,bool),
    P('hsp_overlap','Minimum HSP anchor overlap','Fraction of the anchor query region that a segment must cover in the exploratory MSA.',.5,float,0,1),
    P('frame','Query frame for translated analyses','Unset: best HSP frame, otherwise +1. Does not restrict the six frames searched by BLASTX.',None,int,choices=(-3,-2,-1,1,2,3)),
    P('sframe','Reference frame for Needle/Water','Unset: best HSP frame, otherwise +1. Applies to TBLASTN/TBLASTX.',None,int,choices=(-3,-2,-1,1,2,3)),
    P('gapopen','EMBOSS gap opening','Needle/Water penalty only; depends on the selected matrix.',10,float,0),
    P('gapextend','EMBOSS gap extension','Needle/Water penalty only.',.5,float,0),
    P('matrix','EMBOSS matrix','Unset: EDNAFULL for DNA, EBLOSUM62 for proteins (EMBOSS defaults).','',str),
    P('free_ends','Free Needle ends','AlignView penalizes sequence ends by default; native EMBOSS leaves them free. No effect on Water.',False,bool),
    P('needle_endopen','Needle terminal gap opening','Unset: same as AlignView gapopen; ignored when sequence ends are free.',None,float,0),
    P('needle_endextend','Needle terminal gap extension','Unset: same as AlignView gapextend; ignored when sequence ends are free.',None,float,0),
    P('fasttree_model','FastTree evolutionary model','auto: JC for DNA, JTT for proteins. GTR: DNA; JTT/WAG/LG: proteins.','auto',str,choices=('auto','JC','GTR','JTT','WAG','LG')),
    P('fasttree_gamma','Gamma20 reassessment','After CAT optimization, rescore and rescale branches under Gamma20; does not repeat topology inference.',False,bool),
    P('fasttree_boot','Local-support resampling','SH-like support, not a tree bootstrap. 0: disable local support.',None,int,0),
    P('fasttree_seed','FastTree random seed','Set the random seed for stochastic calculations.',None,int,1),
]
BY_NAME={p.name:p for p in PARAMETERS}
FILTERS=('hits','identity','coverage','bitscore','preset')


def add_arguments(parser):
    groups={}
    for p in PARAMETERS:
        if p.flag in parser._option_string_actions:continue
        label=next((label for prefix,label in (
            ('blast_','BLAST parameters'),('diamond_','DIAMOND parameters'),
            ('mmseqs_','MMseqs2 parameters'),('mafft_','MAFFT parameters'),
            ('muscle_','MUSCLE 5 parameters'),('clustalo_','Clustal Omega parameters'),
            ('fasttree_','FastTree parameters'),('needle_','Needle parameters')) if p.name.startswith(prefix)),
            'Other scientific parameters')
        if label not in groups:groups[label]=parser.add_argument_group(label)
        kwargs=dict(default=p.default,help=p.label+': '+p.effect)
        if p.kind is bool:kwargs['action']=argparse.BooleanOptionalAction
        else:
            kwargs['type']=p.kind
            if p.choices:kwargs['choices']=p.choices
        action=groups[label].add_argument(p.flag,**kwargs)
        # Python 3.10's BooleanOptionalAction appends a default to help; newer
        # releases do not. Keep our canonical help identical for every version.
        action.help=kwargs['help']
    profiles=parser.add_argument_group('Reusable calculation profiles')
    profiles.add_argument('--load-profile',metavar='JSON',help='Load a calculation profile; explicit CLI options take priority')
    profiles.add_argument('--write-profile',metavar='JSON',help='Write the calculation profile and exit without running an engine')


def validate(args,*,context=None,mode=None):
    for p in PARAMETERS:
        try:p.check(getattr(args,p.name,p.default))
        except ValueError as exc:raise ValueError(f'{p.flag}: {exc}') from exc
    if context=='msa' and args.aligner=='mafft' and args.mafft=='auto' and args.mafft_maxiterate is not None:
        raise ValueError('--mafft-maxiterate requires --mafft global, local or fast: --auto selects its own strategy.')
    if context=='vsearch' and args.vsearch_min_length>args.vsearch_max_length:
        raise ValueError('VSEARCH: minimum length exceeds maximum length.')
    if context=='msa' and args.aligner=='clustalo' and args.clustalo_full_iter and not args.clustalo_iterations:
        raise ValueError('--clustalo-full-iter requires --clustalo-iterations > 0.')
    if context=='search' and mode in ('blastn','blastp','blastx','tblastn','tblastx') and args.task:
        tasks={'blastn':('blastn','blastn-short','megablast','dc-megablast'),
               'blastp':('blastp','blastp-short','blastp-fast'),'blastx':('blastx','blastx-fast'),
               'tblastn':('tblastn','tblastn-fast'),'tblastx':()}
        if args.task not in tasks[mode]:raise ValueError(f'Task {args.task!r} incompatible with {mode}. Choice: '+', '.join(tasks[mode]))


def relevant(context,args,mode='',*,strategy=''):
    names=[]
    if context=='search':
        names=['evalue','threads','search_candidates']
        if mode.startswith('diamond-'):
            names+=['sensitivity','diamond_matrix','diamond_gapopen','diamond_gapextend','diamond_masking','diamond_comp_stats','diamond_max_hsps']
            if mode=='diamond-blastx':names+=['diamond_min_orf']
        elif mode=='mmseqs-blastp':names+=['mmseqs_sensitivity','mmseqs_mask','mmseqs_comp_bias']
        elif mode in ('blastn','blastp','blastx','tblastn','tblastx'):
            if mode!='tblastx':names+=['task']
            names+=['blast_word_size','blast_masking','blast_soft_masking','blast_max_hsps']
            if mode in ('blastn','blastx','tblastx'):names+=['blast_strand']
            names+=['blast_reward','blast_penalty'] if mode=='blastn' else ['blast_matrix']
            if mode!='tblastx':names+=['blast_gapopen','blast_gapextend']
            if mode in ('blastp','blastx','tblastn'):names+=['blast_comp_stats']
    elif context=='msa':
        names=['threads']
        engine='mafft' if strategy=='fragments' else args.aligner
        if engine=='mafft':names+=['mafft','mafft_op','mafft_ep','mafft_maxiterate']
        elif engine=='muscle':names+=['muscle_method','muscle_perturb','muscle_perm','muscle_refineiters']
        else:names+=['clustalo_iterations','clustalo_full','clustalo_full_iter']
        if strategy=='hsp':names+=['hsp_overlap']
        elif mode in ('blastx','diamond-blastx'):names+=['frame']
    elif context in ('needle','water'):
        names=['gapopen','gapextend','matrix']
        if context=='needle':names+=['free_ends','needle_endopen','needle_endextend']
        if mode in ('blastx','diamond-blastx','tblastx'):names+=['frame']
        if mode in ('tblastn','tblastx'):names+=['sframe']
    elif context=='tree':names=['fasttree_model','fasttree_gamma','fasttree_boot','fasttree_seed']
    elif context=='translation':names=['genetic_code','reference_genetic_code','orf_starts']
    elif context=='vsearch':
        common=['vsearch_min_length','vsearch_max_length']
        selections={
            'search':['threads','vsearch_id','vsearch_iddef','vsearch_strand','vsearch_query_cov','vsearch_target_cov',
                      'vsearch_maxaccepts','vsearch_maxrejects','vsearch_minwordmatches','vsearch_masking','vsearch_n_mismatch'],
            'cluster':['threads','vsearch_id','vsearch_iddef','vsearch_strand','vsearch_masking','vsearch_n_mismatch','vsearch_sizein'],
            'dereplicate':['vsearch_strand','vsearch_sizein'],
            'chimeras':['vsearch_strand','vsearch_sizein','vsearch_masking','vsearch_abskew','vsearch_minh'],
        }
        names=selections.get(strategy,['threads',*[p.name for p in PARAMETERS if p.name.startswith('vsearch_')]])+common
        names=list(dict.fromkeys(names))
    if (context=='search' and mode in ('blastx','tblastn','tblastx','diamond-blastx')) or (context in ('msa','needle','water') and mode in ('blastx','diamond-blastx','tblastn','tblastx')):
        names+=['genetic_code','reference_genetic_code']
    return [BY_NAME[name] for name in names]


def snapshot(context,args,mode='',strategy=''):
    result={p.name:getattr(args,p.name,p.default) for p in relevant(context,args,mode,strategy=strategy)}
    if context=='search':result['candidate_limit']=args.search_candidates or max(100,args.hits*10)
    if context=='search' and mode in ('blastx','tblastn','tblastx','diamond-blastx'):result['genetic_code']=getattr(args,'genetic_code',1)
    if context=='msa':result.update(aligner='mafft' if strategy=='fragments' else args.aligner,strategy=strategy)
    if context=='vsearch':result['reported_limit']=args.hits
    result['mode']=mode
    return result


def profile(args):
    return dict(schema=1,alignview_version=VERSION,parameters={name:getattr(args,name)
        for name in (*BY_NAME,*FILTERS)})


def write_profile(args,filename):
    path=Path(filename).expanduser()
    # A profile is a new file, so a typo cannot overwrite a project/input.
    with path.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(profile(args),ensure_ascii=False,indent=2)+'\n')
    return path


def load_profile(args,filename,explicit=()):
    try:data=json.loads(Path(filename).expanduser().read_text(encoding='utf-8'))
    except (OSError,json.JSONDecodeError) as exc:raise ValueError(f'Cannot read profile: {exc}') from exc
    if not isinstance(data,dict) or data.get('schema')!=1 or not isinstance(data.get('parameters'),dict):
        raise ValueError('Invalid AlignView profile (expected schema 1).')
    values=data['parameters'];unknown=set(values)-set(BY_NAME)-set(FILTERS)
    if unknown:raise ValueError('Unknown profile parameters: '+', '.join(sorted(unknown)))
    candidate=argparse.Namespace(**vars(args))
    for name,value in values.items():
        if name in BY_NAME:value=BY_NAME[name].check(value)
        elif name=='preset':
            if value=='sensible':value='sensitive'  # Read profiles saved before the English-only interface.
            if value not in ('standard','sensitive','strict'):raise ValueError('Unknown threshold preset.')
        elif isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (name=='hits' and (value<1 or int(value)!=value)) or (name in ('identity','coverage') and value>100):
            raise ValueError('Invalid filter: '+name)
        if name not in explicit:setattr(candidate,name,value)
    validate(candidate)
    vars(args).update(vars(candidate))


def extra_cli(args):
    """New flags needed to reproduce settings in a generated Slurm job."""
    legacy={'threads','evalue','task','sensitivity','aligner','mafft','frame','gapopen','gapextend','matrix','free_ends'}
    result=[]
    for p in PARAMETERS:
        value=getattr(args,p.name,p.default)
        if p.name in legacy or value is None or value==p.default:continue
        result+=([p.flag] if p.kind is bool and value else ['--no-'+p.name.replace('_','-')] if p.kind is bool else [p.flag,str(value)])
    return result
