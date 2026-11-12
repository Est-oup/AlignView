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
            if type(value) is not bool:raise ValueError('attend oui ou non')
        elif self.kind in (int,float):
            if isinstance(value,bool) or not isinstance(value,(int,float)):raise ValueError('attend un nombre')
            if not math.isfinite(value) or (self.kind is int and int(value)!=value):raise ValueError('nombre invalide')
            if self.minimum is not None and value<self.minimum:raise ValueError(f'minimum : {self.minimum:g}')
            if self.maximum is not None and value>self.maximum:raise ValueError(f'maximum : {self.maximum:g}')
            value=self.kind(value)
        elif not isinstance(value,str):raise ValueError('attend du texte')
        if self.choices and value not in self.choices:raise ValueError('choix : '+', '.join(map(str,self.choices)))
        if self.kind is str and not self.choices and value and not re.fullmatch(r'[A-Za-z0-9_.+-]+',value):
            raise ValueError('utiliser un nom de matrice ou de tâche, sans espace')
        return value

    def display(self,args):
        value=getattr(args,self.name,self.default)
        if value is None:
            if self.name=='search_candidates':return f'Auto : {max(100,args.hits*10)} candidats moteur'
            if self.name in ('frame','sframe'):return 'Meilleur HSP, sinon +1 (AlignView)'
            if self.name=='needle_endopen':return f'Comme gapopen : {args.gapopen:g} (AlignView)'
            if self.name=='needle_endextend':return f'Comme gapextend : {args.gapextend:g} (AlignView)'
            if self.name=='mafft_maxiterate':return 'Auto : choix MAFFT ; global/local : 1000 ; fast : 0'
        if value is None or value in ('native',''):return 'Défaut du moteur installé'
        if self.kind is bool:return 'Oui' if value else 'Non'
        return str(value)


def P(name,label,effect,default=None,kind=float,lo=None,hi=None,choices=()):
    return Parameter(name,label,effect,default,kind,lo,hi,tuple(choices))


PARAMETERS = [
    P('genetic_code','Table de traduction QUERY (NCBI)','Code utilisé pour les traductions internes et les recherches traduites de la requête.',1,int,choices=tuple(TABLES)),
    P('reference_genetic_code','Table de traduction des références','Vide : même table que QUERY ; transmis à TBLASTN/TBLASTX et aux traductions des références.',None,int,choices=tuple(TABLES)),
    P('orf_starts','Codons de départ des ORF','atg : ATG uniquement ; table : tous les initiateurs de la table NCBI choisie.','atg',str,choices=('atg','table')),
    P('vsearch_id','Identité minimale VSEARCH','Fraction de 0 à 1, distincte du filtre d’identité locale BLAST.',.97,float,0,1),
    P('vsearch_iddef','Définition d’identité VSEARCH','0 à 4 ; 2 : exclut les gaps terminaux. Lire la documentation VSEARCH.',2,int,choices=(0,1,2,3,4)),
    P('vsearch_strand','Brins VSEARCH','both : recherche et regroupement sur les deux brins ; UCHIME de novo travaille sur le brin plus.','both',str,choices=('plus','both')),
    P('vsearch_query_cov','Couverture minimale de la requête VSEARCH','Fraction de 0 à 1, calculée par VSEARCH.',0.,float,0,1),
    P('vsearch_target_cov','Couverture minimale des références VSEARCH','Fraction de 0 à 1, calculée par VSEARCH.',0.,float,0,1),
    P('vsearch_maxaccepts','Candidats acceptés avant arrêt','0 : pas de plafond de candidats acceptés ; distinct du maximum de références rapportées.',0,int,0),
    P('vsearch_maxrejects','Candidats rejetés avant arrêt','0 : pas d’arrêt sur rejets ; mettre aussi maxaccepts et minwordmatches à 0 pour une recherche exhaustive.',32,int,0),
    P('vsearch_minwordmatches','Mots communs minimum VSEARCH','0 désactive ce préfiltre ; utile aux séquences courtes, avec un coût plus élevé.',10,int,0),
    P('vsearch_masking','Masquage VSEARCH','dust : faible complexité ; none : aucun. Les entrées internes sont normalisées en majuscules.','dust',str,choices=('none','dust')),
    P('vsearch_n_mismatch','Compter les ambiguïtés comme différences','Évite qu’un N soit traité comme une identité avec toute base compatible.',False,bool),
    P('vsearch_sizein','Lire les abondances des en-têtes','Utiliser les annotations ;size=N; existantes. Sinon chaque entrée compte pour une séquence.',False,bool),
    P('vsearch_min_length','Longueur minimale VSEARCH','Les séquences plus courtes sont exclues et comptées dans le résumé.',1,int,1),
    P('vsearch_max_length','Longueur maximale VSEARCH','Les séquences plus longues sont exclues et comptées dans le résumé.',50000,int,1),
    P('vsearch_abskew','Rapport d’abondance UCHIME','Abondance minimale des parents relativement à la chimère candidate.',2.,float,1),
    P('vsearch_minh','Score minimal de chimère UCHIME','Plus grand : détection plus stricte, avec moins de candidats détectés.',.28,float,1e-12),
    P('threads','CPU','Ressources des moteurs qui acceptent des threads ; FastTree standard reste monothread.',4,int,1),
    P('evalue','E-value maximale','Seuil transmis au moteur. Le modifier nécessite une nouvelle recherche.',1e-5,float,1e-300),
    P('search_candidates','Limite de candidats du moteur','Auto = max(100, 10 × références retenues). Ce plafond peut changer les résultats ; distinct des 100 hits affichés.',None,int,1),
    P('task','Tâche BLAST','Change graines et score. Sans option : BLASTN utilise megablast ; autres programmes : tâche native.','',str),
    P('sensitivity','Sensibilité DIAMOND','Recherche plus sensible = davantage de calcul ; ce réglage est propre à DIAMOND.','sensitive',str,choices=('fast','sensitive','more-sensitive','very-sensitive','ultra-sensitive')),
    P('blast_word_size','Taille de graine BLAST','Petites graines : détection plus sensible, calcul plus coûteux ; valeurs autorisées selon la tâche.',None,int,2),
    P('blast_reward','Score de match ADN','Gain pour une base identique ; dépend des coûts de mismatch et de gap.',None,int,1),
    P('blast_penalty','Score de mismatch ADN','Score négatif pour une substitution ; change le score et les statistiques.',None,int,hi=-1),
    P('blast_matrix','Matrice protéique BLAST','Ex. BLOSUM62, BLOSUM45, PAM30. Les coûts de gap doivent être compatibles.',None,str),
    P('blast_gapopen','Ouverture de gap BLAST','Coût propre à BLAST, indépendant de Needle/Water. TBLASTX est non gappé.',None,int,0),
    P('blast_gapextend','Extension de gap BLAST','Coût propre à BLAST. Seules certaines combinaisons matrice/coûts sont acceptées.',None,int,0),
    P('blast_masking','Masquage DUST / SEG','Masque la faible complexité. Le désactiver peut révéler des correspondances peu spécifiques.','native',str,choices=('native','yes','no')),
    P('blast_soft_masking','Masquage souple BLAST','true : masque les graines ; false : masque aussi l’extension. Défaut dépendant du programme.','native',str,choices=('native','true','false')),
    P('blast_strand','Brin de la requête ADN','both : deux brins ; plus/minus : limite la recherche à un brin.','both',str,choices=('both','plus','minus')),
    P('blast_max_hsps','HSP maximum par référence','Limite les segments rapportés, donc potentiellement la couverture calculée ensuite.',None,int,1),
    P('blast_comp_stats','Statistiques de composition BLAST','0 : désactivées ; 1/2/3 : corrections NCBI selon le moteur. Change scores et E-values.','native',str,choices=('native','0','1','2','3')),
    P('diamond_matrix','Matrice DIAMOND','Matrice protéique ; les coûts de gap doivent être compatibles.',None,str),
    P('diamond_gapopen','Ouverture de gap DIAMOND','Coût propre à DIAMOND, indépendant des autres moteurs.',None,int,0),
    P('diamond_gapextend','Extension de gap DIAMOND','Pénalise l’allongement des gaps.',None,int,0),
    P('diamond_masking','Masquage DIAMOND','0 : désactivé ; 1 : tantan ; seg : SEG (DIAMOND ≥ 2.0.12).','native',str,choices=('native','0','1','seg')),
    P('diamond_comp_stats','Composition DIAMOND','0 : désactivée ; 1 : correction compositionnelle ; 3/4 : ajustements de matrice, selon version.','native',str,choices=('native','0','1','3','4')),
    P('diamond_max_hsps','HSP maximum DIAMOND','0 : tous ; plusieurs HSP peuvent améliorer la couverture et augmenter le coût.',None,int,0),
    P('diamond_min_orf','ORF minimum DIAMOND BLASTX','AlignView utilisait déjà 1, contrairement au choix automatique de DIAMOND. 0 : automatique natif.',1,int,0),
    P('mmseqs_sensitivity','Sensibilité MMseqs2 (-s)','De 1 (rapide) à 7.5 (sensible). Indépendante du réglage DIAMOND.',None,float,1,7.5),
    P('mmseqs_mask','Masquage MMseqs2','1 : faible complexité masquée ; 0 : désactivé.','native',str,choices=('native','0','1')),
    P('mmseqs_comp_bias','Correction de composition MMseqs2','0/1 : désactive/active la correction du biais de composition.','native',str,choices=('native','0','1')),
    P('aligner','Moteur MSA','Les algorithmes et leurs paramètres ne sont pas interchangeables.','mafft',str,choices=('mafft','muscle','clustalo')),
    P('mafft','Stratégie MAFFT','auto : choix natif ; global : G-INS-i ; local : L-INS-i ; fast : FFT-NS-2.','auto',str,choices=('auto','global','local','fast')),
    P('mafft_op','Ouverture de gap MAFFT (--op)','Pénalité de gap en alignement groupe à groupe ; indépendante d’EMBOSS.',None,float,0),
    P('mafft_ep','Offset MAFFT (--ep)','Paramètre d’offset de MAFFT ; ne se confond pas avec gapextend d’EMBOSS.',None,float,0),
    P('mafft_maxiterate','Raffinements MAFFT','Global/local : 1000 par défaut AlignView. Pour régler ce nombre, choisir global, local ou fast plutôt que auto.',None,int,0,1000),
    P('muscle_method','Algorithme MUSCLE 5','align : PPP ; super5 : prévu pour des ensembles plus grands. MUSCLE 3 non pris en charge.','align',str,choices=('align','super5')),
    P('muscle_perturb','Graine de perturbation MUSCLE','0 : HMM standard ; une autre graine explore la variabilité de l’alignement.',0,int,0),
    P('muscle_perm','Permutation de l’arbre guide','Change le parcours de l’arbre. Une seule permutation produit un seul MSA.','none',str,choices=('none','abc','acb','bca')),
    P('muscle_refineiters','Raffinements MUSCLE 5','Nombre d’itérations de raffinement ; coût et résultat peuvent changer.',None,int,0),
    P('clustalo_iterations','Itérations Clustal Omega','Cycles combinés arbre guide / HMM. 0 : aucun raffinement itératif.',0,int,0),
    P('clustalo_full','Distances complètes Clustal Omega','Matrice de distances complète plutôt que mBed ; coût quadratique.',False,bool),
    P('clustalo_full_iter','Distances complètes aux itérations','Matrice complète pendant les raffinements ; requiert des itérations.',False,bool),
    P('hsp_overlap','Recouvrement minimal de l’ancre HSP','Fraction de la région requête ancre recouverte pour retenir un segment dans le MSA exploratoire.',.5,float,0,1),
    P('frame','Cadre QUERY pour les analyses traduites','Vide : cadre du meilleur HSP, sinon +1. Ne réduit pas les six cadres recherchés par BLASTX.',None,int,choices=(-3,-2,-1,1,2,3)),
    P('sframe','Cadre de référence pour Needle/Water','Vide : cadre du meilleur HSP, sinon +1. Applicable à TBLASTN/TBLASTX.',None,int,choices=(-3,-2,-1,1,2,3)),
    P('gapopen','Ouverture de gap EMBOSS','Pénalité Needle/Water uniquement ; dépend de la matrice utilisée.',10,float,0),
    P('gapextend','Extension de gap EMBOSS','Pénalité Needle/Water uniquement.',.5,float,0),
    P('matrix','Matrice EMBOSS','Vide : EDNAFULL en ADN, EBLOSUM62 en protéines (défauts EMBOSS).','',str),
    P('free_ends','Extrémités libres Needle','AlignView pénalise les extrémités par défaut ; EMBOSS natif les laisse libres. Sans effet sur Water.',False,bool),
    P('needle_endopen','Ouverture de gap terminal Needle','Vide : même coût que gapopen AlignView ; ignoré si extrémités libres.',None,float,0),
    P('needle_endextend','Extension de gap terminal Needle','Vide : même coût que gapextend AlignView ; ignoré si extrémités libres.',None,float,0),
    P('fasttree_model','Modèle d’évolution FastTree','auto : JC pour ADN, JTT pour protéines. GTR : ADN ; JTT/WAG/LG : protéines.','auto',str,choices=('auto','JC','GTR','JTT','WAG','LG')),
    P('fasttree_gamma','Réévaluation Gamma20','Après optimisation CAT, réévalue et redimensionne les branches sous Gamma20 ; ne refait pas la topologie.',False,bool),
    P('fasttree_boot','Rééchantillonnages du support local','Support SH-like, pas un bootstrap d’arbres. 0 : désactive les supports locaux.',None,int,0),
    P('fasttree_seed','Graine FastTree','Fixe la graine aléatoire pour les calculs stochastiques.',None,int,1),
]
BY_NAME={p.name:p for p in PARAMETERS}
FILTERS=('hits','identity','coverage','bitscore','preset')


def add_arguments(parser):
    groups={}
    for p in PARAMETERS:
        if p.flag in parser._option_string_actions:continue
        label=next((label for prefix,label in (
            ('blast_','Paramètres BLAST'),('diamond_','Paramètres DIAMOND'),
            ('mmseqs_','Paramètres MMseqs2'),('mafft_','Paramètres MAFFT'),
            ('muscle_','Paramètres MUSCLE 5'),('clustalo_','Paramètres Clustal Omega'),
            ('fasttree_','Paramètres FastTree'),('needle_','Paramètres Needle')) if p.name.startswith(prefix)),
            'Autres paramètres scientifiques')
        if label not in groups:groups[label]=parser.add_argument_group(label)
        kwargs=dict(default=p.default,help=p.label+' : '+p.effect)
        if p.kind is bool:kwargs['action']=argparse.BooleanOptionalAction
        else:
            kwargs['type']=p.kind
            if p.choices:kwargs['choices']=p.choices
        groups[label].add_argument(p.flag,**kwargs)
    profiles=parser.add_argument_group('Profils de calcul réutilisables')
    profiles.add_argument('--load-profile',metavar='JSON',help='Charger un profil de calcul ; les options CLI explicites sont prioritaires')
    profiles.add_argument('--write-profile',metavar='JSON',help='Écrire le profil de calcul puis quitter, sans lancer de moteur')


def validate(args,*,context=None,mode=None):
    for p in PARAMETERS:
        try:p.check(getattr(args,p.name,p.default))
        except ValueError as exc:raise ValueError(f'{p.flag} : {exc}') from exc
    if context=='msa' and args.aligner=='mafft' and args.mafft=='auto' and args.mafft_maxiterate is not None:
        raise ValueError('--mafft-maxiterate nécessite --mafft global, local ou fast : --auto choisit lui-même la stratégie.')
    if context=='vsearch' and args.vsearch_min_length>args.vsearch_max_length:
        raise ValueError('VSEARCH : la longueur minimale dépasse la longueur maximale.')
    if context=='msa' and args.aligner=='clustalo' and args.clustalo_full_iter and not args.clustalo_iterations:
        raise ValueError('--clustalo-full-iter nécessite --clustalo-iterations > 0.')
    if context=='search' and mode in ('blastn','blastp','blastx','tblastn','tblastx') and args.task:
        tasks={'blastn':('blastn','blastn-short','megablast','dc-megablast'),
               'blastp':('blastp','blastp-short','blastp-fast'),'blastx':('blastx','blastx-fast'),
               'tblastn':('tblastn','tblastn-fast'),'tblastx':()}
        if args.task not in tasks[mode]:raise ValueError(f'Tâche {args.task!r} incompatible avec {mode}. Choix : '+', '.join(tasks[mode]))


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
    except (OSError,json.JSONDecodeError) as exc:raise ValueError(f'Profil illisible : {exc}') from exc
    if not isinstance(data,dict) or data.get('schema')!=1 or not isinstance(data.get('parameters'),dict):
        raise ValueError('Profil AlignView invalide (schema 1 attendu).')
    values=data['parameters'];unknown=set(values)-set(BY_NAME)-set(FILTERS)
    if unknown:raise ValueError('Paramètres de profil inconnus : '+', '.join(sorted(unknown)))
    candidate=argparse.Namespace(**vars(args))
    for name,value in values.items():
        if name in BY_NAME:value=BY_NAME[name].check(value)
        elif name=='preset':
            if value not in ('standard','sensible','strict'):raise ValueError('Profil de seuil inconnu.')
        elif isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or (name=='hits' and (value<1 or int(value)!=value)) or (name in ('identity','coverage') and value>100):
            raise ValueError('Filtre invalide : '+name)
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
