"""Task-based navigation. Scientific calculations stay in Session and av_core."""
from __future__ import annotations

from contextlib import redirect_stdout
from dataclasses import dataclass, replace
import io
import json
from pathlib import Path
import re
import shutil

from av_core import sequence_stats, sequence_type
from av_defaults import SEARCH_PRESETS
from av_parameters import (relevant, snapshot, load_profile, write_profile, BY_NAME)
from av_start import choose_file, choose_sequence, startup
from av_tools import CommandError
from av_ui import HomeRequested, MenuItem


def item(key,label,description='',*,enabled=True,reason=''):
    return MenuItem(str(key),label,description,enabled=enabled,reason=reason)


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
        if count<2:reason='Il faut au moins deux séquences. Ajouter un fichier dans Données.'
        elif s.mode in ('tblastn','tblastx'):
            reason='Plusieurs cadres de traduction sont possibles. Utiliser les segments détectés dans les options d’alignement.'
        elif s.single and len({sequence_type(r.seq) for r in s.query_list})>1:
            reason='Le fichier mélange des types de séquences. Séparer ADN et protéines.'
        calculated=getattr(s,'search_parameters',{})
        current=snapshot('search',s.a,s.mode)
        stale=bool(calculated) and any(name in calculated and calculated[name]!=value
            for name,value in current.items() if name!='candidate_limit')
        return cls(len(s.query_list),len(s.references),len(s.selected),bool(s.msa_records),
                   bool(s.hsps) or getattr(s,'search_completed',False),not reason,reason,stale)


HOME_ITEMS = (
    ('1','Explorer','Lire les séquences, les correspondances et les alignements'),
    ('2','Analyser','Aligner des séquences ou rechercher des similarités'),
    ('3','Données','Changer les fichiers, la requête ou le projet'),
    ('4','Enregistrer et partager','Conserver le projet, exporter ou créer un rapport'),
    ('5','Réglages et aide','Logiciels, paramètres et explication des méthodes'),
)


class Workbench:
    def __init__(self,session):
        self.s=session
        self.ui=session.ui

    @property
    def state(self):return WorkspaceState.from_session(self.s)

    def refresh_context(self):
        s=self.s;state=self.state
        source=(f'{s.source_query_path.name} · {state.sequence_count} séquences · active : {s.query.id}' if s.single
                else f'Requête : {s.query.id} · {state.reference_count} références')
        results=[]
        if state.search_stale:results.append('Recherche à relancer : réglages modifiés')
        if state.searched:results.append(f'{state.hit_count} correspondances')
        if state.aligned:results.append(f'alignement prêt ({len(s.msa_records[0].seq):,} colonnes)')
        run=getattr(s,'calculation_runs',{}).get('msa',{}).get('parameters',{})
        if run and run!=snapshot('msa',s.a,s.mode,run.get('strategy','')):
            results.append('alignement calculé avec les anciens réglages')
        if not results:results.append('aucun calcul lancé')
        saved='modifications à enregistrer' if getattr(s,'dirty',True) else 'projet enregistré'
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
                self.ui.page((self.s.work/artifact).read_text(errors='replace'),title or 'Résultats')
            self.refresh_context()
            return result
        except KeyboardInterrupt:
            self.ui.notice='Action interrompue. Le projet reste ouvert.'
        except (ValueError,CommandError,OSError) as exc:
            self.ui.error(str(exc))
            self.ui.page(str(exc),'Action non terminée')
        return None

    def run(self):
        self.ui.route=['Accueil']
        self.ui.positions.pop(('Accueil','Projet'),None)
        while True:
            state=self.state
            hint=('L’alignement est prêt : ouvrir Explorer pour le lire.' if state.aligned else
                  'La recherche est terminée : ouvrir Explorer pour examiner les correspondances.' if state.searched else
                  'Commencer par Analyser, ou consulter les séquences dans Explorer.')
            try:
                choice=self.choose('Projet',[item(*entry) for entry in HOME_ITEMS],hint=hint,
                    default='1' if state.aligned or state.searched else '2',
                    shortcuts={'v':'view','s':'save','c':'data','t':'software','h':'help','x':'exit'})
                if choice in (None,'exit'):
                    if self.confirm_exit():return
                    continue
                handlers={'1':self.explore,'2':self.analyze,'3':self.data,'4':self.save,
                          '5':self.settings,'view':self.view_alignment,'save':self.s.save_checkpoint,
                          'data':self.data,'software':self.software,'help':self.help}
                completed=self.run_action(handlers[choice])
                if choice=='2' and completed:self.ui.positions[('Accueil','Projet')]='1'
                if self.s.pending:return
            except HomeRequested:
                self.ui.positions[('Accueil','Projet')]='1' if self.state.aligned or self.state.searched else '2'
                self.ui.notice='Retour à l’accueil. Les données et résultats sont conservés.'
            except KeyboardInterrupt:
                self.ui.notice='Action annulée. Le projet reste ouvert.'
            except (ValueError,CommandError,OSError) as exc:
                self.ui.error(str(exc));self.ui.page(str(exc),'Action non terminée')

    def explore(self):
        with self.ui.section('Explorer'):
            while True:
                state=self.state
                ready=state.aligned or (self.s.single and len(self.s.query_list)>=2 and len({len(r.seq) for r in self.s.query_list})==1)
                choice=self.choose('Explorer les résultats',[
                    item('1','Lire l’alignement','Ouvrir la visionneuse et retrouver sa position',enabled=ready,
                         reason='Aucun alignement disponible. Choisir Analyser → Créer un alignement.'),
                    item('2','Correspondances et couvertures',f'{state.hit_count} références retenues par la recherche',enabled=state.searched,
                         reason='Lancer une recherche dans Analyser avec un fichier de références.'),
                    item('3','Consulter les séquences',f'{state.sequence_count} séquences chargées, identifiants et contenu'),
                    item('4','Qualité des séquences','Longueur, composition et caractères ambigus'),
                    item('5','Statistiques de l’alignement','Conservation, distances et arbre exploratoire',enabled=state.aligned,
                         reason='Créer ou ouvrir un alignement avant de calculer ses statistiques.'),
                ],default='1' if ready else '2' if state.searched else '3')
                if choice is None:return
                actions={'1':self.view_alignment,'2':self.search_results,'3':self.sequences,
                         '4':self.s.show_stats,'5':self.alignment_statistics}
                self.run_action(actions[choice])

    def view_alignment(self):
        state=self.state
        if not state.aligned and not (self.s.single and len(self.s.query_list)>=2 and len({len(r.seq) for r in self.s.query_list})==1):
            self.ui.page('Aucun alignement n’est disponible.\n\nOuvrir Analyser → Créer un alignement.','Ouvrir un alignement')
            return
        self.s.show_msa_window()

    def sequences(self):
        with self.ui.section('Séquences'):
            while True:
                records=self.s.query_list if self.s.single else [self.s.query,*self.s.references]
                selected=choose_sequence(self.ui,records,'Consulter une séquence')
                if selected is None:return
                stats=sequence_stats(selected)
                text=[selected.id,selected.description,'',f'Type : {stats["kind"]}',
                      f'Longueur : {len(selected.seq):,} caractères','',
                      *[f'{i+1:>8}  {selected.seq[i:i+80]}' for i in range(0,len(selected.seq),80)]]
                self.ui.page('\n'.join(text),'Séquence : '+selected.id)

    def hits(self):
        with self.ui.section('Correspondances'):
            while True:
                rows=self.s.coverage_rows()
                if not rows:
                    self.ui.page('Aucune correspondance retenue.\n\nRéduire les filtres dans Réglages → Recherche,\nou relancer une recherche plus sensible dans Analyser.','Recherche terminée sans correspondance')
                    return
                options=[item(str(i),r['subject'],f'{r["best_identity"]:.1f}% identité locale · {r["query_coverage"]:.1f}% de la requête · {r["hsps"]} segments') for i,r in enumerate(rows,1)]
                choice=self.choose('Correspondances',options,f'{len(rows)} références retenues · limite : {self.s.a.hits}',
                    hint='Choisir une référence pour lire ses segments alignés. / recherche son identifiant.')
                if choice is None:return
                selected=rows[int(choice)-1]
                output=io.StringIO()
                with redirect_stdout(output):
                    print(f'Référence : {selected["subject"]}\nIdentité : {selected["best_identity"]:.1f}% sur le meilleur segment local')
                    print(f'Couverture de la requête : {selected["query_coverage"]:.1f}% (union des segments)')
                    for hsp in self.s.selected[selected['subject']]:self.s.print_alignment(hsp)
                self.ui.page(output.getvalue(),'Segments alignés : '+selected['subject'])

    def search_results(self):
        with self.ui.section('Résultats de recherche'):
            while True:
                count=len(self.s.selected)
                choice=self.choose('Correspondances et couvertures',[
                    item('1','Choisir une correspondance',f'{count} références ; rechercher et lire leurs segments'),
                    item('2','Couverture de la requête','Régions détectées ou sans correspondance',enabled=bool(count),reason='Aucune référence retenue après filtrage.'),
                    item('3','Couverture des références','Régions détectées sur chaque séquence de référence',enabled=bool(count),reason='Aucune référence retenue après filtrage.'),
                    item('4','Identité et couverture','Comparer la similarité locale et la part de la requête couverte',enabled=bool(count),reason='Aucune référence retenue après filtrage.'),
                    item('5','Qualité des segments','Substitutions, gaps et cadres de lecture',enabled=bool(count),reason='Aucune référence retenue après filtrage.'),
                ])
                if choice is None:return
                if choice=='1':self.hits()
                elif choice=='3':
                    output=io.StringIO()
                    with redirect_stdout(output):self.s.show_coverage(by_ref=True)
                    self.ui.page(output.getvalue(),'Couverture des références')
                else:self.run_action({'2':self.s.show_coverage,'4':self.s.identity_coverage,'5':self.s.hsp_diagnostics}[choice])

    def alignment_statistics(self):
        with self.ui.section('Statistiques de l’alignement'):
            while True:
                count=len(self.s.msa_records or [])
                choice=self.choose('Analyser l’alignement',[
                    item('1','Conservation et consensus','Repérer les colonnes stables, variables et ambiguës'),
                    item('2','Distances entre séquences','Distances p, en excluant les gaps de chaque paire'),
                    item('3','Arbre exploratoire','Calculer un arbre approximatif avec FastTree',enabled=count>=3,
                         reason='Un arbre nécessite au moins trois séquences alignées.'),
                ])
                if choice is None:return
                if choice=='3':self.tree()
                else:self.run_action({'1':self.s.msa_analysis,'2':self.s.msa_distances}[choice],changes=True)

    def analyze(self):
        with self.ui.section('Analyser'):
            while True:
                state=self.state
                choice=self.choose('Que souhaitez-vous calculer ?',[
                    item('1','Créer un alignement','Comparer les positions de plusieurs séquences',
                         enabled=state.can_align or (bool(self.s.selected) and self.s.mode in ('tblastn','tblastx')),
                         reason=state.align_reason),
                    item('2','Rechercher des similarités','Comparer la requête à des séquences de référence',enabled=not self.s.single,
                         reason='Ajouter un fichier de références dans Données pour lancer une recherche.'),
                    item('3','Comparer deux séquences','Alignement global ou local avec Needle / Water',enabled=not self.s.single,
                         reason='Ajouter des références dans Données pour comparer deux séquences.'),
                    item('4','Autres analyses','Motifs, traduction, ORF et profils de composition'),
                ])
                if choice is None:return
                action={'1':self.align,'2':self.search,'3':self.pairwise,'4':self.other_analyses}[choice]
                if self.run_action(action):return

    def align(self):
        with self.ui.section('Alignement'):
            strategy=('all' if self.s.single or self.s.mode.startswith('global-') else 'hits' if self.s.selected else 'all')
            if self.s.mode in ('tblastn','tblastx'):strategy='hsp'
            labels={'all':'Toutes les séquences','hits':'Requête et références retenues','fragments':'Ajouter une séquence partielle','hsp':'Segments détectés seulement'}
            while True:
                scope=labels[strategy]
                engine='mafft' if strategy=='fragments' else self.s.a.aligner
                choice=self.choose('Créer un alignement',[
                    item('1','Lancer l’alignement',f'{scope} · {engine.upper()} · {self.s.a.threads} CPU'),
                    item('2','Choisir le moteur',f'Actuellement : {engine.upper()}',enabled=strategy!='fragments',
                         reason='L’ajout de fragments utilise MAFFT en deux étapes. Changer les séquences à aligner pour choisir un autre moteur.'),
                    item('3','Choisir les séquences',scope,enabled=not self.s.single,
                         reason='Dans un fichier unique, toutes ses séquences sont alignées.'),
                    item('4','Paramètres de calcul',self.parameter_summary('msa',strategy=strategy)),
                    item('5','Voir les commandes et la préparation','Inspecter les arguments avant de lancer le moteur'),
                    item('6','Profils de calcul','Charger ou enregistrer ses réglages pour les réutiliser'),
                ],hint=('Fragments : références alignées avec la stratégie choisie, puis ajout MAFFT --addfragments --auto.' if strategy=='fragments' else
                    'Les paramètres de calcul changent l’alignement. Le résultat sera ouvert dans la visionneuse.'))
                if choice is None:return False
                if choice=='2':self.choose_aligner();continue
                if choice=='3':
                    picked=self.alignment_scope()
                    if picked:strategy=picked
                    continue
                if choice=='4':self.edit_parameters('msa',strategy=strategy);continue
                if choice=='5':self.preview(lambda:self.s.msa_run(strategy,preview=True));continue
                if choice=='6':self.profiles();continue
                program='mafft' if strategy=='fragments' else self.s.a.aligner
                aliases=('muscle','muscle5') if program=='muscle' else (program,)
                if not self.ensure_engine(aliases,program):continue
                result=self.run_action(lambda:self.s.msa_run(strategy),changes=True)
                return bool(result)

    def alignment_scope(self):
        state=self.state
        return self.choose('Séquences à aligner',[
            item('hits','Références retenues','La requête et les séquences ayant une correspondance',enabled=bool(self.s.selected) and state.can_align,
                 reason='Nécessite une recherche avec des références retenues et un cadre compatible.'),
            item('all','Toutes les références','Inclut aussi les références sans correspondance',enabled=state.can_align,reason=state.align_reason),
            item('fragments','Ajouter une séquence partielle','Ajout de la requête à un alignement des références avec MAFFT',enabled=state.can_align,reason=state.align_reason),
            item('hsp','Uniquement les segments détectés','Comparaison exploratoire des régions locales, sans les extrémités',enabled=bool(self.s.selected),
                 reason='Lancer une recherche avant de comparer ses segments.'),
        ])

    def choose_aligner(self):
        choice=self.choose('Moteur d’alignement',[
            item('mafft','MAFFT','Choix polyvalent ; stratégie automatique par défaut'),
            item('muscle','MUSCLE 5','Alternative pour ADN et protéines'),
            item('clustalo','Clustal Omega','Alternative pour les alignements multiples'),
        ],f'Moteur actuel : {self.s.a.aligner.upper()}')
        if choice:
            self.s.a.aligner=choice;self.s.save_settings();self.ui.notice='Moteur sélectionné : '+choice.upper()

    def recommended_mode(self):
        from av_start import suggest_modes
        from collections import Counter
        qtype=sequence_type(self.s.query.seq)
        rtype=Counter(sequence_type(r.seq) for r in self.s.references).most_common(1)[0][0]
        mode,_,_=suggest_modes(qtype,rtype)
        if mode in ('blastp','blastx') and not self.s.tools.has(mode) and self.s.tools.has('diamond'):
            mode='diamond-'+mode
        return mode

    def choose_search_method(self):
        from alignview import MODES
        from av_start import suggest_modes
        from collections import Counter
        qtype=sequence_type(self.s.query.seq)
        rtype=Counter(sequence_type(r.seq) for r in self.s.references).most_common(1)[0][0]
        default,_,_=suggest_modes(qtype,rtype)
        allowed={'blastn':('blastn','tblastx'),'blastp':('blastp','diamond-blastp','mmseqs-blastp'),
                 'blastx':('blastx','diamond-blastx'),'tblastn':('tblastn',)}.get(default,())
        return self.choose('Méthode de recherche',[item(mode,MODES[mode][0],MODES[mode][1],
            enabled=not allowed or mode in allowed,reason=f'Incompatible avec les types détectés : {qtype} → {rtype}.')
            for mode in MODES if not mode.startswith('global-')],f'{qtype} → {rtype} · recommandation : {default or "à vérifier"}')

    def search(self):
        with self.ui.section('Recherche'):
            mode=self.s.mode if self.s.mode and not self.s.mode.startswith('global-') else self.recommended_mode()
            if not mode:mode=self.choose_search_method()
            if not mode:return False
            while True:
                choice=self.choose('Rechercher des similarités',[
                    item('1','Lancer la recherche',f'{mode.upper()} · {self.s.a.preset} · jusqu’à {self.s.a.hits} références'),
                    item('2','Choisir la méthode',f'Actuellement : {mode.upper()}'),
                    item('3','Régler la recherche',f'E-value ≤ {self.s.a.evalue:g} · {self.s.a.threads} CPU'),
                    item('4','Paramètres du moteur',self.parameter_summary('search',mode)),
                    item('5','Voir les commandes','Indexation et recherche ; aucun moteur lancé par l’aperçu'),
                    item('6','Profils de calcul','Réutiliser ou partager ses paramètres'),
                ],hint='Calcul : graines, masquage, score. Filtres : sélection après calcul. Les deux peuvent changer les hits retenus.')
                if choice is None:return False
                if choice=='2':
                    selected=self.choose_search_method()
                    if selected:mode=selected
                    continue
                if choice=='3':self.search_settings(mode);continue
                if choice=='4':self.edit_parameters('search',mode);continue
                if choice=='5':self.preview(lambda:self.s.search_preview(mode));continue
                if choice=='6':self.profiles();continue
                program='diamond' if mode.startswith('diamond-') else 'mmseqs' if mode=='mmseqs-blastp' else mode
                if not self.ensure_engine((program,),'mmseqs' if program=='mmseqs' else 'diamond' if program=='diamond' else 'blast'):continue
                if program not in ('diamond','mmseqs') and not self.ensure_engine(('makeblastdb',),'blast'):continue
                previous=(self.s.mode,self.s.a.mode);previous_hsps=self.s.hsps
                try:
                    self.s.mode=self.s.a.mode=mode
                    self.s.check_type();self.s.search();self.s.dirty=True
                except BaseException:
                    if self.s.hsps is previous_hsps:self.s.mode,self.s.a.mode=previous
                    else:self.s.dirty=True
                    raise
                self.ui.notice=f'Recherche terminée : {len(self.s.selected)} références retenues. Ouvrir Explorer.'
                return True

    def ensure_engine(self,aliases,profile):
        if any(self.s.tools.available(command) for command in aliases):return True
        from av_setup import setup_tools
        with self.ui.section('Logiciel manquant'):
            while True:
                choice=self.choose('Installer '+profile.upper(),[
                    item('1','Installer dans l’environnement AlignView','Utiliser Conda, Mamba ou Micromamba'),
                    item('2','Voir les autres possibilités','Module HPC, chemin explicite et installation manuelle'),
                ],hint='Commande nécessaire : '+(' ou '.join(aliases))+'. Aucun calcul n’a été lancé.')
                if choice is None:return False
                if choice=='2':
                    self.ui.page('Module HPC : charger le module fourni par le cluster.\n\nChemin local : --tools-dir /chemin/bin\nExécutable précis : --executable '+aliases[0]+'=/chemin/'+aliases[0]+
                        '\n\nInstallation : alignview --setup '+profile+' --install-tools','Installer le moteur nécessaire')
                    continue
                setup_tools(self.ui,profile)
                if self.ui.yesno('Lancer cette installation ?',False):
                    setup_tools(self.ui,profile,install=True);self.s.tools.refresh_managed()
                    return any(self.s.tools.available(command) for command in aliases)

    def pairwise(self):
        with self.ui.section('Deux séquences'):
            choice=self.choose('Comment comparer les deux séquences ?',[
                item('needle','De bout en bout','Alignement global avec Needle ; inclut les extrémités'),
                item('water','Le meilleur segment local','Alignement local optimal avec Water'),
            ])
            if choice and self.ensure_engine((choice,),'emboss'):
                ref,_=self.s.ref_select()
                if ref is None:return
                while True:
                    action=self.choose('Préparer '+choice.upper(),[
                        item('1','Lancer la comparaison',f'{self.s.query.id} / {ref.id}'),
                        item('2','Paramètres de calcul',self.parameter_summary(choice)),
                        item('3','Voir les commandes','Matrice, pénalités, extrémités et traduction'),
                        item('4','Profils de calcul','Charger ou enregistrer les paramètres'),
                    ],hint='Needle : extrémités pénalisées par défaut AlignView ; option extrémités libres pour les relâcher.' if choice=='needle' else 'Water optimise un segment local ; les pénalités d’extrémités ne s’appliquent pas.')
                    if action is None:return
                    if action=='2':self.edit_parameters(choice)
                    elif action=='3':self.preview(lambda:self.s.run_pairwise(choice,reference=ref,preview=True))
                    elif action=='4':self.profiles()
                    else:self.run_action(lambda:self.s.run_pairwise(choice,reference=ref),changes=True);return

    def tree(self):
        with self.ui.section('Arbre exploratoire'):
            while True:
                choice=self.choose('Préparer FastTree',[
                    item('1','Calculer l’arbre',self.parameter_summary('tree')),
                    item('2','Paramètres de calcul','Modèle, Gamma20, support local et graine'),
                    item('3','Voir la commande','Inspecter les options FastTree avant le calcul'),
                    item('4','Profils de calcul','Charger ou enregistrer les paramètres'),
                ],hint='Arbre approximatif ; le support SH-like ne constitue pas un bootstrap d’arbres.')
                if choice is None:return
                if choice=='2':self.edit_parameters('tree')
                elif choice=='3':self.preview(lambda:self.s.make_tree(preview=True))
                elif choice=='4':self.profiles()
                else:
                    if not self.ensure_engine(('FastTree','fasttree'),'fasttree'):continue
                    self.run_action(self.s.make_tree,changes=True,artifact='tree_fasttree.nwk',title='Arbre exploratoire (Newick)');return

    def other_analyses(self):
        with self.ui.section('Autres analyses'):
            while True:
                nucleic=sequence_type(self.s.query.seq)=='ADN/ARN'
                choice=self.choose('Analyses de séquences',[
                    item('1','Rechercher un motif','Motif IUPAC et positions sur les deux brins',enabled=nucleic,reason='Cette analyse de motifs IUPAC nécessite une séquence ADN/ARN.'),
                    item('2','Traduction et régions codantes','Traduire, rechercher des ORF ou extraire une région',enabled=nucleic,reason='La traduction nécessite une séquence ADN/ARN.'),
                    item('3','Composition et profils','GC, entropie et fréquences de mots'),
                    item('4','Dotplot','Visualiser des correspondances de mots entre deux séquences',enabled=len(self.s.query_list)>1 or not self.s.single,
                         reason='Charger au moins deux séquences pour les comparer.'),
                ])
                if choice is None:return
                actions={'1':lambda:self.run_action(self.s.motif,changes=True,artifact='motifs.tsv',title='Positions des motifs'),
                         '2':self.coding_tools,'3':self.composition_tools,
                         '4':lambda:self.run_action(self.s.make_dotplot,changes=True)}
                self.run_action(actions[choice])

    def coding_tools(self):
        choice=self.choose('Traduction et régions',[
            item('1','Traduire les six cadres','Code génétique standard ; export des séquences protéiques'),
            item('2','Rechercher des ORF simples','Segments ATG → STOP ; analyse exploratoire'),
            item('3','Extraire une région','Choisir ses coordonnées et son orientation'),
        ])
        if choice:self.run_action({'1':self.s.translate_frames,'2':self.s.find_orfs,'3':self.s.extract_region}[choice],changes=True)

    def composition_tools(self):
        choice=self.choose('Composition',[
            item('1','Qualité et composition globale','Longueur, GC et caractères ambigus'),
            item('2','Profil le long de la séquence','GC, complexité ou gaps selon les fenêtres'),
            item('3','Fréquences de mots','Compter les k-mers exacts'),
        ])
        if choice:self.run_action({'1':self.s.show_stats,'2':self.s.profiles,'3':self.s.composition}[choice],changes=choice!='1')

    def data(self):
        with self.ui.section('Données'):
            while True:
                can_change=len(self.s.query_list)>1
                choice=self.choose('Données et projets',[
                    item('1','Choisir la séquence active' if self.s.single else 'Changer la requête active',
                         f'{len(self.s.query_list)} séquences ; les analyses ciblées utilisent la séquence active',enabled=can_change,
                         reason='Le fichier chargé ne contient qu’une seule séquence.'),
                    item('2','Ouvrir d’autres fichiers','Choisir un alignement, des séquences ou des références'),
                    item('3','Revenir à une analyse précédente',f'{len(self.s.history)} analyses dans l’historique',enabled=bool(self.s.history),
                         reason='L’historique se remplit lorsque les fichiers ou la requête changent.'),
                    item('4','Ouvrir un projet enregistré','Retrouver ses fichiers, paramètres et calculs'),
                    item('5','Voir les données chargées','Fichiers, requête, paramètres et résultats disponibles'),
                ],hint='Changer de données conserve l’analyse courante dans l’historique de cette session.')
                if choice is None:return
                if choice=='1':
                    record=choose_sequence(self.ui,self.s.query_list,'Choisir la séquence active' if self.s.single else 'Choisir la requête active',self.s.query.id)
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
                    selected=self.choose('Analyses précédentes',choices)
                    if selected:self.s.pending=('resume',self.s.history[int(selected)-1]);return
                elif choice=='4':
                    folder=self.ui.ask('Dossier du projet (q pour annuler)')
                    if folder and folder.casefold()!='q':
                        from av_project import read_project
                        root,_=read_project(folder);self.s.pending=('resume',root);return
                elif choice=='5':self.run_action(self.s.show_dashboard)

    def save(self):
        with self.ui.section('Enregistrer'):
            while True:
                choice=self.choose('Enregistrer et partager',[
                    item('1','Enregistrer le projet','Tout conserver pour reprendre le travail sur un autre poste'),
                    item('2','Créer un rapport HTML','Un fichier à consulter et partager dans un navigateur'),
                    item('3','Exporter les résultats','Choisir les tables, alignements et figures à conserver'),
                    item('4','Préparer un calcul sur HPC','Script Slurm inclus dans le projet à enregistrer ; sans soumission'),
                ],hint='Chemins relatifs depuis '+str(Path.cwd())+' · Projet : reprendre · HTML : partager · Exports : réutiliser')
                if choice is None:return
                self.run_action({'1':self.s.save_checkpoint,'2':self.s.open_report,'3':self.export,
                                 '4':self.s.slurm_script}[choice])

    def export(self):
        with self.ui.section('Exporter'):
            choice=self.choose('Que souhaitez-vous exporter ?',[
                item('standard','Résultats essentiels','FASTA, qualité, correspondances, couvertures et provenance'),
                item('publication','Tables, alignements et figures','Ajoute CSV, GFF3, NEXUS/PHYLIP, SVG et HTML'),
                item('complete','Toutes les données disponibles','Ajoute JSON/JSONL, régions sans correspondance et SHA256'),
            ])
            if choice is None:return
            destination=self.ui.ask('Dossier destination (q pour annuler)',f'alignview_{re.sub("[^A-Za-z0-9_.-]","_",self.s.query.id)[:35]}_exports')
            if destination.casefold()=='q':return
            target=Path(destination).expanduser().resolve()
            if target.exists():raise ValueError('Le dossier existe déjà : '+str(target))
            self.s.a.export_profile=choice
            self.s.export_summary();target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copytree(self.s.work/'exports',target)
            self.ui.notice='Résultats exportés dans '+str(target)

    def settings(self):
        with self.ui.section('Réglages'):
            while True:
                choice=self.choose('Réglages et aide',[
                    item('1','Logiciels et installation','Voir ce qui est disponible et installer les moteurs manquants'),
                    item('2','Réglages de recherche',f'Profil {self.s.a.preset} · maximum {self.s.a.hits} références'),
                    item('3','Réglages d’alignement',self.s.a.aligner.upper()+' : stratégie, gaps et raffinements'),
                    item('4','Comprendre les résultats','Différence entre recherche, alignement, couverture et identité'),
                    item('5','Calculs effectués et provenance','Paramètres réellement utilisés, versions et commandes'),
                    item('6','Profils de calcul','Charger ou enregistrer des réglages réutilisables'),
                    item('7','Documentation complète','Guides de toutes les fonctions, sorties, HTML, options et dépannage'),
                ])
                if choice is None:return
                self.run_action({'1':self.software,'2':self.search_settings,'3':self.alignment_settings,
                                 '4':self.help,'5':self.provenance,'6':self.profiles,'7':self.documentation}[choice])

    def documentation(self):
        from av_docs import documentation_menu
        documentation_menu(self.ui)

    def search_settings(self,mode=None):
        with self.ui.section('Recherche'):
            while True:
                choice=self.choose('Réglages de recherche',[
                    item('1','Profil de seuil E-value',self.s.a.preset+' ; sensibilité DIAMOND si applicable'),
                    item('2','Nombre maximal de références',str(self.s.a.hits)+' références retenues au maximum'),
                    item('3','E-value maximale',f'{self.s.a.evalue:g} ; plus petit = plus strict'),
                    item('4','Nombre de CPU',str(self.s.a.threads)),
                    item('5','Filtres des résultats','Identité, couverture et bitscore ; sans refaire la recherche'),
                    item('6','Paramètres du moteur','Graines, masquage, score et sensibilité, selon la méthode'),
                ],hint='Les filtres réutilisent les résultats existants. Une modification de recherche s’applique au prochain calcul.')
                if choice is None:return
                if choice=='1':self.choose_preset()
                elif choice in ('2','3','4'):
                    name,label,cast,minimum={'2':('hits','Nombre maximal de références',int,1),
                        '3':('evalue','E-value maximale',float,1e-300),'4':('threads','Nombre de CPU',int,1)}[choice]
                    changed=self.s._change_numeric(name,label,cast,minimum)
                    if changed and name=='hits' and self.state.searched:self.s.apply_filters()
                elif choice=='5':self.filters()
                elif choice=='6':
                    method=mode or self.s.mode
                    if method in ('msa','global-dna','global-aa',''):
                        if not self.s.single:method=self.recommended_mode()
                        else:
                            self.ui.page('Choisir une méthode dans Analyser → Rechercher des similarités, après avoir ajouté des références.','Méthode de recherche requise');continue
                    self.edit_parameters('search',method)

    def choose_preset(self):
        choice=self.choose('Profil de seuil E-value',[
            item('standard','Standard','E-value 1e-5 · réglage de départ'),
            item('sensible','Plus permissif','E-value 1e-2 ; DIAMOND very-sensitive si utilisé'),
            item('strict','Plus strict','E-value 1e-10 ; DIAMOND sensitive si utilisé'),
        ],f'La limite reste à {self.s.a.hits} références, quel que soit le profil.')
        if not choice:return False
        old=(self.s.a.evalue,self.s.a.sensitivity)
        self.s.a.preset=choice
        for name,value in SEARCH_PRESETS[choice].items():setattr(self.s.a,name,value)
        self.s.save_settings();self.ui.notice=f'Profil {choice} appliqué. Limite conservée : {self.s.a.hits} références.'
        return old!=(self.s.a.evalue,self.s.a.sensitivity)

    def filters(self):
        while True:
            choice=self.choose('Filtrer les correspondances',[
                item('1','Identité locale minimale',f'{self.s.a.identity:g}% par segment aligné'),
                item('2','Couverture minimale de la requête',f'{self.s.a.coverage:g}% par union des segments'),
                item('3','Bitscore minimal',str(self.s.a.bitscore)),
                item('4','Réinitialiser les filtres','Revenir à 0 pour identité, couverture et bitscore'),
            ],hint='Les seuils filtrent les correspondances déjà calculées. Une sélection différente invalide l’ancien alignement.')
            if choice is None:return
            if choice=='4':
                self.s.a.identity=self.s.a.coverage=self.s.a.bitscore=0
                changed=True
            else:
                name,label,maximum={'1':('identity','Identité minimale (%)',100),
                    '2':('coverage','Couverture minimale (%)',100),'3':('bitscore','Bitscore minimal',None)}[choice]
                changed=self.s._change_numeric(name,label,float,0,maximum)
            if changed:
                if self.state.searched:self.s.apply_filters()
                else:self.s.save_settings()
                self.ui.notice=f'Filtres appliqués : {len(self.s.selected)} références retenues.'

    def parameter_summary(self,context,mode=None,strategy=''):
        fields=relevant(context,self.s.a,mode or self.s.mode,strategy=strategy)
        custom=[p.label+' : '+p.display(self.s.a) for p in fields if getattr(self.s.a,p.name)!=p.default]
        if custom:return ' · '.join(custom[:2])+(' · autres réglages modifiés' if len(custom)>2 else '')
        if context=='msa' and (strategy=='fragments' or self.s.a.aligner=='mafft'):
            return 'MAFFT '+self.s.a.mafft+' · gaps natifs · raffinements 1000 en global/local'
        return 'Valeurs de départ AlignView et défauts du moteur ; ouvrir pour les consulter'

    def alignment_settings(self):
        while True:
            choice=self.choose('Réglages d’alignement',[
                item('1','Choisir le moteur',self.s.a.aligner.upper()),
                item('2','Paramètres de calcul',self.parameter_summary('msa')),
                item('3','Profils de calcul','Charger ou enregistrer les réglages'),
            ])
            if choice is None:return
            if choice=='1':self.choose_aligner()
            elif choice=='2':self.edit_parameters('msa')
            else:self.profiles()

    def edit_parameters(self,context,mode=None,strategy=''):
        mode=mode or self.s.mode
        engine=('MAFFT' if strategy=='fragments' else self.s.a.aligner.upper()) if context=='msa' else mode.upper() if context=='search' else 'FASTTREE' if context=='tree' else context.upper()
        with self.ui.section('Paramètres '+engine):
            while True:
                fields=relevant(context,self.s.a,mode,strategy=strategy)
                choice=self.choose('Paramètres de calcul · '+engine,
                    [item(str(i),p.label,p.display(self.s.a)+' · '+p.effect) for i,p in enumerate(fields,1)]+
                    [item('r','Rétablir les valeurs de départ','Réinitialiser uniquement les paramètres de ce calcul')],
                    hint=('Traduction : code génétique standard (1) uniquement. ' if mode in ('blastx','diamond-blastx','tblastn','tblastx') else '')+
                        'F1 : effet du paramètre. / : rechercher. Défaut du moteur = option non transmise, selon sa version.')
                if choice is None:return
                if choice=='r':
                    for p in fields:setattr(self.s.a,p.name,p.default)
                    self.s.save_settings();self.ui.notice='Valeurs de départ rétablies pour '+engine;continue
                p=fields[int(choice)-1]
                old=getattr(self.s.a,p.name)
                if p.choices or p.kind is bool:
                    values=p.choices if p.choices else (True,False)
                    picked=self.choose(p.label,[item(str(i),('Oui' if value else 'Non') if p.kind is bool else str(value),
                        'Valeur de départ' if value==p.default else p.effect) for i,value in enumerate(values,1)]+
                        [item('r','Valeur de départ','Revenir à '+str(p.default if p.default is not None else 'automatique'))],p.effect)
                    if picked is None:continue
                    value=p.default if picked=='r' else values[int(picked)-1]
                else:
                    answer=self.ui.ask(p.label,str(old) if old is not None else 'défaut',
                        p.effect+'\nEntrée : conserver · défaut : réinitialiser · q : annuler')
                    if answer.casefold() in ('q','annuler'):continue
                    try:
                        value=p.default if answer.casefold() in ('défaut','defaut','native','auto') else p.kind(answer)
                        p.check(value)
                    except ValueError as exc:
                        self.ui.page(str(exc)+'\n\nLe paramètre est conservé.','Valeur invalide');continue
                setattr(self.s.a,p.name,value);self.s.save_settings()
                self.ui.notice=p.label+' : '+p.display(self.s.a)+' ; appliqué au prochain calcul.'

    def preview(self,prepare):
        result=self.run_action(prepare)
        if result:self.ui.page(result,'Commandes et préparation du calcul')

    def profiles(self):
        while True:
            choice=self.choose('Profils de calcul',[
                item('1','Charger un profil JSON','Appliquer des paramètres validés ; les calculs précédents gardent leur provenance'),
                item('2','Enregistrer les réglages actuels','Créer un JSON réutilisable par --load-profile sur Linux ou HPC'),
                item('3','Consulter le profil actuel','Voir les réglages de tous les moteurs et les filtres'),
            ],hint='Un profil contient des paramètres, sans séquences ni chemins d’exécutables. Charger ne relance aucun moteur.')
            if choice is None:return
            if choice=='3':
                from av_parameters import profile
                self.ui.page(json.dumps(profile(self.s.a),ensure_ascii=False,indent=2),'Profil actuel');continue
            filename=self.ui.ask('Fichier JSON (q pour annuler)', 'alignview-profile.json' if choice=='2' else '')
            if not filename or filename.casefold()=='q':continue
            try:
                if choice=='2':write_profile(self.s.a,filename);self.ui.notice='Profil enregistré : '+filename
                else:
                    old=tuple(getattr(self.s.a,n) for n in ('hits','identity','coverage','bitscore'))
                    load_profile(self.s.a,filename)
                    if self.state.searched and old!=tuple(getattr(self.s.a,n) for n in ('hits','identity','coverage','bitscore')):self.s.apply_filters()
                    else:self.s.save_settings()
                    self.ui.notice='Profil chargé. Les nouveaux paramètres s’appliqueront au prochain calcul.'
            except (OSError,ValueError) as exc:self.ui.page(str(exc),'Profil non appliqué')

    def provenance(self):
        runs=getattr(self.s,'calculation_runs',{})
        text=['Les réglages actuels peuvent différer de ceux ayant produit les résultats.',
              'Chaque calcul réussi conserve ses paramètres et les versions/commandes réellement exécutées.',
              'null / native : option non transmise ; consulter la version du moteur et sa documentation.',
              'Les recherches traduites et traductions internes utilisent le code génétique standard (1).',
              '',json.dumps(runs,ensure_ascii=False,indent=2) if runs else 'Aucune provenance détaillée de calcul (projet ancien ou alignement importé).']
        if not runs and self.s.tools.history:text+=['',json.dumps(self.s.tools.history,ensure_ascii=False,indent=2)]
        self.ui.page('\n'.join(text),'Calculs effectués et provenance')

    def software(self):
        from av_setup import software_hub
        with self.ui.section('Logiciels'):software_hub(self.ui,self.s.tools)

    def help(self):
        self.ui.page('Un FASTA contient des séquences ; il ne constitue pas nécessairement un alignement.\n\n'
            'ALIGNEMENT MULTIPLE\nComparer les positions de séquences déjà supposées homologues.\n'
            'MAFFT, MUSCLE et Clustal Omega construisent cet alignement.\n\n'
            'RECHERCHE DE SIMILARITÉS\nComparer une requête à des références. BLAST et DIAMOND trouvent des\n'
            'segments localement similaires. Une correspondance peut comporter plusieurs segments (HSP).\n\n'
            'IDENTITÉ\nPart des caractères identiques dans un segment local.\n'
            'Elle ne mesure pas l’identité de la séquence entière.\n\n'
            'COUVERTURE\nPart de la requête touchée par l’union des segments retenus.\n'
            'Sans correspondance détectée ne signifie pas sans homologie.\n\n'
            'ENREGISTRER\nProjet : fichiers et calculs pour reprendre.\nRapport HTML : présentation à partager.\n'
            'Exports : données pour d’autres logiciels.\n\n'
            'NAVIGATION\nÉchap ou flèche gauche : retour au menu précédent.\n'
            'F2 : retour direct à l’accueil. F1 : détails du choix sélectionné.\n'
            '/ : rechercher dans une liste. Numéro puis Entrée : ouvrir un choix.','Comprendre AlignView')

    def confirm_exit(self):
        if not getattr(self.s,'dirty',True) and getattr(self.s,'last_saved',None) and Path(self.s.last_saved).is_dir():
            self.s.exit_decision='saved';return True
        choice=self.choose('Quitter AlignView ?',[
            item('1','Enregistrer et quitter','Choisir une destination et conserver le projet complet'),
            item('2','Quitter sans enregistrer','Abandonner les modifications ; les copies déjà enregistrées restent conservées'),
            item('3','Continuer le travail','Revenir au projet'),
        ])
        if choice=='1':
            if self.s.save_checkpoint():self.s.exit_decision='saved';return True
        elif choice=='2':self.s.exit_decision='discard';return True
        return False
