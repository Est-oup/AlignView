# Audit et comparaison — AlignView 0.1.0

Audit réalisé le 10 novembre 2026 à partir du prototype local v6. Les modifications sont dans le nouveau dossier `alignview/` ; la base v6 est conservée.

La [revue critique de navigation](UX_REVIEW.md) décrit la révision suivante : cinq rubriques, parcours guidés, prérequis, sauvegarde explicite et maximum de 100 références indépendant des profils.

La [revue des calculs](CALCULATION_REVIEW.md) détaille les paramètres auparavant cachés ou ignorés, leur correction, la comparaison des interfaces et les essais avec les huit moteurs réels. Le [guide des paramètres](PARAMETERS.md) décrit leur utilisation.

## Positionnement

Le parcours proposé est : ouvrir des séquences → choisir une question → calculer → explorer → conserver un projet portable. La valeur recherchée est la continuité de ce parcours sur Linux et SSH, avec des moteurs remplaçables et une provenance explicite.

Comparaison fonctionnelle des documentations officielles consultées ; les concurrents n’ont pas été chronométrés ou exécutés sur cette machine.

| Outil | Points utiles pour cette refonte | Choix dans AlignView |
| --- | --- | --- |
| [Jalview](https://www.jalview.org/help/html/features/overview.html) | Vue globale navigable, indication de la région visible, accès par clic | Barre globale et clic vers une région dans le terminal ; colonnes et curseur explicitement séparés |
| [AliView](https://ormbunkar.se/aliview/) | Navigation, zoom, repères de consensus, choix de moteurs d’alignement, formats d’échange | Accès rapide aux colonnes et séquences, différences au consensus, moteurs multiples et exports ouverts |
| [alv](https://github.com/arvestad/alv) | Visualisation terminal, couleurs biologiques, aperçu compact et plusieurs formats d’entrée | Visionneuse embarquée sans dépendance obligatoire ; ALV reste une alternative facultative |
| [alan](https://github.com/mpdunne/alan) | Consultation d’alignements dans un terminal Linux avec les outils Unix | Pager intégré avec recherche et défilement dans les deux axes, indépendant de la configuration de `less` |

Les logiciels de bureau conservent des fonctions absentes ici, notamment l’édition manuelle riche et des vues structurales. AlignView vise l’exploration et l’orchestration dans un terminal. Ce tableau ne démontre ni une originalité scientifique ni une supériorité de performance.

## Défauts traités

| Constat dans la base | Correction et vérification |
| --- | --- |
| Pieds de page plus larges que le terminal et hauteur variant avec le message d’état | Rendu borné en cellules Unicode, largeur maximale de `colonnes − 1`, ligne de message réservée, test à plusieurs dimensions et benchmark avec v6 |
| Décodage limité des séquences de touches dans les menus | Lecteur commun UTF-8, CSI/SS3, modificateurs et souris ; tests dans de vrais PTY |
| Déplacement latéral lent, pages peu explicites | Pas de dix colonnes, pages horizontales, pages de séquences distinctes, molette et vue globale cliquable |
| Position perdue en rouvrant le MSA | Conservation de la visionneuse dans la session |
| Fichiers et séquences coupés à des limites arbitraires | Suppression des plafonds de 35 dossiers / 100 fichiers, 120 requêtes et 40 références ; filtre sur l’identifiant complet |
| Sorties textuelles longues tronquées ou difficiles à consulter | Pager vertical/horizontal avec recherche ; HSP et fichiers accessibles en entier |
| Changer de requête ou de données nécessite de quitter | Parcours « Données et projets », sauvegarde temporaire du contexte précédent et historique reprenable ; test du parcours complet dans un PTY |
| Pas de reprise portable de la session calculée | Manifeste de projet, entrées copiées, SHA256, MSA actif, paramètres, HSP, historique ; reprise sans recalcul et sans chemins d’origine |
| Risque de copie partielle à la sauvegarde | Préparation dans un dossier voisin puis renommage ; refus d’écraser ; tests d’erreur de copie |
| Entrées modifiées ou supprimées pendant une session | Copie initiale immuable ; test de suppression du fichier source avant sauvegarde ; hash calculé une seule fois pour les réglages |
| MSA, arbre et rapports périmés après certains changements | Invalidation des résultats dérivés lorsque la sélection des HSP change, et de l’arbre lorsque le MSA est remplacé |
| Erreur générique en cas d’outil absent | Message par outil, plan Conda/Mamba, installation dédiée, diagnostic JSON et découverte du préfixe sans activation |
| Alias de moteurs / modules difficiles à distinguer | Chemins individuels ; échec d’un nom ne masque pas l’autre alias ; test `muscle` / `muscle5` |
| Sous-processus pouvant survivre à une interruption | Groupe de processus séparé, terminaison puis arrêt forcé si nécessaire ; provenance de l’interruption |
| FASTA gzip non lisibles dans le parcours standard | Lecture gzip et entrées FASTA normalisées pour les moteurs ; test d’intégration |
| Changement de page HTML réinitialisant la position verticale | Conservation des décalages ; molette/pavé tactile horizontal, début/fin et tailles de pages ; DOM limité à 100 séquences par page |
| Scripts Slurm utilisant des entrées temporaires | Chemins relatifs vers les entrées du projet sauvegardé et dossier de sortie de relance distinct |
| Numérotation et identité des prototypes dans le produit | Nom AlignView, version publique 0.1.0 partagée par CLI, HTML, exports et distribution |

## Mesures locales reproductibles

Commande :

```bash
python3 scripts/benchmark.py --rows 100 --columns 10000 --repeats 100 \
  --baseline ../AlignView_Pro_v6 --output docs/benchmark-results.json
```

Jeu synthétique de 100 séquences × 10 000 colonnes, terminal fixé à 80×24, message d’état actif. Les [résultats bruts](benchmark-results.json) conservent plateforme, Python et temps.

| Mesure de cadrage | Prototype v6 | AlignView 0.1.0 |
| --- | --- | --- |
| Largeur maximale d’une ligne visible | 139 cellules | 79 cellules |
| Lignes physiques avec retours automatiques | 25 | 23 |

La nouvelle version élimine le dépassement dans ce scénario. La vérification des cellules Unicode ajoute du coût au rendu ; elle ne rend pas le calcul biologique plus rapide. Les temps mesurés, de l’ordre de la milliseconde par image pour ce scénario, sont locaux et ne constituent pas une garantie pour d’autres machines ou jeux de données.

## Validation et portée

Les tests automatisés couvrent calculs, coordonnées traduites, exports, projets, plans d’installation, longues listes, navigation et changements de requête. Les tests d’acceptation complètent les contrôles des touches avec les parcours complets. Les simulations isolent les formats et commandes ; elles ne valident pas le comportement biologique des moteurs réels.

Le rapport autonome a également été exécuté dans Chrome sans fenêtre : pages horizontales, pavé tactile, maintien de la position verticale, dernière colonne, pages de séquences, filtre au-delà de la première page et jeu sans colonne variable. Le script reproductible est `scripts/check_report_browser.py`.

Les [essais réels locaux](real-tools-results.json) ont réussi pour MAFFT, BLASTN, DIAMOND, MUSCLE 5, Clustal Omega, EMBOSS Needle, FastTree et MMseqs2, avec vérification des sorties et création de projets. Un [deuxième passage](tuned-tools-results.json) utilise des paramètres scientifiques modifiés. Le script `scripts/validate_real_tools.py` et le workflow manuel `real-tools.yml` permettent de reproduire ces contrôles. La matrice Python de GitHub est préparée ; aucun résultat de CI distante n’est revendiqué avant exécution dans le dépôt créé.

Le wheel, le sdist et le zipapp doivent être reconstruits après tout changement. La procédure est dans [RELEASING.md](RELEASING.md). La recette Bioconda reste un modèle tant que le dépôt public et le SHA256 de la release ne sont pas connus. Son organisation suit les [directives Bioconda](https://bioconda.github.io/contributor/guidelines.html) ; l’installation des outils suit les [canaux recommandés](https://bioconda.github.io/#with-conda).

## Limites à conserver visibles

- Les FASTA et le MSA restent chargés en mémoire ; pas de promesse de taille illimitée ou d’indexation de fichiers géants.
- L’import interactif reste FASTA/FASTA gzip. NEXUS, PHYLIP, Stockholm et Clustal sont des formats d’export.
- Une recherche compare une requête active à des références ; pas de planificateur all-vs-all.
- Le rapport est plafonné à 1 000 séquences de MSA et 2 000 HSP, avec mention de cette limite. Les exports applicatifs sont complets.
- Les alignements complets TBLASTN/TBLASTX nécessitent encore une sélection biologique explicite des régions/cadres ; aucune traduction globale arbitraire n’est ajoutée.
- ORF ATG→STOP et FastTree restent exploratoires ; les distances p ne sont pas un modèle évolutif corrigé.
- Les modules HTML historiques sont maintenus. Une refonte ultérieure de leur construction devra préserver les calculs et formats déjà vérifiés.
- Les versions finales de tous les moteurs, le cluster cible, macOS et les logiciels tiers recevant les exports restent à valider avant de présenter une version 1.0 stabilisée.
