# Application, sorties et sauvegardes

## Le dossier à conserver pour l’application

Le dossier **`AlignView/`** contient les sources, le lanceur, les exemples, les guides, les scripts de fabrication, les tests et les distributions préparées. Il peut être déplacé ou copié entier. Il ne dépend d’aucun fichier Python situé dans son dossier parent.

```text
AlignView/
  alignview.sh          lanceur principal
  alignview.py          entrée Python
  av_*.py               fonctions de l’application
  alignview.pyz         application autonome en un fichier, documentation comprise
  README.md             point d’entrée des guides
  docs/                 guides, captures et rapports de validation
  examples/             petits FASTA d’exemple
  scripts/              construction et vérification
  tests/                tests et outils simulés
  packaging/            préparation Bioconda
  dist/                 wheel, sources et archive portable
  archives/             anciennes versions locales, facultatives et inutilisées
  tools/                moteurs Conda locaux, uniquement si installés ici
```

Les anciennes versions sont rangées dans `archives/` lorsqu’elles sont présentes. Elles ne sont ni chargées par l’application courante, ni incluses dans les distributions. Les jeux de données et les résultats voisins restent indépendants. Le raccourci `alignview.sh` éventuellement présent à côté d’`AlignView/` pointe vers son lanceur ; il est facultatif et inutile lorsqu’on copie le dossier.

Python et les moteurs scientifiques sont des programmes externes. Copier les sources d’AlignView suffit à déplacer **l’application**, mais ne copie pas une installation système de BLAST, MAFFT ou FastTree. Voir [installation et moteurs locaux](INSTALLATION.md).

## Quelle base pour les chemins ?

Tous les chemins relatifs d’entrée et de sortie partent du **dossier courant du terminal**, appelé aussi dossier de lancement. `pwd` l’affiche. Le lanceur localise ses sources sans changer ce dossier.

```bash
cd /travail/mon_analyse
/logiciels/AlignView/alignview.sh seqs.fasta --msa view --report --save --out resultats/lecture_01 --no-menu
```

Ici, l’entrée est `/travail/mon_analyse/seqs.fasta`, le projet est `/travail/mon_analyse/resultats/lecture_01`, et le rapport est dans ce projet. Rien n’est enregistré dans `/logiciels/AlignView` par cette commande. Si l’on lance depuis l’intérieur de l’application, les sorties relatives se trouveront naturellement dans ce dossier : choisir un autre dossier de travail pour les séparer.

Un chemin absolu, par exemple `--out /scratch/moi/projet_01`, utilise directement cette destination. Les espaces dans un chemin doivent être entourés de guillemets dans le shell : `--out "mes resultats/projet_01"`. Dans une saisie de l’interface, coller le chemin directement.

## Ce que chaque option conserve

| Commande ou action | Écriture durable | Reprise des calculs possible ? |
| --- | --- | --- |
| Exploration seule, sans sauvegarde | Aucune | Non après fermeture |
| `--report` seul | Rapport préparé dans la session temporaire | Non après fermeture sans autre action |
| `--save` | Projet `alignview_YYYYMMDD_HHMMSS/` dans le dossier courant | Oui |
| `--save --out dossier` | Projet dans `dossier/` | Oui |
| `--out dossier` seul | Aucune sauvegarde automatique | Non après fermeture |
| `--html-out rapport.html` | Un HTML autonome | Consultation uniquement |
| `--export-dir exports` | Les exports du profil choisi | Réutilisation externe, pas un projet complet |
| `--export-profile publication` sans destination ni sauvegarde | Prépare un choix de fichiers ; rien n’est conservé seul | Non |

Ces options peuvent se combiner. `--html-out` reconstruit le rapport et n’a pas besoin de `--report`. Le profil publication/complet construit aussi le HTML lorsqu’un export ou une sauvegarde est effectivement demandé.

## Un projet sauvegardé

Le dossier comprend notamment :

```text
mon_projet/
  session.json               état nécessaire à la reprise
  inputs/                    copies des FASTA utilisés
  alignment_active.fasta     alignement conservé, si disponible
  parameters.json            réglages actuels et historique des calculs
  commands.json              commandes et versions des moteurs
  logs/                      journaux des calculs
  tree_fasttree.nwk          arbre, si calculé
  report.html                rapport, si généré
  figures/                   figures SVG, si générées
  exports/                   fichiers du profil choisi
  history/                   analyses antérieures, si présentes
```

D’autres fichiers intermédiaires et résultats peuvent être présents. Pour transférer un projet, copier **tout le dossier**, pas seulement `session.json`. Les entrées ont leurs SHA256 ; si elles sont modifiées après sauvegarde, la reprise le signale. Déplacer ou supprimer les fichiers d’origine n’empêche pas de reprendre le projet.

```bash
alignview --resume /chemin/vers/mon_projet
```

La reprise copie le projet dans une nouvelle session temporaire et conserve les résultats existants. Les moteurs ne sont pas relancés par défaut. Enregistrer ensuite crée une nouvelle copie ; le dossier original n’est pas écrasé.

## Reprendre pour exporter ou recalculer

```bash
# Reconstruire le rapport avec la version actuelle de l’application.
alignview --resume mon_projet --html-out actualise.html --no-menu

# Préparer de nouveaux exports, sans refaire de calcul scientifique.
alignview --resume mon_projet --export-profile complete --export-dir exports_02 --no-menu

# Refaire explicitement le MSA puis l’arbre et conserver une nouvelle copie.
alignview --resume mon_projet --msa align --tree --report --save --out projet_02 --no-menu
```

Les actions `--msa`, `--tree`, `--pairwise`, `--dotplot` et `--motif` ne sont lancées à la reprise que si elles figurent explicitement dans la commande. Des options de paramètres seules, comme `--mafft-op`, préparent les réglages futurs ; elles ne refont pas un résultat existant. Changer la requête ou la méthode de recherche se fait dans **Données** ou **Analyser**.

## Écrasement, erreurs et temporaires

Une sauvegarde prépare un dossier voisin puis le renomme ; elle ne fusionne pas avec un projet existant. En batch, une destination existante fait échouer la commande. En interactif, le logiciel propose un autre nom ou demande confirmation lorsqu’un HTML existant est choisi.

Les calculs temporaires utilisent le dossier temporaire de Python, généralement `/tmp`. Pour de grosses analyses sur un cluster, on peut choisir avant de lancer :

```bash
mkdir -p /scratch/mon_compte/tmp
export TMPDIR=/scratch/mon_compte/tmp
alignview sequences.fasta --msa align --save --out /scratch/mon_compte/projet_01 --no-menu
```

Une interruption ou erreur après des calculs peut créer un projet de récupération `alignview_recovery_…/project`. Le chemin est alors affiché : le copier ou le reprendre avant le nettoyage des temporaires du système. Un arrêt brutal ou une panne ne garantit pas cette récupération. L’historique dans l’application reste lui aussi temporaire tant qu’il n’a pas été enregistré.
