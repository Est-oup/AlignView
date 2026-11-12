# Premiers pas

## 1. Lancer depuis le dossier où travailler

L’application est dans **un seul dossier `AlignView/`**. Pour commencer, Python 3.10+ suffit. Choisir son dossier de travail, puis lancer l’application par son chemin :

```bash
cd /chemin/vers/mes_donnees
/chemin/vers/AlignView/alignview.sh --demo
```

Les chemins ci-dessus sont à remplacer par les tiens. Le lanceur conserve le dossier de travail ; il ne se déplace pas dans celui de l’application. Depuis l’intérieur d’`AlignView/`, `./alignview.sh --demo` fonctionne aussi.

La démonstration contient quatre séquences déjà alignées. Elle n’installe rien et ne lance aucun moteur scientifique. Dans **Explorer → Lire l’alignement**, déplacer le curseur avec les flèches. `?` affiche l’aide ; Échap revient au menu et F2 à l’accueil.

## 2. Choisir son intention

Lancer sans `--demo` pour sélectionner ses fichiers :

```bash
/chemin/vers/AlignView/alignview.sh
```

| Choix de départ | Fichier à fournir | Étape suivante |
| --- | --- | --- |
| Ouvrir un alignement | FASTA déjà aligné | Explorer → Lire l’alignement |
| Aligner des séquences | FASTA brut contenant au moins deux séquences | Analyser → Créer un alignement → Lancer |
| Comparer à des références | Un FASTA de requêtes puis un FASTA de références | Analyser → Rechercher des similarités |
| Reprendre un projet | Dossier contenant `session.json` | Retrouver les résultats dans Explorer |

Le sélecteur montre les fichiers FASTA et FASTA gzip. `/` filtre les noms ; **Saisir un chemin** permet de coller un chemin complet. Dans une longue liste, PgUp/PgDn accélère le déplacement. Une séquence se reconnaît à son identifiant placé juste après `>` dans le FASTA.

Deux séquences de même longueur ne sont pas forcément alignées. Choisir **Ouvrir un alignement** seulement si le fichier constitue déjà un alignement ; sinon demander sa création.

## 3. Savoir ce qui est prêt

L’accueil garde cinq entrées : **Explorer**, **Analyser**, **Données**, **Enregistrer et partager**, **Réglages et aide**. Le contexte indique les données actives, les résultats disponibles et l’état de sauvegarde.

Un choix indisponible indique son prérequis. Par exemple, un arbre nécessite au moins trois séquences alignées et le logiciel FastTree. Avant un calcul, le moteur et les réglages sont visibles ; **Voir les commandes** permet de vérifier ce qui serait exécuté.

Si un moteur manque, AlignView propose son installation ou un chemin d’exécutable. La lecture d’un alignement existant reste possible sans ce moteur. Voir [installer les outils](INSTALLATION.md).

## 4. Conserver ce qui compte

Les calculs se font dans une session temporaire. Dans **Enregistrer et partager**, choisir :

| Besoin | Action | Résultat conservé |
| --- | --- | --- |
| Reprendre plus tard | Enregistrer le projet | Dossier portable avec entrées et résultats |
| Partager une exploration | Créer un rapport HTML | Un fichier `.html` autonome |
| Réutiliser les données ailleurs | Exporter les résultats | Dossier de tables, alignements et figures selon le profil |

Un chemin comme `resultats/essai_01` est relatif au **dossier de lancement**, pas au dossier de l’application. Choisir une destination libre : un ancien projet n’est pas écrasé. Le HTML à lui seul ne permet pas de reprendre les calculs dans AlignView.

Fermer sans enregistrer supprime les temporaires. L’interface demande quoi faire si le travail a changé. Une copie de projet déjà enregistrée reste conservée.

## 5. Obtenir un arbre et sa figure

Ouvrir ou créer un alignement, puis **Explorer → Statistiques de l’alignement → Arbre exploratoire → Calculer l’arbre**. Ensuite créer le HTML. Son bouton **Voir l’arbre** ouvre la figure ; **Télécharger SVG** enregistre une image vectorielle.

La figure donne les identifiants, les branches, les supports internes et une échelle lorsque les longueurs sont présentes. La racine du dessin est arbitraire. Les supports locaux FastTree sont des SH-like, pas des bootstraps. Voir [interpréter le rapport](HTML.md).

## Des commandes à adapter

Les exemples suivants utilisent la commande `alignview` installée. Sans installation, remplacer `alignview` par `/chemin/vers/AlignView/alignview.sh`. Les entrées sont dans le dossier de travail et les destinations doivent être libres.

```bash
# Lire un alignement existant, puis garder seulement le HTML.
alignview aligne.fasta --msa view --html-out rapport.html --no-menu

# Aligner des séquences et garder un projet complet avec HTML.
alignview sequences.fasta --msa align --report --save --out resultats/msa_01 --no-menu

# Comparer une requête ADN à des références protéiques.
alignview requete.fna references.faa --mode blastx --msa hits --report --save --out resultats/recherche_01 --no-menu

# Calculer l’arbre d’un alignement existant et conserver le projet.
alignview aligne.fasta --msa view --tree --report --save --out resultats/arbre_01 --no-menu

# Reconstruire uniquement le HTML d’un projet ; aucun moteur relancé.
alignview --resume resultats/arbre_01 --html-out rapport_actualise.html --no-menu
```

`--no-menu` exécute puis termine. `--save` conserve le projet ; `--out` seul ne le fait pas. `--html-out` conserve le HTML sans sauvegarder un projet. Ces différences sont détaillées dans [Dossiers et sauvegardes](FILES.md).

Pour continuer : [toutes les fonctions](USAGE.md), [les paramètres](PARAMETERS.md), [le dépannage](TROUBLESHOOTING.md).

## Essayer VSEARCH avec les exemples fournis

Depuis le dossier de l’application : `./alignview.sh examples/example_vsearch.fna --vsearch cluster --html-out vsearch_demo.html --no-menu`. Pour choisir une requête et la comparer : `./alignview.sh examples/example_vsearch.fna examples/example_vsearch_refs.fna --vsearch search`, puis choisir l’identifiant dans la liste. Ces données sont synthétiques et servent à découvrir le parcours.
