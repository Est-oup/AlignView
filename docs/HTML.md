# Lire, créer et partager le rapport HTML

## Créer un fichier qui reste disponible

Dans l’application : **Enregistrer et partager → Créer un rapport HTML**, puis choisir le nom du fichier. Le HTML est copié à cet endroit et reste consultable après fermeture de la session.

En ligne de commande :

```bash
# Un seul HTML à partir d’un alignement existant.
alignview aligne.fasta --msa view --html-out rapport.html --no-menu

# HTML et tous les résultats dans un projet.
alignview sequences.fasta --msa align --report --save --out mon_projet --no-menu

# Refaire le HTML d’un projet, sans relancer les moteurs.
alignview --resume mon_projet --html-out actualise.html --no-menu

# Tables et HTML dans un dossier d’export, sans projet de reprise.
alignview --resume mon_projet --export-profile publication --export-dir a_partager --no-menu
```

`--html-out` suffit à demander la génération. `--report` seul prépare un fichier temporaire ; il faut le conserver avec `--save`, `--html-out` ou un export publication/complet. Les destinations relatives partent du dossier de lancement. Un fichier existant est refusé en batch.

Ouvrir le fichier dans un navigateur avec JavaScript activé. Le rapport contient ses données, son style, ses scripts et ses SVG ; il ne contacte pas de serveur et ne dépend pas d’un CDN. On peut copier ou envoyer le fichier seul. Il contient des séquences et des identifiants : partager avec les personnes auxquelles ces données sont destinées.

## Navigation et lecture

Le résumé propose **Par où commencer ?** et des accès aux données présentes. Les vues sans données ne sont pas proposées. Le bouton **Résumé** et l’historique du navigateur permettent de revenir à une vue précédente.

| Vue | Ce qu’elle montre | Comment la lire |
| --- | --- | --- |
| Résumé | Données, compteurs, qualité et couverture | Vérifier le mode et le nombre de références retenues |
| Correspondances | Références et HSP | Filtrer ou trier, puis choisir une référence et un HSP |
| Séquences | Longueurs, GC, ambiguïtés, profils | Le GC concerne les séquences nucléotidiques |
| Alignement | Visionneuse MSA et métriques | Déplacer les colonnes et parcourir les pages de séquences |
| Dotplot | Points de mots identiques | Une diagonale est descriptive ; ce n’est pas un nouvel alignement |
| Figures | Atlas SVG conditionnel | Filtrer par nom et par domaine, puis déplier une figure |
| Exporter | Tables et séquences embarquées | Téléchargements locaux, limités aux données intégrées au HTML |
| Paramètres et aide → Paramètres | Réglages actuels, calculs effectués et commandes | Distinguer les réglages préparés des paramètres historiques |
| Paramètres et aide → Comprendre | Méthodes et limites | Lire avant une interprétation biologique |

## Lire un alignement long

Les commandes de la visionneuse déplacent les colonnes et les pages de séquences séparément. Le défilement vertical reste conservé lorsqu’on change de page horizontale. Le pavé tactile horizontal fonctionne aussi. Un filtre d’identifiant permet d’atteindre une séquence située au-delà de la première page.

Une page comporte au plus 100 séquences. Le rapport intègre au plus 1 000 séquences de MSA. Le filtre et les téléchargements portent sur les données intégrées ; ils ne chargent pas les lignes absentes du fichier. La vue des sites variables ne modifie pas le MSA.

## L’arbre

Calculer d’abord l’arbre dans **Explorer → Statistiques de l’alignement → Arbre exploratoire**, ou ajouter `--tree` à une commande qui ouvre/crée un MSA. Il faut au moins trois séquences alignées et FastTree pour ce calcul.

Le HTML propose alors **Voir l’arbre**, qui ouvre **Figures → Phylogénie** et déplie la figure. Pour un arbre très grand, faire défiler son cadre. **Télécharger SVG** conserve la figure entière, y compris les parties qui dépassent du cadre à l’écran. Le projet contient aussi `figures/45_arbre_fasttree.svg` et `tree_fasttree.nwk`.

L’horizontale représente les longueurs de branches en **substitutions par site**. La verticale organise les étiquettes et ne constitue pas une échelle biologique. La racine du dessin est arbitraire : elle ne désigne pas un ancêtre inféré ni un groupe externe. Les étiquettes internes du Newick sont conservées ; les supports FastTree locaux SH-like ne sont pas des bootstraps.

Si les longueurs sont absentes ou toutes nulles, le dessin est explicitement un **cladogramme**, sans échelle de distance. Un Newick mal formé produit un message d’indisponibilité ; une ancienne figure n’est pas réutilisée. Un projet ancien contenant déjà un Newick obtient sa figure en reconstruisant le HTML, sans réexécuter FastTree.

## Quelles figures sont possibles ?

L’atlas dépend des données. Une figure absente ne signifie pas automatiquement une erreur.

| Domaine | Exemples de figures | Prérequis |
| --- | --- | --- |
| Séquences | Composition, longueurs, GC, motifs, ORF | Séquences et analyse correspondante lorsqu’elle est requise |
| HSP et coordonnées | Longueurs, nombres, étendues, pistes par référence | Recherche ayant retenu des HSP |
| Couvertures | Union des segments, couverture par référence, couverture cumulée | HSP et longueurs disponibles |
| Qualité des HSP | Identité, substitutions, gaps, bitscores, E-values | HSP retenus |
| Traduction | Cadres, contrôle du rapport 3 nt / aa | Recherche traduite et cadres disponibles |
| MSA | Conservation, entropie, gaps, occupation, consensus, distances p | Alignement disponible ; certaines matrices sont limitées à 45 séquences |
| Phylogénie | Arbre FastTree | Newick calculé et valide |

Le catalogue `figures_catalog.tsv` associe chaque SVG à sa catégorie, sa description et ses unités. Les histogrammes de HSP décrivent des segments locaux : un HSP n’est pas une séquence entière et plusieurs segments d’une même référence peuvent se chevaucher.

## Téléchargements, impression et limites

Chaque figure propose un **SVG vectoriel**, utile pour modifier les dimensions sans pixellisation. **Imprimer / PDF** utilise le navigateur ; les figures se déplient à l’impression. Pour un arbre très grand, le SVG autonome reste le moyen de conserver l’ensemble lisiblement.

Les téléchargements du navigateur utilisent seulement le contenu du HTML : au plus **2 000 HSP**, **301 statistiques de séquences** et **1 000 séquences de MSA**. Pour les tableaux complets, utiliser les exports depuis l’application. Le HTML est un instantané : modifier les paramètres dans AlignView ne modifie pas automatiquement un fichier HTML déjà enregistré. Le reconstruire après une nouvelle analyse.

Le rapport affiche des résultats, mais n’exécute ni MAFFT, ni BLAST, ni FastTree. Pour changer un calcul, revenir dans AlignView. Pour reprendre le travail, conserver aussi le dossier de projet.
