# Utiliser AlignView

Pour une première découverte, commencer par [les premiers pas](QUICKSTART.md). Pour la conservation et les chemins, lire [Dossiers et sauvegardes](FILES.md). Dans les commandes ci-dessous, `alignview` peut être remplacé par `/chemin/vers/AlignView/alignview.sh`.

## Trouver le bon parcours

Au démarrage, choisir le résultat attendu :

- **Ouvrir un alignement** : FASTA déjà aligné. Aucun moteur n’est lancé.
- **Aligner des séquences** : FASTA brut. Le calcul se lance ensuite depuis Analyser.
- **Comparer à des références** : fichier de requêtes, puis fichier de références.
- **Reprendre un projet** : réouvrir des données et calculs enregistrés.
- **Essayer un exemple** : découvrir la visionneuse sans installer de moteur.

Une longueur égale entre séquences ne prouve pas qu’un fichier est aligné. « Ouvrir un alignement » exprime cette connaissance de l’utilisateur. Si les longueurs diffèrent, le fichier reste chargé et AlignView propose de créer un alignement.

## Se repérer

Le fil en haut indique sa place : `Accueil › Analyser › Alignement`. Le contexte indique la séquence active, les références, les résultats et la sauvegarde. Un choix grisé explique son prérequis ; le sélectionner ne lance pas de calcul inadapté.

Échap et ← reviennent d’un niveau. F2 revient à l’accueil, même depuis la visionneuse et les résultats. F1 affiche les détails d’un choix. Ctrl+C annule une saisie ou une action ; le projet reste ouvert. `q` à l’accueil propose l’enregistrement si le travail a changé.

## Lire les résultats

**Explorer → Lire l’alignement** ouvre le MSA. Flèches pour se déplacer, `[` / `]` pour un écran horizontal, PgUp/PgDn pour les séquences, `g` pour une colonne ou un pourcentage, `s` pour un identifiant et `/` pour un motif. `v` masque les colonnes non variables, `d` souligne les différences. `?` donne toutes les commandes. Les coordonnées restent celles du MSA avec ses gaps.

**Explorer → Correspondances et couvertures** regroupe les références, les régions détectées et les diagnostics. Choisir une référence ouvre ses segments alignés. `/` filtre les identifiants.

**Explorer → Statistiques de l’alignement** donne conservation, consensus, distances et arbre exploratoire. Dans la matrice, les indices renvoient aux identifiants affichés. Les pages larges se déplacent avec ←/→ ; les pages longues avec PgUp/PgDn.

Après **Arbre exploratoire → Calculer l’arbre**, AlignView crée le Newick et une figure vectorielle `figures/45_arbre_fasttree.svg`. Créer ensuite le rapport HTML depuis **Enregistrer et partager** : **Voir l’arbre** ouvre directement la figure, également accessible dans **Figures → Phylogénie**. **Télécharger SVG** conserve la figure seule ; le HTML contient l’arbre et reste autonome.

Les longueurs horizontales et la barre d’échelle représentent des substitutions par site. Les identifiants et les supports internes du Newick sont conservés ; les supports locaux FastTree (SH-like) ne sont pas des bootstraps. La racine du dessin est arbitraire, sans enracinement biologique. En l’absence de longueurs utilisables, la figure indique un cladogramme sans échelle de distance. Un projet repris contenant déjà `tree_fasttree.nwk` obtient sa figure en reconstruisant le HTML, sans relancer FastTree.

```bash
alignview sequences_alignees.fasta --msa view --tree --report --save --out projet_arbre --no-menu
```

## Choisir une séquence ou des données

**Données → Choisir la séquence active** sélectionne la cible des motifs, ORF, traductions, extractions et profils. Dans un FASTA unique, cette sélection conserve le fichier complet pour l’alignement multiple.

Avec deux FASTA, **Changer la requête active** ouvre un autre contexte de comparaison. Les analyses précédentes restent dans **Données → Revenir à une analyse précédente**. Cet historique est temporaire jusqu’à l’enregistrement.

**Ouvrir d’autres fichiers** propose les mêmes parcours qu’au démarrage. **Ouvrir un projet enregistré** reprend les calculs sans les chemins d’origine. Les entrées sont vérifiées par SHA256.

## Régler la recherche

La limite initiale est **100 références retenues**. Les profils ne changent pas ce nombre :

| Profil | E-value maximale | Sensibilité DIAMOND |
| --- | --- | --- |
| Standard | 1e-5 | sensitive |
| Plus sensible | 1e-2 | very-sensitive |
| Plus strict | 1e-10 | sensitive |

Identité, couverture et bitscore filtrent les résultats déjà calculés. Si la sélection change, l’alignement dérivé est invalidé. E-value et paramètres du moteur s’appliquent à la prochaine recherche ; le contexte indique qu’une relance est nécessaire. Les paramètres du calcul effectué restent distincts dans la provenance.

Dans chaque écran de préparation, **Paramètres de calcul / Paramètres du moteur** permet de modifier les contrôles scientifiques applicables. **Voir les commandes** inspecte le calcul sans lancer de moteur. **Profils de calcul** charge ou enregistre un JSON réutilisable. Les réglages avancés et leurs effets sont décrits dans le [guide des paramètres](PARAMETERS.md).

Les profils sont statistiques, pas des seuils taxonomiques universels. Les projets existants conservent leur propre limite.

## Préparer un moteur manquant

Avant un calcul, AlignView vérifie le moteur nécessaire. S’il manque, l’écran propose son installation ou les possibilités module HPC / chemin local. **Réglages → Logiciels et installation** donne le diagnostic complet. Le moteur reste visible dans l’aperçu du calcul et modifiable dans les réglages.

## Conserver le travail

Dans **Enregistrer et partager** :

- **Projet** : garder entrées, paramètres, calculs, alignements et historique pour reprendre.
- **Rapport HTML** : partager une exploration autonome sans serveur ni réseau.
- **Exports** : réutiliser les données dans R, Python ou d’autres logiciels.
- **Script HPC** : préparer un Slurm ; aucune soumission automatique.

Un enregistrement réussi affiche « projet enregistré ». Les modifications suivantes réactivent « modifications à enregistrer ». Fermer un projet enregistré ne redemande pas une seconde copie. Les copies déjà enregistrées sont conservées si les nouvelles modifications sont abandonnées.

Le HTML propose des points d’entrée « Par où commencer ? ». Les vues vides sont masquées ; Paramètres et Comprendre sont regroupés. Les figures se déplient à la demande. Résumé et l’historique du navigateur permettent de retrouver une vue.

## Automatiser

```bash
alignview sequences.fasta --msa align --report --no-menu --save --out msa_projet
alignview query.fna refs.faa --mode blastx --msa hits --report --no-menu --save --out recherche_projet
alignview query.fna refs.faa --mode blastx --query-id ma_requete --hits 200 --no-menu --save --out projet
```

Les options explicites et le batch restent utilisables. Les valeurs par défaut sont communes à la CLI, aux menus et à la sélection de références.

## Limites scientifiques

La référence fonctionnelle ci-dessous complète les parcours d’introduction. Les limites restent applicables à chaque analyse.

La couverture est l’union des segments ; l’identité est locale à un segment. BLASTX conserve nt/aa/aa, TBLASTN aa/nt/aa, TBLASTX nt/nt/aa. Les dotplots traduits ont des axes en aa.

Les ORF ATG→STOP et FastTree restent exploratoires. TBLASTN/TBLASTX nécessite une sélection biologique des régions et cadres avant un alignement intégral. La recherche cible une requête active, pas un calcul all-vs-all.

## Choisir le mode de recherche

Le premier FASTA fournit la requête active et le second les références. Si le premier contient plusieurs séquences, une seule est recherchée à la fois ; en batch, fournir `--query-id` avec le premier mot de son en-tête. Changer de requête dans Données ouvre un autre contexte de comparaison.

| Requête | Références | Mode | Coordonnées Q / référence / HSP |
| --- | --- | --- | --- |
| ADN/ARN | ADN/ARN | `blastn` | nt / nt / nt |
| Protéine | Protéines | `blastp`, `diamond-blastp`, `mmseqs-blastp` | aa / aa / aa |
| ADN/ARN | Protéines | `blastx`, `diamond-blastx` | nt / aa / aa |
| Protéine | ADN/ARN | `tblastn` | aa / nt / aa |
| ADN/ARN | ADN/ARN, comparaison traduite | `tblastx` | nt / nt / aa |

Les modes `global-dna` et `global-aa` chargent deux FASTA homogènes sans recherche locale initiale ; les comparaisons globales se lancent dans Analyser. Un FASTA unique sert à l’exploration et au MSA et ne nécessite pas de second fichier artificiel.

L’écran de recherche propose **Choisir la méthode**, **Régler la recherche**, **Paramètres du moteur**, **Voir les commandes** et **Profils de calcul** avant lancement. BLAST peut nécessiter `makeblastdb` en plus de l’exécutable de recherche ; l’index est préparé dans la session temporaire.

## Choisir les séquences du MSA

Un fichier unique utilise toutes ses séquences pour le MSA, même lorsqu’une autre séquence active est choisie. Avec deux fichiers :

| Portée | Contenu | Usage |
| --- | --- | --- |
| Références retenues / `--msa hits` | Requête et références retenues après recherche | Comparer les homologues détectés |
| Toutes les références / `--msa all` | Requête et toutes les références | Comparer un ensemble déjà choisi biologiquement |
| Ajouter une séquence partielle / `--msa fragments` | Références alignées puis ajout de la requête | Requête courte ou fragment ; deux étapes MAFFT explicites |
| Segments détectés / `--msa hsp` | Segments recouvrant suffisamment la région ancre | Explorer des domaines/régions locales sans les extrémités |
| Ouvrir / `--msa view` | Un FASTA unique déjà aligné | Lire sans recalcul |

`--msa align` est la forme pratique pour un FASTA unique brut. Avec deux FASTA en recherche, `align`, `full` et `global` sont des alias historiques du MSA des références retenues ; `all` demande explicitement toutes les références. Dans les modes globaux, toutes les références sont utilisées.

Les HSP peuvent ne pas représenter exactement le même domaine ; régler leur seuil d’ancre avec `--hsp-overlap` et conserver une interprétation exploratoire. Les fragments imposent MAFFT en deux étapes, ce qui est annoncé même lorsqu’un autre moteur MSA était précédemment sélectionné.

En BLASTX, une requête entière peut être traduite dans un cadre choisi pour un MSA protéique : vérifier la préparation et `--frame`. Pour TBLASTN/TBLASTX, utiliser les segments détectés ou préparer biologiquement les protéines/régions appropriées. Un nouveau MSA invalide l’arbre et ses vues dérivées ; les réglages futurs seuls ne remplacent pas les paramètres historiques d’un résultat conservé.

## Toutes les touches de la visionneuse

| Action | Commande |
| --- | --- |
| Colonne / séquence | ← → / ↑ ↓ ou h l / k j |
| Dix colonnes | Maj+← / Maj+→ |
| Un écran de colonnes | Ctrl+← / Ctrl+→ ou `[` / `]` |
| Une page de séquences | PgUp / PgDn |
| Première / dernière colonne | Home / End |
| Première / dernière séquence | Ctrl+Home / Ctrl+End |
| Colonne ou pourcentage | `g`, puis `500` ou `75%` |
| Séquence par identifiant | `s` |
| Motif exact / motif du consensus | `/`, puis motif ou `c:motif` |
| Occurrence suivante / précédente | `n` / `N` |
| Sites variables / différences au consensus | `v` / `d` |
| Rapport HTML / export du HTML | `w` / `e` |
| Aide / retour / accueil | `?` / q ou Échap / F2 |

La molette déplace les séquences, Maj+molette les colonnes ; cliquer dans la vue globale atteint une région. Les coordonnées restent celles du MSA. Masquer les sites invariants ne retire pas de colonnes du fichier original.

## Needle et Water : comparer deux séquences

**Analyser → Comparer deux séquences** propose **De bout en bout** (Needle) ou **Le meilleur segment local** (Water). Choisir une référence, vérifier les paramètres et la commande, puis lancer.

En batch : `--pairwise global` ou `--pairwise local`, avec `--ref N` pour l’indice de référence. Ces calculs utilisent EMBOSS. Ils ne remplacent ni une recherche BLAST, ni un MSA. Needle pénalise les extrémités par défaut dans AlignView ; `--free-ends` les relâche. Les coûts terminaux et internes se règlent séparément.

`--gapopen`, `--gapextend` et `--matrix` sont des paramètres EMBOSS. Ils ne règlent ni MAFFT, ni BLAST, ni DIAMOND. Dans un mode traduit, vérifier les cadres et la traduction complète annoncés par l’aperçu.

## Motifs, traductions, ORF, profils et dotplots

Les analyses interactives ciblées de **Analyser → Autres analyses** utilisent la **séquence active**, choisie dans Données.

| Action | Saisie | Résultat de session |
| --- | --- | --- |
| Rechercher un motif | Motif IUPAC, ex. `ATG` ou `NNNGTT` | Positions sur les deux brins, `motifs.tsv` |
| Traduire les six cadres | Code génétique standard | Six protéines, `six_frames.faa` |
| Rechercher des ORF simples | Longueur minimale en aa | ORF complètes ATG → STOP, `orfs.tsv`, `orfs.faa` |
| Extraire une région | Début et fin inclusifs, complément inverse éventuel | `query_DEBUT_FIN.fasta` |
| Qualité et composition globale | Sans moteur externe | Longueurs, GC pour ADN/ARN, ambiguïtés, entropie |
| Profil le long de la séquence | `gc`, `entropy` ou `gap`, fenêtre et pas | `profile_METRIQUE.tsv` et `.svg` |
| Fréquences de mots | k de 1 à 8 | Comptages et fréquences, `kmers_K.tsv` |
| Dotplot | Seconde séquence et longueur de mot k | `dotplot.json` et `dotplot.svg` |

Les extractions sont 1-based inclusives. Le complément inverse concerne les acides nucléiques. Les ORF restent des détecteurs ATG→STOP simples. Le dotplot affiche des mots exactement identiques ; un petit k augmente généralement les points et les répétitions, et le logiciel peut limiter le nombre de points pour préserver l’affichage.

Les aperçus de motifs, ORF, profils et k-mers peuvent être plus courts que les fichiers. Pour conserver toutes les analyses internes et auxiliaires, enregistrer le **projet complet** : un profil d’export ne copie pas nécessairement chaque fichier de la session. En batch, `--motif` analyse les séquences d’un FASTA unique, ou la requête active en mode comparaison.

## Statistiques et provenance

**Conservation et consensus** calcule les statistiques par colonne et le consensus. **Distances entre séquences** calcule les distances p pour au plus 250 séquences ; les gaps sont exclus par paire, donc toutes les paires ne sont pas forcément comparées sur les mêmes positions.

FastTree requiert au moins trois séquences alignées. Choisir un modèle adapté à l’alphabet, puis lire le Newick et la figure. [Le guide du HTML](HTML.md) explique son échelle, sa racine arbitraire et ses supports locaux.

**Calculs effectués et provenance** relie les résultats réussis aux paramètres, commandes et versions qui les ont produits. Les journaux stderr distincts permettent de relire plusieurs tentatives. Les paramètres actuels peuvent différer ; une relance reste nécessaire pour produire un résultat avec ces nouveaux paramètres.

## Documentation et calcul sur HPC

La documentation est accessible au démarrage et dans **Réglages et aide → Documentation complète**. `alignview --docs` affiche l’index hors ligne ; `--docs cli` donne les options complètes.

**Enregistrer et partager → Préparer un calcul sur HPC** crée `alignview_job.sbatch`. Enregistrer le projet pour conserver les entrées auxquelles il se réfère. Vérifier ressources Slurm, chemins et moteurs sur le cluster. Aucune soumission n’est automatique ; la relance produit ses résultats dans un dossier distinct du projet sauvegardé.

## VSEARCH, langue et codes de traduction

**Analyser → Analyses VSEARCH** donne quatre parcours nucléotidiques indépendants : recherche globale sur la requête active, déduplication, regroupement et chimères de novo. **Explorer → Résultats VSEARCH** permet de relire les résumés, aperçus et paramètres ; le HTML ajoute une vue VSEARCH. Voir [le guide complet](VSEARCH.md).

La sélection de requête est proposée avant toute comparaison interactive multirequête, y compris le parcours guidé. En batch, `--query-id` reste obligatoire.

**Réglages et aide → Langue / Language** change la langue de présentation. `--lang en` ou `ALIGNVIEW_LANG=en` active l’anglais depuis le lancement. Les identifiants et données scientifiques ne sont pas traduits. Les [guides anglais](en/README.md) sont aussi embarqués.

Le choix de table se trouve dans **Traduction et régions codantes → Choisir la table de traduction**, et dans les paramètres des calculs traduits. `--genetic-code` règle QUERY et `--reference-genetic-code` les références. Voir [les limites et exemples](GENETIC_CODES.md).
