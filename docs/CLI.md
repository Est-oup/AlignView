# Référence de la ligne de commande

Cette page est générée depuis le parseur réel par `python3 scripts/generate_cli_docs.py`. Ne pas modifier ses tables à la main.

Dans les exemples, `alignview` désigne la commande installée. Sans installation, remplacer ce mot par `/chemin/vers/AlignView/alignview.sh`.

## Choisir une commande

```bash
alignview --demo
alignview --docs demarrage
alignview aligne.fasta --msa view
alignview sequences.fasta --msa align --report --save --out mon_projet --no-menu
alignview requetes.fna references.faa --mode blastx --query-id mon_identifiant --save --out recherche --no-menu
alignview aligne.fasta --msa view --tree --html-out arbre.html --no-menu
alignview --resume mon_projet --html-out rapport_actualise.html --no-menu
alignview --resume mon_projet --msa align --tree --report --save --out nouveau_calcul --no-menu
```

## Règles à retenir

- Les chemins relatifs partent du dossier depuis lequel la commande est lancée ; le lanceur ne change pas de dossier.
- `--save --out DOSSIER` conserve un projet portable. `--out` seul ne déclenche pas la sauvegarde.
- `--html-out FICHIER.html` conserve seulement le HTML autonome. `--export-dir DOSSIER` conserve les exports du profil choisi.
- `--report` prépare le rapport en mémoire de session ; ajouter une destination ou une sauvegarde pour le conserver.
- Un fichier ou dossier de destination existant est refusé en batch ; choisir un nouveau nom.
- `--resume PROJET` ne relance aucun moteur par défaut. Seules les actions explicitement ajoutées (`--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif`) sont exécutées.
- En reprise, les paramètres explicites deviennent les réglages futurs. Les paramètres historiques des résultats déjà calculés restent conservés. `--report` et `--html-out` reconstruisent le HTML sans relancer de moteur.
- Un profil JSON fournit les valeurs de départ ; les options explicites sont prioritaires. Un profil de seuil (`--preset`) est distinct d’un profil de fichiers (`--export-profile`).

Les modes, les stratégies et leurs limites sont expliqués dans [les parcours](USAGE.md), [les paramètres](PARAMETERS.md) et [le guide des fichiers](FILES.md).

## Aide générale

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `-h / --help` | — | show this help message and exit |

## Entrées et reprise

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `query` | — | FASTA de la requête, ou FASTA unique à explorer |
| `--query-id` | '' | Identifiant FASTA à utiliser comme requête quand le premier fichier en contient plusieurs |
| `references` | — | FASTA des références |
| `--mode` | — | Programme de recherche ou alignement Valeurs : msa, vsearch, blastn, blastp, blastx, tblastn, tblastx, diamond-blastp, diamond-blastx, mmseqs-blastp, global-dna, global-aa. |
| `--resume` | — | Reprendre un projet sauvegardé sans relancer ses calculs |
| `--demo` | False | Explorer un exemple embarqué sans outil externe |

## Recherche et filtres

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--task` | '' | Tâche NCBI compatible avec le programme |
| `--evalue` | 1e-05 | E-value maximale (1e-5 par défaut) |
| `--identity` | 0 | Identité locale minimale % (0) |
| `--coverage` | 0 | Couverture QUERY par union des HSP % (0) |
| `--bitscore` | 0 | Bitscore local minimal (0) |
| `--hits` | 100 | Nombre maximal de références retenues (100 par défaut) |
| `--preset` | 'standard' | Profil de recherche : standard, sensible ou strict ; les options explicites sont prioritaires Valeurs : standard, sensible, strict. |
| `--threads` | 4 | Nombre de threads (4) |
| `--sensitivity` | 'sensitive' |  Valeurs : fast, sensitive, more-sensitive, very-sensitive, ultra-sensitive. |

## Alignements et analyses

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--msa` | 'none' |  Valeurs : none, view, align, full, hits, all, global, fragments, hsp. |
| `--mafft` | 'auto' |  Valeurs : auto, global, local, fast. |
| `--aligner` | 'mafft' | Moteur MSA : MAFFT, MUSCLE 5, Clustal Omega Valeurs : mafft, muscle, clustalo. |
| `--tree` | False | Inférer un arbre exploratoire FastTree à partir du MSA calculé |
| `--pairwise` | 'none' |  Valeurs : none, global, local. |
| `--ref` | 1 | Indice d’une référence (1) |
| `--frame` | — | Forcer un cadre de traduction Valeurs : -3, -2, -1, 1, 2, 3. |
| `--gapopen` | 10 |  |
| `--gapextend` | 0.5 |  |
| `--matrix` | '' | Matrice de substitution EMBOSS |
| `--free-ends` | False | Ne pas pénaliser les gaps terminaux Needle |
| `--dotplot` | False | Générer dotplot QUERY / meilleure référence |
| `--kmer` | 9 | Longueur des mots du dotplot (9) |
| `--motif` | '' | Motif IUPAC sur la requête |
| `--vsearch` | 'none' | Analyse nucléotidique VSEARCH : recherche, déduplication, regroupement ou chimères Valeurs : none, search, dereplicate, cluster, chimeras. |

## Rapports, exports et sauvegardes

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--report` | False | Préparer le HTML et ses SVG ; conserver avec --save, --export-dir ou --html-out |
| `--html-out` | '' | Enregistrer un HTML autonome à cet endroit, même sans --save ; refuse d’écraser |
| `--no-open` | False | Ne pas ouvrir automatiquement le HTML dans le navigateur |
| `--out` | '' | Dossier du projet avec --save ; chemin relatif au dossier de lancement, créé à la fin |
| `--save` | False | Sauvegarde explicite à la fin (batch compris) |
| `--export-profile` | 'standard' | standard / publication / complete : fichiers préparés en cas de sauvegarde Valeurs : standard, publication, complete. |
| `--export-dir` | '' | Exporter explicitement les fichiers tabulaires dans ce dossier, même sans --save |

## Moteurs scientifiques et installation

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--executable` | [] | Chemin exact vers un outil externe, répétable, ex. blastx=/opt/ncbi/bin/blastx |
| `--tools-dir` | [] | Ajouter des exécutables locaux au PATH d’AlignView (répétable) |
| `--module-policy` | 'auto' | auto : modules HPC si commande absente ; never : PATH et chemins seulement Valeurs : auto, never. |
| `--diagnose` | False | Diagnostic PATH/modules puis sortie |
| `--setup` | — | Préparer l’installation de moteurs via Conda/Mamba Valeurs : alv, blast, clustalo, diamond, emboss, fasttree, full, mafft, minimal, mmseqs, muscle, vsearch. |
| `--install-tools` | False | Exécuter le plan --setup dans un environnement dédié |
| `--tools-prefix` | — | Préfixe Conda dédié (ou ALIGNVIEW_TOOLS_PREFIX) |
| `--diagnose-json` | False | Diagnostic machine-readable JSON, sans chargement de modules |

## Interface et documentation

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--list-genetic-codes` | False | Afficher les tables de traduction NCBI disponibles, puis quitter |
| `--no-menu` | False | Sans menu interactif |
| `--simple-ui` | False | Menus numérotés classiques |
| `--no-color` | False | Sans couleurs ANSI |
| `--list-tools` | False | Lister les fonctionnalités |
| `--docs` | — | Consulter un guide hors ligne, sans ouvrir de données (défaut : index) Valeurs : index, demarrage, parcours, fichiers, installation, parametres, html, exports, cli, depannage, glossaire, vsearch, translation, developpement, publication, audit, getting-started, quickstart, usage, files, parameters, troubleshooting, glossary, development, releasing, genetics. |
| `--version` | — | show program's version number and exit |
| `--lang` | 'fr' | Langue de l’interface, du HTML et des guides (ou ALIGNVIEW_LANG) Valeurs : fr, en. |

## Autres paramètres scientifiques

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--genetic-code` | 1 | Table de traduction QUERY (NCBI) : Code utilisé pour les traductions internes et les recherches traduites de la requête. Valeurs : 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| `--reference-genetic-code` | — | Table de traduction des références : Vide : même table que QUERY ; transmis à TBLASTN/TBLASTX et aux traductions des références. Valeurs : 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33. |
| `--orf-starts` | 'atg' | Codons de départ des ORF : atg : ATG uniquement ; table : tous les initiateurs de la table NCBI choisie. Valeurs : atg, table. |
| `--vsearch-id` | 0.97 | Identité minimale VSEARCH : Fraction de 0 à 1, distincte du filtre d’identité locale BLAST. |
| `--vsearch-iddef` | 2 | Définition d’identité VSEARCH : 0 à 4 ; 2 : exclut les gaps terminaux. Lire la documentation VSEARCH. Valeurs : 0, 1, 2, 3, 4. |
| `--vsearch-strand` | 'both' | Brins VSEARCH : both : recherche et regroupement sur les deux brins ; UCHIME de novo travaille sur le brin plus. Valeurs : plus, both. |
| `--vsearch-query-cov` | 0.0 | Couverture minimale de la requête VSEARCH : Fraction de 0 à 1, calculée par VSEARCH. |
| `--vsearch-target-cov` | 0.0 | Couverture minimale des références VSEARCH : Fraction de 0 à 1, calculée par VSEARCH. |
| `--vsearch-maxaccepts` | 0 | Candidats acceptés avant arrêt : 0 : pas de plafond de candidats acceptés ; distinct du maximum de références rapportées. |
| `--vsearch-maxrejects` | 32 | Candidats rejetés avant arrêt : 0 : pas d’arrêt sur rejets ; mettre aussi maxaccepts et minwordmatches à 0 pour une recherche exhaustive. |
| `--vsearch-minwordmatches` | 10 | Mots communs minimum VSEARCH : 0 désactive ce préfiltre ; utile aux séquences courtes, avec un coût plus élevé. |
| `--vsearch-masking` | 'dust' | Masquage VSEARCH : dust : faible complexité ; none : aucun. Les entrées internes sont normalisées en majuscules. Valeurs : none, dust. |
| `--vsearch-n-mismatch / --no-vsearch-n-mismatch` | False | Compter les ambiguïtés comme différences : Évite qu’un N soit traité comme une identité avec toute base compatible. |
| `--vsearch-sizein / --no-vsearch-sizein` | False | Lire les abondances des en-têtes : Utiliser les annotations ;size=N; existantes. Sinon chaque entrée compte pour une séquence. |
| `--vsearch-min-length` | 1 | Longueur minimale VSEARCH : Les séquences plus courtes sont exclues et comptées dans le résumé. |
| `--vsearch-max-length` | 50000 | Longueur maximale VSEARCH : Les séquences plus longues sont exclues et comptées dans le résumé. |
| `--vsearch-abskew` | 2.0 | Rapport d’abondance UCHIME : Abondance minimale des parents relativement à la chimère candidate. |
| `--vsearch-minh` | 0.28 | Score minimal de chimère UCHIME : Plus grand : détection plus stricte, avec moins de candidats détectés. |
| `--search-candidates` | — | Limite de candidats du moteur : Auto = max(100, 10 × références retenues). Ce plafond peut changer les résultats ; distinct des 100 hits affichés. |
| `--hsp-overlap` | 0.5 | Recouvrement minimal de l’ancre HSP : Fraction de la région requête ancre recouverte pour retenir un segment dans le MSA exploratoire. |
| `--sframe` | — | Cadre de référence pour Needle/Water : Vide : cadre du meilleur HSP, sinon +1. Applicable à TBLASTN/TBLASTX. Valeurs : -3, -2, -1, 1, 2, 3. |

## Paramètres BLAST

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--blast-word-size` | — | Taille de graine BLAST : Petites graines : détection plus sensible, calcul plus coûteux ; valeurs autorisées selon la tâche. |
| `--blast-reward` | — | Score de match ADN : Gain pour une base identique ; dépend des coûts de mismatch et de gap. |
| `--blast-penalty` | — | Score de mismatch ADN : Score négatif pour une substitution ; change le score et les statistiques. |
| `--blast-matrix` | — | Matrice protéique BLAST : Ex. BLOSUM62, BLOSUM45, PAM30. Les coûts de gap doivent être compatibles. |
| `--blast-gapopen` | — | Ouverture de gap BLAST : Coût propre à BLAST, indépendant de Needle/Water. TBLASTX est non gappé. |
| `--blast-gapextend` | — | Extension de gap BLAST : Coût propre à BLAST. Seules certaines combinaisons matrice/coûts sont acceptées. |
| `--blast-masking` | 'native' | Masquage DUST / SEG : Masque la faible complexité. Le désactiver peut révéler des correspondances peu spécifiques. Valeurs : native, yes, no. |
| `--blast-soft-masking` | 'native' | Masquage souple BLAST : true : masque les graines ; false : masque aussi l’extension. Défaut dépendant du programme. Valeurs : native, true, false. |
| `--blast-strand` | 'both' | Brin de la requête ADN : both : deux brins ; plus/minus : limite la recherche à un brin. Valeurs : both, plus, minus. |
| `--blast-max-hsps` | — | HSP maximum par référence : Limite les segments rapportés, donc potentiellement la couverture calculée ensuite. |
| `--blast-comp-stats` | 'native' | Statistiques de composition BLAST : 0 : désactivées ; 1/2/3 : corrections NCBI selon le moteur. Change scores et E-values. Valeurs : native, 0, 1, 2, 3. |

## Paramètres DIAMOND

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--diamond-matrix` | — | Matrice DIAMOND : Matrice protéique ; les coûts de gap doivent être compatibles. |
| `--diamond-gapopen` | — | Ouverture de gap DIAMOND : Coût propre à DIAMOND, indépendant des autres moteurs. |
| `--diamond-gapextend` | — | Extension de gap DIAMOND : Pénalise l’allongement des gaps. |
| `--diamond-masking` | 'native' | Masquage DIAMOND : 0 : désactivé ; 1 : tantan ; seg : SEG (DIAMOND ≥ 2.0.12). Valeurs : native, 0, 1, seg. |
| `--diamond-comp-stats` | 'native' | Composition DIAMOND : 0 : désactivée ; 1 : correction compositionnelle ; 3/4 : ajustements de matrice, selon version. Valeurs : native, 0, 1, 3, 4. |
| `--diamond-max-hsps` | — | HSP maximum DIAMOND : 0 : tous ; plusieurs HSP peuvent améliorer la couverture et augmenter le coût. |
| `--diamond-min-orf` | 1 | ORF minimum DIAMOND BLASTX : AlignView utilisait déjà 1, contrairement au choix automatique de DIAMOND. 0 : automatique natif. |

## Paramètres MMseqs2

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--mmseqs-sensitivity` | — | Sensibilité MMseqs2 (-s) : De 1 (rapide) à 7.5 (sensible). Indépendante du réglage DIAMOND. |
| `--mmseqs-mask` | 'native' | Masquage MMseqs2 : 1 : faible complexité masquée ; 0 : désactivé. Valeurs : native, 0, 1. |
| `--mmseqs-comp-bias` | 'native' | Correction de composition MMseqs2 : 0/1 : désactive/active la correction du biais de composition. Valeurs : native, 0, 1. |

## Paramètres MAFFT

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--mafft-op` | — | Ouverture de gap MAFFT (--op) : Pénalité de gap en alignement groupe à groupe ; indépendante d’EMBOSS. |
| `--mafft-ep` | — | Offset MAFFT (--ep) : Paramètre d’offset de MAFFT ; ne se confond pas avec gapextend d’EMBOSS. |
| `--mafft-maxiterate` | — | Raffinements MAFFT : Global/local : 1000 par défaut AlignView. Pour régler ce nombre, choisir global, local ou fast plutôt que auto. |

## Paramètres MUSCLE 5

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--muscle-method` | 'align' | Algorithme MUSCLE 5 : align : PPP ; super5 : prévu pour des ensembles plus grands. MUSCLE 3 non pris en charge. Valeurs : align, super5. |
| `--muscle-perturb` | 0 | Graine de perturbation MUSCLE : 0 : HMM standard ; une autre graine explore la variabilité de l’alignement. |
| `--muscle-perm` | 'none' | Permutation de l’arbre guide : Change le parcours de l’arbre. Une seule permutation produit un seul MSA. Valeurs : none, abc, acb, bca. |
| `--muscle-refineiters` | — | Raffinements MUSCLE 5 : Nombre d’itérations de raffinement ; coût et résultat peuvent changer. |

## Paramètres Clustal Omega

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--clustalo-iterations` | 0 | Itérations Clustal Omega : Cycles combinés arbre guide / HMM. 0 : aucun raffinement itératif. |
| `--clustalo-full / --no-clustalo-full` | False | Distances complètes Clustal Omega : Matrice de distances complète plutôt que mBed ; coût quadratique. |
| `--clustalo-full-iter / --no-clustalo-full-iter` | False | Distances complètes aux itérations : Matrice complète pendant les raffinements ; requiert des itérations. |

## Paramètres Needle

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--needle-endopen` | — | Ouverture de gap terminal Needle : Vide : même coût que gapopen AlignView ; ignoré si extrémités libres. |
| `--needle-endextend` | — | Extension de gap terminal Needle : Vide : même coût que gapextend AlignView ; ignoré si extrémités libres. |

## Paramètres FastTree

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--fasttree-model` | 'auto' | Modèle d’évolution FastTree : auto : JC pour ADN, JTT pour protéines. GTR : ADN ; JTT/WAG/LG : protéines. Valeurs : auto, JC, GTR, JTT, WAG, LG. |
| `--fasttree-gamma / --no-fasttree-gamma` | False | Réévaluation Gamma20 : Après optimisation CAT, réévalue et redimensionne les branches sous Gamma20 ; ne refait pas la topologie. |
| `--fasttree-boot` | — | Rééchantillonnages du support local : Support SH-like, pas un bootstrap d’arbres. 0 : désactive les supports locaux. |
| `--fasttree-seed` | — | Graine FastTree : Fixe la graine aléatoire pour les calculs stochastiques. |

## Profils de calcul réutilisables

| Option | Défaut | Rôle et valeurs |
| --- | --- | --- |
| `--load-profile` | — | Charger un profil de calcul ; les options CLI explicites sont prioritaires |
| `--write-profile` | — | Écrire le profil de calcul puis quitter, sans lancer de moteur |
