"""English presentation catalog. Named placeholders preserve user data."""

_TEXT = '''
Accueil|Home
Projet|Project
Explorer|Explore
Analyser|Analyze
Données|Data
Enregistrer et partager|Save and share
Réglages et aide|Settings and help
Lire les séquences, les correspondances et les alignements|Browse sequences, matches and alignments
Aligner des séquences ou rechercher des similarités|Align sequences or search for similarities
Changer les fichiers, la requête ou le projet|Change files, the active query or the project
Conserver le projet, exporter ou créer un rapport|Save the project, export results or create a report
Logiciels, paramètres et explication des méthodes|Tools, parameters and explanations
Que souhaitez-vous faire ?|What would you like to do?
Que souhaitez-vous calculer ?|What would you like to calculate?
Démarrer|Start
Ouvrir un alignement|Open an alignment
Lire un FASTA déjà aligné, sans refaire les calculs|Read an aligned FASTA without recalculating it
Aligner des séquences|Align sequences
Choisir un FASTA de séquences à aligner|Choose a FASTA file containing sequences to align
Comparer à des références|Compare against references
Choisir une requête et un fichier de références|Choose a query file and a reference file
Reprendre un projet|Resume a project
Retrouver des données et des résultats enregistrés|Restore saved data and results
Essayer un exemple|Try an example
Découvrir la visionneuse sans installer de moteur|Try the alignment viewer without installing an engine
Préparer les logiciels|Set up tools
Installer ou retrouver les moteurs sur Linux / HPC|Install or locate engines on Linux / HPC
Lire la documentation|Read the manual
Premiers pas, fonctions, sauvegardes et dépannage hors ligne|Getting started, features, saving and offline troubleshooting
Choisir le résultat attendu. Les méthodes et paramètres restent accessibles ensuite.|Choose the result you need. Methods and parameters remain available afterwards.
Choisir la requête active|Choose the active query
Choisir la séquence active|Choose the active sequence
Changer la requête active|Change the active query
Requête sélectionnée : {id} · {count} caractères|Selected query: {id} · {count} characters
{count} séquences · analyses ciblées sur une requête|{count} sequences · analyses focus on one query
{count} caractères · {kind}|{count} characters · {kind}
ADN/ARN|DNA/RNA
protéine|protein
inconnu|unknown
Séquences|Sequences
Références|References
Requête|Query
Alignement|Alignment
Programme|Program
Seuils|Thresholds
Recherche|Search
Sauvegarde|Saving
Dossier de lancement|Launch directory
aucune — un FASTA|none — single FASTA
non généré|not generated
généré|generated
modifications à enregistrer|unsaved changes
projet enregistré|project saved
aucun calcul lancé|no calculation started
Recherche à relancer : réglages modifiés|Search needs rerunning: parameters changed
{count} correspondances|{count} matches
alignement prêt ({count} colonnes)|alignment ready ({count} columns)
alignement calculé avec les anciens réglages|alignment uses the previous parameters
Retour à l’accueil. Les données et résultats sont conservés.|Back to Home. Data and results are retained.
Action annulée. Le projet reste ouvert.|Action cancelled. The project remains open.
Action non terminée|Action not completed
Explorer les résultats|Explore results
Lire l’alignement|Read the alignment
Ouvrir la visionneuse et retrouver sa position|Open the viewer at its previous position
Aucun alignement disponible. Choisir Analyser → Créer un alignement.|No alignment available. Choose Analyze → Create an alignment.
Correspondances et couvertures|Matches and coverage
{count} références retenues par la recherche|{count} references retained by the search
Lancer une recherche dans Analyser avec un fichier de références.|Run a search in Analyze using a reference file.
Consulter les séquences|Browse sequences
{count} séquences chargées, identifiants et contenu|{count} sequences loaded, identifiers and content
Qualité des séquences|Sequence quality
Longueur, composition et caractères ambigus|Length, composition and ambiguous symbols
Statistiques de l’alignement|Alignment statistics
Conservation, distances et arbre exploratoire|Conservation, distances and exploratory tree
Créer ou ouvrir un alignement avant de calculer ses statistiques.|Create or open an alignment before calculating its statistics.
Consulter une séquence|Browse a sequence
Séquence : {id}|Sequence: {id}
Type|Type
Longueur|Length
Correspondances|Matches
Correspondances trouvées|Matches found
{count} références retenues · limite : {limit}|{count} retained references · limit: {limit}
Choisir une référence pour lire ses segments alignés. / recherche son identifiant.|Choose a reference to inspect its aligned segments. / searches its identifier.
Segments alignés : {id}|Aligned segments: {id}
Résultats de recherche|Search results
Choisir une correspondance|Choose a match
Couverture de la requête|Query coverage
Régions détectées ou sans correspondance|Detected and unmatched regions
Couverture des références|Reference coverage
Régions détectées sur chaque séquence de référence|Detected regions on each reference sequence
Identité et couverture|Identity and coverage
Comparer la similarité locale et la part de la requête couverte|Compare local similarity and the fraction of the query covered
Qualité des segments|Segment quality
Substitutions, gaps et cadres de lecture|Mismatches, gaps and reading frames
Aucune référence retenue après filtrage.|No references retained after filtering.
Analyser l’alignement|Analyze the alignment
Conservation et consensus|Conservation and consensus
Repérer les colonnes stables, variables et ambiguës|Identify conserved, variable and ambiguous columns
Distances entre séquences|Sequence distances
Distances p, en excluant les gaps de chaque paire|Pairwise p-distances, excluding gaps in each pair
Arbre exploratoire|Exploratory tree
Calculer un arbre approximatif avec FastTree|Infer an approximate tree using FastTree
Un arbre nécessite au moins trois séquences alignées.|A tree requires at least three aligned sequences.
Créer un alignement|Create an alignment
Comparer les positions de plusieurs séquences|Compare positions across several sequences
Rechercher des similarités|Search for similarities
Comparer la requête à des séquences de référence|Compare the query against reference sequences
Ajouter un fichier de références dans Données pour lancer une recherche.|Add a reference file in Data to run a search.
Comparer deux séquences|Compare two sequences
Alignement global ou local avec Needle / Water|Global or local alignment using Needle / Water
Ajouter des références dans Données pour comparer deux séquences.|Add references in Data to compare two sequences.
Autres analyses|Other analyses
Motifs, traduction, ORF et profils de composition|Motifs, translation, ORFs and composition profiles
Toutes les séquences|All sequences
Requête et références retenues|Query and retained references
Ajouter une séquence partielle|Add a partial sequence
Segments détectés seulement|Detected segments only
Lancer l’alignement|Run the alignment
Choisir le moteur|Choose the engine
Actuellement : {engine}|Currently: {engine}
Choisir les séquences|Choose sequences
Dans un fichier unique, toutes ses séquences sont alignées.|All sequences in a single FASTA are aligned.
Paramètres de calcul|Calculation parameters
Voir les commandes et la préparation|Inspect commands and input preparation
Inspecter les arguments avant de lancer le moteur|Inspect arguments before running the engine
Profils de calcul|Calculation profiles
Réutiliser ou partager ses paramètres|Reuse or share parameters
Charger ou enregistrer ses réglages pour les réutiliser|Load or save reusable settings
Les paramètres de calcul changent l’alignement. Le résultat sera ouvert dans la visionneuse.|Calculation parameters affect the alignment. The result opens in the viewer.
Séquences à aligner|Sequences to align
Références retenues|Retained references
La requête et les séquences ayant une correspondance|The query and matching reference sequences
Toutes les références|All references
Inclut aussi les références sans correspondance|Also include references without a match
Ajout de la requête à un alignement des références avec MAFFT|Add the query to a reference alignment using MAFFT
Uniquement les segments détectés|Detected segments only
Comparaison exploratoire des régions locales, sans les extrémités|Exploratory comparison of local regions, excluding flanks
Lancer une recherche avant de comparer ses segments.|Run a search before comparing its segments.
Choix polyvalent ; stratégie automatique par défaut|General-purpose choice; automatic strategy by default
Alternative pour ADN et protéines|Alternative for DNA and proteins
Alternative pour les alignements multiples|Alternative for multiple alignments
Moteur actuel : {engine}|Current engine: {engine}
Méthode de recherche|Search method
recommandation|recommendation
Lancer la recherche|Run the search
Choisir la méthode|Choose the method
Régler la recherche|Adjust the search
Paramètres du moteur|Engine parameters
Voir les commandes|Inspect commands
Indexation et recherche ; aucun moteur lancé par l’aperçu|Indexing and search; preview runs no analysis engine
Logiciel manquant|Missing tool
Installer dans l’environnement AlignView|Install in the AlignView environment
Utiliser Conda, Mamba ou Micromamba|Use Conda, Mamba or Micromamba
Voir les autres possibilités|See other options
Module HPC, chemin explicite et installation manuelle|HPC module, explicit path or manual installation
Lancer cette installation ?|Run this installation?
Logiciels et installation|Tools and installation
Deux séquences|Two sequences
Comment comparer les deux séquences ?|How should the two sequences be compared?
De bout en bout|End to end
Alignement global avec Needle ; inclut les extrémités|Global alignment using Needle, including sequence ends
Le meilleur segment local|Best local segment
Alignement local optimal avec Water|Optimal local alignment using Water
Lancer la comparaison|Run the comparison
Matrice, pénalités, extrémités et traduction|Matrix, penalties, sequence ends and translation
Charger ou enregistrer les paramètres|Load or save parameters
Calculer l’arbre|Infer the tree
Modèle, Gamma20, support local et graine|Model, Gamma20, local support and random seed
Voir la commande|Inspect the command
Inspecter les options FastTree avant le calcul|Inspect FastTree options before inference
Arbre approximatif ; le support SH-like ne constitue pas un bootstrap d’arbres.|Approximate tree; SH-like support is not a tree bootstrap.
Analyses de séquences|Sequence analyses
Rechercher un motif|Find a motif
Motif IUPAC et positions sur les deux brins|IUPAC motif and positions on both strands
Cette analyse de motifs IUPAC nécessite une séquence ADN/ARN.|IUPAC motif analysis requires a DNA/RNA sequence.
Traduction et régions codantes|Translation and coding regions
Traduire, rechercher des ORF ou extraire une région|Translate, find ORFs or extract a region
La traduction nécessite une séquence ADN/ARN.|Translation requires a DNA/RNA sequence.
Composition et profils|Composition and profiles
GC, entropie et fréquences de mots|GC, entropy and k-mer frequencies
Visualiser des correspondances de mots entre deux séquences|Visualize identical words shared by two sequences
Charger au moins deux séquences pour les comparer.|Load at least two sequences to compare them.
Traduction et régions|Translation and regions
Traduire les six cadres|Translate all six frames
Code génétique standard ; export des séquences protéiques|Selected NCBI code; export protein sequences
Rechercher des ORF simples|Find simple ORFs
Segments ATG → STOP ; analyse exploratoire|Start-to-stop segments; exploratory analysis
Extraire une région|Extract a region
Choisir ses coordonnées et son orientation|Choose coordinates and orientation
Choisir la table de traduction|Choose the translation table
Codons de départ des ORF|ORF start codons
Table de traduction QUERY (NCBI)|Query translation table (NCBI)
Table de traduction des références|Reference translation table
Code NCBI pour la traduction|NCBI translation code
Même table que QUERY|Use the query table
Hériter de la table de la requête|Inherit the query translation table
Composition|Composition
Qualité et composition globale|Overall quality and composition
Longueur, GC et caractères ambigus|Length, GC and ambiguous symbols
Profil le long de la séquence|Sequence profile
GC, complexité ou gaps selon les fenêtres|Sliding-window GC, complexity or gaps
Fréquences de mots|K-mer frequencies
Compter les k-mers exacts|Count exact k-mers
Données et projets|Data and projects
Ouvrir d’autres fichiers|Open other files
Choisir un alignement, des séquences ou des références|Choose an alignment, sequences or references
Revenir à une analyse précédente|Return to a previous analysis
Ouvrir un projet enregistré|Open a saved project
Retrouver ses fichiers, paramètres et calculs|Restore files, parameters and calculations
Voir les données chargées|Inspect loaded data
Fichiers, requête, paramètres et résultats disponibles|Files, query, parameters and available results
Le fichier chargé ne contient qu’une seule séquence.|The loaded file contains only one sequence.
Analyses précédentes|Previous analyses
Enregistrer|Save
Enregistrer le projet|Save the project
Tout conserver pour reprendre le travail sur un autre poste|Keep everything to resume work on another machine
Créer un rapport HTML|Create an HTML report
Un fichier à consulter et partager dans un navigateur|A file to view and share in a browser
Exporter les résultats|Export results
Choisir les tables, alignements et figures à conserver|Choose tables, alignments and figures to keep
Préparer un calcul sur HPC|Prepare an HPC run
Script Slurm inclus dans le projet à enregistrer ; sans soumission|Slurm script saved with the project; no automatic submission
Exporter|Export
Que souhaitez-vous exporter ?|What would you like to export?
Résultats essentiels|Essential results
FASTA, qualité, correspondances, couvertures et provenance|FASTA, quality, matches, coverage and provenance
Tables, alignements et figures|Tables, alignments and figures
Ajoute CSV, GFF3, NEXUS/PHYLIP, SVG et HTML|Adds CSV, GFF3, NEXUS/PHYLIP, SVG and HTML
Toutes les données disponibles|All available data
Ajoute JSON/JSONL, régions sans correspondance et SHA256|Adds JSON/JSONL, unmatched regions and SHA256
Dossier destination (q pour annuler)|Destination directory (q to cancel)
Dossier du projet (q pour annuler)|Project directory (q to cancel)
Dossier du projet (q = annuler)|Project directory (q to cancel)
Réglages|Settings
Réglages de recherche|Search settings
Réglages d’alignement|Alignment settings
Voir ce qui est disponible et installer les moteurs manquants|Inspect available tools and install missing engines
Profil de seuil E-value|E-value threshold preset
Nombre maximal de références|Maximum number of references
E-value maximale|Maximum E-value
Nombre de CPU|CPU count
Filtres des résultats|Result filters
Identité, couverture et bitscore ; sans refaire la recherche|Identity, coverage and bitscore without rerunning the search
Graines, masquage, score et sensibilité, selon la méthode|Seeds, masking, scores and sensitivity, depending on the method
Standard|Standard
Plus permissif|More permissive
Plus strict|Stricter
Comprendre les résultats|Understand results
Différence entre recherche, alignement, couverture et identité|Differences between search, alignment, coverage and identity
Calculs effectués et provenance|Completed calculations and provenance
Paramètres réellement utilisés, versions et commandes|Actual parameters, versions and commands
Charger ou enregistrer des réglages réutilisables|Load or save reusable settings
Documentation complète|Complete manual
Guides de toutes les fonctions, sorties, HTML, options et dépannage|Guides to features, outputs, HTML, options and troubleshooting
Documentation|Documentation
Choisir un guide|Choose a guide
Tous les guides, classés par besoin|All guides, organized by task
Premiers pas|Getting started
Un exemple et les premiers résultats, étape par étape|An example and your first results, step by step
Utiliser toutes les fonctions|Use every feature
Menus, calculs et exploration des résultats|Menus, analyses and result exploration
Dossiers et sauvegardes|Directories and saving
Application, sorties, temporaires et reprise|Application, outputs, temporary files and resuming
Installer et déplacer AlignView|Install and move AlignView
Python, moteurs scientifiques et Linux / HPC|Python, scientific engines and Linux / HPC
Comprendre les paramètres|Understand parameters
Effets, défauts et profils par moteur|Effects, defaults and profiles for each engine
Lire et partager le rapport HTML|Read and share the HTML report
Figures, arbre, filtres et téléchargements|Figures, tree, filters and downloads
Réutiliser les résultats|Reuse results
Formats, coordonnées et profils d’export|Formats, coordinates and export profiles
Toutes les options de commande|All command-line options
Référence générée à partir du parseur réel|Reference generated from the actual parser
Résoudre un problème|Troubleshooting
Démarrage, moteurs, résultats et récupération|Startup, engines, results and recovery
Comprendre les termes|Understand the terminology
Développer et tester|Develop and test
Architecture, tests et fabrication des distributions|Architecture, tests and distribution builds
Préparer une publication|Prepare a release
Archives, validations et limites de la release|Archives, validation and release limitations
Portée des validations|Scope of validation
Mesures, moteurs réels et limites connues|Measurements, actual engines and known limitations
Réglages actuels|Current settings
Valeur de départ|Initial value
Valeurs de départ|Initial values
Valeur invalide|Invalid value
Rétablir les valeurs de départ|Restore initial values
Réinitialiser uniquement les paramètres de ce calcul|Reset only this calculation's parameters
Le paramètre est conservé.|The parameter is unchanged.
Oui|Yes
Non|No
défaut|default
automatique|automatic
Retour|Back
Choix (q : retour)|Choice (q: back)
Choix inconnu ; saisir le numéro affiché.|Unknown choice; enter the displayed number.
Entrée pour continuer…|Press Enter to continue…
↑↓ choisir · Entrée ouvrir · ←/Échap retour · F2 accueil|↑↓ select · Enter open · ←/Esc back · F2 home
↑↓ choisir · Entrée ouvrir · q quitter|↑↓ select · Enter open · q quit
{position}/{count} · / rechercher · PgUp/PgDn pages · F1 détails|{position}/{count} · / search · PgUp/PgDn pages · F1 details
Lignes {start}–{end}/{count} · décalage horizontal {offset}|Lines {start}–{end}/{count} · horizontal offset {offset}
↑↓ / PgUp/PgDn défiler · ←→ latéral · / chercher · n suivant|↑↓ / PgUp/PgDn scroll · ←→ horizontal · / search · n next
Échap / q retour · F2 accueil|Esc / q back · F2 home
Recherche : {text}▏ · Entrée ouvre le résultat|Search: {text}▏ · Enter opens the result
Choix : {text} · Entrée pour ouvrir|Choice: {text} · Enter to open
Aucun résultat pour cette recherche.|No results for this search.
Aucun résultat. Échap pour effacer la recherche.|No results. Esc clears the search.
Aucun résultat.|No results.
indisponible|unavailable
Résultats|Results
Fichiers du projet|Project files
Outil|Tool
Exécutable|Executable
introuvable|not found
chemin explicite|explicit path
environnement AlignView|AlignView environment
Commandes préparées à partir des mêmes arguments que l’exécution.|Commands prepared from the same arguments used for execution.
Aucun moteur n’est lancé par cet aperçu.|This preview runs no analysis engine.
Les options absentes gardent les défauts de la version installée.|Omitted options keep the installed version's defaults.
Les modules HPC peuvent résoudre le chemin au moment du lancement.|HPC modules may resolve the executable path at run time.
Dossier de travail|Working directory
Journal stderr|Standard error log
Rapport construit · fichier HTML autonome · aucune connexion réseau nécessaire.|Report created · standalone HTML file · no network connection required.
HTML temporaire préparé ; sauvegarde avec --save / demande à la fermeture.|Temporary HTML prepared; retain it with --save, --html-out or saving on exit.
HTML indépendant créé : {path}|Standalone HTML saved: {path}
Résultats enregistrés : {path}|Results saved: {path}
Projet enregistré : {path}|Project saved: {path}
Projet repris : {id} · {count} HSP · aucun calcul relancé|Project resumed: {id} · {count} HSPs · no calculations rerun
Fin : résultats temporaires effacés (--save pour les conserver).|Finished: temporary results removed (use --save to retain them).
Ouvrir report.html dans le navigateur, sans connexion nécessaire.|Open report.html in a browser; no network connection is required.
Session interrompue.|Session interrupted.
Voir l’arbre|View the tree
Télécharger SVG|Download SVG
Afficher la figure|Show figure
Figures|Figures
Atlas des figures|Figure atlas
Filtrer les figures|Filter figures
Tous les domaines|All categories
Domaine|Category
Phylogénie|Phylogeny
Séquences & qualité|Sequences and quality
HSP et coordonnées|HSPs and coordinates
Qualité des HSP|HSP quality
Couvertures|Coverage
Traduction|Translation
Résumé|Overview
Vue générale|Overview
Par où commencer ?|Where should I start?
Ce rapport présente les données enregistrées. Choisir une vue ci-dessous ; les nouveaux calculs se lancent dans l’application AlignView.|This report presents saved results. Choose a view below; new analyses run in the AlignView application.
Paramètres|Parameters
Comprendre|Understand
Paramètres et aide|Parameters and help
Exporter les données|Export data
Lire l’alignement|Read the alignment
Examiner les correspondances|Inspect matches
Consulter les séquences|Browse sequences
Imprimer / PDF|Print / PDF
Arbre exploratoire · FastTree|Exploratory tree · FastTree
Arbre exploratoire FastTree|Exploratory FastTree tree
{count} séquences · Racine de dessin arbitraire, sans enracinement biologique|{count} sequences · Arbitrary drawing root; no biological rooting
Longueurs horizontales : substitutions par site.|Horizontal lengths: substitutions per site.
Cladogramme : longueurs absentes ou toutes nulles ; distances non représentées.|Cladogram: missing or all-zero branch lengths; distances are not represented.
Étiquettes internes : supports locaux FastTree (SH-like), pas des bootstraps.|Internal labels: local FastTree support (SH-like), not bootstrap values.
Interprétation exploratoire ; dépend de la qualité du MSA et du modèle choisi.|Exploratory interpretation; depends on MSA quality and the selected model.
Branches et étiquettes du fichier Newick. Racine de dessin arbitraire, sans enracinement biologique. Les supports FastTree sont locaux (SH-like), pas des bootstraps.|Branches and labels from the Newick file. Arbitrary drawing root; no biological rooting. FastTree support values are local (SH-like), not bootstraps.
Branche|Branch
non précisée|not specified
Analyses VSEARCH|VSEARCH analyses
Résultats VSEARCH|VSEARCH results
Recherche globale, déduplication, groupes et chimères|Global search, dereplication, clustering and chimeras
Recherche globale, groupes, séquences uniques et chimères|Global search, clusters, unique sequences and chimeras
Lancer une analyse VSEARCH dans Analyser.|Run a VSEARCH analysis in Analyze.
VSEARCH nécessite des séquences nucléotidiques.|VSEARCH requires nucleotide sequences.
Quelle analyse VSEARCH ?|Which VSEARCH analysis?
Recherche globale contre les références|Global search against references
Comparer la requête active aux références nucléotidiques|Compare the active query against nucleotide references
Dédupliquer les séquences|Dereplicate sequences
Regrouper les séquences identiques et conserver leurs abondances|Merge identical sequences and retain their abundances
Regrouper les séquences|Cluster sequences
Créer des groupes par identité et exporter leurs représentants|Create identity-based clusters and export their representatives
Détecter les chimères de novo|Detect chimeras de novo
Déduplication puis UCHIME ; les abondances sont importantes|Dereplication followed by UCHIME; abundances matter
Charger des références pour une recherche globale.|Load references for a global search.
Préparer VSEARCH : {action}|Prepare VSEARCH: {action}
Lancer l’analyse|Run the analysis
Paramètres VSEARCH|VSEARCH parameters
Identité, brins, couvertures, abondances et limites|Identity, strands, coverage, abundances and limits
Moteur, chemins, normalisation et paramètres réels|Engine, paths, normalization and actual parameters
Vérifier VSEARCH et sa version|Inspect VSEARCH and its version
La recherche utilise la requête active : {id}|Search uses the active query: {id}
Analyse VSEARCH terminée : {action}|VSEARCH analysis completed: {action}
Voir les analyses VSEARCH|View VSEARCH analyses
Outil et paramètres utilisés|Tool and parameters used
Analyses nucléotidiques. La recherche VSEARCH est globale : ni E-value ni bitscore BLAST ne sont calculés.|Nucleotide analyses. VSEARCH search uses global alignment; BLAST E-values and bitscores are not calculated.
Aperçus limités à 200 lignes par analyse ; fichiers complets dans le sous-dossier vsearch du projet et des exports.|Previews are limited to 200 rows per analysis; full files are in the vsearch subdirectory of the project and exports.
Identité selon iddef ; qcov/tcov VSEARCH en %, distincts d’une union de HSP locaux.|Identity follows iddef; VSEARCH qcov/tcov are percentages and differ from a union of local HSPs.
Les abondances et la composition du jeu influencent UCHIME ; une absence de chimère détectée n’est pas une validation absolue.|Abundances and dataset composition affect UCHIME; no detected chimera is not an absolute validation.
Aucun résultat retenu pour cet aperçu.|No results retained for this preview.
input_sequences|Input sequences
eligible_sequences|Sequences within length limits
excluded_by_length|Sequences excluded by length
matches|Matches
matched_references|Matching references
representatives|Representatives
total_abundance|Total abundance
unique|Unique sequences
chimeras|Chimeras
nonchimeras|Non-chimeras
borderline|Borderline
Identité minimale VSEARCH|Minimum VSEARCH identity
Définition d’identité VSEARCH|VSEARCH identity definition
Brins VSEARCH|VSEARCH strands
Couverture minimale de la requête VSEARCH|Minimum VSEARCH query coverage
Couverture minimale des références VSEARCH|Minimum VSEARCH reference coverage
Candidats acceptés avant arrêt|Accepted candidates before stopping
Candidats rejetés avant arrêt|Rejected candidates before stopping
Mots communs minimum VSEARCH|Minimum shared words for VSEARCH
Masquage VSEARCH|VSEARCH masking
Compter les ambiguïtés comme différences|Count ambiguous symbols as mismatches
Lire les abondances des en-têtes|Read header abundances
Longueur minimale VSEARCH|Minimum VSEARCH sequence length
Longueur maximale VSEARCH|Maximum VSEARCH sequence length
Rapport d’abondance UCHIME|UCHIME abundance skew
Score minimal de chimère UCHIME|Minimum UCHIME chimera score
Commandes et versions|Commands and versions
Enregistrer les modifications|Save changes
Quitter sans enregistrer|Quit without saving
Continuer le travail|Continue working
Enregistrement|Saving
Avant de quitter|Before quitting
Conserver les fichiers de cette session ?|Keep this session's files?
Le HTML existe déjà : {path}. Choisir un autre nom.|HTML already exists: {path}. Choose another name.
Extension .html requise.|The .html extension is required.
Fichier introuvable : {path}|File not found: {path}
Identifiant QUERY introuvable : {id}|Query identifier not found: {id}
FASTA QUERY multiple : préciser --query-id en mode non interactif.|Multiple query sequences: provide --query-id in non-interactive mode.
Table NCBI inconnue : {id}|Unknown NCBI translation table: {id}
Table {id} : arrêts contextuels ; la détection simple d’ORF nécessite un contexte CDS absent ici.|Table {id}: context-dependent stops; simple ORF detection requires CDS context, which is unavailable here.
VSEARCH nécessite des séquences nucléotidiques ; les protéines ne sont pas prises en charge.|VSEARCH requires nucleotide sequences; proteins are not supported.
VSEARCH : aucune séquence dans les limites de longueur choisies.|VSEARCH: no sequences meet the selected length limits.
VSEARCH recherche : charger un FASTA de références.|VSEARCH search: load a reference FASTA.
Créer un MSA avant --tree : ajouter --msa align, ou utiliser les menus.|Create an MSA before --tree: add --msa align, or use the menus.
Cette stratégie MSA nécessite un FASTA de références.|This MSA strategy requires a reference FASTA.
Documentation locale · / rechercher un guide · Échap retour|Offline manual · / find a guide · Esc back
ALIGNEMENTS / EXPLORATION GUIDÉE|ALIGNMENTS / GUIDED EXPLORATION
SESSION / paramètres et données en mémoire|SESSION / data and parameters in memory
ÉTAT DE LA SESSION|SESSION STATUS
RAPPORT / EXPLORATION DANS LE NAVIGATEUR|REPORT / BROWSER EXPLORATION
TRADUCTION SIX CADRES / TABLE NCBI {code}|SIX-FRAME TRANSLATION / NCBI TABLE {code}
{count} séquences · gaps retirés et U→T · {eligible} dans les limites de longueur.|{count} sequences · gaps removed and U→T · {eligible} meet length limits.
VSEARCH : recherche globale|VSEARCH: global search
VSEARCH : regroupement|VSEARCH: clustering
VSEARCH : déduplication|VSEARCH: dereplication
VSEARCH : chimères de novo|VSEARCH: de novo chimera detection
FastTree : arbre approximatif|FastTree: approximate tree
Arbre exploratoire créé : tree_fasttree.nwk (Newick)|Exploratory tree created: tree_fasttree.nwk (Newick)
FastTree fournit une phylogénie approximative ; vérifier modèle, support et qualité du MSA pour publication.|FastTree provides an approximate phylogeny; check the model, support and MSA quality before publication.
Composition de la séquence QUERY|Query sequence composition
Distribution des longueurs|Length distribution
Distribution des longueurs des références|Reference length distribution
Longueur nucléotidique / GC (%)|Nucleotide length / GC (%)
Identité locale / étendue HSP sur QUERY|Local identity / HSP span on the query
Identités des HSP (%)|HSP identities (%)
Longueurs alignées (HSP)|Aligned HSP lengths
Cadres de lecture QUERY (HSP)|Query reading frames (HSP)
Conservation MSA (%)|MSA conservation (%)
Entropie de Shannon / colonne|Shannon entropy per column
Colonnes avec gaps (%)|Columns with gaps (%)
Distances p entre séquences (MSA)|Pairwise p-distances (MSA)
Distribution des taux de gaps par colonne|Distribution of per-column gap rates
Occurrences de motifs par séquence|Motif occurrences per sequence
Longueurs des ORF détectées|Detected ORF lengths
HSP par référence|HSPs per reference
Couverture QUERY, union des HSP (%)|Query coverage, union of HSPs (%)
Couverture dans les références (%)|Reference coverage (%)
Bitscore des HSP|HSP bitscores
Significativité −log₁₀(E-value)|Significance −log₁₀(E-value)
% colonnes avec gap par HSP|Columns with gaps per HSP (%)
Identité locale × bitscore|Local identity × bitscore
Gaps par séquence du MSA (%)|Gaps per MSA sequence (%)
Distribution de l’entropie des colonnes du MSA|Distribution of MSA column entropy
Occupation du MSA|MSA occupancy
Occupancy du MSA|MSA occupancy
Variabilité en fonction des gaps|Variability versus gaps
Composition du consensus MSA|MSA consensus composition
Distribution de la conservation du MSA|Distribution of MSA conservation
Séquences sans gap (%)|Sequences without gaps (%)
Entropie moyenne (bits)|Mean entropy (bits)
Conservation (%)|Conservation (%)
Longueur (nt)|Length (nt)
Longueur (aa)|Length (aa)
Colonnes alignées (aa)|Aligned columns (aa)
Colonnes alignées (nt)|Aligned columns (nt)
Étendue QUERY (%)|Query span (%)
Identité HSP (%)|HSP identity (%)
Acides aminés|Amino acids
Substitutions (colonnes aa)|Mismatches (aa columns)
Substitutions (colonnes nt)|Mismatches (nt columns)
Colonnes avec gap|Columns with gaps
Couverture QUERY (%)|Query coverage (%)
Couverture RÉFÉRENCE (%)|Reference coverage (%)
'''

EN=dict(line.split('|',1) for line in _TEXT.strip().splitlines() if line and '|' in line)

_MORE = '''
Entrées et reprise|Input and resuming
Recherche et filtres|Search and filtering
Alignements et analyses|Alignments and analyses
Rapports, exports et sauvegardes|Reports, exports and saving
Moteurs scientifiques et installation|Scientific engines and installation
Interface et documentation|Interface and documentation
Autres paramètres scientifiques|Other scientific parameters
Profils de calcul réutilisables|Reusable calculation profiles
Paramètres BLAST|BLAST parameters
Paramètres DIAMOND|DIAMOND parameters
Paramètres MMseqs2|MMseqs2 parameters
Paramètres MAFFT|MAFFT parameters
Paramètres MUSCLE 5|MUSCLE 5 parameters
Paramètres Clustal Omega|Clustal Omega parameters
Paramètres FastTree|FastTree parameters
Paramètres Needle|Needle parameters
AlignView — exploration, alignements, visualisation et diagnostic en terminal.|AlignView — terminal sequence exploration, alignment, visualization and diagnostics.
Exemples|Examples
FASTA de la requête, ou FASTA unique à explorer|Query FASTA, or a single FASTA to explore
Identifiant FASTA à utiliser comme requête quand le premier fichier en contient plusieurs|FASTA identifier to use as the query when the first file contains several sequences
FASTA des références|Reference FASTA
Programme de recherche ou alignement|Search or alignment method
Reprendre un projet sauvegardé sans relancer ses calculs|Resume a saved project without rerunning its calculations
Explorer un exemple embarqué sans outil externe|Explore a bundled example without external tools
Tâche NCBI compatible avec le programme|NCBI task compatible with the selected program
E-value maximale (1e-5 par défaut)|Maximum E-value (default: 1e-5)
Identité locale minimale %% (0)|Minimum local identity %% (0)
Couverture QUERY par union des HSP %% (0)|Query coverage from the union of HSPs %% (0)
Bitscore local minimal (0)|Minimum local bitscore (0)
Nombre maximal de références retenues (100 par défaut)|Maximum retained references (default: 100)
Profil de recherche : standard, sensible ou strict ; les options explicites sont prioritaires|Search preset: standard, sensible or strict; explicit options take priority
Nombre de threads (4)|Number of threads (4)
Moteur MSA : MAFFT, MUSCLE 5, Clustal Omega|MSA engine: MAFFT, MUSCLE 5, Clustal Omega
Inférer un arbre exploratoire FastTree à partir du MSA calculé|Infer an exploratory FastTree tree from the MSA
Indice d’une référence (1)|Reference index (1)
Forcer un cadre de traduction|Select a translation frame
Matrice de substitution EMBOSS|EMBOSS substitution matrix
Ne pas pénaliser les gaps terminaux Needle|Do not penalize terminal gaps in Needle
Générer dotplot QUERY / meilleure référence|Generate a query / best-reference dotplot
Longueur des mots du dotplot (9)|Dotplot word length (9)
Motif IUPAC sur la requête|IUPAC motif in the query
Analyse nucléotidique VSEARCH : recherche, déduplication, regroupement ou chimères|VSEARCH nucleotide analysis: search, dereplication, clustering or chimeras
Préparer le HTML et ses SVG ; conserver avec --save, --export-dir ou --html-out|Prepare HTML and SVG figures; retain them using --save, --export-dir or --html-out
Enregistrer un HTML autonome à cet endroit, même sans --save ; refuse d’écraser|Save a standalone HTML here, even without --save; refuses to overwrite
Ne pas ouvrir automatiquement le HTML dans le navigateur|Do not automatically open HTML in the browser
Dossier du projet avec --save ; chemin relatif au dossier de lancement, créé à la fin|Project directory with --save; relative to the launch directory, created at the end
Sauvegarde explicite à la fin (batch compris)|Explicitly save at the end, including batch runs
standard / publication / complete : fichiers préparés en cas de sauvegarde|standard / publication / complete: files prepared when saving
Exporter explicitement les fichiers tabulaires dans ce dossier, même sans --save|Explicitly export files to this directory, even without --save
Chemin exact vers un outil externe, répétable, ex. blastx=/opt/ncbi/bin/blastx|Exact external-tool path, repeatable, e.g. blastx=/opt/ncbi/bin/blastx
Ajouter des exécutables locaux au PATH d’AlignView (répétable)|Add local executable directories to AlignView's PATH (repeatable)
auto : modules HPC si commande absente ; never : PATH et chemins seulement|auto: try HPC modules for missing commands; never: use existing paths and the managed environment
Diagnostic PATH/modules puis sortie|Report available executables and versions, then exit
Préparer l’installation de moteurs via Conda/Mamba|Prepare engine installation using Conda/Mamba
Exécuter le plan --setup dans un environnement dédié|Run the --setup plan in a dedicated environment
Préfixe Conda dédié (ou ALIGNVIEW_TOOLS_PREFIX)|Dedicated Conda prefix (or ALIGNVIEW_TOOLS_PREFIX)
Diagnostic machine-readable JSON, sans chargement de modules|Machine-readable JSON diagnostics without loading modules
Afficher les tables de traduction NCBI disponibles, puis quitter|List available NCBI translation tables and exit
Langue de l’interface, du HTML et des guides (ou ALIGNVIEW_LANG)|Interface, HTML and manual language (or ALIGNVIEW_LANG)
Consulter un guide hors ligne, sans ouvrir de données (défaut : index)|Read an offline guide without opening data (default: index)
Lister les fonctionnalités|List features
Sans menu interactif|Without interactive menus
Menus numérotés classiques|Classic numbered menus
Sans couleurs ANSI|Without ANSI colors
Code utilisé pour les traductions internes et les recherches traduites de la requête.|Code used for internal query translation and translated query searches.
Vide : même table que QUERY ; transmis à TBLASTN/TBLASTX et aux traductions des références.|Unset: use the query table; passed to TBLASTN/TBLASTX and internal reference translation.
atg : ATG uniquement ; table : tous les initiateurs de la table NCBI choisie.|atg: ATG only; table: all initiators from the selected NCBI table.
Fraction de 0 à 1, distincte du filtre d’identité locale BLAST.|Fraction from 0 to 1; distinct from the BLAST local identity filter.
0 à 4 ; 2 : exclut les gaps terminaux. Lire la documentation VSEARCH.|0 to 4; 2 excludes terminal gaps. See the VSEARCH documentation.
both : recherche et regroupement sur les deux brins ; UCHIME de novo travaille sur le brin plus.|both: search and cluster both strands; de novo UCHIME uses the plus strand.
Fraction de 0 à 1, calculée par VSEARCH.|Fraction from 0 to 1, calculated by VSEARCH.
0 : pas de plafond de candidats acceptés ; distinct du maximum de références rapportées.|0: no accepted-candidate limit; distinct from the reported reference limit.
0 : pas d’arrêt sur rejets ; mettre aussi maxaccepts et minwordmatches à 0 pour une recherche exhaustive.|0: do not stop after rejections; also set maxaccepts and minwordmatches to 0 for exhaustive search.
0 désactive ce préfiltre ; utile aux séquences courtes, avec un coût plus élevé.|0 disables this prefilter; useful for short sequences, at greater computational cost.
dust : faible complexité ; none : aucun ; soft : masquage par minuscules.|dust: low-complexity masking; none: no masking; soft: lowercase masking.
Évite qu’un N soit traité comme une identité avec toute base compatible.|Prevent N from counting as a match to every compatible nucleotide.
Utiliser les annotations ;size=N; existantes. Sinon chaque entrée compte pour une séquence.|Use existing ;size=N; annotations. Otherwise each input record counts as one sequence.
Les séquences plus courtes sont exclues et comptées dans le résumé.|Shorter sequences are excluded and counted in the summary.
Les séquences plus longues sont exclues et comptées dans le résumé.|Longer sequences are excluded and counted in the summary.
Abondance minimale des parents relativement à la chimère candidate.|Minimum parent abundance relative to the candidate chimera.
Plus grand : détection plus stricte, avec moins de candidats détectés.|Higher values make detection stricter, with fewer candidates detected.
Ressources des moteurs qui acceptent des threads ; FastTree standard reste monothread.|Resources for engines that support threads; standard FastTree remains single-threaded.
Seuil transmis au moteur. Le modifier nécessite une nouvelle recherche.|Threshold sent to the engine. Changing it requires a new search.
Auto = max(100, 10 × références retenues). Ce plafond peut changer les résultats ; distinct des 100 hits affichés.|Auto = max(100, 10 × retained references). This candidate limit can affect results and differs from the 100 displayed hits.
Change graines et score. Sans option : BLASTN utilise megablast ; autres programmes : tâche native.|Changes seeds and scores. Unset: BLASTN uses megablast; other programs use their native task.
Recherche plus sensible = davantage de calcul ; ce réglage est propre à DIAMOND.|More sensitive search requires more work; this setting is specific to DIAMOND.
Petites graines : détection plus sensible, calcul plus coûteux ; valeurs autorisées selon la tâche.|Smaller seeds increase sensitivity and cost; allowed values depend on the task.
Gain pour une base identique ; dépend des coûts de mismatch et de gap.|Reward for an identical nucleotide; depends on mismatch and gap costs.
Score négatif pour une substitution ; change le score et les statistiques.|Negative mismatch score; affects scoring and statistics.
Ex. BLOSUM62, BLOSUM45, PAM30. Les coûts de gap doivent être compatibles.|Examples: BLOSUM62, BLOSUM45, PAM30. Gap costs must be compatible.
Coût propre à BLAST, indépendant de Needle/Water. TBLASTX est non gappé.|BLAST-specific cost, independent of Needle/Water. TBLASTX is ungapped.
Coût propre à BLAST. Seules certaines combinaisons matrice/coûts sont acceptées.|BLAST-specific cost. Only supported matrix/cost combinations are accepted.
Masque la faible complexité. Le désactiver peut révéler des correspondances peu spécifiques.|Mask low complexity. Disabling masking may reveal nonspecific matches.
true : masque les graines ; false : masque aussi l’extension. Défaut dépendant du programme.|true: mask seeds; false: also mask extension. The default depends on the program.
both : deux brins ; plus/minus : limite la recherche à un brin.|both: both strands; plus/minus: restrict the search to one strand.
Limite les segments rapportés, donc potentiellement la couverture calculée ensuite.|Limit reported segments and potentially the coverage calculated from them.
0 : désactivées ; 1/2/3 : corrections NCBI selon le moteur. Change scores et E-values.|0: disabled; 1/2/3: NCBI corrections depending on the engine. Affects scores and E-values.
Matrice protéique ; les coûts de gap doivent être compatibles.|Protein matrix; gap costs must be compatible.
Coût propre à DIAMOND, indépendant des autres moteurs.|DIAMOND-specific cost, independent of other engines.
Pénalise l’allongement des gaps.|Penalize gap extension.
0 : désactivé ; 1 : tantan ; seg : SEG (DIAMOND ≥ 2.0.12).|0: disabled; 1: tantan; seg: SEG (DIAMOND ≥ 2.0.12).
0 : désactivée ; 1 : correction compositionnelle ; 3/4 : ajustements de matrice, selon version.|0: disabled; 1: composition correction; 3/4: matrix adjustments, depending on the version.
0 : tous ; plusieurs HSP peuvent améliorer la couverture et augmenter le coût.|0: all; multiple HSPs may increase coverage and computational cost.
AlignView utilisait déjà 1, contrairement au choix automatique de DIAMOND. 0 : automatique natif.|AlignView uses 1 rather than DIAMOND's automatic choice. 0: native automatic behavior.
De 1 (rapide) à 7.5 (sensible). Indépendante du réglage DIAMOND.|From 1 (fast) to 7.5 (sensitive). Independent of DIAMOND sensitivity.
1 : faible complexité masquée ; 0 : désactivé.|1: low complexity masked; 0: disabled.
0/1 : désactive/active la correction du biais de composition.|0/1: disable/enable composition-bias correction.
Les algorithmes et leurs paramètres ne sont pas interchangeables.|Algorithms and their parameters are not interchangeable.
auto : choix natif ; global : G-INS-i ; local : L-INS-i ; fast : FFT-NS-2.|auto: native choice; global: G-INS-i; local: L-INS-i; fast: FFT-NS-2.
Pénalité de gap en alignement groupe à groupe ; indépendante d’EMBOSS.|Group-to-group alignment gap penalty; independent of EMBOSS.
Paramètre d’offset de MAFFT ; ne se confond pas avec gapextend d’EMBOSS.|MAFFT offset parameter; different from EMBOSS gapextend.
Global/local : 1000 par défaut AlignView. Pour régler ce nombre, choisir global, local ou fast plutôt que auto.|Global/local: AlignView defaults to 1000 iterations. Choose global, local or fast to set this number explicitly.
align : PPP ; super5 : prévu pour des ensembles plus grands. MUSCLE 3 non pris en charge.|align: PPP; super5: intended for larger datasets. MUSCLE 3 is unsupported.
0 : HMM standard ; une autre graine explore la variabilité de l’alignement.|0: standard HMM; another seed explores alignment variability.
Change le parcours de l’arbre. Une seule permutation produit un seul MSA.|Change the guide-tree traversal. One permutation produces one MSA.
Nombre d’itérations de raffinement ; coût et résultat peuvent changer.|Refinement iteration count; cost and results may change.
Cycles combinés arbre guide / HMM. 0 : aucun raffinement itératif.|Combined guide-tree/HMM cycles. 0: no iterative refinement.
Matrice de distances complète plutôt que mBed ; coût quadratique.|Full distance matrix instead of mBed; quadratic cost.
Matrice complète pendant les raffinements ; requiert des itérations.|Use a full matrix during refinement; requires iterations.
Fraction de la région requête ancre recouverte pour retenir un segment dans le MSA exploratoire.|Fraction of the anchor query region that a segment must cover in the exploratory MSA.
Vide : cadre du meilleur HSP, sinon +1. Ne réduit pas les six cadres recherchés par BLASTX.|Unset: best HSP frame, otherwise +1. Does not restrict the six frames searched by BLASTX.
Vide : cadre du meilleur HSP, sinon +1. Applicable à TBLASTN/TBLASTX.|Unset: best HSP frame, otherwise +1. Applies to TBLASTN/TBLASTX.
Pénalité Needle/Water uniquement ; dépend de la matrice utilisée.|Needle/Water penalty only; depends on the selected matrix.
Pénalité Needle/Water uniquement.|Needle/Water penalty only.
Vide : EDNAFULL en ADN, EBLOSUM62 en protéines (défauts EMBOSS).|Unset: EDNAFULL for DNA, EBLOSUM62 for proteins (EMBOSS defaults).
AlignView pénalise les extrémités par défaut ; EMBOSS natif les laisse libres. Sans effet sur Water.|AlignView penalizes sequence ends by default; native EMBOSS leaves them free. No effect on Water.
Vide : même coût que gapopen AlignView ; ignoré si extrémités libres.|Unset: same as AlignView gapopen; ignored when sequence ends are free.
Vide : même coût que gapextend AlignView ; ignoré si extrémités libres.|Unset: same as AlignView gapextend; ignored when sequence ends are free.
auto : JC pour ADN, JTT pour protéines. GTR : ADN ; JTT/WAG/LG : protéines.|auto: JC for DNA, JTT for proteins. GTR: DNA; JTT/WAG/LG: proteins.
Après optimisation CAT, réévalue et redimensionne les branches sous Gamma20 ; ne refait pas la topologie.|After CAT optimization, rescore and rescale branches under Gamma20; does not repeat topology inference.
Support SH-like, pas un bootstrap d’arbres. 0 : désactive les supports locaux.|SH-like support, not a tree bootstrap. 0: disable local support.
Fixe la graine aléatoire pour les calculs stochastiques.|Set the random seed for stochastic calculations.
Taille de graine BLAST|BLAST seed size
Gain identité BLASTN|BLASTN match reward
Pénalité substitution BLASTN|BLASTN mismatch penalty
Matrice BLAST protéines|BLAST protein matrix
Ouverture de gap BLAST|BLAST gap opening
Extension de gap BLAST|BLAST gap extension
Masquage BLAST|BLAST masking
Soft masking BLAST|BLAST soft masking
Brin BLAST|BLAST strand
HSP maximum par référence|Maximum HSPs per reference
Statistiques de composition BLAST|BLAST composition statistics
Matrice DIAMOND|DIAMOND matrix
Ouverture de gap DIAMOND|DIAMOND gap opening
Extension de gap DIAMOND|DIAMOND gap extension
Masquage DIAMOND|DIAMOND masking
Composition DIAMOND|DIAMOND composition correction
HSP maximum DIAMOND|Maximum DIAMOND HSPs
ORF minimum DIAMOND BLASTX|Minimum DIAMOND BLASTX ORF
Sensibilité MMseqs2 (-s)|MMseqs2 sensitivity (-s)
Masquage MMseqs2|MMseqs2 masking
Correction de composition MMseqs2|MMseqs2 composition correction
Moteur MSA|MSA engine
Stratégie MAFFT|MAFFT strategy
Ouverture de gap MAFFT (--op)|MAFFT gap opening (--op)
Offset MAFFT (--ep)|MAFFT offset (--ep)
Raffinements MAFFT|MAFFT refinement iterations
Algorithme MUSCLE 5|MUSCLE 5 algorithm
Graine de perturbation MUSCLE|MUSCLE perturbation seed
Permutation de l’arbre guide|Guide-tree permutation
Raffinements MUSCLE 5|MUSCLE 5 refinement iterations
Itérations Clustal Omega|Clustal Omega iterations
Distances complètes Clustal Omega|Full Clustal Omega distances
Distances complètes aux itérations|Full distances during iterations
Recouvrement minimal de l’ancre HSP|Minimum HSP anchor overlap
Cadre QUERY pour les analyses traduites|Query frame for translated analyses
Cadre de référence pour Needle/Water|Reference frame for Needle/Water
Ouverture de gap EMBOSS|EMBOSS gap opening
Extension de gap EMBOSS|EMBOSS gap extension
Matrice EMBOSS|EMBOSS matrix
Extrémités libres Needle|Free Needle ends
Ouverture de gap terminal Needle|Needle terminal gap opening
Extension de gap terminal Needle|Needle terminal gap extension
Modèle d’évolution FastTree|FastTree evolutionary model
Réévaluation Gamma20|Gamma20 reassessment
Rééchantillonnages du support local|Local-support resampling
Graine FastTree|FastTree random seed
Plafond de candidats du moteur|Engine candidate limit
Langue / Language|Language / Langue
Français / English · menus, HTML et documentation|French / English · menus, HTML and documentation
Séquences en entrée|Input sequences
Dans les limites de longueur|Within length limits
Exclues par longueur|Excluded by length
Références correspondantes|Matching references
Représentants|Representatives
Abondance totale|Total abundance
Séquences uniques|Unique sequences
Chimères|Chimeras
Non chimériques|Non-chimeras
Indéterminées|Borderline
←→ / h l : colonne · Maj+←→ : 10 colonnes|←→ / h l: column · Shift+←→: 10 columns
Ctrl+←→ ou [ ] : une largeur d’écran|Ctrl+←→ or [ ]: one screen width
↑↓ / k j : séquence · PgUp/PgDn : page de séquences|↑↓ / k j: sequence · PgUp/PgDn: sequence page
Home/End : première/dernière colonne|Home/End: first/last column
Ctrl+Home/End : première/dernière séquence|Ctrl+Home/End: first/last sequence
Molette : séquences · Maj+molette : colonnes|Wheel: sequences · Shift+wheel: columns
Clic dans la vue globale : atteindre une région|Click the overview to jump to a region
g : colonne absolue (1-based) ou pourcentage, ex. 75%|g: absolute column (1-based) or percentage, e.g. 75%
s : atteindre une séquence par son identifiant|s: go to a sequence by identifier
/ : motif exact, c:motif pour le consensus · n / N : suivant/précédent|/: exact motif, c:motif for consensus · n / N: next/previous
v : sites variables · d : différences au consensus|v: variable sites · d: differences from consensus
w : rapport HTML · e : exporter le HTML|w: HTML report · e: export HTML
q/Échap : retour · F2 : accueil|q/Esc: back · F2: home
Coordonnées = colonnes du MSA, gaps compris.|Coordinates are MSA columns, including gaps.
La vue variable masque des colonnes sans modifier le MSA.|The variable-site view hides columns without modifying the MSA.
Appuyer sur une touche pour revenir.|Press any key to return.
AlignView / Aide de l’alignement|AlignView / Alignment help
'''
EN.update(line.split('|',1) for line in _MORE.strip().splitlines() if line and '|' in line)

_STATUS = '''
SESSION / paramètres et données en mémoire|SESSION / parameters and data in memory
{file} · {count} séq. · {length} positions|{file} · {count} sequences · {length} positions
{count} HSP · {refs} références retenues|{count} HSPs · {refs} retained references
{count} séquences × {columns} colonnes · colonnes {start}–{end}|{count} sequences × {columns} columns · columns {start}–{end}
Modules chargés : {modules}|Loaded modules: {modules}
aucun, PATH utilisé|none; using PATH
Dernière copie enregistrée : {path}|Last saved copy: {path}
Chargé|Loaded
Taille|Length
VARIABLES|VARIABLE SITES
Module activé : {module} ({command})|Module loaded: {module} ({command})
Commande absente : {command}|Missing command: {command}
Charger un profil de calcul ; les options CLI explicites sont prioritaires|Load a calculation profile; explicit CLI options take priority
Écrire le profil de calcul puis quitter, sans lancer de moteur|Write the calculation profile and exit without running an engine
Profil de calcul enregistré : {path}|Calculation profile saved: {path}
Exemples :|Examples:
sans correspondance|no match
rechercher|search
annuler|cancel
détails|details
retour|back
Séquence|Sequence
Entropie|Entropy
Observations|Observations
aucune|none
Régions couvertes par la recherche|Regions covered by the search
Données et calculs du projet|Project data and calculations
Qualité des segments alignés|Aligned segment quality
Logiciels disponibles|Available tools
Conservation de l’alignement|Alignment conservation
Distances entre les séquences|Sequence distances
Enregistrer les modifications avant de quitter ?|Save changes before quitting?
Les modifications de cette analyse n’ont pas été enregistrées.|Changes in this analysis have not been saved.
Conserver les modifications|Keep changes
Enregistrer dans un nouveau dossier|Save in a new directory
Abandonner les modifications|Discard changes
Les copies déjà enregistrées restent conservées|Previously saved copies are retained
Revenir au projet|Return to the project
Continuer l’analyse|Continue the analysis
Chemin du FASTA|FASTA path
Saisir un chemin|Enter a path
Coller un chemin complet ou relatif|Paste an absolute or relative path
Dossier parent  ..|Parent directory  ..
Alignement déjà calculé|Existing alignment
Séquences à aligner|Sequences to align
Fichier de requêtes|Query file
Séquences de référence|Reference sequences
Logiciels / installation|Tools / installation
État et versions|Availability and versions
Exécutables, origine et commande de réparation|Executables, source and setup command
Installer le socle conseillé|Install the recommended core
BLAST + DIAMOND + MAFFT dans un environnement dédié|BLAST + DIAMOND + MAFFT in a dedicated environment
Installer un outil précis|Install one tool
Ajouter uniquement ce qui manque|Add only what is missing
Installer tous les moteurs|Install all engines
Rechercher les modules HPC|Look for HPC modules
Charger à la demande selon la politique de la session|Load on demand according to the session policy
Reprendre l’analyse|Resume the analysis
Choisir un outil|Choose a tool
Outils scientifiques / diagnostic|Scientific tools / diagnostics
Installation en cours ; les messages du gestionnaire restent visibles.|Installation in progress; package-manager messages remain visible.
Outils installés. AlignView les détecte sans activer un environnement.|Tools installed. AlignView detects them without activating an environment.
Plan uniquement. Ajouter --install-tools pour exécuter cette commande.|Plan only. Add --install-tools to run this command.
Exécuter cette installation dans le préfixe dédié ?|Run this installation in the dedicated prefix?
Environnement dédié|Dedicated environment
Défaut du moteur|Engine default
Valeurs de départ rétablies pour {engine}|Initial values restored for {engine}
Paramètres de calcul · {engine}|Calculation parameters · {engine}
Paramètres {engine}|Parameters {engine}
Longueur minimale en acides aminés|Minimum length in amino acids
Début (1-indexée)|Start (1-based)
Fin (incluse)|End (inclusive)
Complément inverse ?|Reverse complement?
Métrique : gc, entropy ou gap|Metric: gc, entropy or gap
Fenêtre|Window
Pas|Step
Longueur k du motif exact|Exact-match word length k
Motif IUPAC (ex. ATG, NNNGTT)|IUPAC motif (e.g. ATG, NNNGTT)
Conserver une copie du projet avant de quitter.|Save a copy of the project before quitting.
RAPPORT / UN SEUL FICHIER HTML|REPORT / SINGLE HTML FILE
OUVRIR / RAPPORT AUTONOME|OPEN / STANDALONE REPORT
Fichier HTML de destination|Destination HTML file
Le HTML sera conservé indépendamment de la session (action explicite).|The HTML will be retained independently of the session (explicit export).
Le rapport existe. Remplacer le fichier ?|The report exists. Replace the file?
Copie annulée.|Copy cancelled.
Ouvrir en local|Open locally
Le HTML existe déjà|The HTML already exists
Page d’accueil|Home page
Voir les données|Inspect data
Les fichiers complets TSV/JSON, alignements FASTA, cartes SVG et logs sont disponibles uniquement si la session a été enregistrée.|Complete TSV/JSON files, FASTA alignments, SVG maps and logs are available when the session has been saved.
Un fichier à consulter et partager dans un navigateur|A file to view and share in a browser
Rapport autonome|Standalone report
Aucun appel réseau · données locales|No network requests · local data
Une absence de HSP ne signifie pas une absence d'homologie.|No HSP does not prove the absence of homology.
Résumé de la recherche|Search summary
Qualité de la requête|Query quality
Profils le long de la requête|Profiles along the query
Explorer les références|Explore references
Trier|Sort
Bitscore|Bitscore
Couverture Q|Query coverage
Identité|Identity
Cadres|Frames
Positions couvertes et non couvertes|Covered and uncovered positions
Alignement local, colonne par colonne|Local alignment, column by column
Montrer seulement les blocs contenant des différences|Show only blocks containing differences
Statistiques FASTA|FASTA statistics
Type probable|Predicted type
Sans gaps|Ungapped
Profil GC / entropie|GC / entropy profile
Données et résultats disponibles|Available data and results
Qualité de la première séquence|First-sequence quality
Un alignement est disponible dans la visionneuse.|An alignment is available in the viewer.
Aucun alignement multiple n’a été calculé.|No multiple sequence alignment has been calculated.
À calculer|Not calculated
Fichier|File
Rapport autonome|Standalone report
{count} séquences intégrées au rapport.|{count} sequences included in the report.
Paramètres de cette recherche|Search parameters
Réglages actuels et paramètres des calculs effectués|Current settings and parameters of completed calculations
Commandes externes exécutées|Executed external commands
Reproductibilité|Reproducibility
Exporter les résultats|Export results
Exports depuis le navigateur|Browser exports
Tous les fichiers sont créés localement dans ton navigateur, sans serveur.|All files are created locally in your browser without a server.
Format recommandé selon l’usage|Recommended formats by use
Besoin|Purpose
Pourquoi|Why
Figures pour article|Figures for a paper
Vectoriel, exportable indépendamment depuis l'atlas|Vector graphics, independently exportable from the atlas
Statistiques pour R / Python|Statistics for R / Python
Données tabulaires lisibles|Readable tabular data
Régions alignées / non couvertes|Aligned / uncovered regions
Coordonnées clairement documentées|Explicitly documented coordinates
MSA vers d’autres logiciels|MSA for other software
Compatibilité interopérable après alignement|Interoperable alignment exchange
Rapport à partager|Report to share
HTML autonome|Standalone HTML
Une seule page, sans réseau ni serveur|One page, without network or server
Comprendre ce qui est calculé|Understand what is calculated
Méthodes & limites|Methods and limitations
Alignement global et multiple|Global and multiple alignment
Statistiques MSA et phylogénie|MSA statistics and phylogeny
BLASTX et coordonnées|BLASTX and coordinates
Couverture, chevauchement et zones non détectées|Coverage, overlaps and undetected regions
Exporter SVG|Export SVG
Réglages utilisés|Settings used
'''
EN.update(line.split('|',1) for line in _STATUS.strip().splitlines() if line and '|' in line)

EN.update({
    'E ≤ {evalue} · id ≥ {identity}% · couverture QUERY ≥ {coverage}%':'E ≤ {evalue} · identity ≥ {identity}% · query coverage ≥ {coverage}%',
    'MSA · {count} séquences × {columns} colonnes · colonnes {start}–{end}':'MSA · {count} sequences × {columns} columns · columns {start}–{end}',
    'Table NCBI {code} ; export des séquences protéiques':'NCBI table {code}; export protein sequences',
    'dust : faible complexité ; none : aucun. Les entrées internes sont normalisées en majuscules.':'dust: low complexity; none: no masking. Internal input sequences are normalized to uppercase.',
    'Tables de traduction':'Translation tables',
    'Codes NCBI, cadres, initiateurs et arrêts contextuels':'NCBI codes, frames, initiators and context-dependent stops',
})

_REPORT = '''
{count} observations · {low} à {high}|{count} observations · {low} to {high}
{count} observations|{count} observations
Séquences et alignements · {version}|Sequences and alignments · {version}
Référence|Reference
Unités :|Units:
Unités de cette recherche :|Search units:
← Résumé|← Overview
◐ Thème|◐ Theme
limitées|limited
Le rapport est un instantané des données enregistrées.|The report is a snapshot of saved data.
Analyse exploratoire en fenêtres glissantes. Le GC ne s'applique qu'aux séquences nucléotidiques.|Exploratory sliding-window analysis. GC applies only to nucleotide sequences.
Attention : la couverture est calculée par union des coordonnées des HSP retenus, sans compter plusieurs fois les chevauchements. Identité = proportion dans le HSP, pas identité sur l'ensemble des deux séquences. Une région sans HSP peut néanmoins être homologue.|Coverage is the union of retained HSP coordinates, counting overlaps only once. Identity is the fraction within an HSP, not across both complete sequences. A region without an HSP may still be homologous.
Calcul sur résidus observés, gaps et caractères manquants exclus.|Calculated from observed residues, excluding gaps and missing characters.
Chaque HSP est mesuré séparément, les chevauchements ne sont pas soustraits.|Each HSP is measured separately; overlaps are not subtracted.
Comparaison des séquences alignées ; les gaps ne sont pas comptés comme substitutions.|Comparison of aligned sequences; gaps are not counted as mismatches.
Consensus calculé par colonne ; n’est pas une séquence biologique observée.|Column-wise consensus; not an observed biological sequence.
Conservation calculée sur les résidus non-gap, dépend des séquences retenues.|Conservation uses non-gap residues and depends on the selected sequences.
Coordonnées QUERY/RÉFÉRENCE en {query}/{subject} ; alignement en {aligned}.|Query/reference coordinates in {query}/{subject}; alignment in {aligned}.
Coordonnées sur QUERY dans son alphabet d’origine, jamais des colonnes alignées.|Query coordinates use its native alphabet, never alignment columns.
Couverture calculée dans les coordonnées propres à chaque référence.|Coverage uses each reference's own coordinate system.
Couverture par référence · deux coordonnées indépendantes|Coverage per reference · two independent coordinate systems
Couverture relative des deux séquences|Relative coverage of both sequences
Distance p · gaps exclus · {count} séquences max|p-distance · gaps excluded · up to {count} sequences
Entropie locale (fenêtre ≈ {window} colonnes)|Local entropy (window ≈ {window} columns)
Exporter les paramètres, les checksums et les logs. Les figures seules ne suffisent pas à reproduire la recherche.|Export parameters, checksums and logs. Figures alone cannot reproduce the search.
Figure exploratoire ; interprétation conditionnée par les entrées et paramètres.|Exploratory figure; interpretation depends on inputs and parameters.
HSP et coordonnées.|HSPs and coordinates.
Qualité des HSP.|HSP quality.
Identité locale (%)|Local identity (%)
La couverture QUERY est calculée par union des intervalles couverts par référence. Les intervalles non couverts représentent des positions sans HSP détecté, pas une preuve de non-homologie.|Query coverage is the union of intervals covered for each reference. Uncovered intervals mean no HSP was detected, not proof of non-homology.
La longueur du segment est exprimée dans l’unité propre aux références.|Segment length uses the reference's native unit.
Les blocs search_run et calculation_runs conservent les réglages ayant produit les résultats. Les autres valeurs sont les réglages actuels pour le prochain calcul. null/native signifie : option laissée au moteur installé.|search_run and calculation_runs retain the settings that produced existing results. Other values are current settings for the next calculation. null/native leaves an option to the installed engine.
Les exports ci-dessous utilisent les données intégrées au rapport : elles peuvent être|The downloads below use data embedded in the report, which may be
Les figures sont calculées à partir des données de la session. Les statistiques des HSP représentent des|Figures use the session data. HSP statistics describe
Les gaps peuvent refléter de la biologie, des longueurs différentes ou des artefacts d’alignement.|Gaps can reflect biology, different lengths or alignment artifacts.
Les intervalles sont calculés à partir de l'union des HSP. Les positions d'une requête BLASTX sont nucléotidiques.|Intervals use the union of HSPs. BLASTX query positions are nucleotides.
Les valeurs dépendent de la longueur des alignements locaux ; comparaison descriptive.|Values depend on local-alignment lengths; this comparison is descriptive.
Localisent des HSP. Les scores, E-values et identités concernent les segments détectés. Les E-values dépendent de l’espace de recherche.|Locate HSPs. Scores, E-values and identities refer to detected segments. E-values depend on the search space.
Needle réalise un alignement global pairwise. MAFFT, MUSCLE et Clustal Omega construisent des MSA et peuvent différer. L’alignement n’est pas une recherche d’homologues.|Needle performs global pairwise alignment. MAFFT, MUSCLE and Clustal Omega build MSAs and can differ. Alignment is not a homolog search.
Nombre de colonnes contenant au moins un gap ; en BLASTX, les colonnes représentent des acides aminés.|Columns containing at least one gap; in BLASTX, columns represent amino acids.
Nombre de références avec HSP selon la position QUERY|References with an HSP along query positions
Pistes indépendantes : les positions QUERY et RÉFÉRENCE utilisent leurs unités propres (nt ou aa).|Independent tracks: query and reference positions use their own units (nt or aa).
Plusieurs HSP peuvent se chevaucher : ce nombre ne représente pas des régions indépendantes.|HSPs can overlap: this count does not represent independent regions.
Position exprimée en colonnes du MSA ({unit}), pas en positions natives de QUERY.|Positions are MSA columns ({unit}), not native query positions.
Positions RÉFÉRENCE ({unit})|Reference positions ({unit})
Positions QUERY ({unit})|Query positions ({unit})
Pourcentage de positions QUERY touchées par au moins un HSP de la référence.|Percentage of query positions covered by at least one HSP from the reference.
Pourcentages calculés sur des dénominateurs distincts ; ne pas les convertir en longueurs comparables.|Percentages use different denominators; do not turn them into directly comparable lengths.
QUERY : {query} · RÉFÉRENCE : {subject} · alignement HSP : {aligned} · intervalles 1-based inclusifs|Query: {query} · Reference: {subject} · HSP alignment: {aligned} · 1-based inclusive intervals
QUERY est nucléotidique et les références sont protéiques ; les positions HSP sont donc exprimées dans des unités distinctes. Les BED concernent toujours les coordonnées de leur séquence respective.|The query is nucleotide and references are proteins, so HSP coordinates use different units. BED files always use their respective sequence coordinates.
Références ayant au moins un HSP à la position échantillonnée; ne mesure pas l’abondance.|References with at least one HSP at the sampled position; not an abundance measure.
Régions sans HSP (QUERY, {unit})|Regions without an HSP (query, {unit})
Sites variables, entropie et distances p décrivent le MSA calculé, qui doit être biologiquement homogène. FastTree produit un arbre exploratoire : il ne remplace pas une validation phylogénétique.|Variable sites, entropy and p-distances describe the calculated MSA, which must be biologically coherent. FastTree produces an exploratory tree, not a replacement for phylogenetic validation.
Un dotplot représente les positions partageant un même mot de longueur k. Les répétitions abondantes sont limitées pour préserver la lisibilité. Il ne remplace pas un alignement.|A dotplot shows positions sharing the same word of length k. Frequent repeats are limited for readability. It does not replace an alignment.
Un résultat local n'est pas une identité globale. Deux séquences peuvent partager un domaine très similaire alors que leurs régions non alignées sont longues ou divergentes.|A local result is not whole-sequence identity. Two sequences may share a highly similar domain while having long or divergent unaligned regions.
Une entropie élevée sur des colonnes peu occupées doit être interprétée avec prudence.|High entropy in sparsely occupied columns requires cautious interpretation.
Vert : au moins un HSP · Gris : aucune région détectée par la recherche|Green: at least one HSP · Gray: no region detected by this search
Vert = positions QUERY couvertes · bleu = positions référence couvertes · gris = sans HSP détecté|Green: covered query positions · Blue: covered reference positions · Gray: no detected HSP
aux 2 000 premiers HSP, 301 statistiques de séquences et 1 000 séquences de MSA. Pour les résultats complets, utiliser le centre d’export du logiciel.|to the first 2,000 HSPs, 301 sequence-quality records and 1,000 MSA sequences. Use application exports for complete results.
et non l’identité globale des génomes ou protéines. Les données MSA sont disponibles après alignement.|rather than whole-genome or whole-protein identity. MSA data are available after alignment.
Étendue de chaque HSP / QUERY (%)|Each HSP span / query (%)
Étendue de chaque HSP sur QUERY|Each HSP's query span
Étendue de chaque HSP sur sa référence|Each HSP's reference span
Étendue relative QUERY (%)|Relative query span (%)
Arbre FastTree issu du MSA ; longueurs en substitutions/site lorsqu’elles sont disponibles. Racine de dessin arbitraire, sans enracinement biologique. Supports locaux SH-like, pas des bootstraps.|FastTree tree from the MSA; branch lengths are substitutions/site when available. The drawing root is arbitrary, without biological rooting. Local SH-like supports are not bootstraps.
alignements locaux|local alignments
colonnes {start}–{end}|columns {start}–{end}
{count} séquences × {columns} colonnes|{count} sequences × {columns} columns
'''
EN.update(line.split('|',1) for line in _REPORT.strip().splitlines() if line and '|' in line)

EN.update({
    'Les analyses VSEARCH sont prêtes : ouvrir Explorer pour les consulter.':'VSEARCH analyses are ready: open Explore to inspect them.',
    'L’alignement est prêt : ouvrir Explorer pour le lire.':'The alignment is ready: open Explore to read it.',
    'La recherche est terminée : ouvrir Explorer pour examiner les correspondances.':'The search is complete: open Explore to inspect matches.',
    'Commencer par Analyser, ou consulter les séquences dans Explorer.':'Start with Analyze, or browse sequences in Explore.',
    'Exports scientifiques prêts : {count} fichiers · profil {profile}':'Scientific exports ready: {count} files · profile {profile}',
    'Ajoute EMBOSS, MUSCLE, Clustal Omega, MMseqs2, FastTree, VSEARCH, ALV':'Adds EMBOSS, MUSCLE, Clustal Omega, MMseqs2, FastTree, VSEARCH and ALV',
    'Figure descriptive : consulter sa légende dans le rapport HTML':'Descriptive figure: see its caption in the HTML report',
    'Préciser un FASTA en mode --no-menu':'Provide a FASTA in --no-menu mode',
    '--install-tools nécessite --setup PROFIL':'--install-tools requires --setup PROFILE',
    '--resume ne se combine pas avec des fichiers ou --demo':'--resume cannot be combined with input files or --demo',
    '--html-out attend un nom de fichier avec l’extension .html':'--html-out requires a filename ending in .html',
    '--mode msa uniquement avec un seul FASTA':'--mode msa requires a single FASTA',
    'Ce mode nécessite un FASTA de références.':'This mode requires a reference FASTA.',
    '--msa view nécessite un seul FASTA déjà aligné.':'--msa view requires a single aligned FASTA.',
    'Pour changer de méthode ou de séquence après --resume, utiliser les menus Analyser / Données.':'To change the method or query after --resume, use the Analyze / Data menus.',
    'Utiliser --vsearch sans mode BLAST/DIAMOND ; les analyses ont des métriques différentes.':'Use --vsearch without a BLAST/DIAMOND mode; these analyses use different metrics.',
    'Choisir un fichier — Alignement déjà calculé':'Choose a file — Existing alignment',
    'Choisir un fichier — Séquences à aligner':'Choose a file — Sequences to align',
    'Choisir un fichier — Fichier de requêtes':'Choose a file — Query file',
    'Choisir un fichier — Séquences de référence':'Choose a file — Reference sequences',
    'Dossier : {path}':'Directory: {path}',
    'Tu peux coller un chemin complet ; un dossier ouvre son contenu':'Paste a complete path; a directory opens its contents',
    'FASTA invalide : {error}':'Invalid FASTA: {error}',
    '↑↓ choisir · / rechercher · Entrée ouvrir · p chemin · q annuler':'↑↓ select · / search · Enter open · p path · q cancel',
})
