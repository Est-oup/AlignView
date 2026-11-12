# Changelog

## Documentation, portabilité et analyses — 10 novembre 2026

- Guides complets français/anglais, premiers pas, parcours, fichiers, HTML, exports, dépannage, glossaire, VSEARCH et codes génétiques ; références CLI générées et contrôles de liens.
- Documentation hors ligne dans le terminal, le wheel et le zipapp ; choix de langue par --lang fr/en, variable ALIGNVIEW_LANG ou menu.
- Sélection de la requête avant comparaison multirequête guidée ; liste filtrable et absence de choix implicite de la première séquence.
- Quatre analyses VSEARCH nucléotidiques : recherche globale, déduplication, regroupement et UCHIME de novo ; paramètres contextualisés, fichiers natifs, provenance, reprise et vue HTML dédiées.
- 27 tables NCBI embarquées, codes séparés QUERY/références, initiateurs ORF réglables ; arrêts contextuels signalés et refusés par le détecteur d’ORF simple.
- Moteur, version, origine et chemin visibles avant calcul ; exécutable et version rappelés à l’exécution.
- HTML autonome conservé avec --html-out ; options explicites de reprise et profil d’export respectés sans relance implicite des anciens calculs.
- Lanceur fonctionnel depuis tout dossier et via un lien symbolique ; détection d’un préfixe tools/ local, sorties relatives au dossier de lancement.
- Anciennes versions regroupées sous archives/ ; archive portable avec un unique dossier AlignView/, sans moteurs ni anciennes sorties.

## Figure d’arbre — 10 novembre 2026

- Figure SVG générée automatiquement avec l’arbre FastTree et intégrée au rapport HTML autonome.
- Accès « Voir l’arbre », domaine Phylogénie, téléchargement SVG et inclusion dans les exports publication/complet.
- Longueurs de branches, échelle en substitutions/site, identifiants complets et supports locaux conservés ; racine de dessin explicitement arbitraire.
- Cladogramme signalé lorsque les longueurs sont absentes ou toutes nulles ; aucune figure périmée après remplacement de l’alignement.
- Les projets contenant déjà un Newick gagnent la figure à la reconstruction du rapport, sans nouveau calcul FastTree.

## Révision des paramètres scientifiques — 10 novembre 2026

- Réglages propres à chaque moteur accessibles avant recherche, MSA, Needle/Water et arbre ; aide sur leurs effets et leurs défauts.
- Graines, masquage, scores BLAST ; contrôles DIAMOND ; sensibilité MMseqs2 réelle et limite de candidats transmise.
- Stratégies, gaps et raffinements MSA ; choix MUSCLE 5 et itérations Clustal Omega.
- Stratégie du MSA des HSP respectée, seuil d’ancre réglable, pipeline de fragments MAFFT explicitement annoncé.
- Coûts terminaux Needle séparés ; modèle, Gamma20, support local et graine FastTree.
- Aperçu des commandes sans moteur, profils JSON validés, priorité aux options CLI explicites et paramètres Slurm conservés.
- Provenance par calcul réussi, anciens réglages distingués des prochains, journaux de relance conservés et sorties MSA/arbre validées avant remplacement.
- Comparaison des interfaces Jalview/AliView/UGENE et essais des huit moteurs réels avec valeurs de départ et modifiées.
- Limite de 100 références maintenue et code génétique standard indiqué explicitement.

## Révision de la navigation — 10 novembre 2026

- Cinq rubriques stables ; démarrage par intention et aperçu avant calcul.
- Contexte permanent, séquence active, sauvegarde, prérequis et retour F2.
- Navigation centralisée, menus dupliqués supprimés, descriptions lisibles.
- Touches rapides conservées ; aides imbriquées restaurées ; aucun rafraîchissement inutile au repos.
- 100 références par défaut dans CLI, menus et sélection ; limite indépendante des profils.
- Séquence active choisie aussi dans un FASTA unique, sans réduire l’entrée du MSA.
- Paramètres calculés distincts des réglages futurs ; historique temporaire distinct d’une sauvegarde.
- Rapport guidé, vues vides masquées, navigation du navigateur et figures dépliables.
- Qualité et matrices consultables au-delà des anciens aperçus de 45 / 12 lignes.

## 0.1.0 — préparation du 10 novembre 2026

Première numérotation publique d’**AlignView**, issue des prototypes internes v2–v6.

- Nom simplifié dans le terminal, les rapports, la CLI et les distributions.
- Navigation commune : UTF-8, touches CSI complètes, modificateurs, souris, pages, mémoire de sélection et rendu borné au terminal.
- Visionneuse MSA : déplacements rapides, molette, vue globale cliquable, accès en pourcentage ou par identifiant, motif précédent/suivant et position conservée.
- Sorties longues avec recherche et défilement dans les deux axes.
- Sélecteurs sans troncature arbitraire ; prise en charge des FASTA gzip.
- Changement de requête ou de fichiers avec historique temporaire des analyses.
- Projets portables : entrées, SHA256, reprise sans recalcul, sauvegarde atomique et récupération après erreur.
- Installation Conda/Mamba/Micromamba par outil ou profil, préfixe géré détecté sans activation du shell.
- Diagnostic détaillé ou JSON, aide ciblée pour les moteurs manquants, versions dans la provenance.
- Arrêt des sous-processus scientifiques lors d’une interruption.
- HTML : défilement vertical conservé, pavé tactile horizontal, tailles de pages et pagination des séquences.
- Wheel, sdist, zipapp, exemples embarqués, CI, recette Bioconda et documentation de publication.
- Audit comparatif, benchmark reproductible et validations réelles distinctes des simulations.
