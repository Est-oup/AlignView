# Formats, unités et interprétation des exports

## Choix rapide

| Question | Format conseillé | Détail |
|---|---|---|
| Travailler sous R / Python | TSV (ou CSV) | Lisible avec `readr::read_tsv`, pandas |
| Archiver les paramètres et HSP | JSON / JSONL | JSONL complet HSP, une observation par ligne |
| Visualiser les régions de séquences | BED | Coordonnées 0-based half-open, comme BED standard |
| Annoter des similarités dans un viewer | GFF3 | Features `match` sur QUERY ; coordonnées 1-based inclusives |
| Reprendre un alignement multiple | FASTA, NEXUS, relaxed PHYLIP | Inclut consensus et sites variables quand possible |
| Faire une figure d'article | SVG | Format vectoriel, une figure par fichier |
| Partager une exploration autonome | HTML | Un seul fichier HTML, sans réseau, contenant les figures intégrées |
| Avoir un PDF | Depuis le HTML : Imprimer / PDF | Fonction du navigateur, pas un moteur PDF installé par AlignView |

Depuis l’accueil, « Rapports et exports » ouvre le centre d’export. Le raccourci `s` enregistre un projet portable. Dans la visionneuse MSA, `e` exporte le HTML. En batch : `--export-profile standard|publication|complete --export-dir dossier`, sans `--save` pour ne conserver que les exports.

## Fichiers

- `sequence_qc.tsv/csv/json` : source du FASTA, longueur, ungapped, gaps, GC (nucléotides seulement), %N, entropie.
- `hsps_all.tsv` et `hsps_selected.tsv` : tous les HSP et les HSP après filtres, avec alignements par colonnes et frames.
- `hits_by_reference.tsv` : un résumé par référence, y compris couverture de QUERY et de la référence (union sans double comptage).
- `coverage_regions.tsv` : intervalles couverts ou sans HSP selon QUERY ou référence.
- `hsp_query.bed`, `hsp_reference.bed` : segments des HSP ; bed score de 0 à 1000.
- `hsp_features.gff3` : segments HSP sous forme de `match` GFF3, avec attribut `Target` sur les références.
- `query_without_hsp.bed`, `references_without_hsp.bed` : régions sans HSP détecté, **pas une absence prouvée d'homologie**.
- `msa_alignment.fasta`, `msa_consensus.fasta`, `msa_variable_sites.fasta` : alignement et sous-ensemble des colonnes variables.
- `msa_alignment.nex`, `msa_alignment.phy`, `msa_alignment.sto`, `msa_alignment.aln` : export pour logiciel phylogénétique ; le PHYLIP est *relaxed*, avec tables de correspondance des identifiants (PHYLIP, Stockholm et CLUSTAL).
- `msa_positions.tsv/csv`, `msa_pairwise_distances.tsv/csv` : statistiques par colonne du MSA et distances p entre séquences.
- `tree_fasttree.nwk`, `figures/45_arbre_fasttree.svg` : arbre FastTree et figure vectorielle avec longueurs, supports et échelle. La figure est intégrée au HTML ; le Newick et le SVG sont aussi inclus dans les exports publication/complet. La racine du dessin est arbitraire ; les supports SH-like sont locaux, pas des bootstraps.
- `provenance.json` : versions du projet/Python, plateforme, paramètres, commandes, modules chargés, hashes d'entrée.
- `MANIFEST.tsv` : tailles et empreintes SHA256 des exports ; auto-exclu pour éviter un checksum récursif.

Les fichiers conditionnels ne sont créés **que lorsque les données existent** : pas de NEXUS sans MSA, pas de HSP pour un FASTA unique non comparé.

## Attention à l'unité des coordonnées

| Format / résultat | Origine | Fin | Unité |
|---|---|---|---|
| HSP BLAST/TSV/GFF3 | 1-based | inclusive | QUERY en nucléotides pour BLASTX, en aa pour BLASTP |
| BED | 0-based | exclusive | même unité que la séquence d'origine du côté considéré |
| MSA | 1-based | colonne, gaps inclus | positions d'alignement, pas coordonnées génomiques |

Les distances p sont calculées en ignorant les caractères gaps/manquants dans chaque paire ; les colonnes ainsi comparables varient selon la paire. La couverture est calculée sur l'union des HSP filtrés. Les identités sont celles des **HSP locaux**, pas une identité globale recalculée.

**Rapport HTML** : pour limiter sa taille, il intègre une sélection (2 000 HSP, 301 statistiques de séquences, 1 000 séquences de MSA). La visionneuse affiche des pages de 100 séquences. Les boutons du navigateur exportent uniquement les données intégrées. Pour les résultats complets, exporter depuis l’application.

## Audit d’unités et figures

`coordinate_units.tsv` précise les unités pour le mode employé. `coordinate_audit.tsv` compare les étendues natives aux résidus non-gap des alignements, avec ratio 3 nt/aa pour les séquences traduites. Le rapport HTML autonome, les SVG et `figures_catalog.tsv` sont inclus dans les profils publication et complet, si générés. `MANIFEST.tsv` fournit les hashes SHA256 des fichiers exportés.
