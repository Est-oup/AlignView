# Export formats and coordinates

## Choose a destination and profile

Use **Save and share → Export results**. Batch examples:

```bash
alignview --resume project_01 --export-profile publication --export-dir exports_01 --no-menu
alignview --resume project_01 --export-profile complete --export-dir exports_02 --no-menu
```

Standard exports contain FASTA, sequence quality, local-search tables, coverage, BED and provenance. Publication adds CSV, GFF3, alignment exchange formats, HTML and SVG figures. Complete adds detailed JSON/JSONL, unmatched regions and a SHA256 inventory. Conditional files require the underlying data. VSEARCH output directories are preserved in all profiles when available.

An export directory is not a full resumable project. Save a project to retain every internal/auxiliary analysis and history. Browser exports use the report's embedded subset, not necessarily the complete application data.

## Main files

| File | Purpose |
| --- | --- |
| sequence_qc.tsv/csv/json | Origin, identifier, alphabet, lengths, gaps, GC where appropriate, ambiguity and entropy |
| query_active.fasta / references.fasta | Active query and reference sequences |
| hsps_all.tsv / hsps_selected.tsv | Calculated and filtered local HSPs, with aligned characters and frames |
| hits_by_reference.tsv | Reference-level identity, scores and coverage summaries |
| coverage_regions.tsv | Covered/uncovered coordinate intervals |
| hsp_query.bed / hsp_reference.bed | HSP intervals in native sequence coordinates |
| hsp_features.gff3 | Match features and explicit Target attributes |
| msa_alignment.fasta / msa_consensus.fasta / msa_variable_sites.fasta | MSA, consensus and variable-column subset |
| msa_alignment.nex/.phy/.sto/.aln | NEXUS, relaxed PHYLIP, Stockholm and Clustal exchange |
| msa_phylip_ids.tsv / msa_stockholm_ids.tsv / msa_clustal_ids.tsv | Identifier mappings for format-specific aliases |
| msa_positions.tsv / msa_pairwise_distances.tsv | Per-column metrics and pairwise p-distances |
| tree_fasttree.nwk / figures/45_fasttree_tree.svg | Tree topology/lengths/supports and vector figure |
| vsearch/ | Complete VSEARCH result files, previews, parameters and provenance |
| report.html / figures/ / figures_catalog.tsv | Standalone presentation, individual SVGs and figure descriptions |
| provenance.json | Current settings, historical calculation blocks, engine versions, commands and input hashes |
| MANIFEST.tsv | File sizes and SHA256 values in the complete export; excludes itself |

## Coordinates and interpretation

| Format | Start | End | Unit |
| --- | --- | --- | --- |
| HSP/TSV/GFF3 | 1-based | Inclusive | Native query/reference unit |
| BED | 0-based | Half-open | Native sequence unit |
| MSA | 1-based column | Column, including gaps | Aligned nt or aa |

BLASTX uses nt/aa/aa for query/reference/HSP columns; TBLASTN uses aa/nt/aa; TBLASTX uses nt/nt/aa. `coordinate_units.tsv` and `coordinate_audit.tsv` preserve these distinctions. Reverse-strand coordinates need their strand/frame context.

Coverage is an interval union, not a sum of overlapping HSP lengths. Identity is local to a segment, not a recalculated whole-sequence identity. Pairwise p-distances exclude gaps/missing characters independently per pair.

VSEARCH Matches.tsv uses native VSEARCH fields (query, target, id, alnlen, qilo/qihi, tilo/tihi, qcov/tcov, qstrand, qrow/trow). It contains no invented E-value/bitscore. Clustering/dereplication UC files map records to Representatives; size annotations retain abundances. Chimera files separate reported classifications. See [VSEARCH](VSEARCH.md).
