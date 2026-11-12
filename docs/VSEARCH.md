# VSEARCH nucleotide analyses

Use **Advanced → Additional analyses → VSEARCH analyses**, or select VSEARCH as the nucleotide comparison method. Proteins are rejected before execution. The preparation screen shows the engine, version, source and executable; commands and normalized inputs can be previewed. Install with `--setup vsearch --install-tools` when needed.

## Four workflows

| Action | Input | Main outputs |
| --- | --- | --- |
| search | Active query against the reference FASTA | vsearch/search/Matches.tsv |
| dereplicate | All records of the first FASTA | vsearch/dereplicate/Unique sequences.fasta and dereplication.uc |
| cluster | All records of the first FASTA | vsearch/cluster/centroids.fasta and clusters.uc |
| Chimeras | All first-FASTA records, dereplicated before UCHIME | vsearch/Chimeras/Unique sequences.fasta, Chimeras.fasta, Non-chimeras.fasta, Borderline.fasta and uchime.tsv |

The actions preserve source files and the active MSA. They produce their own native metrics, not BLAST HSPs. Several actions can coexist in one project. Repeating an action replaces that action's result only after successful execution and validation; logs retain separate attempts.

## Search and query selection

Multi-query comparison opens a searchable query selector before calculation, including guided comparison. Use Data to change the active query later. Batch selection is explicit:

```bash
alignview queries.fna refs.fna --vsearch search --query-id chosen_id --vsearch-id 0.97 --report --save --out vsearch_01 --no-menu
```

Search uses global pairwise alignment and a word prefilter. It does not enumerate all motif occurrences inside a long target. The identity threshold is a fraction from 0 to 1. Reported identity is a percentage under iddef (default 2). qcov/tcov follow VSEARCH's coverage definitions, rather than a union of local HSPs. No nucleotide E-value or bitscore is calculated or invented; BLAST E-value/bitscore filters do not apply. `--hits` limits reported references separately from candidate acceptance/rejection limits.

For an intentionally exhaustive small comparison, set all three controls to zero:

```bash
alignview query.fna refs.fna --vsearch search --vsearch-id 0.9 --vsearch-maxaccepts 0 --vsearch-maxrejects 0 --vsearch-minwordmatches 0 --vsearch-masking none --save --out exhaustive_01 --no-menu
```

This can cost much more on large inputs. Lowering identity alone does not remove the word prefilter.

## Dereplication, clustering and counts

```bash
alignview reads.fna --vsearch dereplicate --save --out unique_01 --no-menu
alignview reads.fna --vsearch cluster --vsearch-id 0.97 --report --save --out clusters_01 --no-menu
alignview counted.fna --vsearch cluster --vsearch-sizein --save --out clusters_02 --no-menu
```

Counts are written as ;size=N; annotations. By default, each FASTA record counts once. Enable `--vsearch-sizein` when input annotations already contain counts to sum. A 97% cluster is an algorithmic group, not a species identification.

## De novo chimera screening

```bash
alignview reads.fna --vsearch Chimeras --report --save --out chimeras_01 --no-menu
```

AlignView dereplicates first and retains abundances, then runs UCHIME de novo. Reported Chimeras, non-Chimeras and Borderline sequences are separated. Counts, dataset composition, available parents, masks and thresholds affect detection. Uniform-abundance references do not substitute for abundance-bearing amplicon data. No detected chimera is not a guarantee of no Chimeras. De novo UCHIME works on the plus strand, independently of the strand setting used by other steps.

## Preparation and retained results

Alignment gaps are removed and U is normalized to T; this preparation is announced and recorded. Length bounds are explicit and exclusions are counted. Internal FASTA normalization uses uppercase, so only dust/no masking is offered rather than pretending to preserve lowercase masks. Genetic-code settings do not translate VSEARCH inputs: these remain nucleotides.

The HTML includes summaries, up to 200 preview rows per action, settings, actual commands and versions. Complete TSV/UC/FASTA files are retained in the project's vsearch/ directory and in every application export profile. The workflows were exercised with VSEARCH 2.32.0; external-version compatibility remains subject to engine option support and logs.

Engine references: [global search](https://torognes.github.io/vsearch/commands/vsearch-usearch_global.1.html), [output fields](https://torognes.github.io/vsearch/misc/vsearch-userfields.7.html), [command manual](https://torognes.github.io/vsearch/).
