# Use every feature

## Start and navigate

Choose an aligned FASTA, raw sequences, query/reference files, a saved project or the bundled demo. A FASTA with equal-length records is not proof of alignment. FASTA and gzip FASTA are supported inputs; NEXUS, PHYLIP, Stockholm and Clustal are export formats.

Home keeps five entries: **Explore**, **Analyze**, **Data**, **Save and share**, **Settings and help**. The context identifies the active query, available results and unsaved changes. Esc/Left goes back, F2 returns Home, F1 explains a choice, and `/` filters a long list. Numeric choices are confirmed with Enter. Ctrl+C cancels an input/action where possible while retaining the open project.

## Queries and references

For comparison, the first FASTA contains possible queries and the second contains references. If several queries are present, choose one from the searchable selector before calculation, including in the guided workflow. Batch runs require `--query-id` with the first word after `>` in the chosen header. A search compares one active query at a time.

**Data → Change the active query** opens another comparison context. Previous analyses remain in session history until saving. With a single FASTA, selecting an active sequence changes targeted analyses such as translation, extraction and motifs, while retaining every record for MSA creation.

## Similarity searches

| Query | References | Methods | Query/reference/HSP coordinates |
| --- | --- | --- | --- |
| DNA/RNA | DNA/RNA | BLASTN, or VSEARCH global search | nt/nt/nt |
| Protein | Proteins | BLASTP, DIAMOND BLASTP, MMseqs2 | aa/aa/aa |
| DNA/RNA | Proteins | BLASTX, DIAMOND BLASTX | nt/aa/aa |
| Protein | DNA/RNA | TBLASTN | aa/nt/aa |
| DNA/RNA | DNA/RNA, translated comparison | TBLASTX | nt/nt/aa |

Before a local search, select the method, adjust parameters and inspect commands. The preparation screen shows the engine, version, source and path when available. Missing tools offer installation or explicit executable paths.

A new local search retains up to 100 references by default. Standard/More permissive/Stricter presets adjust E-value and DIAMOND sensitivity, not that maximum. Identity, coverage and bitscore filters reuse calculated HSPs. Changing E-value, seeds, masking or engine parameters requires a new search. Changed result selection can invalidate an MSA and tree derived from it.

**Explore → Matches and coverage** opens references and their HSPs, query coverage, reference coverage, identity/coverage comparisons and segment diagnostics. Identity is local to an HSP; coverage is a union of coordinate intervals. No detected match is not proof of no homology.

## Multiple alignments

A single FASTA aligns all its sequences. With query/reference files, choose the intended scope:

| Scope | Contents |
| --- | --- |
| `--msa hits` | Query and retained references |
| `--msa all` | Query and all references |
| `--msa fragments` | Align references, then add the partial query with MAFFT |
| `--msa hsp` | Local segments sufficiently overlapping the query anchor region |
| `--msa view` | Read one existing aligned FASTA without calculation |

`--msa align` is convenient for a single raw FASTA. With two FASTA files in search modes, `align`, `full` and `global` are historical aliases for retained-reference alignment; `all` explicitly requests every reference. Global comparison modes use all references.

MAFFT, MUSCLE 5 and Clustal Omega have separate parameters. Fragment addition explicitly uses MAFFT in two stages. HSP alignment is exploratory: selected regions may represent different domains. `--hsp-overlap` adjusts anchor overlap.

BLASTX whole-sequence MSA translates the query in the selected frame and genetic code; stop symbols are replaced by X as announced. TBLASTN/TBLASTX require explicit biological region/frame selection; use suitable HSPs or prepared proteins instead of an arbitrary whole-genome translation.

## Viewer controls

| Action | Keys |
| --- | --- |
| Column / sequence | Left/Right / Up/Down, or h/l / k/j |
| Ten columns | Shift+Left/Right |
| Horizontal screen | Ctrl+Left/Right or [ / ] |
| Sequence page | PgUp / PgDn |
| First/last column | Home / End |
| First/last sequence | Ctrl+Home / Ctrl+End |
| Absolute column or percentage | g, then 500 or 75% |
| Sequence identifier | s |
| Exact motif / consensus motif | /, then motif or c:motif |
| Next/previous occurrence | n / N |
| Variable sites / differences from consensus | v / d |
| HTML report / HTML export | w / e |
| Help / back / Home | ? / q or Esc / F2 |

The wheel scrolls sequences; Shift+wheel scrolls columns. Click the overview to jump. Viewer position is retained. Coordinates are MSA columns, including gaps. Hiding invariant sites does not modify the alignment file.

## Pairwise comparison, statistics and trees

**Analyze → Compare two sequences** offers end-to-end Needle or best-local-segment Water. Choose a reference, inspect the matrix and gap costs, then run. Batch options are `--pairwise global|local` and `--ref N`. EMBOSS costs are independent of BLAST, DIAMOND and MAFFT. AlignView penalizes Needle ends by default; `--free-ends` relaxes them. Check frames and genetic codes for translated comparisons.

**Explore → Alignment statistics** provides conservation, consensus, pairwise p-distances and FastTree inference. P-distances exclude gaps per pair and support at most 250 sequences; different pairs can use different comparable positions. FastTree needs at least three aligned sequences. Its supports are local SH-like, not tree bootstraps. The figure's root is only a drawing choice. See [HTML](HTML.md).

## Other analyses

| Action | Input/settings | Temporary result |
| --- | --- | --- |
| Motif search | IUPAC nucleotide motif | Both-strand positions in `motifs.tsv` |
| Six-frame translation | NCBI code and frame orientation | `six_frames.faa` |
| Simple ORFs | Minimum aa length, code and start policy | `orfs.tsv`, `orfs.faa` |
| Region extraction | 1-based inclusive bounds; optional reverse complement | `query_START_END.fasta` |
| Quality and composition | Active/loaded sequences | Length, GC, ambiguity and entropy |
| Sliding profiles | gc, entropy or gap; window and step | `profile_METRIC.tsv/.svg` |
| K-mer frequencies | k from 1 to 8 | `kmers_K.tsv` |
| Dotplot | Second sequence and exact-word length | `dotplot.json/.svg` |

Interactive targeted analyses use the active sequence. Batch `--motif` uses all records of a single FASTA, or the active query in comparison mode. Some screen previews are shorter than complete files. Save the full project to retain auxiliary internal analyses; an export profile does not necessarily copy every temporary file.

Alternative translation tables are explained in [Genetic codes](GENETIC_CODES.md). VSEARCH offers global search, dereplication, clustering and chimera screening in [its own guide](VSEARCH.md); its results are not BLAST HSPs.

## Save, share, resume and reproduce

Use **Save and share** for a portable project, standalone HTML, export directory or Slurm script. Relative paths start from the launch directory. A Slurm script uses saved project inputs and is never submitted automatically. Check resources, executable paths and cluster modules before sbatch.

**Settings and help → Completed calculations and provenance** preserves the parameters, versions, commands and logs that produced a successful result. Current settings can differ. Replacing a successful MSA removes its now-stale derived tree. A failed rerun does not replace the previous valid MSA/tree. See [Files](FILES.md), [Parameters](PARAMETERS.md) and [CLI](CLI.md).
