# Use every feature

## Open a dataset or a task

AlignView opens a file-summary menu before any table or reader. Five explained choices cover sequence inspection, alignment, reference search, other analyses, and files/previous results. A calculation returns to a result-ready summary; Open this result opens its reader. Esc returns to file actions. No fixed scientific route is imposed. Explicit `view` declares that a FASTA is already aligned.

```bash
alignview view aligned.fasta
alignview align raw.fasta
alignview compare queries.fna references.fna
alignview resume project_01
```

FASTA and gzip FASTA are supported imports. NEXUS, PHYLIP, Stockholm and Clustal are export formats. Existing flag-based commands remain supported. `--help` covers ordinary tasks; `--help-all` and [CLI](CLI.md) cover every engine setting.

## Workspace controls

| Action | Control |
| --- | --- |
| Main file actions | Choose a labeled menu entry with arrows and Enter |
| Other analyses | File menu → Other analyses |
| Loaded files / previous results | File menu → Files and previous results |
| Add another dataset | o |
| Preview HTML without a save destination | w |
| Inspect or run the displayed contextual action | Enter |
| Prepare alignment / comparison | a / c |
| Optional save, preview or export | s |
| Search available commands | : or Ctrl+P |
| Filter a list | / |
| Select a reference in Matches | Space |
| Navigate rows and pages | Up/Down, PgUp/PgDn |
| Return from a result page | Esc / q |

The command palette exposes analysis, data/result navigation, input links, output reuse, preview, optional saving/export, advanced settings, tools and help. Unavailable commands explain their prerequisites. Esc returns to the preceding screen; cursor, list filter and reader position are retained. Menus accept mouse clicks. In the file menu, q closes; in a reader, q returns. Specialized table/reader controls are shown only after you open that view.

Switching views, filtering and opening an existing result never launch a scientific engine. Contextual preparation screens show inputs, query, sequence scope and engine before execution. Detailed native controls, command previews and reusable profiles stay in the calculation's advanced options.

## Query and reference evidence

The first comparison FASTA contains candidate queries; the second contains references. Find matching sequences offers a checkbox selection for one or several query records (Space to check, a for all, c to clear, / to find, Enter to confirm). Each selected query runs independently against the same references, with results organized in a query group. Open a query to inspect its own evidence; identity and coverage are never pooled. Noninteractive comparison still requires `--query-id ID` rather than silently taking the first record. The sequence list and Choose active query command expose this choice; Other analyses asks for the target sequence, or uses the highlighted browser record.

BLAST and DIAMOND find local HSPs; MMseqs2 is a protein-search alternative. The proposed method depends on the detected types. VSEARCH global nucleotide search is a separate method with its own identity/coverage definitions and no BLAST E-values or bitscores.

Matches shows reference evidence in a filterable list. Choose a reference to inspect its segments. Displayed best identity belongs to the strongest-scoring HSP; it is not the maximum identity of an arbitrarily short segment or whole-sequence identity. Query/reference coverage is an interval union. A reference can contain several HSPs. No detected match is not proof of non-homology.

Use **Space** to select references for an MSA, then **Align selected references**. This selection changes the next alignment's scope, not the original search or its statistical meaning. The preparation screen identifies the selected inputs. A full-reference alignment remains a distinct choice.

Comparison preparation removes pre-existing gap characters (`-` and `.`), converts nucleotide U to T, and maps unknown `?` positions to N or X for the detected alphabet. Comparison coordinates and coverage denominators use those prepared ungapped records; original source copies remain available for provenance. MSA columns are a separate coordinate system.

The initial retained-reference limit is 100. Standard, sensitive and strict presets adjust search thresholds without changing that result limit. Native engine candidate limits are separate controls and can affect what the search finds. Filters reuse calculated HSPs. Previously completed alignments and trees remain in Results with their original scope. New search settings apply to the next calculation while historical provenance remains associated with existing results.

## Build an alignment

In a single FASTA, Space checks records for alignment. Two or more checks align that exact selection; no checks align the complete dataset. The active sequence for targeted tools does not otherwise reduce MSA scope. With references, choose the scope: all sequences, retained references, selected references, fragments or detected segments. MAFFT, MUSCLE 5 and Clustal Omega have separate settings. Fragment addition explicitly uses MAFFT in two steps.

Whole-sequence BLASTX alignment translates the query using its selected frame and NCBI code; announced stop-to-X handling prepares protein alignment input. TBLASTN/TBLASTX require biological frame/region selection; use appropriate HSPs or prepared proteins rather than arbitrary whole-sequence translation. HSP alignment remains exploratory because local fragments can represent different domains.

## Read and extract alignment evidence

The full-screen viewer supports horizontal and vertical movement, identifier/motif search, variable-site filtering, a chosen trace/reference sequence and explicit aligned-region export. Coordinates remain MSA columns, including gaps. Hiding invariant sites changes the view, not the alignment file.

| Action | Viewer keys |
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
| Variable sites / differences | v / d |
| Use current row as trace / return to consensus | r |
| Filter rows by identifier | f |
| Check/uncheck a row for export | Space |
| Mark region start/end | m at each boundary |
| Clear region and row checks | x |
| Export the selected rows/region | e |
| Save a project snapshot | Ctrl+S |
| Help / back | ? / q or Esc |

The wheel scrolls rows; Shift+wheel scrolls columns. Clicking the overview jumps to a region. The workspace retains view state when you return. Press m at the first and last MSA columns, check rows with Space if needed, and press e to export. The new folder contains `alignment.fasta`, `columns.tsv` and `selection.json`, including selected identifiers, original MSA columns and the source alignment hash. Explicit exports describe rows and columns; they do not silently create an ungapped native-coordinate interval.

## Tree and detailed statistics

Create an exploratory tree uses FastTree on an existing alignment with at least three rows. The Tree view reuses its Newick result; opening it does not rerun FastTree. The SVG preserves branch lengths and complete identifiers. Its drawing root is arbitrary, supports are local SH-like values rather than bootstraps, and missing/all-zero lengths produce an explicitly labeled cladogram.

Additional analyses and settings contains conservation, consensus, sequence quality, detailed segment diagnostics and pairwise p-distances. Distances exclude gaps per pair and support at most 250 sequences; different pairs can use different comparable positions.

## Supplementary analyses

Search **Commands** directly for VSEARCH nucleotide analyses, Compare two sequences, Find a nucleotide motif, Translate all six frames, Find simple ORFs, Extract a query region, Create a dotplot or Choose translation table. **Additional analyses and settings** also contains composition, profiles, detailed scientific settings and HPC preparation.

| Action | Scope/settings | Result |
| --- | --- | --- |
| VSEARCH | Nucleotide global search, dereplication, clustering or chimera screening | Native files under `vsearch/` |
| Needle / Water | Chosen query/reference, separate EMBOSS matrix and gap costs | Global/local pairwise alignment |
| Motif search | IUPAC nucleotide motif, both strands | `motifs.tsv` |
| Six-frame translation | NCBI table and orientation | `six_frames.faa` |
| Simple ORFs | Minimum aa length, table and start policy | `orfs.tsv`, `orfs.faa` |
| Region extraction | 1-based inclusive coordinates and optional reverse complement | Region FASTA |
| Composition and profiles | GC, entropy, gaps, k-mer size, window and step | Tables and SVG profiles |
| Dotplot | Compatible sequence alphabets and exact-word length | `dotplot.json`, SVG |

Targeted interactive analyses use the active sequence. Batch `--motif` uses every record of a single FASTA, or the active query in a comparison. Some previews are bounded; projects retain complete generated files. See [VSEARCH](VSEARCH.md), [Genetic codes](GENETIC_CODES.md), [Parameters](PARAMETERS.md) and [Exports](EXPORTS.md).

## Save, share and reproduce

Completed interactive attempts accumulate in Results without requiring a project directory. Save / Export separately offers a resumable project, temporary HTML preview, durable HTML export and individual exports. Ctrl+S asks for a project directory on the first save. Later modified saves create a new sibling such as `<name>_snapshot_001`; earlier copies are unchanged. Saving an unchanged project creates no duplicate snapshot. FASTA alignment, reference evidence TSV, tree SVG and a full result folder serve different downstream uses. Paths are relative to the caller's working directory. Batch outputs and alignment-region exports refuse an existing destination. Interactive HTML sharing asks explicitly before replacing a report.

Every completed interactive analysis has an independent temporary snapshot. Adding datasets or opening older results retains the catalog; Back restores the previous screen. Save a project to retain the complete catalog after closing. New snapshots include a `project_manifest.json` inventory of file sizes and SHA256 hashes covering inputs, state, results, logs and nested history. It detects accidental changes; it is not authenticated proof against deliberate rewriting of both data and checksums. Older valid projects without the manifest remain readable. Resume alone reruns no engines; explicit CLI actions can request new calculations. A new alignment gets a separate entry. Earlier alignments and their trees remain readable with their original input links. A failed rerun preserves previous completed entries.

Tools and execution history exposes actual executables, versions, commands and logs. Profiles contain validated future settings; each completed calculation preserves its historical settings separately. Slurm preparation creates a script but never submits it. See [Files](FILES.md), [HTML](HTML.md) and [Validation scope](AUDIT.md).
