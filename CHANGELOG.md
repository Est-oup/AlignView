# Changelog

## Unreleased — user journey review and Bioconda preparation

- Make aligned-region inspection explicit, with input provenance, native coordinates, orientation, translation frames and gaps interpreted consistently in the terminal and offline report.
- Label the native bounds of each displayed HSP page, outline its query/target intervals and export the complete selected two-row alignment as FASTA.
- Distinguish no detected segments from matches excluded by filters; enable native VSEARCH TSV export without fabricated local-search metrics.
- Make `--no-menu` and report commands noninteractive on real terminals, including alignment readers, result pagers, save questions and occupied destinations; stop claiming unsaved data is a saved project.
- Prepare noarch Bioconda recipe generation against verified source bytes, installed CLI/manual/project/report smoke tests and a concrete submission handoff; source publication remains separate.
- Run offline report checks in CI and include release reconstruction resources in the source archive.

## Unreleased — clear search inputs and targets

- Rename Find matching references to Find matching sequences and distinguish query sequences from the local target FASTA throughout search preparation.
- Keep query/target filenames and scope visible before running, in Matches and in result summaries, including compact terminal menus.
- Change the target directly from preparation while preserving selected queries, or review full paths and target-sequence previews before Run.
- Keep the original query/target inputs visible across offline report views and in grouped comparison summaries.

## Unreleased — result organization and alignment reading

- Rename completed attempts from their result summary or Commands; keep display names through save/reopen without modifying scientific files or lineage.
- Detach the complete result catalog when reopening a project, so moving or deleting the saved source folder does not lose access to earlier attempts.
- Find exact motifs across alignment gaps in the offline report, with residue and original MSA ranges, next-hit wrapping and an option for literal aligned-text search.
- Keep HTML view exports scoped to the displayed sequence page as well as the displayed columns; replace active local evidence only after a successful VSEARCH comparison.

## Unreleased — selected query comparisons

- Select one or several FASTA query records with checkboxes, all/clear controls and identifier filtering before comparing with references.
- Run a selection against one reference dataset, retaining independent query identity/coverage, parameters, commands and results in a resumable group.
- Preserve successful query results after a failure or cancellation; record per-query statuses and attempt logs.
- Reopen grouped/per-query evidence without recalculation and export an honest group HTML summary. Noninteractive `--query-id` behavior remains unchanged.

## Unreleased — explicit file and result menus

- Loading FASTA or the demo now opens a named file summary and five explained actions; no table or alignment reader opens automatically.
- Completed alignments/searches return to a result-ready summary; opening the actual reader is explicit.
- Group other analyses into small menus and choose the individual sequence before a targeted operation.
- Return to the recognizable file menu after canceled preparation or reading; expose files, previous results, optional exports and help through an explained menu.
- Closing a file without analyses no longer requests a save.

## Unreleased — free sequence exploration

- Open FASTA data directly, add datasets, and choose any analysis without a goal chooser or metro route.
- Keep each completed interactive attempt independently in Results, with frozen inputs and tree-to-alignment links; project saving stays optional.
- Restore view filters, cursors and alignment reading positions on Back; support clicking catalog controls and entries.
- Reuse translations, extracted regions, ORFs and alignments as data without manual export/reimport.
- Preview HTML without a destination prompt; use direct browser navigation and Back/Forward. Durable HTML export remains separate.
- Extend portable projects with catalog metadata and retained snapshots; preserve batch calculation behavior.

## 1.0.0 — persistent analysis workspace

- Added goal-first startup, a persistent workflow metro map, explicit reference-selection stage and Next step navigation in the terminal and offline report. Query/reference roles and required inputs are explained before calculations.

- Replaced menu-bucket navigation with persistent Alignment, Matches, Tree and Data views, contextual primary actions, searchable commands and one Save/share sheet.
- Added explicit view/align/compare/report/resume/tools task commands with concise common help and complete `--help-all`; retained legacy flag forms.
- Added reference selection from evidence for targeted alignments, trace-sequence inspection and explicit aligned-region export.
- Split CLI, session and application orchestration behind the stable entry module; replaced chained report rewrites with one payload/figure/report composition pipeline.
- Reviewed scientific result consistency, project integrity, failure handling, sequence scopes and coordinate semantics.
- Defined a supported release contract and one runtime inventory; deterministic portable builds, version-scoped checksums and actual installed-wheel/zipapp/portable verification.
- Rewrote task guides around the resulting workflow, preserving one English manual and local offline operation.

## English-only application — 12 November 2026

- English is the sole source language for menus, CLI help, prompts, errors, report HTML/JavaScript, SVG labels, comments, docstrings and documentation.
- Removed runtime translation, language selection and duplicate manuals. The single offline manual lives in `docs/`.
- GitHub issue/PR templates, workflow labels, examples and developer guides use English.
- Presentation changes preserve user sequence identifiers, native tool fields, input data and historical calculation settings.
- Figure filenames and the search preset use English names. Rebuilds of older projects regenerate current English reports without rerunning scientific engines.

## Simplified workflows and CI — 12 November 2026

- Direct Home access to Results, Compare with references, Align sequences, Save and share, Files and query, and Advanced.
- Guided preparation and immediate result viewing; diagnostics, detailed statistics, settings, tools, provenance and HPC grouped under Advanced.
- Visible query selection before comparison and direct addition of missing references.
- Individual FASTA, TSV and SVG exports without recalculation or overwriting; read-only sessions close without an unnecessary save prompt.
- Readable terminal tree topology with explicit numerical lengths and standalone SVG output.
- Updated documentation and navigation preview.
- Python 3.10 CLI documentation is independent of argparse's group-title and BooleanOptionalAction help differences. Matrix jobs no longer cancel other versions after a failure.

## Documentation, portability and analyses — 10 November 2026

- Bundled guides covering first use, workflows, files, HTML, exports, troubleshooting, terminology, VSEARCH and genetic codes; generated CLI reference and link checks.
- Offline manual support in source, wheel and zipapp distributions.
- Explicit query selection for guided multi-record comparisons; no silent first-query selection.
- Four native VSEARCH workflows: global search, dereplication, clustering and UCHIME de novo, with contextual settings, native outputs, provenance, resuming and an HTML view.
- 27 bundled NCBI codes, separate query/reference tables and configurable ORF start policies. Simple ORF detection rejects codes requiring CDS context.
- Engine, version, origin and path displayed before calculation; executable/version retained with execution records.
- Persistent standalone HTML through `--html-out`; resuming honors explicit actions and export profiles without replaying saved action flags.
- Launcher supports arbitrary working directories and symbolic links; detects a local `tools/` prefix and preserves caller-relative output paths.
- Historical sources moved under `archives/`; the portable archive contains one application folder without bundled engines or earlier outputs.

## Tree figure — 10 November 2026

- FastTree Newick automatically produces an embedded SVG, with a tree shortcut, phylogeny filter and SVG download.
- Branch lengths, scale, complete identifiers and local support values are retained; the drawing root is explicitly arbitrary.
- Missing or zero lengths produce a labeled cladogram. Replacing an alignment invalidates its old derived tree.

## Scientific parameter review — 10 November 2026

- Typed engine-specific controls shared by CLI, menus, profiles and provenance.
- BLAST, DIAMOND and MMseqs2 search controls; MAFFT, MUSCLE 5 and Clustal Omega alignment controls; separate EMBOSS and FastTree settings.
- Inspectable command preparation and successful-calculation snapshots preserve the distinction between future settings and existing results.
- Staged MSA/tree replacement, unique execution logs and explicit HSP/fragment preparation.
- Recorded baseline and modified-parameter runs on eight real engines; see [Calculation review](docs/CALCULATION_REVIEW.md).

## Navigation review — 10 November 2026

- Task-oriented navigation, visible prerequisites and consistent keyboard behavior.
- Sensitivity presets retain the independent 100-reference result limit.
- Active sequence selection, portable input copies, resumable history and explicit saving.
- Terminal and browser regression checks for navigation, result reading and exports.

## 0.1.0 preparation — 10 November 2026

- Standard-library terminal application with local/translated searches, alignments, sequence analyses and standalone reports.
- Portable projects, scientific exchange formats, engine diagnostics, HPC modules and Slurm preparation.
- Wheel, source archive, portable ZIP, zipapp and checksum build scripts.
