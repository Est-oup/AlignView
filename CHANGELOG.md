# Changelog

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
