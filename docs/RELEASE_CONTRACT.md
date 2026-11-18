# Version 1.0 release contract

This contract defines the supported behavior of AlignView 1.0. It is an engineering acceptance checklist, not a claim that every biological dataset, terminal, cluster or engine version has been validated.

## The core workflow

AlignView loads sequence datasets before asking for analyses. A file-summary menu explains five possible actions before opening any data table or reader. After calculation, a result-ready summary offers an explicit Open this result action and a return to the file menu. Independent completed attempts remain accessible without saving a project. Back retains reading context. Native scientific prerequisites apply without imposing a global route.

The command palette exposes every operation, input links and output reuse. HTML preview requires no export destination; optional project saving and durable exports serve later sessions and external reuse. The CLI retains explicit task verbs and documented flag commands.

## Data and scientific behavior

- Input accepts FASTA and gzip-compressed FASTA. Equal lengths alone are not proof of alignment.
- An interactive comparison can select one or several query records. Each native search uses one query and retains independent metrics; a group links the resulting comparisons without pooling identity or coverage. A multi-record query FASTA cannot silently become its first record in noninteractive comparison.
- Alignments use the displayed sequence scope. A chosen reference subset is distinct from the complete reference database.
- Displayed best identity belongs to the strongest-scoring HSP; maximum segment identity is a separate statistic. Query/reference coverage uses interval unions, counting overlaps once. No detected match is not proof of non-homology.
- Translated searches retain native nucleotide and amino-acid units. MSA columns, native coordinates and aligned amino-acid columns are not interchangeable.
- VSEARCH retains native global-search, dereplication, clustering and chimera results. It does not receive fabricated BLAST E-values or bitscores.
- Genetic codes and ORF starts are explicit biological choices. The simple ORF detector cannot infer context-dependent stops without CDS context.
- FastTree output is exploratory. Its drawing root is arbitrary and its local SH-like supports are not bootstrap values. Missing branch lengths produce a labeled cladogram.

## Calculation and storage behavior

- Opening, viewing or exporting an existing result does not silently rerun a scientific engine. Resuming runs only explicitly requested actions.
- Calculation preparation identifies inputs, active query, sequence scope and engine. Completed runs retain actual versions, arguments and settings separately from future settings.
- A failed rerun cannot replace a valid result with a partial file. Earlier interactive snapshots remain readable when new calculations or filters change the working state; each tree retains its source alignment link.
- Projects include immutable input copies and can resume after the original files move. New manifests verify saved input/state/result/log/history bytes and sizes; older valid projects remain readable. The checks detect accidental edits, not authenticated tamper resistance.
- Relative paths start from the caller's working directory, including when the launcher is invoked through a symbolic link.
- Batch outputs and aligned-region exports refuse existing destinations. Interactive HTML sharing requires explicit confirmation before replacing a report. Repeated changed Ctrl+S saves create new numbered sibling snapshots, while a clean save creates no duplicate. Saving is never a prerequisite for navigation, analysis or preview. Closing offers an explicit choice to retain or discard temporary work.
- HTML is a standalone offline snapshot with inline scripts, styles, data and figures. Sequence identifiers are escaped safely and remain intact.

## Runtime and distribution behavior

AlignView uses Python 3.10 or newer and has no mandatory third-party Python dependencies. Scientific executables are optional and required only for the selected calculation. The application and manual use English; user identifiers, source paths and raw engine output remain original data.

The source folder, installed wheel, zipapp and portable folder must provide the same runtime modules, release version and offline manual. `pyproject.toml` is the runtime inventory. Scientific engines, caches, historical prototypes and previous caller outputs are excluded from release payloads. The portable archive extracts under one `AlignView/` folder.

Every release set contains a wheel, source archive, versioned zipapp, portable ZIP and a checksum manifest whose paths stay within `dist/`. Builds normalize archive order, ownership, permissions and timestamps. Reproducibility comparisons require identical source, `SOURCE_DATE_EPOCH`, Python and build-tool versions.

## Release gates

| Gate | Evidence required |
| --- | --- |
| Scientific/application regressions | Full unit/integration suite on supported Python versions |
| Interactive workflows | Pseudo-terminal scenarios, small-screen navigation and persistent workspace state |
| Offline report | Browser checks for MSA navigation, matches, tree viewing, Unicode identifiers and downloads |
| Native integrations | Recorded representative real-engine runs, distinguished from mock executables |
| CLI/manual consistency | Generated reference check and local documentation-link checks |
| Artifact contents | Declared module inventory, version/metadata agreement, current-source comparison and checksums |
| Installed behavior | Actual pip-installed wheel, isolated zipapp and extracted portable launcher from a separate working directory |

`scripts/build_release.py` runs the local Python suite, CLI/manual checks and distribution checks, then creates the release files. The Python-version matrix, browser checks and real-engine checks are separate gates. A release can be prepared locally without publishing it. Target-cluster validation, representative difficult biological datasets and registry publication remain separate activities. See [Validation scope](AUDIT.md) for the limits of recorded tests.
