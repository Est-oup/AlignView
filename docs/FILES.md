# Application files, outputs and saving

## One application directory

Copy the complete `AlignView/` directory to move the source application. It contains `alignview.sh`, `alignview.py`, the `av_*.py` modules, documentation, examples, scripts, tests and distributions. No Python source in the parent directory is required. `archives/`, when present locally, contains unused historical versions and is excluded from distributions. An optional parent launcher is only a symbolic link to the actual launcher.

The portable ZIP contains one top-level `AlignView/` folder. The `.pyz` is an alternative single-file application including documentation. Python and scientific engines remain external programs. An optional Conda prefix may be installed under `AlignView/tools/`; relocation of a Conda prefix is separate from copying Python sources.

## Relative paths

All relative input and destination paths start from the terminal's current directory (`pwd`). The launcher locates its sources without changing directory.

```bash
cd /work/my_analysis
/software/AlignView/alignview.sh sequences.fasta --msa view --report --save --out results/view_01 --no-menu
```

This reads `/work/my_analysis/sequences.fasta` and saves `/work/my_analysis/results/view_01/`. The report is `report.html` inside that project. It does not write these results into `/software/AlignView`. If you launch from inside the application folder, relative outputs naturally go there. Use absolute paths to choose a different destination, and quote shell paths containing spaces.

## What is retained?

| Option/action | Durable output | Resumable project? |
| --- | --- | --- |
| Exploration without saving | None after closing | No |
| `--report` alone | Temporary HTML only | No |
| `--save` | A dated `alignview_YYYYMMDD_HHMMSS/` directory in the current directory | Yes |
| `--save --out directory` | Project at that destination | Yes |
| `--out directory` alone | No automatic save | No |
| `--html-out report.html` | One standalone HTML | No, presentation only |
| `--export-dir directory` | Files from the selected export profile | No, external reuse |

The options can be combined. Ctrl+S is the project snapshot shortcut: the first save asks for a destination, later modified saves use a fresh `<name>_snapshot_NNN` sibling, and a clean save creates no new copy. Earlier saved snapshots are not overwritten. `--html-out` implies HTML generation and needs no `--report`. Publication/complete exports also generate HTML when saving or exporting is requested. Existing destinations are refused in batch; choose another name for an existing destination.

## Project layout and resuming

`--no-menu` is noninteractive even when launched from a real terminal: alignment output is static, pagers are disabled and no save question is asked. Add `--save` or an explicit export destination to keep outputs. Existing batch destinations produce an error rather than an interactive rename question. The `report` command also writes its HTML and exits without opening an alignment reader.

A project includes `session.json`, input copies under `inputs/`, current parameters, historical calculation parameters, commands, versions and logs. Available alignments, tree, figures, HTML, VSEARCH results, exports and session history are retained too. Copy the entire directory, not only `session.json`. New snapshots also contain `project_manifest.json`, with file sizes and SHA256 hashes for inputs, state, results, logs and nested history. Resuming verifies recorded files and rejects accidental edits. This is integrity checking, not authenticated tamper resistance: someone able to replace both files and checksums can forge a new inventory. Older valid schema-1 projects without this manifest remain readable. Original source paths are no longer required.

```bash
alignview --resume my_project
alignview --resume my_project --html-out updated.html --no-menu
alignview --resume my_project --export-profile complete --export-dir exports_02 --no-menu
alignview --resume my_project --msa align --tree --report --save --out project_02 --no-menu
```

Resuming alone reruns no scientific engine. Only explicit actions such as `--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif` or `--vsearch` run again. Explicit parameter options prepare future settings; they do not rewrite historical results. Use Data / Choose active query and the comparison preparation command to change the query or search method. Saving creates a new copy and does not modify the original project.

## Interactive result catalog

The running application keeps an independent temporary snapshot after each completed analysis. Data lists imported sequence sets; Results lists attempts with input links. Reopening a result copies its stored state into the working session, so subsequent calculations cannot overwrite that entry. No destination prompt is needed. Batch commands keep their existing explicit saving behavior.

Project schema 1 gains optional `workspace_entry`, `result_sources`, `sequence_selection`, `navigation_state`, `comparison_selection` and `comparison_batch` fields. Catalog entries are complete validated snapshots under `history/`; individual result files and original inputs remain portable. The implementation currently copies snapshots, so repeated large analyses use additional scratch space. It does not provide permanent autosave or simultaneous background calculations.

A multi-query comparison stores its selected identifiers and run statuses in `comparison_batch`. Each completed query has an independent catalog snapshot and exact scientific inputs/results. The group refers to its members by stable IDs; no temporary directory paths are needed to find them after saving. Group attempt logs are retained under `attempts/`. Group HTML exports summarize statuses; complete evidence can be exported from each query or retained in the full project.

## Temporary files and recovery

Python's temporary directory, usually `/tmp`, holds calculations. On HPC, set `TMPDIR` to a suitable existing scratch directory before launching. Relative final outputs still use the launch directory, not `TMPDIR`.

After some interrupted or failed calculations, AlignView prints a recovery project path under `alignview_recovery_…/project`. Resume or copy it before system temporary-file cleanup. A forced termination or machine failure cannot guarantee recovery. The complete catalog remains temporary until a project is saved. HTML preview files also disappear when their owning session closes; export a report explicitly to keep it.
