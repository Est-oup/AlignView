# Application files, outputs and saving

## One application directory

Copy the complete `AlignView/` directory to move the source application. It contains `alignview.sh`, `alignview.py`, the `av_*.py` modules, documentation, examples, scripts, tests and distributions. No Python source in the parent directory is required. `archives/`, when present locally, contains unused historical versions and is excluded from distributions. An optional parent launcher is only a symbolic link to the actual launcher.

The portable ZIP contains one top-level `AlignView/` folder. The `.pyz` is an alternative single-file application including documentation. Python and scientific engines remain external programs. An optional Conda prefix may be installed under `AlignView/tools/`; relocation of a Conda prefix is separate from copying Python sources.

## Relative paths

All relative input and destination paths start from the terminal's current directory (`pwd`). The launcher locates its sources without changing directory.

```bash
cd /work/my_analysis
/software/AlignView/alignview.sh sequences.fasta --lang en --msa view --report --save --out results/view_01 --no-menu
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

The options can be combined. `--html-out` implies HTML generation and needs no `--report`. Publication/complete exports also generate HTML when saving or exporting is requested. Existing destinations are refused in batch; the interactive HTML action asks before replacing an existing file.

## Project layout and resuming

A project includes `session.json`, input copies under `inputs/`, current parameters, historical calculation parameters, commands, versions and logs. Available alignments, tree, figures, HTML, VSEARCH results, exports and session history are retained too. Copy the entire directory, not only `session.json`. Input SHA256 hashes detect modifications; original source paths are no longer required.

```bash
alignview --resume my_project --lang en
alignview --resume my_project --lang en --html-out updated.html --no-menu
alignview --resume my_project --lang en --export-profile complete --export-dir exports_02 --no-menu
alignview --resume my_project --lang en --msa align --tree --report --save --out project_02 --no-menu
```

Resuming alone reruns no scientific engine. Only explicit actions such as `--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif` or `--vsearch` run again. Explicit parameter options prepare future settings; they do not rewrite historical results. Use the Data/Analyze menus to change the query or search method. Saving creates a new copy and does not modify the original project.

## Temporary files and recovery

Python's temporary directory, usually `/tmp`, holds calculations. On HPC, set `TMPDIR` to a suitable existing scratch directory before launching. Relative final outputs still use the launch directory, not `TMPDIR`.

After some interrupted or failed calculations, AlignView prints a recovery project path under `alignview_recovery_…/project`. Resume or copy it before system temporary-file cleanup. A forced termination or machine failure cannot guarantee recovery. Session history remains temporary until a project is saved.
