# AlignView

A Linux/SSH terminal application for inspecting sequences, comparing queries with references, building alignments and sharing standalone HTML/SVG results. Python 3.10 or newer is required; there are no mandatory third-party Python packages. The application, source code and documentation use English.

## Start

Keep the complete application in one `AlignView/` folder. Launch it from your data directory:

```bash
/path/to/AlignView/alignview.sh --demo
/path/to/AlignView/alignview.sh
/path/to/AlignView/alignview.sh --docs quickstart
```

The launcher preserves your working directory. Relative inputs and outputs start there. The demo opens an existing alignment and requires no scientific engine. Press `q` to reach Home.

## Choose a task

| Home action | Purpose |
| --- | --- |
| Results | Read available alignments, matches, trees and VSEARCH results |
| Compare with references | Choose the active query, inspect settings and run a search |
| Align sequences | Check the inputs and engine, then create an MSA |
| Save and share | Save a project, standalone HTML, FASTA, TSV, SVG or result folder |
| Files and query | Switch queries, inspect sequences, open files and restore history |
| Advanced | Additional analyses, detailed statistics, settings, tools, provenance, HPC and help |

Opening an alignment launches the viewer immediately. Guided calculations open a preparation screen before execution. Search results show local identity and query coverage; detailed segment diagnostics remain in Advanced. Arrow keys select, Enter opens, Esc goes back, F1 explains a choice and F2 returns Home. Merely reading a file does not require saving on exit.

![Home screen](docs/navigation-preview.svg)

## Save and share

| Need | Option |
| --- | --- |
| Resume calculations later | `--save --out project_01` |
| Keep one standalone report | `--html-out report.html` |
| Reuse exported tables and figures | `--export-profile publication --export-dir results_01` |

`--out` alone does not save. `--report` alone prepares temporary HTML. Existing output destinations are refused. The HTML contains its scripts, styles, data and SVG figures; it works offline without a server. A saved project includes immutable input copies, result files, settings and provenance.

```bash
# Read an existing alignment and keep its report.
./alignview.sh aligned.fasta --msa view --html-out report.html --no-menu

# Align sequences and retain a portable project.
./alignview.sh sequences.fasta --msa align --report --save --out project_01 --no-menu

# Choose one query from a multi-record FASTA.
./alignview.sh queries.fna references.faa --mode blastx --query-id query_2 \
  --save --out comparison_01 --no-menu

# Run a native VSEARCH global search.
./alignview.sh queries.fna references.fna --vsearch search --query-id query_2 \
  --report --save --out vsearch_01 --no-menu

# Rebuild HTML from a project without rerunning its engines.
./alignview.sh --resume project_01 --html-out updated.html --no-menu
```

## Scientific tools

Supported tools include BLAST, DIAMOND, MMseqs2, MAFFT, MUSCLE 5, Clustal Omega, EMBOSS Needle/Water, FastTree and VSEARCH. They are optional and needed only for their respective calculations. Engine names, versions and paths are visible before calculation; actual commands and historical settings remain in provenance.

```bash
./alignview.sh --diagnose-json
./alignview.sh --setup mafft
./alignview.sh --setup mafft --install-tools
```

Tool installations can live in the application’s `tools/` environment. HPC modules, executable directories and explicit executable paths are also supported. See [Installation](docs/INSTALLATION.md).

VSEARCH provides nucleotide global search, dereplication, clustering and de novo chimera screening. Its native results remain distinct from BLAST HSPs and do not have invented E-values or bitscores. [Translation tables](docs/GENETIC_CODES.md) include 27 offline NCBI codes with separate query/reference settings and configurable ORF starts.

FastTree requires at least three aligned sequences. The terminal displays its topology; **View the tree** in HTML opens the SVG with identifiers, branch lengths and local supports. The drawing root is arbitrary, and local SH-like support is not bootstrap support.

## Documentation and development

[Manual index](docs/README.md) · [Getting started](docs/QUICKSTART.md) · [Usage](docs/USAGE.md) · [Files](docs/FILES.md) · [Parameters](docs/PARAMETERS.md) · [HTML](docs/HTML.md) · [VSEARCH](docs/VSEARCH.md) · [Exports](docs/EXPORTS.md) · [CLI](docs/CLI.md) · [Troubleshooting](docs/TROUBLESHOOTING.md).

The manual is bundled with source, wheel, portable ZIP and zipapp distributions. Run `./alignview.sh --docs` before loading any data, or choose **Advanced → Help → Complete manual** in an open project.

```bash
python3 -m unittest discover -s tests -v
python3 scripts/generate_cli_docs.py --check
python3 scripts/check_docs.py
```

See [Contributing](CONTRIBUTING.md), [Development](docs/DEVELOPMENT.md), [Calculation review](docs/CALCULATION_REVIEW.md) and [Validation scope](docs/AUDIT.md). Local engine tests demonstrate integration behavior, not a ranking of biological accuracy. Large FASTA inputs remain in memory; the report is an offline snapshot rather than a calculation server.
