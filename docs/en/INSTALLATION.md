# Install and move AlignView

## Portable folder

Copy the complete `AlignView/` folder and run `/path/to/AlignView/alignview.sh --lang en --demo`. Python 3.10+ and Linux are the primary target. No third-party Python package is needed for the application. The launcher keeps the caller's current directory. If executable permissions were lost, use `bash /path/to/AlignView/alignview.sh` or restore the executable bit.

The portable ZIP extracts to one `AlignView/` folder. `python3 /path/to/alignview.pyz --lang en --demo` is the single-file alternative, including offline documentation. Neither distribution embeds scientific binaries.

For an installed command, create a Python environment and install `dist/alignview-0.1.0-py3-none-any.whl` with pip. Installing the source folder is also supported with the build dependencies available. The optional `scripts/install_user.sh` copies the zipapp into `~/.local/bin/alignview`; skip it if you prefer to retain only the portable folder.

## Scientific engines

| Calculation | Engine | Setup profile |
| --- | --- | --- |
| FASTA/MSA reading, internal analyses, SVG and HTML | None | None |
| NCBI local/translated searches | BLAST+ | blast |
| Protein and translated-query search | DIAMOND | diamond |
| Protein search alternative | MMseqs2 | mmseqs |
| MSA and fragment addition | MAFFT | mafft |
| Alternative MSA | MUSCLE 5 / Clustal Omega | muscle / clustalo |
| Global/local pairwise alignment | EMBOSS | emboss |
| Exploratory tree | FastTree | fasttree |
| Nucleotide global search, dereplication, clustering, chimeras | VSEARCH | vsearch |
| Optional external viewer | ALV | alv |

`minimal` includes BLAST, DIAMOND and MAFFT. `full` adds all other supported engines, including VSEARCH. The internal viewer does not need ALV. MUSCLE 3 is unsupported.

```bash
./alignview.sh --lang en --diagnose
./alignview.sh --lang en --setup vsearch
./alignview.sh --lang en --setup vsearch --install-tools
./alignview.sh --lang en --setup full --install-tools
```

Micromamba, Mamba or Conda must already be available. Without `--install-tools`, the command only shows a plan. The default dedicated prefix is `$XDG_DATA_HOME/alignview/tools` or `~/.local/share/alignview/tools`. It is detected without activating a shell environment.

## Keep engines under the application folder

From inside `AlignView/`:

```bash
./alignview.sh --setup full --tools-prefix "$PWD/tools" --install-tools
./alignview.sh --diagnose
```

The launcher detects `tools/conda-meta` on later runs. An explicit `ALIGNVIEW_TOOLS_PREFIX` or `--tools-prefix` takes priority. Direct Python, wheel and zipapp runs can select the same prefix explicitly. Conda binaries and scripts may contain installation-specific paths: copying their directory does not guarantee relocation to another machine/path. Check diagnostics and recreate that prefix when needed.

Existing tools can be selected using repeated `--tools-dir /path/bin` or `--executable mafft=/path/mafft`. Resolution order is explicit executable, tool directories, existing PATH, dedicated prefix, then on-demand HPC modules when enabled.

## HPC and platforms

Use `--module-policy auto` to try modules for missing programs or `never` to use existing paths and loaded modules only. Normal diagnostics do not load new modules. The tools menu offers an explicit module search. Set `TMPDIR` to suitable scratch space and choose the final output path explicitly. Standard FastTree remains single-threaded.

Slurm scripts are prepared but never submitted automatically. Check partition, resources and engine availability before sbatch. Transfer standalone HTML to your local machine for browser viewing; no temporary web server is required. Linux/SSH is the validated target; Windows requires a Linux environment such as WSL, and macOS/other targets need separate validation.
