# Getting started

## Launch from your data directory

Keep the complete application in one `AlignView/` folder. Python 3.10+ is enough to use an existing alignment and the demo; no Python package installation or scientific engine is required.

```bash
cd /path/to/my_data
/path/to/AlignView/alignview.sh --demo
```

The launcher preserves your directory. Relative inputs and outputs start there. From inside the application folder, `./alignview.sh` also works.

## Load first, choose analyses later

Run `alignview sequences.fasta` to see a summary of what was loaded and a small menu of actions. Without arguments, choose **Load sequences** or **Open a saved project**. **d** opens the demo, **?** the manual and **t** the tools. Loading data does not start a calculation or ask for a project directory.

| Task | Direct command | What happens |
| --- | --- | --- |
| Inspect existing aligned sequences | `alignview view aligned.fasta` | Opens the alignment without an engine |
| Align raw sequences | `alignview align sequences.fasta` | Shows the calculation plan before running |
| Find matching sequences | `alignview compare queries.fna references.faa` | Chooses the active query, then prepares a search |
| Make one shareable report | `alignview report aligned.fasta` | Writes a standalone HTML file and exits |
| Continue a saved analysis | `alignview resume my_project` | Opens its workspace without replaying engines |
| Inspect installed engines | `alignview tools` | Prints availability, paths and versions |

Use `--html-out study.html` to choose a report name. `report` also accepts a saved project, rebuilding its HTML without a scientific calculation. Existing flag-style commands remain supported.

Equal-length records are not automatically treated as an alignment. Choose Look at the sequences → Read an existing alignment only when the FASTA is already aligned.

## Start from the file menu

The first screen names the file and shows the sequence count, detected alphabet and ungapped length range. No analysis has run. Five choices explain what you can do:

1. **Look at the sequences:** choose a sequence browser or explicitly read an existing alignment.
2. **Align sequences:** prepare an alignment of the loaded or checked records.
3. **Find matching sequences:** check one or several query records and select the target FASTA to search against. The Query/Queries and Target lines show the filenames and scope before Run. Target dataset changes the target; Review search inputs shows full paths and IDs. Space checks rows, a selects all, / finds a record and Enter confirms. Each query keeps independent evidence; multi-query results open as a query group.
4. **Other analyses:** a few groups for sequence checks, motifs/regions, translation/ORFs, alignment statistics/trees, pairwise comparison and nucleotide collections.
5. **Files and previous results:** choose/add files, reopen analyses, save/export or open settings/help.

Use Up/Down and Enter, a displayed number followed by Enter, or click an item. Esc cancels a preparation screen and returns to the file menu. F2 also returns from nested screens to the file menu. q closes from the file menu; it returns from a reader.

After Run, the result-ready screen identifies what was produced and offers **Open this result**, **Try different settings**, **Back to file actions**, **Previous results** and **Save or export**. An alignment reader does not open before you choose Open this result. Reading position is retained when you reopen it.

In the sequence browser, Space checks records and a prepares their alignment. No checks means all records. Targeted tools explicitly ask which sequence to use when several are loaded, or use the highlighted browser record. The command palette remains available with : or Ctrl+P. See [Exploration examples](WORKFLOWS.md).

## Keep what matters

Every completed interactive analysis stays in Previous results until closing, without an explicit save. Use Save / Export only to retain a resumable workspace or reusable files after closing. Ctrl+S saves a project snapshot: the first save asks for a directory, later changed saves use a new numbered sibling, and unchanged saves create no duplicate. Closing after unsaved analyses offers Continue working, Close without saving, or Save a copy and close. Calculations remain temporary until explicitly retained.

```bash
# Build an MSA without interactive menus and keep the project.
alignview align sequences.fasta --no-menu --save --out results/msa_01

# Search with one explicit query.
alignview compare queries.fna references.faa --query-id query_2 \
  --mode blastx --no-menu --save --out results/search_01

# Infer a tree from an existing alignment.
alignview view aligned.fasta --tree --no-menu --save --out results/tree_01

# Share an existing project without rerunning its engines.
alignview report results/tree_01 --html-out study.html
```

`--save` retains a complete project. `--out` alone does not save. `--html-out` retains one independent report. `--report` alone prepares temporary HTML. Existing destinations are refused; choose another name.

For installed engines, run `alignview tools`; `--setup PROFILE` shows an installation plan and `--install-tools` explicitly executes it. VSEARCH, alternative translation tables, diagnostics and HPC actions are available through commands. See [Usage](USAGE.md), [Files](FILES.md), [Installation](INSTALLATION.md) and [Troubleshooting](TROUBLESHOOTING.md).

## Try a small comparison

The bundled `examples/demo_query.fna` and `examples/demo_references.faa` form one coherent synthetic comparison: a coding nucleotide query, an exact protein, a related variant and an unrelated protein. They demonstrate the workflow and are not biological reference data. BLAST+ is required for this calculation.

```bash
# Run from the application folder, or use absolute example paths.
alignview compare examples/demo_query.fna examples/demo_references.faa \
  --mode blastx --query-id synthetic_query
```

Review the query and engine plan, run the search, inspect the strongest segment and coverage in Matches, select references with Space, then create the selected-reference alignment. With at least three aligned rows, create a tree. Save the project or share its HTML when ready. Each completed calculation remains independently accessible in Results; other sequence analyses can be run whenever useful.
