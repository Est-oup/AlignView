# Getting started

## Launch from your working directory

The complete application lives inside **one `AlignView/` folder**. Python 3.10+ is enough to explore the demo; no scientific engine or Python package installation is required.

```bash
cd /path/to/my_data
/path/to/AlignView/alignview.sh --demo
/path/to/AlignView/alignview.sh
```

Replace the example paths. The launcher keeps your current directory. Relative inputs and outputs belong to that directory, rather than the application directory. From inside the application folder, `./alignview.sh` also works.

## Choose your starting point

| Start option | Input | Next action |
| --- | --- | --- |
| Open an alignment | An aligned FASTA | Results → Read the alignment |
| Align sequences | A FASTA containing at least two sequences | Check the preparation screen, then Run |
| Compare against references | Query FASTA followed by reference FASTA | Select the active query, check the preparation screen, then Run |
| Resume a project | A directory containing `session.json` | Results |
| Try an example | No file needed | Discover the viewer without engines |

A multi-sequence query file opens a searchable query selector before comparison. In batch, select its identifier explicitly with `--query-id`. The identifier is the first word after `>` in the FASTA header.

The six Home actions are **Results**, **Compare with references**, **Align sequences**, **Save and share**, **Files and query**, and **Advanced**. Unavailable actions explain their prerequisites. Esc/Left goes back, F2 returns Home, F1 describes a selected option, and `/` searches long lists.

## Keep the results

Calculations use a temporary session. Use **Save the project** to resume later, **Create an HTML report** to share a standalone presentation, or **Export results** to work in other software. Relative destinations start from your launch directory. Existing destinations are refused in batch.

```bash
# View an existing alignment and retain only HTML.
alignview aligned.fasta --msa view --html-out report.html --no-menu

# Align sequences and save a portable project with HTML.
alignview sequences.fasta --msa align --report --save --out results/msa_01 --no-menu

# Compare a nucleotide query against proteins.
alignview query.fna references.faa --mode blastx --msa hits --report --save --out results/search_01 --no-menu

# Infer a tree from an existing alignment.
alignview aligned.fasta --msa view --tree --report --save --out results/tree_01 --no-menu

# Rebuild HTML without rerunning scientific calculations.
alignview --resume results/tree_01 --html-out updated.html --no-menu
```

These commands use the installed `alignview` name; replace it with the folder launcher's absolute path if needed. `--save` retains a project. `--out` alone does not save. `--html-out` retains one HTML independently of project saving. `--report` alone prepares temporary HTML.

For VSEARCH, choose **Advanced → Additional analyses → VSEARCH analyses**. For alternative genetic codes, choose **Advanced → Additional analyses → Sequences, motifs and translation → Translation and coding regions → Choose the translation table**. Read [Usage](USAGE.md), [Files](FILES.md) and [Troubleshooting](TROUBLESHOOTING.md) before larger runs.

## Bundled VSEARCH examples

From inside the application folder: `./alignview.sh examples/example_vsearch.fna --vsearch cluster --html-out vsearch_demo.html --no-menu`. To choose a query before comparison: `./alignview.sh examples/example_vsearch.fna examples/example_vsearch_refs.fna --vsearch search`. These synthetic records demonstrate the workflow and are not biological references.

Detailed statistics and segment quality are under **Advanced → Checks and statistics**. Extra analyses, scientific settings, tools, HPC, help also live under **Advanced**. The demo opens directly in its viewer; press q to reach Home. Reading an existing alignment does not ask you to save on exit.
