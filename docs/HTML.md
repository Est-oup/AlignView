# Read and share HTML reports


The main menu keeps Overview, Matches, Alignment, VSEARCH and Export when available. **Advanced** groups Sequences, Dotplot, Figures, Parameters and Understanding. **View tree** stays directly accessible from the overview even though the full atlas is in Advanced. Composition diagnostics no longer dominate the first view.

## Create a persistent report

Use **Save and share → Create an HTML report** and choose a destination. From the CLI:

```bash
alignview aligned.fasta --msa view --html-out report.html --no-menu
alignview sequences.fasta --msa align --report --save --out project_01 --no-menu
alignview --resume project_01 --html-out updated.html --no-menu
```

`--html-out` requests generation itself. `--report` alone prepares temporary HTML; add saving or a publication/complete export to retain it. Relative destinations start from the launch directory. Existing batch destinations are refused.

The file contains its data, scripts, styles and SVG figures. Open it in a browser with JavaScript enabled. It needs neither a server nor Internet access and remains available after AlignView closes. Share the one file with the intended recipients of its sequence data. It is an instant snapshot, not a project that can resume calculations.

## Available views

The overview offers entry points for the available results. Empty views are hidden. Use Overview and browser history to navigate.

| View | Contents |
| --- | --- |
| Overview | Available results, counts and search coverage |
| Matches | Reference search, sorting, HSP selection and local alignments |
| Advanced → Sequences | Lengths, nucleotide GC, ambiguity and profiles |
| Alignment | Paged MSA viewer and metrics |
| Advanced → Dotplot | Exact shared words, when calculated |
| Advanced → Figures | Expandable SVG atlas, text/category filtering and downloads |
| VSEARCH | Independent summaries, table previews and engine provenance, when calculated |
| Export | Browser-generated downloads of embedded data |
| Parameters and help | Current/historical parameters, actual commands and method explanations |

MSA horizontal paging preserves vertical position. Sequence pages contain at most 100 rows; filtering can reach records beyond the first page. The HTML embeds at most 1000 MSA sequences. Hiding invariant columns changes only the view.

## Tree figure

Infer a tree from at least three aligned sequences using FastTree, then rebuild the report. **View the tree** opens and expands its figure under Figures/Phylogeny. Scroll the frame for a large tree. **Download SVG** retains the complete vector figure, including parts outside the current frame.

Horizontal lengths and the scale represent substitutions per site. Vertical position organizes labels. The drawing root is arbitrary and does not imply an inferred ancestor or outgroup. Internal Newick labels are retained; FastTree SH-like supports are local, not tree bootstraps. Missing or all-zero lengths use an explicitly labelled cladogram without a distance scale.

A malformed Newick displays an unavailable-figure notice rather than reusing a stale SVG. Rebuilding an older project with a valid Newick adds the figure without rerunning FastTree. Project files are `tree_fasttree.nwk` and `figures/45_fasttree_tree.svg`.

## Figures and downloads

Figures are conditional: sequence composition/quality; HSP scores, spans and gaps; query/reference coverage; translated coordinate/frame checks; MSA conservation, entropy, occupancy, consensus and distances; and the inferred tree. Some distance heatmaps are limited to 45 sequences. `figures_catalog.tsv` documents the generated SVGs and units.

The browser's downloads use only embedded data: at most 2000 HSPs, 301 sequence-quality records, 1000 MSA sequences and 200 preview rows per VSEARCH analysis. Application exports retain their complete result files. VSEARCH identity/coverage are global-alignment metrics and are not BLAST HSP unions.

Print/PDF uses your browser and expands figures. An independent SVG is preferable for a very large tree. The report never runs BLAST, MAFFT, FastTree or VSEARCH. Return to the application for new calculations, and rebuild the HTML after changes.
