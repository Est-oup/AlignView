# Read and share the offline report

The HTML report is the browser companion to the terminal workspace. It opens the most useful available result directly: alignment, matches, tree or native VSEARCH results. The report is a snapshot; scientific calculations still run explicitly in AlignView.

## Create one durable file

```bash
alignview report aligned.fasta --html-out study.html
alignview report my_project --html-out updated.html
```

The first command explicitly treats the FASTA as aligned. The second reads an existing project and rebuilds its report without replaying scientific engines. Legacy commands such as `alignview aligned.fasta --msa view --html-out study.html --no-menu` remain supported.

Press w in the terminal for an immediate temporary preview without choosing a destination. Save / Export → Export HTML report creates a durable copy. Previewing is independent of project saving. `--html-out` requests generation itself; `--report` alone prepares temporary HTML and needs saving or exporting to keep it. Relative destinations start in the launch directory. Batch report commands refuse an existing name; interactive HTML sharing asks explicitly before replacing it.

The file contains its data, JavaScript, styles and SVG figures. Open it in a browser with JavaScript enabled; it needs no server, CDN or Internet connection. Share the complete file with the intended recipients of its sequence data. A report cannot resume calculations like a project can.

## Main views

| View | Purpose |
| --- | --- |
| Alignment | Navigate rows/columns, choose a trace sequence, inspect a residue and export the displayed data |
| Matches | Filter/sort reference evidence, inspect covered positions and HSPs, select references |
| Tree | Read the complete inferred SVG and download it |
| VSEARCH | Inspect independent native summaries, bounded previews and executed settings |
| Export | Download exactly the included or selected data |
| Overview | Check available results and report scope |

Unavailable views are hidden. Browser history and direct view links retain navigation. Sequence data, figures, dotplots and Methods & tools are supplementary views, keeping detailed diagnostics out of ordinary result reading.

## Alignment and reference evidence

Comparison reports keep a **query → target** input panel above the views. It names the query identifier and source FASTA, the target FASTA, target-sequence count and method. Expand Input file paths to distinguish identically named files from different directories. Group summaries name the selected query IDs and the same target dataset. These labels describe the saved comparison, including when an older result is reopened; absent names in legacy data are identified as not recorded.

The alignment viewer keeps horizontal paging independent of vertical position. Search/filter records, choose a trace/reference row and show differences or variable sites. Click a residue to see both its MSA column and ungapped sequence position. These coordinates are explicitly different; hiding a column does not modify the stored alignment.

Exact motif search uses the selected trace (or consensus) and ignores alignment gaps (`-` and `.`) by default. A motif can therefore span a gap, and the match status gives both its ungapped residue range and its original MSA column range. Unknown symbols such as `?` still occupy a residue position; ambiguity codes are matched literally. Turn off **Ignore alignment gaps** to search the aligned text instead. Next match wraps to the first occurrence and explains the wrap. A match hidden by variable-column filtering reveals the complete columns explicitly.

Matches retains the full reference-level evidence table while HSP previews can be bounded. Filters and sorting change the displayed table, not the original calculation. Selected references can constrain the alignment view to those identifiers when the alignment contains them. A local HSP identity is not whole-sequence identity, and coverage counts the union of native intervals.

Select a target and segment to inspect its aligned letters. The whole-segment description gives query/target bounds, orientation and reading frames. Each displayed 140-column page also gives the native query and target positions of those particular letters: gaps do not advance a sequence, reverse strands decrease positions, and translated residues advance three nucleotides on the translated side. Gap-only pages have no residue interval. Inconsistent legacy coordinates are labeled unavailable rather than guessed. The selected segment is outlined on both coverage tracks.

**Download this segment** exports both complete aligned rows, including portions outside the displayed page. FASTA headers retain query/target roles, original identifiers, native ranges, units and frames. It is a local segment alignment, not a full-length alignment or a pooled alignment of several queries.

Browser exports distinguish the current reference view, selected references, filtered full alignment and currently displayed alignment rows and columns. The view export uses the current sequence page; the filtered full-alignment export keeps every filtered row. A column-map TSV records displayed MSA coordinates, especially useful when variable-site filtering removes intervening columns. Filtering and export do not run an engine.

## Tree figure

Create a tree from at least three aligned sequences with FastTree, then rebuild the report. Tree opens the SVG directly. Download SVG retains the complete vector figure, including content outside the visible scroll frame. Figures also includes the tree under Phylogeny.

Horizontal lengths and the scale represent substitutions per site where available. Vertical layout organizes labels. The drawing root is arbitrary rather than an inferred ancestor or outgroup. Internal Newick labels are retained; FastTree SH-like supports are local values, not bootstrap support. Missing or all-zero branch lengths produce a labeled cladogram without a distance scale.

A malformed Newick produces an unavailable-figure notice rather than a stale tree. Rebuilding an older valid project generates the figure without rerunning FastTree. Project files are `tree_fasttree.nwk` and `figures/45_fasttree_tree.svg`.

## Scope, figures and reproducibility

The report explicitly records how much data it embeds. Current bounds are 1000 MSA rows, 2000 HSPs, 5000 sequence-statistic rows and 200 preview rows per VSEARCH analysis. Full MSA metrics and reference-level coverage are computed from the original supplied results. Browser downloads can only export included data; complete generated native files remain in the project/application exports.

Figures are composed once with a catalog of descriptions and units. They cover input composition, local-segment evidence, coverage, translated-coordinate diagnostics, MSA conservation/entropy/occupancy/distances and the tree when applicable. Open a figure only when needed and download its SVG for reuse.

The automatic p-distance figure requires 2–45 aligned sequences and a calculation budget of at most 20,000,000 (`sequence count² × alignment columns`). Larger comparisons omit that figure; the application still offers an explicit distance calculation for up to 250 sequences. These uncorrected distances exclude gaps separately for each pair and are not evolutionary model estimates.

Methods & tools separates future settings from historical calculation settings and records actual engine versions, arguments, input hashes and report limits. A saved report does not update when its source project changes; rebuild it after a new calculation. See [Files](FILES.md), [Exports](EXPORTS.md) and [Validation scope](AUDIT.md).

## Free navigation

The sidebar opens any included view directly. Back and Forward return through visited views; browser history also works. Reading positions, reference selection and filters remain in place. No next step, workflow map or save requirement is imposed. From Matches, Compare selection in alignment opens the included alignment for those identifiers; this changes only the view.

Each exported report describes the currently opened result context. To inspect another historical attempt, open it in the terminal Results list and preview or export its report. The HTML cannot launch scientific engines.
