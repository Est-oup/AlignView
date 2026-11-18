# Design review for AlignView 1.0

This records the earlier v1.0 design review. The current interaction replaces its four-tab navigation with Data / Results, independent analysis snapshots and Back navigation. See [Free exploration](WORKFLOWS.md) and [Interface review](UX_REVIEW.md) for implemented behavior; the source comparisons below are historical design context.

The redesign starts from the research task: choose a query, examine reference evidence, align an intentional sequence set, inspect positions and an exploratory tree, then retain a reproducible result. This review compares documented interaction patterns, not software quality rankings or measured biological accuracy.

## Patterns reviewed in primary documentation

| Source pattern | Decision in AlignView | Deliberate scope |
| --- | --- | --- |
| AliView documents direct alignment navigation, differences against a chosen trace sequence, selected-region FASTA export and configurable external aligners. [AliView manual](https://www.ormbunkar.se/aliview/) | Keep position reading central; add a trace row, explicit selection/export and inspectable engine preparation. | A terminal analysis workspace, without claiming AliView's complete manual-editing or indexed-file feature set. |
| Jalview's overview provides a navigable representation of the alignment and its visible region. [Overview Window](https://www.jalview.org/help/html/features/overview.html) | Preserve an alignment overview and cursor location while moving through rows and columns. | Navigation changes the view, not the underlying sequence data or calculation state. |
| Jalview reference views distinguish the selected sequence's residue positions from alignment columns and highlight differences against that reference. [Reference Sequences](https://www.jalview.org/help/html/calculations/referenceseq.html) | Make trace selection explicit and expose MSA/ungapped coordinates when inspecting a residue. | Gap-containing alignment columns are never relabeled as native nucleotide or amino-acid positions. |
| Jalview documents independent views over shared alignment data and links tree-leaf selection back to alignment sequences. [Multiple Views](https://www.jalview.org/help/html/features/multipleViews.html), [Tree Viewer](https://www.jalview.org/help/html/calculations/treeviewer.html) | Use one persistent project state with Alignment, Matches, Tree and Data views; link reference and tree selections to alignment inspection. | No view switch silently recalculates an engine result or changes statistical meaning. |
| NCBI defines separate tabular fields for local aligned characters, identity, scores, frames, per-HSP coverage and per-subject coverage. [BLAST+ output fields](https://www.ncbi.nlm.nih.gov/sites/books/NBK279684/) | Show reference-level coverage beside a clearly identified strongest HSP, and retain raw segments and native units for inspection. | A high-identity short HSP is not whole-sequence identity; overlapping intervals are not summed twice. |

These manuals establish that the documented patterns exist. They do not establish that AlignView matches another program's speed, memory behavior, accessibility or breadth. No comparative user study or biological accuracy benchmark was performed for this redesign.

## Why the structure changed

A menu hierarchy organizes features, but a researcher usually keeps several related questions active: which reference is relevant, which positions differ, whether the chosen sequence set aligns sensibly, and whether the tree is consistent with that alignment. Re-entering separate Explorer, Analyze and Data menus interrupts that chain.

The v1.0 workspace therefore holds the dataset and view state together. The active tab supplies its next useful action. A command palette makes less frequent capabilities searchable rather than requiring another navigation taxonomy. Reference checkboxes express an analysis scope directly; they do not pretend to rerank or delete the original search.

Trace selection makes the visual comparison answerable: differences are relative to a named sequence or an explicitly labeled consensus. Region exports include original MSA coordinates and selected identifiers, so an extracted view can be understood outside the application. A hidden column is not an invented contiguous native interval.

Save/share separates three intentions: continue calculations, share an offline exploration, or reuse specific files. Executed settings and commands remain distinct from settings prepared for a future run. Immutable source copies, validated result replacement and a saved-file manifest support reproducibility without presenting ordinary checksums as authenticated proof.

## Practical tradeoffs

The target is Linux/SSH work with standard Python, optional native engines and an offline browser report. This favors keyboard operation, explicit calculation plans, readable provenance and portable projects. Mature desktop editors remain useful when the task needs extensive manual alignment editing, richer GUI annotation/structure integration or different import formats. AlignView does not claim feature parity with them.

The standalone report is a bounded snapshot, not a live analysis service. The terminal/project retain complete generated data; browser limits are explicit. Detailed diagnostics remain available without becoming the first screen. The implementation preserves validated scientific helpers where appropriate while replacing the navigation model, report composition and release boundaries.

See [Usage](USAGE.md), [Calculation review](CALCULATION_REVIEW.md), [Interface review](UX_REVIEW.md) and the [release contract](RELEASE_CONTRACT.md) for the resulting behavior and verification scope.
