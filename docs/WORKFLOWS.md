# File menu and exploration examples

Open a FASTA to see what was loaded: filename, sequence count, detected type and ungapped lengths. The file menu then explains five possible actions. Loading alone does not open a table/viewer, calculate an alignment or ask for a save directory.

## Choose from the file menu

| Choice | What happens |
| --- | --- |
| Look at the sequences | Choose a sequence browser or explicitly read a file you know is already aligned |
| Align sequences | Review the sequence scope, engine and settings, then Run or cancel |
| Find matching sequences | Choose query sequences, select the target FASTA and review the search direction |
| Other analyses | Choose an explained group, then the actual operation and sequence scope |
| Files and previous results | Switch/add files, reopen results, save/export, or open settings/help |

No scientific route is mandatory. A tree needs a suitable alignment; a reference search is not required before aligning a compatible sequence set. Requirements are explained on unavailable choices.

Use arrows and Enter, a displayed number followed by Enter, or a mouse click. Esc returns from a menu/reader. F2 returns to file actions. q closes from the file menu and returns from other readers. Pressing Esc at the initial file menu leaves it open; it does not start a task or discard work.

## Example: inspect first, then align

1. Open `alignview sequences.fna`. Read the summary; no calculation has run.
2. Choose Look at the sequences → Browse sequences. Read identifiers and lengths. Esc returns to the file menu.
3. Choose Align sequences. Review the inputs and engine. Esc cancels this preparation without running it; choosing Run starts the calculation.
4. The Alignment ready screen reports the completed alignment. Choose Open this result to enter the reader.
5. Esc/q returns to the result summary. Choose Back to file actions to return to the same loaded file.
6. Run another analysis whenever useful. Previous results retains each completed attempt, including an earlier alignment and its associated tree.

In the sequence browser, Space checks rows for alignment. Two or more checks select the exact scope; no checks means the complete file. In the full-screen alignment reader, n/N still navigates motif occurrences; these are reader controls rather than global steps.

## Example: reference evidence

Choose Find matching sequences. For a multi-record FASTA, a checkbox list selects the queries:

- Space or a row click checks/unchecks a sequence.
- a selects the complete FASTA; c clears the selection.
- / filters identifiers/descriptions. Previously checked sequences remain selected when hidden by a filter.
- Enter confirms a nonempty selection; Esc cancels without starting a search.

The **query** is what you want to find matches for. The **target dataset** is the local FASTA containing sequences to search against. These are roles in this comparison; the target does not have to be a curated reference collection.

Select the target from loaded files or another FASTA. Choosing the query file itself searches within that same dataset. The preparation screen keeps **Query/Queries** and **Target** on separate lines, including on small terminals: original filenames, selected-query count or identifier, and total target-sequence count. Run names the method, for example Run BLASTN search.

**Sequences to search** reopens the query checkbox list. **Target dataset** changes the target while keeping the query selection; canceling keeps the previous target. **Review search inputs** shows full paths, exact selected IDs and a preview of target IDs, lengths and types. Use it to distinguish files with the same name in different directories. No search runs until Run is selected. Selections for one comparison must use the same nucleotide/protein alphabet; compare mixed alphabets in separate runs. A one-record input needs no checkbox screen.

One query produces its ordinary Sequence matches ready result. Several queries produce Sequence comparisons ready: each selected record is searched independently against the same reference dataset. Open this result lists queries and their statuses; choose a completed query, then open its Matches view. Identity, coverage, thresholds and retained references belong to that query, and are never pooled across the selection. VSEARCH keeps its own global metrics.

A failed query does not remove successful ones. Interrupting a group keeps completed results and marks the current/remaining queries canceled/not run. Previous results shows the group rather than duplicating all its query entries. Saving/reopening preserves both the group and each query's native evidence. The group HTML is a summary; open an individual query to preview/export its complete report.

Within a query's Matches table, Space checks references and a prepares their alignment. The query and target filenames remain above the table. Opening an older result restores its own target, even after another comparison has used a different FASTA. Searches and alignments remain separate calculations. Later filters/checks do not remove completed attempts. Noninteractive CLI comparison still uses one explicit `--query-id ID`; the interactive checkbox flow does not silently change batch CLI behavior.

Choose a matching target and press Enter (**View aligned regions**) to see which regions align, in which direction, and how their letters differ. The inspector identifies both original input files and sequence IDs, lists native query/target ranges and frames, then shows each local alignment with identities, substitutions and gaps. Query and target positions are 1-based inclusive in the prepared ungapped inputs; translated searches distinguish nucleotide positions from aligned amino acids. The offline report also lets you select a segment, page through it and download its two aligned rows.

A completed search with no detected segments is different from a search whose segments were all excluded by filters. The Matches view explains which occurred. Relaxing display filters cannot recover a segment the engine never returned; a different search setting requires an explicit new Run.

## Other analyses and reused data

Other analyses offers small groups rather than the complete tool list. Translation and motif/region analyses ask for an individual sequence when the file has several. Selecting a group alone runs nothing. Tree/statistics choices explain that an existing alignment is required.

Translations, extracted regions, ORFs and alignments can be reused through Commands → Use result as sequence data. The derived file opens at the same explanatory file menu. Input links preserve the producing result. No export/reimport is required.

## Keep work only when needed

Completed results are retained automatically in temporary session storage. Previous results opens a named result summary without rerunning an engine. Adding a file also opens its summary/menu and leaves earlier files accessible.

Choose **Rename this result** on a result-ready screen (or in Commands) to give an attempt a recognizable name. Names appear in Previous results and are retained when saving the project. Renaming changes the display name only: identifiers, inputs, calculation settings and scientific files remain intact. Enter keeps the current name; q cancels.

Saving is optional during exploration. Save or export retains a project for later sessions or creates reusable files. A reader's w shortcut previews temporary HTML without a destination prompt. Durable HTML export is separate. Closing a file that has no analyses needs no save prompt; after unsaved calculations, closing offers Continue working, Close without saving, or Save a copy and close.

The offline HTML remains a snapshot with direct view navigation and browser Back/Forward controls. New calculations run in the terminal.

Reopening a project copies its current data and result catalog into session storage. Moving or removing the saved source folder afterwards does not prevent browsing earlier results or saving a new project copy. Changes still need an explicit save to survive closing the session.
