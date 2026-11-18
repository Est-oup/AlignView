# Interface review

The current interaction starts with an explanatory file menu and independent completed results. It replaces both the goal-first route and the automatic sequence-table arrival. This describes implemented behavior, not measured user satisfaction.

## Core behavior

Opening a bare FASTA shows its filename, count, type and ungapped lengths, with five explained actions. It opens no table or reader automatically. Other analyses uses small groups, and sequence operations explicitly select their target. Calculation preparation remains separate from Run.

Completed calculations return to a result-ready summary. Open this result is an explicit reader action; Back to file actions restores the recognizable menu. Files and previous results provides file switching, completed attempts, optional exports and settings. The visible menus work with arrows/Enter, numbers or clicks; shortcuts are optional.

Reference comparison adds a checkbox query selector with all/clear/filter controls. A single Run can compare several selected records with the same references, grouped by query. A failure or interruption preserves completed members. The group and per-query evidence survive project saving; the group HTML explicitly describes a summary rather than pooled evidence.

Find matching sequences names the two input roles explicitly: queries are searched against a local target FASTA. Separate query/target lines remain visible in compact preparation menus, result summaries and Matches. Preparation offers a target change without losing query selection and a full-path/identifier review. Offline reports retain the same input direction above their result views.

Every completed interactive operation retains a separate temporary snapshot. Opening an older attempt does not remove newer ones. Tree input links point to the alignment used. Translations, ORFs, extracted regions and alignments can be reused as sequence data without export/reimport.

Result-ready screens offer Rename this result. Display names persist with the project while stable result identifiers and scientific snapshots retain their original content. Reopened projects detach both current inputs and earlier results from the saved source folder.

Esc returns to the previous screen, with per-object list and alignment-reader state retained. The catalog header and rows also accept mouse clicks. Nested calculation settings retain F2 access to the workspace. The result summary offers Try different settings for supported result kinds, preserving the earlier attempt.

## Saving and browser reading

Project saving is optional until the user wants to keep work after closing. A temporary HTML preview opens without asking for a destination. Durable HTML export is separate. The report uses direct view navigation and Back/Forward controls; it does not prescribe a next step or pretend to run scientific engines.

## Verification and limits

File-menu acceptance checks exercise loading without automatic output, explained choices on small screens, canceled preparation, explicit target selection and explicit reader opening. Exploration tests also exercise independent alignments, tree lineage, selected-record scope, failed/canceled preparation, filter changes, save/reopen, navigation/reader state, highlighted-sequence actions and preview without destination prompts. PTY checks exercise the actual terminal entry point; browser checks exercise history and linked scientific views.

Catalog snapshots currently copy complete working states, which uses additional scratch storage. Completed operations are retained; failed attempts are reported by the existing error UI. Native calculations remain foreground operations. Archive/Undo, a live browser computation service and concurrent job scheduling are not implemented. Observe representative users and large datasets before claiming usability or performance improvements.
