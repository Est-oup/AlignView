# Interface review

Reviewed on 12 November 2026. This describes workflow and implementation choices; it does not claim measured user satisfaction or comparative usability superiority.

## Problems addressed

| Problem | Current behavior |
| --- | --- |
| Results, calculations and diagnostics mixed in the same menus | Direct Home actions; additional work stays under Advanced |
| Opening an alignment required another result-selection step | The viewer opens immediately |
| Guided calculations returned Home before preparation | The selected task opens its preparation screen directly |
| Multi-query comparison could hide the active query | Explicit initial choice, visible query name and a query selector before calculation |
| Missing references led to a disabled search | Compare with references opens the reference picker |
| Segment diagnostics cluttered ordinary result reading | Detailed quality, coverage, conservation and distances are under Advanced |
| Changing a sensitivity preset altered the result limit | The reference limit remains independent, initially 100 |
| Simply reading an existing file requested saving | Read-only sessions close directly; new results/settings/history remain protected |
| Export choices exposed implementation profiles first | Choose FASTA alignment, reference TSV, tree SVG or a result folder |
| HTML repeated composition charts and empty search counters | Overview focuses on available results; detailed charts stay under Advanced |
| Runtime translation could leave mixed-language output | English is written directly in code, reports and the single bundled manual |

## Structure

Home exposes Results, Compare with references, Align sequences, Save and share, Files and query, and Advanced. Results follow available data. Preparation screens keep the launch action visible; engine-specific controls, command previews and profiles are in Advanced options.

Advanced groups additional analyses, checks and statistics, calculation settings, tools and provenance, HPC and help. Its implementation is separated into `av_advanced.py`; ordinary navigation remains in `av_navigation.py`.

Project, HTML and individual exports have distinct purposes: resume, share and reuse. Output paths remain relative to the launch directory and existing destinations are refused. Session history becomes portable when the project is saved.

The HTML report keeps reading and export views in its main navigation. Sequence diagnostics, dotplots, the complete atlas, parameters and method help are in Advanced. The overview retains a direct tree shortcut.

## Verification

Pseudo-terminal tests exercise the demo, viewer, F2 return, nested settings, file selection, reference addition, query changes, history, calculation, saving and closing. The six Home actions fit an 80×24 terminal. Browser checks exercise offline navigation, filtering, horizontal and vertical alignment navigation, expandable figures, printing, tree viewing, SVG download and VSEARCH.

These checks establish observable behavior. They do not replace feedback on real biological datasets, terminal configurations or the target HPC environment. See [Calculation review](CALCULATION_REVIEW.md) and [Validation scope](AUDIT.md).
