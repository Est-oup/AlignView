# User journey and distribution review

Review date: 18 November 2026. The review follows the work of a new user: open data, choose what to search against, inspect where and how sequences align, try another analysis, keep useful outputs and reopen them. It includes actual terminal/browser interactions and installed-distribution checks rather than relying only on menu text.

## Findings and implemented changes

| User task | Problem found | Resulting behavior |
| --- | --- | --- |
| Open sequence matches | The primary action did not clearly say that aligned regions could be inspected. | Enter explicitly offers View aligned regions. The inspector identifies query/target files and IDs, region ranges, direction and frames before displaying aligned letters. |
| Read native positions | The terminal called native positions local HSP coordinates, retained untranslated labels and advanced through dot gaps. | Positions are labeled 1-based inclusive in prepared ungapped inputs. Both gap symbols preserve sequence position; reverse orientation and translated codon steps are respected. |
| Page through a long local alignment | The report showed whole-segment bounds above a partial window, making it easy to assign the wrong positions to displayed residues. | Each page labels its own query/target bounds, while whole-segment bounds, direction and frames remain distinct. Gap-only windows have no invented residue range; inconsistent legacy coordinates are unavailable. |
| Reuse an interesting segment | Users had to retrieve a general HSP table rather than download the segment being read. | Download this segment writes both complete aligned FASTA rows with roles, IDs, ranges, units and frames. Query/target coverage tracks outline the selected segment. |
| Understand an empty search | No engine matches, all matches excluded by filters and stale global-search evidence could look alike. | Completed empty searches explain whether no segment was detected or filters removed all results. Older VSEARCH results do not replace a newer local result. |
| Export global matches | The direct TSV action required local HSPs even when a native VSEARCH search existed. | Global search exports retain native fields and empty tables. A current local search retains priority over older global results. |
| Run a command without interaction | `--no-menu` on a real terminal could still enter the alignment reader, a result pager or a save prompt. The report command was affected too. | Batch/HTML commands finish without those interactive surfaces. Existing batch destinations fail explicitly. Fresh temporary inputs are no longer described as a saved project. |
| Install and reproduce the application | The Bioconda template omitted noarch entry points and had only superficial CLI tests; the source archive omitted some reconstruction resources. | Recipe preparation verifies source bytes and SHA256, defines the noarch entry point and tests installed imports, offline guides, demo, saving, resuming and HTML. Source packaging includes reconstruction requirements and CI files. |

## Reproduce the acceptance scenarios

| Scenario | Evidence and command |
| --- | --- |
| Load without running a tool; cancel preparation; return from readers | `tests/test_file_menu.py`, `tests/test_navigation.py` and real pseudo-terminal journeys |
| Query selection, target changes, exact saved target and grouped results | `tests/test_compare_selection.py`, `tests/test_comparison_orientation.py` |
| Local search, open aligned regions, read inputs and return to Matches | `tests/test_user_journeys.py` |
| Gap/reverse/translated native ranges and true batch behavior on a terminal | `tests/test_user_contract.py` |
| Independent attempts, preserved descendants, renaming, source-independent reopening | `tests/test_exploration.py`, `tests/test_result_names.py`, `tests/test_project_portability.py` |
| Browser paging, residue inspection, coverage, selected segments, downloads and exact view exports | `python3 scripts/check_report_browser.py` |
| Real scientific engines, separately from mocks | `python3 scripts/validate_real_tools.py` and `python3 scripts/validate_vsearch.py` |
| Wheel, source archive, isolated zipapp and extracted portable application | `python3 scripts/build_release.py` and `python3 scripts/verify_release.py` |
| Recipe source identity, no fictitious hash, no overwrite and unsafe archive rejection | `tests/test_bioconda.py` and `python3 scripts/prepare_bioconda.py --help` |

The regular suite uses deterministic mock tools to isolate application behavior. Native integration scripts record the installed executable, version, settings and outcome. Record missing engines as missing rather than treating mock runs as biological validation. Release/recipe receipts belong beside the generated artifacts so their tested source archive is identifiable.

## Bioconda handoff

The local source version is 1.0.0. At review time the public upstream release was still v0.1.0, so that older public archive cannot be used as the source for this candidate. See the [preparation guide](../packaging/bioconda/README.md), [submission handoff](../packaging/bioconda/SUBMISSION.md) and [release procedure](RELEASING.md).

The generator distinguishes a buildable local-file candidate from a recipe whose public HTTPS archive has been downloaded and verified. A successful local build does not validate a future release URL. Publication of the reviewed source, maintainer agreement, final public-source regeneration and the official Bioconda pipeline remain separate steps. No upload, public PR or publication is performed by the preparation scripts.

## Practical limits

The application remains a foreground terminal workbench with a bounded offline report. Multi-query searches retain separate evidence, not one pooled many-to-many alignment. Local HSPs do not establish full-length homology; an explicit MSA or pairwise calculation answers a different question. Catalog snapshots copy working files and therefore need scratch space. These tests do not replace observing representative users, large real datasets, difficult biological cases or target-cluster and macOS validation.
