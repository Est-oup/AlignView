# Validation scope and limitations

AlignView is tested for local-search coordinates, translated units, MSA handling, file exports, input integrity, portable projects, terminal keys and complete interactive workflows. Simulated executables test application behavior independently of installed tools. Real-engine checks are separate and record the versions and parameters they actually exercised.

## Evidence

| Record | Scope |
| --- | --- |
| [Benchmark results](benchmark-results.json) | Recorded local performance checks |
| [Real tools](real-tools-results.json) | Baseline integrations with installed scientific engines |
| [Modified settings](tuned-tools-results.json) | Engine commands with changed calculation parameters |
| [VSEARCH](vsearch-tools-results.json) | Native search, dereplication, clustering and chimera workflows |
| [Calculation review](CALCULATION_REVIEW.md) | Engine settings, command construction and scientific limits |
| [Interface review](UX_REVIEW.md) | Navigation, prerequisites, saving and result readability |

Historical logs describe their original runs. They are not evidence that every later code change, alternative translation table, parameter combination or remote CI job passed. Prepared workflow files do not imply publication or remote CI success.

## Scientific checks

Local HSP identity remains distinct from whole-sequence identity. Coverage is an interval union. Translated searches retain native nucleotide/amino-acid coordinate units. MSA p-distances exclude gaps per pair and can use different comparable positions for different pairs.

The tree SVG preserves Newick branches and labels, explicitly uses an arbitrary drawing root, and identifies a cladogram when distances cannot be drawn. FastTree local SH-like supports are not bootstrap values.

VSEARCH preserves native global metrics rather than assigning BLAST E-values or bitscores. NCBI tables are bundled with source/date/hash information. The simple ORF detector does not infer context-dependent stop codons without CDS context.

## Application checks

Guided comparisons require an active-query choice for multi-record FASTA input. Preparation screens show inputs and engines before launch. Reopening results does not rerun their engines. Saved projects include immutable input copies and preserve historical calculation settings. Resuming executes only explicit actions. English is used directly in source, menus, HTML and the bundled manual; input identifiers and raw tool logs remain unchanged.

## Limits

FASTA/MSA inputs remain in memory. Interactive import accepts FASTA/gzip FASTA, not every exported alignment format. A local search uses one active query rather than a scheduled all-versus-all analysis. HTML includes bounded previews; project and application exports retain full generated files. Trees, ORFs, motifs, clusters and chimera screening remain exploratory. Genetic codes require biological selection. Large difficult datasets, all engine releases, other platforms and the target HPC environment require their own validation.
