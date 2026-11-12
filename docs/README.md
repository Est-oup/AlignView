# AlignView manual

AlignView is a Linux/SSH terminal workbench for sequence exploration, similarity searches, alignments and shareable reports. Start with the result you need, choose the input, inspect the engine and parameters, run the calculation, then save or export.

## Find a guide

| Task | Guide | Offline command |
| --- | --- | --- |
| First run and first results | [Getting started](QUICKSTART.md) | `alignview --docs quickstart` |
| Every analysis and menu | [Usage](USAGE.md) | `alignview --docs usage` |
| Understand directories and saving | [Files](FILES.md) | `alignview --docs files` |
| Install or move the application | [Installation](INSTALLATION.md) | `alignview --docs installation` |
| Adjust scientific settings | [Parameters](PARAMETERS.md) | `alignview --docs parameters` |
| Use VSEARCH | [VSEARCH](VSEARCH.md) | `alignview --docs vsearch` |
| Choose translation tables | [Genetic codes](GENETIC_CODES.md) | `alignview --docs translation` |
| Read and share HTML and figures | [HTML report](HTML.md) | `alignview --docs html` |
| Reuse outputs | [Exports](EXPORTS.md) | `alignview --docs exports` |
| All available flags | [CLI reference](CLI.md) | `alignview --docs cli` |
| Resolve a problem | [Troubleshooting](TROUBLESHOOTING.md) | `alignview --docs troubleshooting` |
| Understand terminology | [Glossary](GLOSSARY.md) | `alignview --docs glossary` |

In commands, `alignview` means the installed command. Without installation, replace it with `/path/to/AlignView/alignview.sh`. The manual is bundled with the portable folder, wheel and zipapp. It is available before loading data, and under **Advanced → Help → Complete manual**.

## English interface and user data

All menus, prompts, error messages, reports and guides use English directly. There is no language selector or translation catalog. User-supplied sequence identifiers, input file names and raw engine output retain their original content. Saved projects preserve their scientific data and historical calculation settings.

See [Validation scope](AUDIT.md), [Calculation review](CALCULATION_REVIEW.md) and [Interface review](UX_REVIEW.md) for implementation choices and known limits.
