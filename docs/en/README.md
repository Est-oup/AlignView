# AlignView manual

AlignView is a Linux/SSH terminal workbench for sequence exploration, similarity searches, alignments and shareable reports. Start with the result you need, choose the input, inspect the engine and parameters, run the calculation, then save or export.

## Find a guide

| Task | Guide | Offline command |
| --- | --- | --- |
| First run and first results | [Getting started](QUICKSTART.md) | `alignview --lang en --docs quickstart` |
| Every analysis and menu | [Usage](USAGE.md) | `alignview --lang en --docs usage` |
| Understand directories and saving | [Files](FILES.md) | `alignview --lang en --docs files` |
| Install or move the application | [Installation](INSTALLATION.md) | `alignview --lang en --docs installation` |
| Adjust scientific settings | [Parameters](PARAMETERS.md) | `alignview --lang en --docs parameters` |
| Use VSEARCH | [VSEARCH](VSEARCH.md) | `alignview --lang en --docs vsearch` |
| Choose translation tables | [Genetic codes](GENETIC_CODES.md) | `alignview --lang en --docs translation` |
| Read and share HTML and figures | [HTML report](HTML.md) | `alignview --lang en --docs html` |
| Reuse outputs | [Exports](EXPORTS.md) | `alignview --lang en --docs exports` |
| All available flags | [CLI reference](CLI.md) | `alignview --lang en --docs cli` |
| Resolve a problem | [Troubleshooting](TROUBLESHOOTING.md) | `alignview --lang en --docs troubleshooting` |
| Understand terminology | [Glossary](GLOSSARY.md) | `alignview --lang en --docs glossary` |

In commands, `alignview` means the installed command. Without installation, replace it with `/path/to/AlignView/alignview.sh`. The manual is bundled with the portable folder, wheel and zipapp. It is available before loading data, and under **Settings and help → Complete manual**.

## Language and data

Use `--lang en`, select **Settings and help → Language / Langue**, or set `ALIGNVIEW_LANG=en`. French remains available with `--lang fr`. Menus, prompts, reports and guides use the selected language. Sequence identifiers, file names, command arguments and raw engine logs retain their original content. Saved projects retain their scientific settings; opening a project does not force its previous interface language.

A FASTA is not necessarily aligned. A local search is not an MSA. A project can be resumed; a standalone HTML is a presentation. Current settings describe the next calculation; historical settings describe existing results.

[Development](DEVELOPMENT.md) · [Releasing](RELEASING.md) · [Validation and limitations](AUDIT.md).
