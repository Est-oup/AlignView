# Contributing

Run `python3 -m unittest discover -s tests -v`. Ordinary tests do not download scientific engines. See [Development](docs/DEVELOPMENT.md) for real-engine validation and distribution builds.

For a bug report, provide a minimal command, AlignView and Python versions, terminal details, `alignview --diagnose-json` output and a small synthetic FASTA. State whether the run was local, over SSH or under Slurm.

Changes to coordinates or statistics need a known expected result. Verify navigation at 80×24 with long lists and alignments. Preserve batch operation and execution without mandatory Python dependencies.

Use English for source identifiers, comments, docstrings, messages, reports, documentation and templates. Write presentation text directly in English; there is no translation catalog or language selector. Preserve user-supplied sequence identifiers, source file names and raw engine output.

The offline manual is in `docs/`. After changing the CLI, run `python3 scripts/generate_cli_docs.py`, then `python3 scripts/check_docs.py`. Each scientific option must be connected to the actual engine command and calculation provenance. Rebuild distributions after changing source or bundled guides.
