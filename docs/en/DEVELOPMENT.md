# Development and verification

The application uses the Python standard library. alignview.py handles CLI/session orchestration; av_navigation.py handles task-based menus; av_core.py handles sequence/HSP models and calculations; av_tools.py handles external processes and provenance. av_parameters.py defines typed settings shared by CLI, menus, profiles and command builders. av_project.py handles portable project snapshots.

av_vsearch.py builds independent VSEARCH workflows and reports their native metrics. av_genetics.py vendors NCBI genetic codes. av_tree.py parses Newick and draws SVG without a dependency. av_ui.py, av_terminal.py and av_msa_tui.py implement terminal interaction. av_i18n.py and av_locale_en.py translate presentation while preserving sequence identifiers and command data. av_docs.py reads offline Markdown from sources or packaged resources.

```bash
python3 -m unittest discover -s tests -v
python3 scripts/generate_cli_docs.py --check
python3 scripts/check_docs.py
python3 scripts/check_report_browser.py
python3 scripts/validate_real_tools.py --require mafft blastn diamond
```

Tests cover biological coordinates, exports, query selection, real pseudo-terminals, resuming, storage safety, command construction, genetic-code consistency, VSEARCH outputs and localization. Mock executables validate orchestration, not biological engine accuracy. Optional real-tool/browser checks are separate and need their external dependencies.

To add a parameter, connect typed validation, relevant menu selection, CLI help, command construction, profiles and successful-calculation provenance. Regenerate the CLI references after changing the parser. Maintain both language guides and the presentation catalog. Preserve raw identifiers, arguments, result fields and engine logs.

Build distributions using the [release procedure](RELEASING.md). The wheel and zipapp must contain both language manuals and all runtime modules. The portable archive must extract under one AlignView/ directory and run from a different current directory. Archives of old prototypes, tools prefixes, generated scientific results and development caches are excluded from that archive.
