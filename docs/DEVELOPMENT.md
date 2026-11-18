# Development and verification

AlignView uses the Python standard library. Its architecture separates task entry, session state, workspace interaction, scientific operations, reporting and distribution. This separation keeps a UI change from changing an engine command implicitly.

## Runtime boundaries

| Area | Responsibility |
| --- | --- |
| `alignview.py` | Stable import and command-entry compatibility |
| `av_cli.py` | Argument parsing, explicit tasks and CLI validation |
| `av_session.py` | Loaded inputs, calculation state and result orchestration |
| `av_app.py` | Application startup, project switching and shutdown |
| `av_workspace.py`, `av_navigation.py` | Persistent workspace, tab state and contextual actions |
| `av_core.py`, `av_units.py` | Biological records, algorithms and explicit coordinate units |
| `av_tools.py`, `av_parameters.py` | Engine execution, validated settings and command provenance |
| `av_library.py` | Temporary dataset/result catalog with immutable snapshots and source links |
| `av_project.py` | Portable immutable project snapshots and integrity checks |
| `av_vsearch.py`, `av_genetics.py` | Native VSEARCH workflows and bundled NCBI codes |
| `av_report.py`, `av_figures.py` | Offline report composition and scientific figures |
| `av_ui.py`, `av_terminal.py`, `av_msa_tui.py` | Terminal rendering, keys and alignment viewing |
| `av_docs.py` | One offline English manual from source or packaged resources |

Older report-module names are compatibility boundaries, not additional rendering pipelines. The runtime module list in `pyproject.toml` is the packaging source of truth; adding a module requires adding it there. Release builders reject imports of undeclared `av_*` modules.

## Verify behavior

```bash
python3 -m unittest discover -s tests -v
python3 scripts/generate_cli_docs.py --check
python3 scripts/check_docs.py
python3 scripts/check_report_browser.py
python3 scripts/validate_real_tools.py --require mafft blastn diamond
python3 scripts/validate_vsearch.py
```

Tests cover coordinate/metric semantics, selected inputs, file exports, project integrity, query selection, pseudo-terminals, state invalidation, command fidelity and provenance. Mock executables validate orchestration rather than biological engine accuracy. Real-engine and browser checks are separate because they require external programs.

## Make a change

Write source identifiers, comments, presentation text and documentation directly in English. Preserve user identifiers, paths, result fields and raw engine logs. New scientific settings need typed validation, CLI/menu access, actual command construction, profile behavior and a successful-calculation snapshot. A scientific behavior change needs a known expected result, including coordinate units and selected sequence scope.

Keep ordinary work in the persistent workspace and expose detailed controls through contextual commands. Check small terminals, long lists and long alignments. Opening a result must not recalculate it or invent a save requirement. Calculation failure must preserve the previous valid result.

Regenerate CLI references when the parser changes. Update the single manual and examples with actual interaction names. Follow the [release procedure](RELEASING.md) to test an installed wheel, isolated zipapp and extracted portable folder. The [release contract](RELEASE_CONTRACT.md) defines the supported v1.0 behavior and limits.
