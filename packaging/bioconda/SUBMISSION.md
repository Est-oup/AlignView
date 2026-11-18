# Submission handoff

This file is a local handoff for the prospective recipe maintainer. It is not evidence that publication or acceptance has occurred.

## Remaining release requirements

- Review the application changes and choose/finalize the release version. The current source declares `1.0.0`.
- Publish the exact reviewed source at a stable versioned HTTPS archive URL. As checked on 18 November 2026, the existing upstream public release was `v0.1.0`, which predates the current source.
- Regenerate the public-source recipe using the actual published archive SHA256. The generator refuses stale runtime/manual/metadata and mismatched downloaded bytes.
- Confirm who agrees to maintain the recipe. The upstream owner `Est-oup` is a candidate, not an automatically assigned volunteer.
- Run the official Bioconda lint/build and clean-container tests for that exact public-source recipe. A successful local-file build does not verify the future public release URL.
- Verify naming/channel availability again immediately before opening the PR.

## Suggested PR title

```text
Add alignview 1.0.0
```

## Suggested PR description

```text
Adds AlignView, a nucleotide/protein sequence evidence and alignment workspace
with a terminal interface and standalone offline HTML reports.

The pure-Python package uses a noarch recipe with an explicit console entry point.
Scientific engines remain optional so existing alignment inspection and reporting
work without installing every supported analysis backend. Users can add BLAST,
DIAMOND, MAFFT and other engines for the calculations they need.

Source: <verified public versioned archive URL>
SHA256: <actual archive SHA256>
Maintainer: <consenting GitHub handle>

Validation actually completed:
- <conda-build version, platform and exact recipe/source tested>
- <official bioconda-utils lint result>
- <clean-container build/test result, or explicitly outstanding>
- Recipe smoke tests cover imports, CLI, packaged offline guides, the bundled
  demo, saving an installed project and rebuilding its standalone report.

Upstream validation:
- <record the application's Python/terminal/browser checks separately>
```

Replace every bracketed field with recorded facts before using this description. Keep local-candidate source checks distinct from a public-source build and distinguish mock integration checks from real scientific-engine runs.
