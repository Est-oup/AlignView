# Prepare and verify a release

`av_version.VERSION` is the application version. The wheel, source archive, zipapp, portable ZIP and checksum manifest must all agree with it. A version number identifies the supported contract; it does not establish biological accuracy or universal platform compatibility. See the [1.0 release contract](RELEASE_CONTRACT.md).

## Build environment

Use a Python environment containing the build tools. The application itself does not depend on them.

```bash
python3 -m pip install 'setuptools>=77' wheel build twine
python3 scripts/build_release.py
```

The command runs the full test suite, checks generated CLI documentation and manual links, builds the wheel/source archive from a fresh temporary source copy, normalizes archive metadata, builds both portable forms, checks distribution metadata and checksums, and exercises the installed artifacts. It publishes nothing and does not install scientific engines.

When the full suite has already passed on the exact current sources, `--skip-tests` avoids repeating it. `--contents-only` explicitly omits execution of the installed artifacts; it is not sufficient for a final release.

## Result files

For a version such as `1.0.0`, the complete release set is:

```text
dist/
  alignview-1.0.0-py3-none-any.whl
  alignview-1.0.0.tar.gz
  alignview-1.0.0.pyz
  AlignView-1.0.0-portable.zip
  SHA256SUMS
```

A matching `alignview.pyz` is also written at the application root for local use. Every checksum path stays within `dist/`; older versions in that directory are not added to the current manifest.

```bash
cd dist
sha256sum --check SHA256SUMS
```

## Independent verification

```bash
python3 scripts/verify_release.py
python3 scripts/check_report_browser.py
python3 scripts/validate_real_tools.py --require mafft blastn diamond needle clustalo FastTree
python3 scripts/validate_vsearch.py
```

The artifact verifier checks the declared module inventory, rejects caches/engines/prototypes, compares packaged runtime files and manuals with current sources, checks embedded versions and validates checksums. It uses pip with `--no-index --no-deps` to install the actual wheel in a temporary target. Application smoke checks run outside the checkout; zipapp imports use Python's isolated mode. Every installation must preserve the caller directory, generate standalone HTML, save a project in a path containing spaces, resume it and rebuild its report without a scientific engine. Both task verbs and compatible legacy flag commands are exercised.

Browser and real-engine checks require their external tools. Record versions, commands, settings and results, and distinguish these checks from mock integration tests. Validate the intended SSH terminals, scheduler environment, export consumers and representative biological datasets separately.

## Reproducibility

The builders sort archive members and normalize timestamps, modes and owners. `SOURCE_DATE_EPOCH` controls archive time; its stable default is 2000-01-01. Compare artifacts only with unchanged source and identical Python/build-tool versions. Source archives and wheels are built from a clean temporary copy so stale `build/lib` files cannot leak into a release.

```bash
SOURCE_DATE_EPOCH=946684800 python3 scripts/build_release.py --skip-tests
```

## Publication

The GitHub artifact workflow builds and uploads files for review; it does not publish to PyPI, create a GitHub release or submit a Bioconda recipe. A tag must match `v` followed by the application version. Finalize repository metadata, authorship/contact and registry availability before publication. Scientific engine binaries are never bundled with the application.

For Bioconda, use the [preparation guide](../packaging/bioconda/README.md) and [submission handoff](../packaging/bioconda/SUBMISSION.md). `scripts/prepare_bioconda.py --local` creates a concrete recipe with the real local archive SHA256 for build testing. Submission mode downloads a stable versioned HTTPS source and checks its independently supplied checksum and runtime/manual contents against these sources. A local-file recipe is not ready to submit, even after a successful conda-build. Publish the reviewed source first, regenerate against that public archive, and validate the final recipe using Bioconda's lint/build pipeline.
