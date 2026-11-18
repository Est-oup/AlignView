# Prepare a Bioconda submission

The recipe is prepared locally; AlignView has not been submitted or published to Bioconda by these scripts. `meta.yaml.template` is a generator input, not a recipe to copy directly into `bioconda-recipes`.

The candidate packages the Python application as `noarch: python`, including its `alignview` console entry point and offline manual. BLAST, DIAMOND, MAFFT and the other scientific engines remain optional. Existing alignments, sequence inspection, the demonstration and reports work without them. To enable the usual calculations, install `alignview blast diamond mafft` together; the existing `environment-full.yml` describes the broader scientific tool set. Engine binaries are never bundled into the Python package.

## Current public source status

Checked on 18 November 2026: the upstream repository exists at [Est-oup/AlignView](https://github.com/Est-oup/AlignView), but its latest public release/tag is `v0.1.0`. The application in this checkout is `1.0.0`. Public PyPI, Bioconda and conda-forge package queries for `alignview` returned HTTP 404. This is a point-in-time check, not a name reservation.

A public `v0.1.0` archive cannot represent these `1.0.0` sources. Publish the reviewed versioned source archive before generating the public-source recipe. Do not substitute a moving branch, a guessed checksum or an older release's checksum. The prospective maintainer must confirm that they accept maintenance; the repository owner's handle is not evidence of consent.

## Build a local candidate

Build the current release first; see [release instructions](../../docs/RELEASING.md). From the application root:

```bash
python3 scripts/build_release.py
python3 scripts/prepare_bioconda.py \
  --local --archive dist/alignview-1.0.0.tar.gz \
  --maintainer Est-oup \
  --output dist/bioconda-local/recipes/alignview
```

Replace the maintainer argument with the real consenting GitHub maintainer. The example `Est-oup` is the current upstream owner's handle and remains a proposed maintainer until confirmed.

The generator inspects the tarball without extracting or executing it, verifies its runtime modules, manual and package metadata against the current checkout, computes the actual SHA256, and creates:

- `dist/bioconda-local/recipes/alignview/meta.yaml`: a buildable recipe with a `file://` source URL, explicitly unsuitable for submission.
- `dist/bioconda-local/recipes/alignview-preparation.json`: source hash, checked files, candidate/public-source status and outstanding validation. Build/lint fields start as false; generating a recipe does not run either check.

The destination must be new. Rebuild the source archive after changing runtime files, manuals or package metadata. A stale archive is rejected.

## Generate the public-source recipe

Upload the final source tarball as an immutable versioned release asset, record its SHA256, and verify the download. The tag, archive metadata and `av_version.VERSION` must agree. Then:

```bash
python3 scripts/prepare_bioconda.py \
  --source-url https://github.com/Est-oup/AlignView/releases/download/v1.0.0/alignview-1.0.0.tar.gz \
  --sha256 ACTUAL_PUBLIC_ARCHIVE_SHA256 \
  --maintainer CONSENTING_GITHUB_HANDLE \
  --output dist/bioconda-public/recipes/alignview
```

These are example publication coordinates, not a claim that the URL currently exists. The generator downloads the HTTPS URL, checks the supplied hash, verifies the archive against the current runtime/manual/metadata, and only then writes the recipe. Passing `--archive dist/alignview-1.0.0.tar.gz` additionally requires the public archive to be byte-identical to the local release artifact. Without this option, a differently compressed/tag-generated source archive is accepted only if its checked file contents are identical.

## Validate and submit

Use a separate temporary Conda environment containing `conda-build` and `bioconda-utils`, with strict channel priority and `conda-forge` before `bioconda`.

A basic local build (candidate or public-source recipe) is:

```bash
conda build --no-anaconda-upload --override-channels \
  -c conda-forge -c bioconda \
  --croot /tmp/alignview-conda-build \
  dist/bioconda-local/recipes/alignview
```

The recipe checks module imports, the console entry point, offline guides, demonstration, project saving and report rebuilding from an installed package. It uses only runtime dependencies and creates its own tiny demonstration data. Upstream unit tests and mock executables are not recipe tests.

After the public-source build passes, copy only the generated `recipes/alignview/` directory into a clean branch of your fork of `bioconda-recipes`. The preparation JSON is review evidence and does not belong in the recipe directory. From that checkout run:

```bash
bioconda-utils lint recipes config.yml --packages alignview
bioconda-utils build recipes config.yml --packages alignview --docker --mulled-build-and-test
```

If Docker is unavailable, omit the two Docker/mulled flags for the local build and keep the clean-container test as outstanding for Bioconda CI. Confirm the current channel name remains available. Review [SUBMISSION.md](SUBMISSION.md), then open a pull request manually. The preparation script never pushes, opens a PR, uploads a package or publishes a release.

## Official requirements consulted

- [Bioconda recipe guidelines](https://bioconda.github.io/contributor/guidelines.html): stable source archive and checksum, bioinformatics scope, pure-Python noarch packages and installed CLI tests.
- [Bioconda lint rules](https://bioconda.github.io/contributor/linting.html): source archives, matching recipe/package names and runtime version exports.
- [Bioconda local validation](https://bioconda.github.io/contributor/building-locally.html): official lint/build commands and clean-container checks.
- [Bioconda contribution workflow](https://bioconda.github.io/contributor/workflow.html): recipe review and publication happen through a PR.
- [Conda-forge noarch guidance](https://conda-forge.org/docs/how-to/basics/noarch/): explicit console entry points and Python-minimum constraints.
