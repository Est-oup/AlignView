# Prepare a release

The prepared local version is 0.1.0. Packaging does not publish to GitHub, PyPI or Bioconda. Validate the actual engines, target terminals, cluster and biological data separately.

```bash
python3 -m unittest discover -s tests -v
python3 scripts/generate_cli_docs.py --check
python3 scripts/check_docs.py
python3 -m build
python3 -m twine check dist/*.whl dist/*.tar.gz
python3 scripts/build_pyz.py
python3 scripts/build_portable.py
python3 scripts/checksums.py
```

Outputs are the wheel, source archive, portable ZIP, zipapp and SHA256SUMS. Run the wheel and zipapp from a different directory; verify offline French/English guides, real input paths, independent output locations and the tree/VSEARCH report features. Extract the portable ZIP into a fresh directory and launch its AlignView/alignview.sh. Check checksums after every rebuild.

Before publication, finalize repository URL, authorship/contact, version and registry availability. The Bioconda recipe remains a template until release URL, checksum and maintainer information are real. The application package does not contain scientific engine binaries.

A stable 1.0 designation needs validation across the intended machines, SSH terminals, engines, cluster workflows, exported formats and representative biological datasets. Existing mocks and a version number do not establish these properties.
