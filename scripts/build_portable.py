#!/usr/bin/env python3
"""Build a reproducible portable application under one AlignView/ directory."""
from pathlib import Path
import argparse
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_version import VERSION
from release_files import portable_files, write_member


def build_portable(destination=None):
    zipapp = ROOT / 'alignview.pyz'
    if not zipapp.is_file():
        raise ValueError('Build alignview.pyz first: python3 scripts/build_pyz.py')
    files = [*portable_files(ROOT), zipapp]
    target = Path(destination) if destination else ROOT / 'dist' / f'AlignView-{VERSION}-portable.zip'
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(target, 'w') as archive:
        for path in sorted(files, key=lambda value: value.relative_to(ROOT).as_posix()):
            relative = path.relative_to(ROOT).as_posix()
            executable = path.suffix in ('.sh', '.pyz') or bool(path.stat().st_mode & 0o111)
            write_member(archive, 'AlignView/' + relative, path.read_bytes(), executable=executable)
    print(target)
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Write the portable ZIP to this path')
    build_portable(parser.parse_args().output)
