#!/usr/bin/env python3
"""Hash only this version's complete release set, all inside dist/."""
import hashlib
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_version import VERSION


def release_artifacts(root=ROOT):
    folder = Path(root) / 'dist'
    names = (f'alignview-{VERSION}-py3-none-any.whl', f'alignview-{VERSION}.tar.gz',
             f'AlignView-{VERSION}-portable.zip', f'alignview-{VERSION}.pyz')
    paths = [folder / name for name in names]
    missing = [path.name for path in paths if not path.is_file()]
    if missing:
        raise ValueError('Build the complete release before checksums: ' + ', '.join(missing))
    return sorted(paths, key=lambda path: path.name)


def checksums():
    paths = release_artifacts()
    output = ROOT / 'dist/SHA256SUMS'
    output.write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n'
                              for path in paths), encoding='utf-8')
    print(output)
    return output


if __name__ == '__main__':
    checksums()
