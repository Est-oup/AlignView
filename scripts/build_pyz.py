#!/usr/bin/env python3
"""Build a reproducible zipapp from the declared runtime inventory."""
from pathlib import Path
import argparse
import shutil
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from av_version import VERSION
from release_files import manual_files, runtime_files, write_member


def build_pyz(destination=None):
    target = Path(destination) if destination else ROOT / 'alignview.pyz'
    target.parent.mkdir(parents=True, exist_ok=True)
    resources = [(path.name, path.read_bytes()) for path in runtime_files(ROOT)]
    resources.extend(('alignview_docs/' + path.relative_to(ROOT / 'docs').as_posix(), path.read_bytes())
                     for path in manual_files(ROOT))
    resources.append(('__main__.py', b'from alignview import main\nraise SystemExit(main())\n'))
    with target.open('wb') as stream:
        stream.write(b'#!/usr/bin/env python3\n')
        with zipfile.ZipFile(stream, 'w') as archive:
            for name, data in sorted(resources):
                write_member(archive, name, data)
    target.chmod(0o755)
    if destination is None:
        release_target = ROOT / 'dist' / f'alignview-{VERSION}.pyz'
        release_target.parent.mkdir(exist_ok=True)
        shutil.copyfile(target, release_target)
        release_target.chmod(0o755)
    print(target)
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Build to this file without changing the default release outputs')
    build_pyz(parser.parse_args().output)
