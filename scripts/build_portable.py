#!/usr/bin/env python3
"""Bundle a complete application under one AlignView/ directory."""
from pathlib import Path
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_version import VERSION


def build_portable():
    if not (ROOT/'alignview.pyz').is_file():
        raise ValueError('Construire d’abord alignview.pyz : python3 scripts/build_pyz.py')
    files=[p for p in ROOT.iterdir() if p.is_file() and (
        p.suffix in ('.py','.pyz','.sh','.md','.toml','.txt','.yml') or p.name in ('LICENSE','MANIFEST.in'))]
    for directory in ('docs','examples','scripts','tests','packaging','.github'):
        files.extend(p for p in (ROOT/directory).rglob('*') if p.is_file()
                     and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.pyo'))
    destination=ROOT/'dist'/f'AlignView-{VERSION}-portable.zip'
    destination.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(destination,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            archive.write(path,'AlignView/'+path.relative_to(ROOT).as_posix())
    print(destination)
    return destination


if __name__=='__main__':build_portable()
