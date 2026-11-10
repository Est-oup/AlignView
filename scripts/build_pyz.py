#!/usr/bin/env python3
"""Créer un executable portable en un fichier à partir des modules Python."""
from pathlib import Path
import tempfile,zipapp,shutil
root=Path(__file__).resolve().parents[1]
files=['alignview.py',*[p.name for p in sorted(root.glob('av_*.py'))]]
with tempfile.TemporaryDirectory(prefix='alignview_pyz_') as tmp:
    temp=Path(tmp)
    for filename in files: shutil.copyfile(root/filename,temp/filename)
    manual=temp/'alignview_docs';manual.mkdir()
    shutil.copyfile(root/'docs/__init__.py',manual/'__init__.py')
    for document in sorted((root/'docs').glob('*.md')):
        shutil.copyfile(document,manual/document.name)
    (temp/'__main__.py').write_text('from alignview import main\nraise SystemExit(main())\n')
    target=root/'alignview.pyz'
    zipapp.create_archive(temp,target=target,interpreter='/usr/bin/env python3',compressed=True)
    target.chmod(0o755)
    print(target)
