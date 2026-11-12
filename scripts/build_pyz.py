#!/usr/bin/env python3
"""Build a single-file portable executable from the Python modules."""
from pathlib import Path
import tempfile,zipapp,shutil
root=Path(__file__).resolve().parents[1]
files=['alignview.py',*[p.name for p in sorted(root.glob('av_*.py'))]]
with tempfile.TemporaryDirectory(prefix='alignview_pyz_') as tmp:
    temp=Path(tmp)
    for filename in files: shutil.copyfile(root/filename,temp/filename)
    manual=temp/'alignview_docs';manual.mkdir()
    shutil.copyfile(root/'docs/__init__.py',manual/'__init__.py')
    for document in sorted((root/'docs').rglob('*.md')):
        target_document=manual/document.relative_to(root/'docs')
        target_document.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(document,target_document)
    (temp/'__main__.py').write_text('from alignview import main\nraise SystemExit(main())\n')
    target=root/'alignview.pyz'
    zipapp.create_archive(temp,target=target,interpreter='/usr/bin/env python3',compressed=True)
    target.chmod(0o755)
    print(target)
