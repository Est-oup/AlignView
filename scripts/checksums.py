#!/usr/bin/env python3
"""Create SHA256SUMS beside the release archives, with portable relative paths."""
import hashlib
from pathlib import Path

root=Path(__file__).resolve().parents[1]
folder=root/'dist';folder.mkdir(exist_ok=True)
files=sorted(folder.glob('*.whl'))+sorted(folder.glob('*.tar.gz'))+sorted(folder.glob('*.zip'))
if (root/'alignview.pyz').is_file():files.append(root/'alignview.pyz')
entries=[]
for path in files:
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    name=path.name if path.parent==folder else '../'+path.name
    entries.append(f'{digest}  {name}')
(folder/'SHA256SUMS').write_text('\n'.join(entries)+'\n')
print(folder/'SHA256SUMS')
