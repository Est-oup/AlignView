#!/usr/bin/env python3
"""Check current manual links, code fences and bundled language resources."""
from pathlib import Path
import re
import sys
from urllib.parse import unquote

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_docs import TOPICS


def main():
    problems=[]
    files=[*ROOT.glob('*.md'),*(ROOT/'docs').rglob('*.md'),*(ROOT/'packaging').rglob('*.md')]
    for path in files:
        text=path.read_text(encoding='utf-8')
        if sum(line.lstrip().startswith('```') for line in text.splitlines())%2:
            problems.append(str(path.relative_to(ROOT))+': unclosed code block')
        for link in re.findall(r'!?\[[^\]]*\]\(([^)\n]+)\)',text):
            target=link.strip().split('#',1)[0]
            if not target or re.match(r'[a-zA-Z][a-zA-Z0-9+.-]*:',target):continue
            target=unquote(target.strip('<>'))
            if not (path.parent/target).exists():
                problems.append(str(path.relative_to(ROOT))+': missing '+target)
    for filename,_,_ in TOPICS.values():
        for relative in (filename,'en/'+filename):
            if not (ROOT/'docs'/relative).is_file():problems.append('Missing bundled guide: '+relative)
    if problems:raise SystemExit('\n'.join(problems))
    print(f'{len(files)} Markdown files: local links, fences and French/English guides passed')


if __name__=='__main__':main()
