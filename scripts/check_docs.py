#!/usr/bin/env python3
"""Check current manual links, code fences and the bundled English manual."""
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
        if not (ROOT/'docs'/filename).is_file():problems.append('Missing bundled guide: '+filename)
    for obsolete in ('README_EN.md', 'av_i18n.py', 'av_locale_en.py', 'docs/en'):
        if (ROOT/obsolete).exists():problems.append('Obsolete duplicate-language resource: '+obsolete)
    if problems:raise SystemExit('\n'.join(problems))
    print(f'{len(files)} Markdown files: local links, fences and English guides passed')


if __name__=='__main__':main()
