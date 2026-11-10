#!/usr/bin/env python3
"""Generate documentation previews from the actual terminal frames, offline."""
from pathlib import Path
import html
import os
import re
import sys
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,parse_args
from av_msa_tui import MSAViewer
from av_navigation import HOME_ITEMS,Workbench,item
from av_terminal import fit
from av_ui import UI


def svg_frame(lines,destination,columns=100):
    width=columns*8.6+28;height=len(lines)*22+36
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Aperçu généré à partir de l’écran AlignView">',
         f'<rect width="{width}" height="{height}" rx="12" fill="#0d1729"/>']
    for i,raw in enumerate(lines):
        line=re.sub(r'\x1b\[[0-9;]*m','',fit(raw,columns-1))
        selected=line.startswith(' ›')
        if selected:svg.append(f'<rect x="8" y="{11+i*22}" width="{width-16}" height="22" rx="4" fill="#19354a"/>')
        fill='#70dfc0' if selected or i==0 else '#d5e0f4' if re.match(r'\s+\d\s',line) else '#a5b6d0'
        svg.append(f'<text x="14" y="{28+i*22}" fill="{fill}" font-family="DejaVu Sans Mono,monospace" font-size="14" xml:space="preserve">{html.escape(line)}</text>')
    svg.append('</svg>');destination.write_text('\n'.join(svg),encoding='utf-8')
    return lines


ui=UI(False)
session=Session(parse_args([str(ROOT/'examples/example_msa.fasta'),'--msa','view','--no-menu']),ui)
try:
    session.load_msa();workbench=Workbench(session);workbench.refresh_context()
    with patch('shutil.get_terminal_size',return_value=os.terminal_size((100,28))):
        home=ui._menu_frame('Projet',[item(*entry) for entry in HOME_ITEMS],hint='L’alignement est prêt : ouvrir Explorer pour le lire.')
        svg_frame(home.lines,ROOT/'docs/navigation-preview.svg')
        svg_frame(MSAViewer(session.msa_records,False)._frame(),ROOT/'docs/terminal-preview.svg')
finally:session.close()
print(ROOT/'docs/navigation-preview.svg')
