#!/usr/bin/env python3
"""Generate documentation previews from the actual terminal frames, offline."""
from pathlib import Path
import html
import os
import re
import sys
import xml.etree.ElementTree as ET
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import Session,parse_args
from av_msa_tui import MSAViewer
from av_workspace import Workspace
from av_app import _register_inputs
from av_terminal import fit
from av_ui import UI


def svg_frame(lines,destination,columns=100):
    width=columns*8.6+28;height=len(lines)*22+36
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Preview generated from the AlignView screen">',
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
    session.library_root=Path(session.temporary.name)/'catalog'
    _register_inputs(session)
    workspace=Workspace(session)
    with patch('shutil.get_terminal_size',return_value=os.terminal_size((100,28))):
        workspace.switch('Home')
        data_frame=Path(session.temporary.name)/'data.svg'
        svg_frame(workspace.render(100,28),data_frame)
        svg_frame(workspace.render(100,28),ROOT/'docs/navigation-preview.svg')
        session.load_msa()
        workspace.switch('ResultHome')
        result_frame=Path(session.temporary.name)/'result.svg'
        svg_frame(workspace.render(100,28),result_frame)
        ns='{http://www.w3.org/2000/svg}'
        ET.register_namespace('',ns[1:-1])
        frames=[ET.parse(path).getroot() for path in (data_frame,result_frame)]
        height=32+24+sum(int(frame.get('height')) for frame in frames)
        preview=ET.Element(ns+'svg',dict(width='920',height=str(height),viewBox=f'0 0 920 {height}',role='img',**{'aria-label':'Actual AlignView terminal frames: file actions and explicit result opening'}))
        ET.SubElement(preview,ns+'rect',dict(width='920',height=str(height),fill='#f0f4f6'))
        position=16
        for frame in frames:
            group=ET.SubElement(preview,ns+'g',dict(transform=f'translate(16,{position})'))
            group.extend(list(frame))
            position+=int(frame.get('height'))+24
        ET.ElementTree(preview).write(ROOT/'docs/free-exploration-preview.svg',encoding='unicode')
        viewer=MSAViewer(session.msa_records,False)
        viewer.workflow_lines=[]
        svg_frame(viewer._frame(),ROOT/'docs/terminal-preview.svg')
finally:session.close()
print(ROOT/'docs/navigation-preview.svg')
