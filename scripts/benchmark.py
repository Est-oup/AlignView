#!/usr/bin/env python3
"""Repeatable synthetic viewport benchmark; optional comparison with the local v6."""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import re
import statistics
import sys
import time
import unicodedata
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import Seq
from av_msa_tui import MSAViewer
from av_version import VERSION


def measure(viewer_class,records,repeats):
    start=time.perf_counter();viewer=viewer_class(records,False)
    viewer.message='Position atteinte ; navigation active.'
    initialization=time.perf_counter()-start
    timings=[];lines=0;maximum=0;physical=0
    with patch('shutil.get_terminal_size',return_value=os.terminal_size((80,24))):
        for i in range(repeats):
            viewer.col=(i*83)%len(records[0].seq);viewer.row=i%len(records)
            output=io.StringIO();start=time.perf_counter()
            with contextlib.redirect_stdout(output):viewer._render()
            timings.append((time.perf_counter()-start)*1000)
            frame=re.sub(r'\x1b\[[0-9;?]*[A-Za-z]','',output.getvalue()).splitlines()
            widths=[sum(0 if unicodedata.combining(c) else 2 if unicodedata.east_asian_width(c) in ('W','F') else 1 for c in line) for line in frame]
            lines=max(lines,len(frame));maximum=max(maximum,max(widths,default=0))
            physical=max(physical,sum(max(1,(width+79)//80) for width in widths))
    return dict(initialization_s=round(initialization,6),
                render_median_ms=round(statistics.median(timings),4),
                render_p95_ms=round(sorted(timings)[int(.95*(len(timings)-1))],4),
                maximum_frame_lines=lines,maximum_visible_line_cells=maximum,
                maximum_terminal_rows_including_wrapping=physical,
                terminal_rows=24,terminal_columns=80,repeats=repeats)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows',type=int,default=100)
    parser.add_argument('--columns',type=int,default=10000)
    parser.add_argument('--repeats',type=int,default=100)
    parser.add_argument('--baseline',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if min(args.rows,args.columns,args.repeats)<1:parser.error('Use positive dimensions and repeats.')
    template=('ACGT'*((args.columns+3)//4))[:args.columns]
    records=[]
    for i in range(args.rows):
        chars=list(template)
        for j in range(i%17,args.columns,97):chars[j]='ACGT-'[(i+j)%5]
        records.append(Seq('seq_'+str(i),'seq_'+str(i),''.join(chars)))
    result=dict(version=VERSION,python=platform.python_version(),platform=platform.platform(),
                synthetic_rows=args.rows,synthetic_columns=args.columns,
                note='Local synthetic timings, no biological accuracy or competitor performance claim.',
                alignview=measure(MSAViewer,records,args.repeats))
    if args.baseline:
        spec=importlib.util.spec_from_file_location('baseline_viewer',args.baseline/'av_msa_tui.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        result['baseline_v6']=measure(module.MSAViewer,records,args.repeats)
    encoded=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8')
    print(encoded,end='')


if __name__=='__main__':main()
