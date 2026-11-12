"""Guided start screen and FASTA picker; no third-party packages required."""
from __future__ import annotations
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from av_core import read_fasta, sequence_type
from av_ui import CYAN, GREEN, GREY, YELLOW, HomeRequested


@dataclass(frozen=True)
class StartRequest:
    intent: str
    query: Path | None = None
    references: Path | None = None

FASTA_EXTENSIONS={'.fa','.faa','.fna','.fasta','.fas','.fsa','.aln','.afa','.fnt'}

def is_fasta_name(path):
    path=Path(path)
    if path.suffix.lower()=='.gz':path=Path(path.stem)
    return path.suffix.lower() in FASTA_EXTENSIONS

def file_description(path:Path):
    try:
        if path.is_dir():
            return 'Directory · Enter to open'
        return f'FASTA  ·  {path.stat().st_size/1024:.0f} Kio'
    except OSError:return 'Cannot access this location'

def brief(path:Path,ui):
    seqs=read_fasta(path)
    kinds={sequence_type(s.seq) for s in seqs}
    kind=next(iter(kinds)) if len(kinds)==1 else 'mixte / ambigu'
    length=sum(len(s.seq) for s in seqs)
    ui.success(f'{path.name}  ·  {len(seqs)} seq.  ·  {length:,} characters  ·  {kind}')
    return seqs

def choose_file(ui, label, start=None):
    """Choose with arrow keys or a typed path; supports directories and parent navigation."""
    cwd=Path(start or Path.cwd()).expanduser().resolve()
    if not cwd.is_dir():cwd=cwd.parent
    while True:
        options=[('p','Enter / paste a path', 'Absolute or relative path · paste or type','INPUT'),
                 ('u','Parent directory  ..','Go to the parent directory','NAVIGATION')]
        folders=[]; files=[]
        try:
            for path in cwd.iterdir():
                if path.name.startswith('.'):continue
                if path.is_dir():folders.append(path)
                elif path.is_file() and is_fasta_name(path):files.append(path)
        except (PermissionError,OSError) as e:ui.warn(str(e))
        folders.sort(key=lambda p:p.name.casefold());files.sort(key=lambda p:p.name.casefold())
        visible=folders+files
        mapping={}
        for i,path in enumerate(visible,1):
            key=str(i);mapping[key]=path
            options.append((key,('📁 ' if path.is_dir() else '◻ ')+path.name,file_description(path),'DIRECTORIES' if path.is_dir() else 'SEQUENCES'))
        selected=ui.menu('Choose a file — '+label,options,f'Directory: {cwd}',
                         '↑↓ select · / search · Enter open · p path · q cancel')
        if selected is None:return None
        if selected=='u':cwd=cwd.parent;continue
        if selected=='p':
            raw=ui.ask('FASTA path',str(cwd)+'/','Paste a complete path; a directory opens its contents')
            raw=raw.strip().strip('\"').strip("'")
            if not raw:return None
            path=Path(os.path.expanduser(raw)).resolve() if Path(os.path.expanduser(raw)).is_absolute() else (cwd/raw).resolve()
        else:path=mapping[selected]
        if path.is_dir():cwd=path;continue
        try:
            brief(path,ui)
            return path
        except (OSError,ValueError) as exc:ui.warn(f'Invalid FASTA: {exc}')

def startup(ui):
    """Ask for an intended result before exposing engine or search options."""
    return_home=len(ui.route)>1
    while True:
        try:
            with ui.section('Start'):
                choice=ui.menu('What would you like to do?',[
                    ('1','Open an alignment','Read an aligned FASTA without recalculating it',''),
                    ('2','Align sequences','Choose a FASTA file containing sequences to align',''),
                    ('3','Compare with references','Choose a query file and a reference file',''),
                    ('4','Resume a project','Restore saved data and results',''),
                    ('d','Try an example','Try the alignment viewer without installing an engine',''),
                    ('t','Set up tools','Install or locate engines on Linux / HPC',''),
                    ('h','Read the manual','Getting started, features, saving and offline troubleshooting',''),
                ],hint='Choose the result you need. Methods and parameters remain available afterwards.')
                if choice is None:return None
                if choice=='d':return StartRequest('demo')
                if choice=='h':
                    from av_docs import documentation_menu
                    documentation_menu(ui);continue
                if choice=='4':
                    folder=ui.ask('Project directory (q to cancel)')
                    if folder and folder.casefold()!='q':
                        from av_project import read_project
                        root,_=read_project(folder)
                        return StartRequest('resume',root)
                    continue
                if choice=='t':
                    import tempfile
                    from av_tools import Toolchain
                    from av_setup import software_hub
                    with tempfile.TemporaryDirectory(prefix='alignview_setup_') as tmp:
                        software_hub(ui,Toolchain(Path(tmp),ui))
                    continue
                label={'1':'Existing alignment','2':'Sequences to align','3':'Query file'}[choice]
                first=choose_file(ui,label)
                if first is None:continue
                if choice in ('1','2'):return StartRequest('view' if choice=='1' else 'align',first)
                second=choose_file(ui,'Reference sequences',first.parent)
                if second:return StartRequest('compare',first,second)
        except HomeRequested:
            if return_home:raise
        except (OSError,ValueError) as exc:
            ui.error(str(exc));ui.page(str(exc),'Ouverture impossible')


def choose_sequence(ui, records, title='CHOOSE A SEQUENCE', current=''):
    items=[(str(i),r.id,f'{len(r.seq):,} characters · {sequence_type(r.seq)}'+(' · active' if r.id==current else ''),'SEQUENCES') for i,r in enumerate(records,1)]
    choice=ui.menu(title,items,f'{len(records):,} sequences · / search by identifier')
    return records[int(choice)-1] if choice is not None else None

MODE_MATRIX={
 ('DNA/RNA','DNA/RNA'):('blastn','global-dna','BLASTN'),
 ('protein','protein'):('blastp','global-aa','BLASTP'),
 ('DNA/RNA','protein'):('blastx',None,'BLASTX'),
 ('protein','DNA/RNA'):('tblastn',None,'TBLASTN'),
}

def suggest_modes(qtype,rtype):
    return MODE_MATRIX.get((qtype,rtype),(None,None,None))
