"""Data-first launcher and FASTA picker; no third-party packages required."""
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
    origin: str = ''

FASTA_EXTENSIONS={'.fa','.faa','.fna','.fasta','.fas','.fsa','.aln','.afa','.fnt'}

def is_fasta_name(path):
    path=Path(path)
    if path.suffix.lower()=='.gz':path=Path(path.stem)
    return path.suffix.lower() in FASTA_EXTENSIONS

def file_description(path:Path):
    try:
        if path.is_dir():
            return 'Directory · Enter to open'
        return f'FASTA  ·  {path.stat().st_size/1024:.0f} KiB'
    except OSError:return 'Cannot access this location'

def brief(path:Path,ui):
    seqs=read_fasta(path)
    kinds={sequence_type(s.seq) for s in seqs}
    kind=next(iter(kinds)) if len(kinds)==1 else 'mixed / ambiguous'
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
    """Load data without deciding a scientific workflow in advance."""
    return_home=ui.route and ui.route[0]=='Workspace'
    while True:
        try:
            with ui.section('Start'):
                ui.workflow=()
                choice=ui.menu('Open a workspace',[
                    ('1','Load sequences','Open a FASTA, inspect its contents and choose any analysis',''),
                    ('2','Open a saved project','Continue with its datasets and previous results',''),
                ],subtitle='Explore freely. Saving is optional; data and results stay available until you close.',hint='d Demo · ? Manual · t Tools · Esc Back',shortcuts={'d':'demo','?':'help','h':'help','t':'tools'})
                if choice is None:return None
                if choice=='demo':return StartRequest('demo')
                if choice=='help':
                    from av_docs import documentation_menu
                    documentation_menu(ui);continue
                if choice=='2':
                    folder=ui.ask('Project directory (q to cancel)')
                    if folder and folder.casefold()!='q':
                        from av_project import read_project
                        root,_=read_project(folder);return StartRequest('resume',root)
                    continue
                if choice=='tools':
                    import tempfile
                    from av_tools import Toolchain
                    from av_setup import software_hub
                    with tempfile.TemporaryDirectory(prefix='alignview_setup_') as tmp:
                        software_hub(ui,Toolchain(Path(tmp),ui))
                    continue
                first=choose_file(ui,'Sequences — nucleotide or protein FASTA')
                if first is None:continue
                return StartRequest('load',first)
        except HomeRequested:
            if return_home:raise
        except (OSError,ValueError) as exc:
            ui.error(str(exc));ui.page(str(exc),'Cannot open this dataset')


def choose_sequence(ui, records, title='CHOOSE A SEQUENCE', current=''):
    items=[(str(i),r.id,f'{len(r.seq):,} characters · {sequence_type(r.seq)}'+(' · active' if r.id==current else ''),'SEQUENCES') for i,r in enumerate(records,1)]
    choice=ui.menu(title,items,f'{len(records):,} sequences · / search by identifier')
    return records[int(choice)-1] if choice is not None else None


def choose_sequences(ui,records,selected=(),title='Choose sequences to compare'):
    """Select query records without changing the caller's committed selection."""
    from av_ui import menu_item
    records=list(records);known={record.id for record in records}
    checked=set(selected)&known
    if len(records)==1:return records
    def confirmed():
        chosen=[record for record in records if record.id in checked]
        if not chosen:return None,'Select at least one sequence.'
        if len({sequence_type(record.seq) for record in chosen})>1:
            return None,'Choose one sequence type per comparison: nucleotides or proteins.'
        return chosen,''
    if ui.simple:
        ui.heading(title)
        for line in ui.context[:2]:ui.info(line)
        for i,record in enumerate(records,1):ui.info(f'{i:>4}  {record.id} · {len(record.seq):,} characters · {sequence_type(record.seq)}')
        while True:
            default=','.join(str(i) for i,record in enumerate(records,1) if record.id in checked)
            answer=ui.ask('Sequences: numbers/ranges, all, or q to cancel',default,'Example: 1,3-5 selects four sequences.')
            if answer.casefold()=='q':return None
            try:
                if answer.casefold()=='all':checked=set(known)
                else:
                    indexes=set()
                    for token in answer.replace(' ',',').split(','):
                        if not token:continue
                        if '-' in token:
                            start,end=map(int,token.split('-',1))
                            if not 1<=start<=end<=len(records):raise ValueError
                            indexes.update(range(start,end+1))
                        else:indexes.add(int(token))
                    if any(not 1<=i<=len(records) for i in indexes):raise ValueError
                    checked={records[i-1].id for i in indexes}
            except ValueError:ui.warn('Use valid sequence numbers, ranges such as 2-4, or all.');continue
            chosen,notice=confirmed()
            if chosen:return chosen
            ui.warn(notice)
    from av_terminal import screen,read_key,draw
    import re
    active=0;phrase='';typing=False;notice='';previous=None;generation=0
    with ui.section('Query selection'),screen(mouse=True):
        while True:
            signature=(ui.size(),active,phrase,typing,notice,generation)
            if signature!=previous:
                items=[menu_item(str(i),('[x] ' if record.id in checked else '[ ] ')+record.id,
                    f'{len(record.seq):,} characters · {sequence_type(record.seq)} · {record.description}')
                    for i,record in enumerate(records,1)]
                frame=ui._menu_frame(title,items,subtitle=f'{len(checked)} of {len(records)} sequences selected',
                    hint='Space checks a sequence · Enter confirms · a All · c Clear · / Find · Esc Cancel',
                    active=active,phrase=phrase,typing=typing,notice=notice)
                frame.lines=[('  Space check · Enter confirm · a All · c Clear · Esc cancel' if '↑↓ select · Enter open' in line else line.replace(' · F1 details','').replace('Enter opens the result','Enter applies the filter')) for line in frame.lines]
                draw(frame.lines);active=frame.active
                previous=(ui.size(),active,phrase,typing,notice,generation)
            key=read_key()
            if key=='tick':continue
            if key=='eof':return None
            if key=='F2':raise HomeRequested
            if key=='esc':
                if typing or phrase:typing=False;phrase='';active=0;continue
                return None
            if typing:
                if key in ('\r','\n'):typing=False
                elif key in ('\x7f','\b'):phrase=phrase[:-1];active=0
                elif len(key)==1 and key.isprintable():phrase+=key;active=0
                continue
            if key in ('q','Q'):return None
            if key=='/':typing=True;phrase='';active=0;continue
            if key in ('a','A'):checked=set(known);generation+=1;notice='All sequences selected.';continue
            if key in ('c','C'):checked.clear();generation+=1;notice='Selection cleared.';continue
            if key in ('\r','\n'):
                chosen,notice=confirmed()
                if chosen:return chosen
                continue
            click=re.fullmatch(r'mouse:0;\d+;(\d+)M',key)
            if click and int(click[1]) in frame.hit_rows:active=frame.hit_rows[int(click[1])];key=' '
            if key==' ' and frame.entries:
                record=records[int(frame.entries[active].key)-1]
                if record.id in checked:checked.remove(record.id)
                else:checked.add(record.id)
                generation+=1;notice='';continue
            delta=1 if key in ('down','j') else -1 if key in ('up','k') else frame.page_size if key=='pagedown' else -frame.page_size if key=='pageup' else 0
            if key.startswith('mouse:64;'):delta=-3
            elif key.startswith('mouse:65;'):delta=3
            active=max(0,min(len(frame.entries)-1,active+delta))
            if key=='home':active=0
            elif key=='end':active=max(0,len(frame.entries)-1)

MODE_MATRIX={
 ('DNA/RNA','DNA/RNA'):('blastn','global-dna','BLASTN'),
 ('protein','protein'):('blastp','global-aa','BLASTP'),
 ('DNA/RNA','protein'):('blastx',None,'BLASTX'),
 ('protein','DNA/RNA'):('tblastn',None,'TBLASTN'),
}

def suggest_modes(qtype,rtype):
    return MODE_MATRIX.get((qtype,rtype),(None,None,None))
