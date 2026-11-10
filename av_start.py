"""Accueil guidé et sélecteur de FASTA ; aucun paquet tiers requis."""
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
            return 'Dossier · Entrée pour ouvrir'
        return f'FASTA  ·  {path.stat().st_size/1024:.0f} Kio'
    except OSError:return 'Accès impossible'

def brief(path:Path,ui):
    seqs=read_fasta(path)
    kinds={sequence_type(s.seq) for s in seqs}
    kind=next(iter(kinds)) if len(kinds)==1 else 'mixte / ambigu'
    length=sum(len(s.seq) for s in seqs)
    ui.success(f'{path.name}  ·  {len(seqs)} séq.  ·  {length:,} caractères  ·  {kind}')
    return seqs

def choose_file(ui, label, start=None):
    """Sélection via flèches ou chemin tapé ; ouverture de dossier et retour parent."""
    cwd=Path(start or Path.cwd()).expanduser().resolve()
    if not cwd.is_dir():cwd=cwd.parent
    while True:
        options=[('p','Saisir / coller un chemin', 'Chemin absolu ou relatif · à coller ou saisir','SAISIE'),
                 ('u','Dossier parent  ..','Remonter au dossier précédent','NAVIGATION')]
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
            options.append((key,('📁 ' if path.is_dir() else '◻ ')+path.name,file_description(path),'DOSSIERS' if path.is_dir() else 'SÉQUENCES'))
        selected=ui.menu('Choisir un fichier — '+label,options,f'Dossier : {cwd}',
                         '↑↓ choisir · / rechercher · Entrée ouvrir · p chemin · q annuler')
        if selected is None:return None
        if selected=='u':cwd=cwd.parent;continue
        if selected=='p':
            raw=ui.ask('Chemin du FASTA',str(cwd)+'/','Tu peux coller un chemin complet ; un dossier ouvre son contenu')
            raw=raw.strip().strip('\"').strip("'")
            if not raw:return None
            path=Path(os.path.expanduser(raw)).resolve() if Path(os.path.expanduser(raw)).is_absolute() else (cwd/raw).resolve()
        else:path=mapping[selected]
        if path.is_dir():cwd=path;continue
        try:
            brief(path,ui)
            return path
        except (OSError,ValueError) as exc:ui.warn(f'FASTA invalide : {exc}')

def startup(ui):
    """Ask for an intended result before exposing engine or search options."""
    return_home=len(ui.route)>1
    while True:
        try:
            with ui.section('Démarrer'):
                choice=ui.menu('Que souhaitez-vous faire ?',[
                    ('1','Ouvrir un alignement','Lire un FASTA déjà aligné, sans refaire les calculs',''),
                    ('2','Aligner des séquences','Choisir un FASTA de séquences à aligner',''),
                    ('3','Comparer à des références','Choisir une requête et un fichier de références',''),
                    ('4','Reprendre un projet','Retrouver des données et des résultats enregistrés',''),
                    ('d','Essayer un exemple','Découvrir la visionneuse sans installer de moteur',''),
                    ('t','Préparer les logiciels','Installer ou retrouver les moteurs sur Linux / HPC',''),
                    ('h','Lire la documentation','Premiers pas, fonctions, sauvegardes et dépannage hors ligne',''),
                ],hint='Choisir le résultat attendu. Les méthodes et paramètres restent accessibles ensuite.')
                if choice is None:return None
                if choice=='d':return StartRequest('demo')
                if choice=='h':
                    from av_docs import documentation_menu
                    documentation_menu(ui);continue
                if choice=='4':
                    folder=ui.ask('Dossier du projet (q pour annuler)')
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
                label={'1':'Alignement déjà calculé','2':'Séquences à aligner','3':'Fichier de requêtes'}[choice]
                first=choose_file(ui,label)
                if first is None:continue
                if choice in ('1','2'):return StartRequest('view' if choice=='1' else 'align',first)
                second=choose_file(ui,'Séquences de référence',first.parent)
                if second:return StartRequest('compare',first,second)
        except HomeRequested:
            if return_home:raise
        except (OSError,ValueError) as exc:
            ui.error(str(exc));ui.page(str(exc),'Ouverture impossible')


def choose_sequence(ui, records, title='CHOISIR UNE SÉQUENCE', current=''):
    items=[(str(i),r.id,f'{len(r.seq):,} caractères · {sequence_type(r.seq)}'+(' · active' if r.id==current else ''),'SÉQUENCES') for i,r in enumerate(records,1)]
    choice=ui.menu(title,items,f'{len(records):,} séquences · / rechercher par identifiant')
    return records[int(choice)-1] if choice is not None else None

MODE_MATRIX={
 ('ADN/ARN','ADN/ARN'):('blastn','global-dna','BLASTN'),
 ('protéine','protéine'):('blastp','global-aa','BLASTP'),
 ('ADN/ARN','protéine'):('blastx',None,'BLASTX'),
 ('protéine','ADN/ARN'):('tblastn',None,'TBLASTN'),
}

def suggest_modes(qtype,rtype):
    return MODE_MATRIX.get((qtype,rtype),(None,None,None))
