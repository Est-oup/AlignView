"""Consistent terminal menus, context, breadcrumbs and readable result pages."""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
import os
import re
import shutil
import sys
import textwrap

from av_terminal import draw, fit, read_key, screen
from av_version import VERSION
from av_i18n import translate as tr, language

RESET='\033[0m'; BOLD='\033[1m'; DIM='\033[2m'
BLUE='\033[38;5;75m'; CYAN='\033[38;5;44m'; GREEN='\033[38;5;78m'
YELLOW='\033[38;5;221m'; RED='\033[38;5;203m'; FG='\033[38;5;255m'
BG='\033[48;5;23m'; GREY='\033[38;5;109m'; HIGHLIGHT='\033[48;5;24m'


class HomeRequested(Exception):
    """Unwind nested menus to the workspace, without changing the session."""


@dataclass(frozen=True)
class MenuItem:
    key: str
    label: str
    description: str = ''
    category: str = ''
    enabled: bool = True
    reason: str = ''


@dataclass
class MenuFrame:
    lines: list[str]
    entries: list[MenuItem]
    active: int
    hit_rows: dict[int, int]
    page_size: int


RESULT_TITLES = {
    'show_hits': 'Correspondances trouvées', 'show_coverage': 'Régions couvertes par la recherche',
    'show_dashboard': 'Données et calculs du projet', 'show_stats': 'Qualité des séquences',
    'guide': 'Comprendre les méthodes', 'hsp_diagnostics': 'Qualité des segments alignés',
    'identity_coverage': 'Identité et couverture', 'list_outputs': 'Fichiers du projet',
    'diagnostics': 'Logiciels disponibles', 'msa_analysis': 'Conservation de l’alignement',
    'msa_distances': 'Distances entre les séquences',
}


class UI:
    def __init__(self, color=True, simple=False):
        self.tty=sys.stdin.isatty() and sys.stdout.isatty()
        self.color=color and self.tty and not os.environ.get('NO_COLOR') and os.environ.get('TERM')!='dumb'
        self.simple=simple or not self.tty
        self.positions={}
        self.route=['Accueil']
        self.context=()
        self.notice=''
        self.status=''; self.messages=[]

    def c(self, style, text):
        return style+tr(text)+RESET if self.color else tr(text)

    def size(self):
        size=shutil.get_terminal_size((104,32))
        return size.columns,size.lines

    @contextmanager
    def section(self, label):
        self.route.append(label)
        try:
            yield
        finally:
            self.route.pop()

    def banner(self, subtitle=''):
        print(self.c(CYAN+BOLD,'  AlignView')+'  '+self.c(GREY,VERSION))
        if subtitle:print('  '+tr(subtitle))

    def heading(self, text):
        print('\n  '+self.c(BLUE+BOLD,text))

    def info(self, text):print('  '+tr(text))
    def success(self, text):
        self.notice=tr(text)
        print('  '+self.c(GREEN,'✓')+' '+tr(text))
    def warn(self, text):
        self.notice=tr(text)
        print('  '+self.c(YELLOW,'!')+' '+tr(text),file=sys.stderr)
    def error(self, text):
        self.notice=tr(text)
        print('  '+self.c(RED,'✕')+' '+tr(text),file=sys.stderr)
    def stage(self, text):print('\n  '+self.c(CYAN,text),flush=True)

    def ask(self, label, default='', hint=''):
        if hint:self.wrapped(hint)
        try:
            value=input('  '+tr(label)+(f' [{tr(default)}]' if default!='' else '')+' : ').strip()
        except EOFError:
            raise KeyboardInterrupt from None
        return value if value else default

    def yesno(self, label, default=False):
        return self.ask(label+(' (y/n)' if language()=='en' else ' (o/n)'),('y' if language()=='en' else 'o') if default else 'n').casefold() in ('o','oui','y','yes')

    def pause(self):
        if not self.tty:return
        try:input('\n  '+tr('Entrée pour continuer…')+' ')
        except (EOFError,KeyboardInterrupt):pass

    def wrapped(self, text, prefix='  ', width=None):
        for line in textwrap.wrap(tr(text),max(10,(width or self.size()[0]-2)-len(prefix))):print(prefix+line)

    def menu(self, title, items, subtitle='', footer='', *, hint='', default=None, shortcuts=None):
        entries=[item if isinstance(item,MenuItem) else MenuItem(*map(str,item)) for item in items]
        entries=[replace(item,label=tr(item.label),description=tr(item.description),category=tr(item.category),reason=tr(item.reason)) for item in entries]
        title,subtitle,footer,hint=map(tr,(title,subtitle,footer,hint))
        if self.simple:return self._menu_simple(title,entries,subtitle,hint,shortcuts or {})
        return self._menu_keys(title,entries,subtitle,footer,hint,default,shortcuts or {})

    def _menu_simple(self, title, items, subtitle, hint, shortcuts):
        self.banner();print('  '+' › '.join(self.route));self.heading(title)
        for line in self.context:self.info(line)
        if subtitle:self.wrapped(subtitle)
        if hint:self.wrapped(hint)
        for item in items:
            self.info(f'{item.key:>2}  {item.label}')
            self.wrapped(item.description if item.enabled else 'Indisponible : '+item.reason,prefix='      ')
        while True:
            value=self.ask('Choix (q : retour)', 'q').casefold()
            if value in ('q','0'):return None
            if value=='accueil':raise HomeRequested
            if value in shortcuts:return shortcuts[value]
            exact=next((item for item in items if item.key.casefold()==value),None)
            matches=[item for item in items if value in item.label.casefold()]
            selected=exact or (matches[0] if len(matches)==1 else None)
            if selected:
                if selected.enabled:return selected.key
                self.warn(selected.reason)
            else:self.warn('Choix inconnu ; saisir le numéro affiché.')

    def _menu_frame(self,title,items,subtitle='',hint='',active=0,phrase='',typing=False,digits='',notice=''):
        width,height=self.size()
        entries=[item for item in items if phrase.casefold() in ' '.join((item.label,item.description,item.reason,item.key)).casefold()]
        active=max(0,min(active,len(entries)-1))
        compact=height<19
        rich=len(items)<=7 and height>=22
        header=[self.c(CYAN+BOLD,'  AlignView')+'  '+self.c(GREY,VERSION),
                '  '+self.c(GREY,' › '.join(self.route)), '  '+self.c(BOLD,title)]
        if not compact:
            header.extend('  '+self.c(GREY,line) for line in self.context[:2])
            if subtitle and hint:
                header.extend('  '+fit(value,max(10,width-4)) for value in (subtitle,hint))
            else:
                header.extend('  '+line for line in textwrap.wrap(hint or subtitle,max(10,width-4))[:2])
            header.append('')
        selected=entries[active] if entries else None
        detail=notice or (('' if rich else selected.description) if selected and selected.enabled else selected.reason if selected else 'Aucun résultat. Échap pour effacer la recherche.')
        details=textwrap.wrap(detail,max(10,width-4))[:2]
        details+=['']*(2-len(details))
        if typing:control='Recherche : '+phrase+'▏ · Entrée ouvre le résultat'
        elif digits:control='Choix : '+digits+' · Entrée pour ouvrir'
        else:control=('↑↓ choisir · Entrée ouvrir · q quitter' if self.route==['Accueil'] and title==tr('Projet')
                      else '↑↓ choisir · Entrée ouvrir · ←/Échap retour · F2 accueil')
        footer=['  '+self.c(YELLOW if notice or (selected and not selected.enabled) else GREY,line) for line in details]
        footer.extend(['', '  '+control, ''])
        if compact:footer=footer[:1]+['  '+control,'']
        unit=2 if rich else 1
        count=max(1,(height-1-len(header)-len(footer))//unit)
        start=min(max(0,active-count//2),max(0,len(entries)-count))
        frame=list(header);hit_rows={}
        for index in range(start,min(len(entries),start+count)):
            item=entries[index]
            label=f' {"›" if index==active else " "} {item.key:>2}  {item.label}'+('  [indisponible]' if not item.enabled else '')
            hit_rows[len(frame)+1]=index
            frame.append(self.c(HIGHLIGHT+FG+BOLD,label) if index==active else self.c(GREY,label) if not item.enabled else label)
            if rich:
                hit_rows[len(frame)+1]=index
                frame.append('       '+self.c(GREY,item.description if item.enabled else item.reason))
        if not entries:frame.append('  Aucun résultat pour cette recherche.')
        footer[-1]='  '+self.c(GREY,f'{active+1 if entries else 0}/{len(entries)} · / rechercher · PgUp/PgDn pages · F1 détails')
        frame.extend(footer)
        return MenuFrame(frame,entries,active,hit_rows,count)

    def _menu_keys(self,title,items,subtitle,footer,hint,default,shortcuts):
        position_key=(*self.route,title)
        saved=self.positions.get(position_key,default)
        active=next((i for i,item in enumerate(items) if item.key==saved),0)
        phrase='';digits='';typing=False;notice=self.notice;self.notice=''
        previous=None
        with screen(mouse=True):
            while True:
                signature=(self.size(),active,phrase,typing,digits,notice)
                if signature!=previous:
                    frame=self._menu_frame(title,items,subtitle,hint,active,phrase,typing,digits,notice)
                    active=frame.active;draw(frame.lines)
                    previous=(self.size(),active,phrase,typing,digits,notice)
                key=read_key()
                if key=='tick':continue
                if key=='eof':raise KeyboardInterrupt
                if key=='F2':
                    if len(self.route)>1:raise HomeRequested
                    continue
                if key=='F1':
                    selected=frame.entries[active] if frame.entries else None
                    # Nested alternate-screen contexts would erase the current frame on exit.
                    # The outer loop redraws the full menu after this details page.
                    self.page((selected.label+'\n\n'+selected.description+('\n\n'+selected.reason if selected.reason else '') if selected else 'Aucun résultat.')+'\n\n'+(hint or subtitle)+'\n\n'+footer,'Détails')
                    previous=None
                    continue
                if key=='esc':
                    if phrase or typing:phrase='';typing=False;active=0
                    elif digits:digits=''
                    elif self.route==['Accueil'] and title=='Projet':continue
                    else:return None
                    notice='';continue
                if typing:
                    if key in ('\r','\n'):typing=False
                    elif key in ('\x7f','\b'):phrase=phrase[:-1];active=0;continue
                    elif len(key)==1 and key.isprintable():phrase+=key;active=0;continue
                    else:continue
                else:
                    if key=='left' and self.route==['Accueil'] and title=='Projet':continue
                    if key in ('q','Q','left'):return None
                    if key=='/':typing=True;phrase='';active=0;notice='';continue
                    if key in ('up','k'):active=max(0,active-1);digits=''
                    elif key in ('down','j'):active=min(len(frame.entries)-1,active+1);digits=''
                    elif key in ('pageup','ctrl-up'):active=max(0,active-frame.page_size)
                    elif key in ('pagedown','ctrl-down'):active=min(len(frame.entries)-1,active+frame.page_size)
                    elif key=='home':active=0
                    elif key=='end':active=max(0,len(frame.entries)-1)
                    elif key in ('\x7f','\b'):digits=digits[:-1]
                    elif key.startswith('mouse:'):
                        match=re.fullmatch(r'mouse:(\d+);\d+;(\d+)([Mm])',key)
                        if match:
                            button,y,event=match.groups()
                            if button=='64':active=max(0,active-3)
                            elif button=='65':active=min(len(frame.entries)-1,active+3)
                            elif button=='0' and event=='M' and int(y) in frame.hit_rows:
                                active=frame.hit_rows[int(y)];key='\r';digits=''
                    elif key.isdigit():
                        digits+=key
                        active=next((i for i,item in enumerate(frame.entries) if item.key==digits),active)
                    elif key in shortcuts:return shortcuts[key]
                    elif len(key)==1:
                        selected=next((item for item in items if item.key.casefold()==key.casefold()),None)
                        if selected:
                            if selected.enabled:self.positions[position_key]=selected.key;return selected.key
                            notice=selected.reason
                            continue
                if key in ('\r','\n','right'):
                    selected=next((item for item in items if item.key==digits),None) if digits else frame.entries[active] if frame.entries else None
                    digits=''
                    if selected:
                        if selected.enabled:self.positions[position_key]=selected.key;return selected.key
                        notice=selected.reason
                    else:notice='Aucun choix correspondant.'
                else:notice=''

    def present(self,action):
        name=getattr(action,'__name__','')
        if self.simple or name not in RESULT_TITLES:return action()
        if name in ('msa_analysis','msa_distances') and not action.__self__.load_msa():return
        import contextlib,io
        output=io.StringIO()
        with contextlib.redirect_stdout(output):result=action()
        self.page(output.getvalue(),RESULT_TITLES[name])
        return result

    def page(self,text,title='Résultats'):
        text,title=tr(text),tr(title)
        if self.simple:print(text);return
        content=[re.sub(r'\x1b\[[0-9;]*m','',line) for line in str(text).expandtabs(4).splitlines()] or ['Aucun résultat.']
        row=col=0;phrase='';typing=False
        maximum=max(map(len,content));previous=None
        with screen(mouse=True):
            while True:
                width,height=self.size();count=max(1,height-6)
                row=max(0,min(row,max(0,len(content)-count)))
                col=max(0,min(col,max(0,maximum-max(1,width-1))))
                signature=(width,height,row,col,phrase,typing)
                if signature!=previous:
                    draw([' '+' › '.join(self.route), ' '+self.c(BOLD,title),
                      f' Lignes {row+1}–{min(len(content),row+count)}/{len(content)} · décalage horizontal {col}',
                      *[line[col:] for line in content[row:row+count]],
                      (' Recherche : '+phrase+'▏') if typing else ' ↑↓ / PgUp/PgDn défiler · ←→ latéral · / chercher · n suivant',
                          ' Échap / q retour · F2 accueil'])
                    previous=signature
                key=read_key()
                if key=='tick':continue
                if key=='F2':raise HomeRequested
                if key=='eof':raise KeyboardInterrupt
                if typing:
                    if key=='esc':typing=False;continue
                    if key in ('\r','\n'):typing=False;key='n'
                    elif key in ('\x7f','\b'):phrase=phrase[:-1];continue
                    elif len(key)==1 and key.isprintable():phrase+=key;continue
                    else:continue
                if key in ('q','Q','esc'):return
                if key in ('down','j'):row+=1
                elif key in ('up','k'):row-=1
                elif key in ('pagedown',' '):row+=count
                elif key in ('pageup','b'):row-=count
                elif key=='right':col+=8
                elif key=='left':col-=8
                elif key in ('home','ctrl-home'):row=col=0
                elif key in ('end','ctrl-end'):row=len(content)-count
                elif key=='/':typing=True;phrase=''
                elif key.startswith('mouse:64;'):row-=3
                elif key.startswith('mouse:65;'):row+=3
                if key=='n' and phrase:
                    order=list(range(row+1,len(content)))+list(range(row+1))
                    found=next((i for i in order if phrase.casefold() in content[i].casefold()),None)
                    if found is not None:row=found
