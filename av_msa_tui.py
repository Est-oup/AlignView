"""Interactive full-screen MSA viewer using standard Python, compatible with SSH terminals.

Read-only, with no dependency on less/alv/curses. MSA coordinates are 1-based,
including gaps. The variable-site view hides invariant columns without modifying the MSA.
"""
from __future__ import annotations
import collections
import os
import select
import shutil
import sys
import termios
import tty

RESET='\x1b[0m'; BOLD='\x1b[1m'; DIM='\x1b[2m'
BLUE='\x1b[38;5;75m'; CYAN='\x1b[38;5;45m'; GREEN='\x1b[38;5;78m'
YELLOW='\x1b[38;5;221m'; RED='\x1b[38;5;203m'; GREY='\x1b[38;5;244m'
WHITE='\x1b[38;5;255m'; ACTIVE='\x1b[48;5;24m'; COLUMN='\x1b[48;5;238m'


def column_summary(seqs):
    """Non-gap consensus, variable sites, gap percentage and conservation, without implicit biological processing."""
    n=len(seqs); L=len(seqs[0].seq)
    consensus=[]; variable=[]; gaps=[]; conservation=[]
    for j in range(L):
        col=[s.seq[j].upper() for s in seqs]
        nongap=[x for x in col if x not in '-.?']
        c=collections.Counter(nongap)
        consensus.append(c.most_common(1)[0][0] if c else '-')
        variable.append(len(c)>1)
        gaps.append((n-len(nongap))/n)
        conservation.append((max(c.values())/len(nongap)) if nongap else 0)
    return ''.join(consensus), variable, gaps, conservation


def msa_static(seqs, width=88, limit=35, start=0, color=False):

    """Readable fallback when no TTY is available."""
    if not seqs:return 'No sequences'
    ncol=len(seqs[0].seq)
    start=max(0,min(max(0,ncol-1),start));stop=min(ncol,start+max(12,width))
    consensus,variable,_,_=column_summary(seqs)
    out=[f'MSA · {len(seqs)} sequences × {ncol} columns · columns {start+1}–{stop}',
         ' ' * 22+f'{start+1}'+(' ' * max(1,stop-start-8))+str(stop)]
    for s in seqs[:limit]:out.append(f'{s.id[:20]:20}  {s.seq[start:stop]}')
    out.append(f'{"CONSENSUS":20}  {consensus[start:stop]}')
    out.append(f'{str("VARIABLE SITES"):20}  '+''.join('^' if v else '·' for v in variable[start:stop]))
    if len(seqs)>limit:out.append(f'… {len(seqs)-limit} additional sequences')
    return '\n'.join(out)


class MSAViewer:
    """A bounded viewport with independent fast horizontal and vertical navigation."""
    def __init__(self,seqs,color=True):
        if not seqs or not seqs[0].seq or len({len(s.seq) for s in seqs})!=1:
            raise ValueError('The viewer requires a nonempty aligned FASTA with equal sequence lengths.')
        self.seqs=seqs;self.length=len(seqs[0].seq)
        self.consensus,self.variable,self.gaps,self.conservation=column_summary(seqs)
        alphabet=set(''.join(x.seq for x in seqs[:3]).upper())
        self.is_nucleic=alphabet.issubset(set('ACGTUNRYKMSWBDHV-.?'))
        self.visible_all=range(self.length)
        self.visible_vars=[i for i,x in enumerate(self.variable) if x]
        self.only_variable=False;self.difference=False;self.color=color and not os.getenv('NO_COLOR')
        self.col=0;self.row=0;self.xoff=0;self.yoff=0
        self.message='';self.search='';self.search_hits=[];self.result=None

    @property
    def cols(self):return self.visible_vars if self.only_variable and self.visible_vars else self.visible_all

    def _c(self,style,s):return style+s+RESET if self.color else s

    def _geom(self):
        w,h=shutil.get_terminal_size((105,35))
        name_width=max(8,min(28,w//4))
        inner=max(1,w-name_width-5)
        height=max(1,h-16)
        return w,h,name_width,inner,height

    def _clamp(self):
        cols=self.cols
        self.col=max(0,min(len(cols)-1,self.col));self.row=max(0,min(len(self.seqs)-1,self.row))
        _,_,_,width,height=self._geom()
        if self.col<self.xoff:self.xoff=self.col
        if self.col>=self.xoff+width:self.xoff=self.col-width+1
        if self.row<self.yoff:self.yoff=self.row
        if self.row>=self.yoff+height:self.yoff=self.row-height+1
        self.xoff=max(0,min(self.xoff,max(0,len(cols)-width)))
        self.yoff=max(0,min(self.yoff,max(0,len(self.seqs)-height)))

    def _paint(self,char,j,selected=False,cons=False):
        if not self.color:return char
        u=char.upper()
        if char in '-.':base=GREY
        elif u in ('*','X','?'):base=RED
        elif self.is_nucleic:base={'A':GREEN,'T':YELLOW,'U':YELLOW,'G':BLUE,'C':CYAN}.get(u,WHITE)
        else:base=GREEN if u in 'AILMFVWY' else BLUE if u in 'KRH' else RED if u in 'DE' else CYAN if u in 'NQST' else WHITE
        if self.difference and not cons and u not in ('-','.') and u!=self.consensus[j]:base=RED+BOLD
        return (ACTIVE+BOLD if selected else '')+base+char+RESET

    def _overview(self,width):
        size=max(1,width-14)
        start=int(self.xoff/max(1,len(self.cols))*size)
        stop=max(start+1,int(min(len(self.cols),self.xoff+self._geom()[3])/len(self.cols)*size))
        return ' Overview ['+''.join('█' if start<=i<stop else '─' for i in range(size))+']'

    def _frame(self):
        self._clamp();w,h,nw,width,height=self._geom();cols=self.cols
        if w<42 or h<18:
            return [' AlignView · MSA',f' Terminal {w}×{h}: resize to at least 42×18.', ' q back · navigation retained on resize']
        indices=cols[self.xoff:self.xoff+width];col=cols[self.col]
        left=lambda s:' '+s.ljust(nw)[:nw]+' │ '
        lines=[self._c(CYAN+BOLD,' AlignView / Alignment'),
               f' {len(self.seqs)} sequences × {self.length:,} columns · view {indices[0]+1:,}–{indices[-1]+1:,} · cursor {col+1:,}',
               self._c(CYAN,self._overview(w-2)),
               self._c(GREY,left('Coord. ×10')+''.join(str((j+1)//10%10) if (j+1)%10==0 else ' ' for j in indices)),
               self._c(GREY,left('Coord. units')+''.join(str((j+1)%10) for j in indices)),
               self._c(YELLOW,left('CURSOR')+''.join('▲' if i==self.col else ' ' for i in range(self.xoff,self.xoff+len(indices))))]
        for row in range(self.yoff,min(len(self.seqs),self.yoff+height)):
            record=self.seqs[row];selected=row==self.row
            name=('›' if selected else ' ')+record.id[:nw-2]
            lines.append(self._c(CYAN if selected else GREY,left(name))+''.join(self._paint(record.seq[j],j,j==col) for j in indices))
        lines.extend(['']*max(0,height-(len(self.seqs)-self.yoff)))
        lines.extend([self._c(GREY,' '+'─'*(w-2)),
                      self._c(GREEN,left('CONSENSUS'))+''.join(self._paint(self.consensus[j],j,j==col,True) for j in indices),
                      self._c(YELLOW,left('VARIABLE SITES')+''.join('▲' if self.variable[j] else '·' for j in indices)),
                      self._c(BLUE,left('CONSERVATION')+''.join('█' if self.conservation[j]>=.95 else '▓' if self.conservation[j]>=.75 else '▒' if self.conservation[j]>=.5 else '░' for j in indices)),
                      f' Row {self.row+1}/{len(self.seqs)} · {self.seqs[self.row].id} · {"variable sites only" if self.only_variable else "all columns"} · gaps {self.gaps[col]*100:.1f}%',
                      ' ←→ columns · ↑↓ sequences · [ ] screen · PgUp/PgDn rows',
                      ' g position · s sequence · / motif · v variable sites · d differences',
                      ' ? all commands · Esc/q back · F2 Home',
                      self._c(YELLOW,' '+self.message)])
        return lines

    def _render(self):
        from av_terminal import draw
        draw(self._frame())

    def _prompt(self,label):
        fd=sys.stdin.fileno()
        try:
            termios.tcsetattr(fd,termios.TCSADRAIN,self._original)
            sys.stdout.write('\x1b[?1000l\x1b[?1006l\x1b[?25h\r\n '+label+' › ');sys.stdout.flush()
            return input().strip()
        except (EOFError,KeyboardInterrupt):return ''
        finally:
            tty.setcbreak(fd,when=termios.TCSANOW)
            sys.stdout.write('\x1b[?1000h\x1b[?1006h\x1b[?25l');sys.stdout.flush()

    def _search(self,phrase):
        if not phrase:return
        motif=phrase.upper().replace(' ','')
        seq=self.consensus if motif.startswith('C:') else self.seqs[self.row].seq.upper()
        if motif.startswith('C:'):motif=motif[2:]
        self.search=phrase;self.search_hits=[];at=0
        while motif and at<len(seq):
            found=seq.find(motif,at)
            if found<0:break
            self.search_hits.append(found);at=found+1
            if len(self.search_hits)>=20000:break
        if self.search_hits:
            current=self.cols[self.col]
            target=next((j for j in self.search_hits if j>=current),self.search_hits[0])
            # An exact motif hit must stay visible, even if its column is invariant.
            self.only_variable=False
            self._goto_genomic(target)
            self.message=f'{len(self.search_hits)} occurrence(s) of « {phrase} » · n next / N previous'
        else:self.message=f'No occurrences of « {phrase} » in the current row.'

    def _goto_genomic(self,j):
        from bisect import bisect_left
        self.col=min(bisect_left(self.cols,j),len(self.cols)-1)
        self.xoff=max(0,self.col-self._geom()[3]//2)

    def _help(self):
        from av_terminal import draw,read_key
        draw([' AlignView / Alignment help','',
              ' ←→ / h l: column · Shift+←→: 10 columns',
              ' Ctrl+←→ or [ ]: one screen width',
              ' ↑↓ / k j: sequence · PgUp/PgDn: sequence page',
              ' Home/End: first/last column',
              ' Ctrl+Home/End: first/last sequence',
              ' Wheel: sequences · Shift+wheel: columns',
              ' Click the overview to jump to a region',
              ' g: absolute column (1-based) or percentage, e.g. 75%',
              ' s: go to a sequence by identifier',
              ' /: exact motif, c:motif for consensus · n / N: next/previous',
              ' v: variable sites · d: differences from consensus',
              ' w: HTML report · e: export HTML',
              ' q/Esc: back · F2: home','',
              ' Coordinates are MSA columns, including gaps.',
              ' The variable-site view hides columns without modifying the MSA.',
              ' Press any key to return.'])
        key=read_key()
        while key=='tick':key=read_key()
        return key

    def _read_key(self):
        from av_terminal import read_key
        return read_key()

    def _handle_key(self,key):
        """Apply navigation without terminal I/O, also used by regression tests."""
        import re
        width,height=self._geom()[3:]
        self.message=''
        if key in ('right','l'):self.col+=1
        elif key in ('left','h'):self.col-=1
        elif key in ('shift-right','L'):self.col+=10
        elif key in ('shift-left','H'):self.col-=10
        elif key in ('ctrl-right',']',' '):
            self.col+=max(1,width-5);self.xoff+=max(1,width-5)
        elif key in ('ctrl-left','['):
            self.col-=max(1,width-5);self.xoff-=max(1,width-5)
        elif key in ('down','j'):self.row+=1
        elif key in ('up','k'):self.row-=1
        elif key in ('pagedown','ctrl-down','J'):self.row+=height;self.yoff+=height
        elif key in ('pageup','ctrl-up','K'):self.row-=height;self.yoff-=height
        elif key=='home':self.col=0
        elif key=='end':self.col=len(self.cols)-1
        elif key=='ctrl-home':self.row=0
        elif key=='ctrl-end':self.row=len(self.seqs)-1
        elif key=='v':
            if not self.visible_vars:self.message='No variable sites in the MSA.'
            else:
                old=self.cols[self.col];self.only_variable=not self.only_variable;self._goto_genomic(old)
        elif key=='d':self.difference=not self.difference
        elif key in ('n','N'):
            if not self.search_hits:self.message='Search for a motif with / first.'
            else:
                current=self.cols[self.col]
                target=(next((i for i in self.search_hits if i>current),self.search_hits[0]) if key=='n'
                        else next((i for i in reversed(self.search_hits) if i<current),self.search_hits[-1]))
                self.only_variable=False;self._goto_genomic(target)
        elif key.startswith('mouse:'):
            match=re.fullmatch(r'mouse:(\d+);(\d+);(\d+)([Mm])',key)
            if match:
                button,x,y,event=match.groups();button=int(button);x=int(x);y=int(y)
                if button in (64,65):self.row+=-3 if button==64 else 3
                elif button in (68,69,66,67):
                    step=max(5,width//5)*(-1 if button in (68,66) else 1)
                    self.col+=step;self.xoff+=step
                elif button==0 and event=='M':
                    if y==3:
                        ratio=max(0,min(1,(x-15)/max(1,self._geom()[0]-17)))
                        self.col=round(ratio*(len(self.cols)-1));self.xoff=max(0,self.col-width//2)
                    elif 7<=y<7+height:
                        self.row=self.yoff+y-7
                        self.col=self.xoff+max(0,min(width-1,x-self._geom()[2]-5))
        self._clamp()

    def run(self):
        self.result=None
        if not (sys.stdin.isatty() and sys.stdout.isatty()):
            print(msa_static(self.seqs));return None
        from av_terminal import screen
        with screen(mouse=True) as original:
            self._original=original
            previous=None
            while True:
                signature=(self._geom(),self.col,self.row,self.xoff,self.yoff,self.only_variable,self.difference,self.message)
                if signature!=previous:
                    self._render()
                    previous=(self._geom(),self.col,self.row,self.xoff,self.yoff,self.only_variable,self.difference,self.message)
                key=self._read_key()
                if key=='tick':continue
                if key=='F2':self.result='home';break
                if key in ('q','Q','esc','eof'):break
                if key in ('w','W'):self.result='web';break
                if key in ('e','E'):self.result='export';break
                if key in ('?','F1'):
                    if self._help()=='F2':self.result='home';break
                    previous=None;continue
                if key in ('g','G'):
                    previous=None
                    value=self._prompt('Go to MSA column (1–'+str(self.length)+', or %)')
                    try:
                        j=round(float(value[:-1])/100*(self.length-1)) if value.endswith('%') else int(value)-1
                        if not 0<=j<self.length:raise ValueError
                        self._goto_genomic(j);self.message='Colonne MSA '+str(self.cols[self.col]+1)
                    except ValueError:self.message='Invalid column number.'
                elif key=='s':
                    previous=None
                    value=self._prompt('Sequence identifier (or part of it)')
                    if value:
                        matches=[i for i,r in enumerate(self.seqs) if value.casefold() in r.id.casefold()]
                        if matches:
                            self.row=next((i for i in matches if i>self.row),matches[0]);self.message=self.seqs[self.row].id
                            self.search_hits=[]
                        else:self.message='Sequence not found.'
                elif key=='/':
                    previous=None
                    self._search(self._prompt('Motif (c: for consensus)'))
                else:self._handle_key(key)
                self._clamp()
        return self.result
