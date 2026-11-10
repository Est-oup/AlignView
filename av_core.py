"""Algorithmes et modèles biologiques AlignView, bibliothèque standard uniquement."""
from __future__ import annotations
import collections
import csv
import hashlib
import math
import re
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable
from av_defaults import DEFAULT_HITS

NUCLEOTIDES = set('ACGTUNRYKMSWBDHV-?.')
DNA = set('ACGTU')
IUPAC = {'A':'A','C':'C','G':'G','T':'T','U':'[TU]','R':'[AG]','Y':'[CT]', 'S':'[GC]','W':'[AT]',
         'K':'[GT]','M':'[AC]','B':'[CGT]','D':'[AGT]','H':'[ACT]','V':'[ACG]','N':'[ACGT]'}
COMPLEMENT = str.maketrans('ACGTURYKMSWBDHVNacgturykmswbdhvn', 'TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn')
CODON = dict(zip(
  [a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],
  list('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')
))

@dataclass
class Seq:
    id: str
    description: str
    seq: str

@dataclass
class HSP:
    query: str
    subject: str
    identity: float
    align_length: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qframe: int = 0
    sframe: int = 0
    qseq: str = ''
    sseq: str = ''

    @property
    def qspan(self): return tuple(sorted((self.qstart, self.qend)))

    @property
    def sspan(self): return tuple(sorted((self.sstart, self.send)))

    @property
    def gap_columns(self): return sum(a == '-' or b == '-' for a,b in zip(self.qseq,self.sseq))

    @property
    def mismatches(self): return sum(a != b and a != '-' and b != '-' for a,b in zip(self.qseq,self.sseq))

    @property
    def matches(self): return sum(a == b and a != '-' for a,b in zip(self.qseq,self.sseq))

    def record(self):
        d=asdict(self)
        d['matches'],d['mismatches'],d['gap_columns']=self.matches,self.mismatches,self.gap_columns
        return d


def fasta_iter(path: Path | str) -> Iterable[Seq]:
    name=None; description=''; parts=[]
    import gzip
    opener=gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path,mode='rt',encoding='utf-8-sig') as handle:
        for line_no,line in enumerate(handle,1):
            line=line.strip()
            if not line or line.startswith(';'): continue
            if line.startswith('>'):
                if name is not None: yield Seq(name,description,''.join(parts).upper())
                description=line[1:].strip()
                if not description: raise ValueError(f'En-tête FASTA vide (ligne {line_no})')
                name=description.split()[0]; parts=[]
            elif name is None:
                raise ValueError(f'Contenu avant le premier en-tête FASTA (ligne {line_no})')
            else:
                if re.search(r'\s',line): line=''.join(line.split())
                parts.append(line)
    if name is not None: yield Seq(name,description,''.join(parts).upper())


def read_fasta(path: Path | str, max_sequences: int | None = None) -> list[Seq]:
    seqs=[]; seen=set()
    for record in fasta_iter(path):
        if record.id in seen: raise ValueError(f'ID FASTA en double : {record.id}')
        if not record.seq: raise ValueError(f'Séquence vide : {record.id}')
        if re.search(r'[^A-Z*?\-.]', record.seq): raise ValueError(f'Caractères invalides dans {record.id}')
        seqs.append(record);seen.add(record.id)
        if max_sequences and len(seqs)>max_sequences: raise ValueError(f'Plus de {max_sequences} séquences ; utiliser un sous-ensemble')
    if not seqs: raise ValueError(f'Aucune séquence FASTA : {path}')
    return seqs


def write_fasta(seqs: Iterable[Seq], path: Path, width: int = 80):
    with open(path,'w',encoding='utf-8') as f:
        for r in seqs:
            f.write('>'+r.description+'\n')
            for i in range(0,len(r.seq),width): f.write(r.seq[i:i+width]+'\n')


def sequence_type(seq: str) -> str:
    a=set(seq.upper())-set('-.?')
    if not a: return 'inconnu'
    if a <= NUCLEOTIDES: return 'ADN/ARN'
    if a <= set('ABCDEFGHIKLMNPQRSTVWXYZJUO*'): return 'protéine'
    return 'inconnu'


def sequence_stats(r: Seq) -> dict:
    s=r.seq; valid=[c for c in s if c not in '-.?']; total=len(s); n=len(valid)
    gc=sum(x in 'GC' for x in valid)
    count=collections.Counter(valid)
    probs=[v/n for v in count.values()] if n else []
    return dict(id=r.id,length=total,ungapped=n,gaps=total-n,kind=sequence_type(s),
                gc_pct=round(100*gc/n,3) if n and sequence_type(s)=='ADN/ARN' else None,
                n_pct=round(100*count['N']/n,3) if n else 0,
                entropy=round(-sum(p*math.log2(p) for p in probs),4),
                ambiguous=sum(c not in 'ACGTU' for c in valid) if sequence_type(s)=='ADN/ARN' else None)


def merge_intervals(intervals: Iterable[tuple[int,int]], end: int | None = None) -> list[tuple[int,int]]:
    clean=[]
    for p,q in intervals:
        a,b=sorted((int(p),int(q)))
        if end is not None: a=max(1,a);b=min(end,b)
        if b>=a: clean.append((a,b))
    merged=[]
    for a,b in sorted(clean):
        if merged and a<=merged[-1][1]+1: merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
        else: merged.append((a,b))
    return merged


def interval_length(intervals): return sum(b-a+1 for a,b in intervals)


def missing_intervals(intervals, length):
    if length<=0: return []
    missing=[]; pos=1
    for a,b in merge_intervals(intervals,length):
        if a>pos: missing.append((pos,a-1))
        pos=b+1
    if pos<=length: missing.append((pos,length))
    return missing


def parse_blast_tab(path, engine='blast') -> list[HSP]:
    """NCBI: 15 colonnes ; DIAMOND: 13, sframe absent. Chaque ligne = un HSP."""
    result=[]
    with open(path,encoding='utf-8') as f:
        for index,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            a=line.rstrip('\n').split('\t')
            expected=15 if engine=='blast' else 13
            if len(a)!=expected: raise ValueError(f'HSP ligne {index} : {len(a)} champs au lieu de {expected}')
            try:
                if engine=='blast':
                    result.append(HSP(a[0],a[1],float(a[2]),int(a[3]),int(a[4]),int(a[5]),int(a[6]),int(a[7]),float(a[8]),float(a[9]),int(float(a[11])),int(float(a[12])),a[13],a[14]))
                else:
                    result.append(HSP(a[0],a[1],float(a[2]),int(a[3]),int(a[4]),int(a[5]),int(a[6]),int(a[7]),float(a[8]),float(a[9]),int(float(a[10])),0,a[11],a[12]))
            except (ValueError,IndexError) as e: raise ValueError(f'Parsing HSP ligne {index} : {e}') from e
    return sorted(result,key=lambda h:(-h.bitscore,h.evalue,h.subject))


def parse_mmseqs_tab(path) -> list[HSP]:
    """MMseqs easy-search, pident est en pourcentage, positions 1-based."""
    rows=[]
    with open(path,encoding='utf-8') as stream:
        for index,line in enumerate(stream,1):
            if not line.strip() or line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=12:raise ValueError(f'MMseqs2 ligne {index} : 12 champs attendus, reçus {len(fields)}')
            try:
                q,t,pid,alen,qs,qe,ts,te,e,bits,qaln,taln=fields
                rows.append(HSP(q,t,float(pid),int(alen),int(qs),int(qe),int(ts),int(te),float(e),float(bits),0,0,qaln,taln))
            except ValueError as exc:raise ValueError(f'MMseqs2 ligne {index} : {exc}') from exc
    return sorted(rows,key=lambda h:(-h.bitscore,h.evalue,h.subject))

def group_hsps(hsps: list[HSP]) -> dict[str,list[HSP]]:
    by=collections.defaultdict(list)
    for h in hsps: by[h.subject].append(h)
    return dict(by)


def select_hits(hsps,query_len,identity=0,coverage=0,max_hits=DEFAULT_HITS,min_bits=0):
    filtered=[h for h in hsps if h.identity>=identity and h.bitscore>=min_bits]
    groups=group_hsps(filtered); selected={}
    for name,rows in groups.items():
        spans=merge_intervals([h.qspan for h in rows],query_len)
        cov=interval_length(spans)/query_len*100 if query_len else 0
        if cov>=coverage: selected[name]=rows
    keys=sorted(selected,key=lambda k: (-selected[k][0].bitscore,selected[k][0].evalue,k))[:max_hits]
    return {k:selected[k] for k in keys}


def coverage_data(hsps, query_len, ref_lengths:dict):
    out=[]
    for subj,grp in group_hsps(hsps).items():
        qs=merge_intervals([h.qspan for h in grp],query_len)
        slen=ref_lengths.get(subj,0)
        ss=merge_intervals([h.sspan for h in grp],slen) if slen else []
        out.append(dict(subject=subj,hsps=len(grp),query_segments=qs,
            query_uncovered=missing_intervals(qs,query_len),query_coverage=interval_length(qs)/max(1,query_len)*100,
            reference_length=slen,reference_segments=ss,reference_uncovered=missing_intervals(ss,slen),
            reference_coverage=interval_length(ss)/slen*100 if slen else None,
            best_identity=grp[0].identity,best_bitscore=grp[0].bitscore,best_evalue=grp[0].evalue,
            frames=sorted(set(h.qframe for h in grp if h.qframe))))
    return out


def reverse_complement(seq): return seq.translate(COMPLEMENT)[::-1].replace('U','T')


def translate_nt(seq, frame=1):
    if frame not in (1,2,3,-1,-2,-3): raise ValueError('Cadre de traduction invalide')
    seq=seq.upper().replace('U','T').replace('-','').replace('.','')
    if frame<0: seq=reverse_complement(seq)
    offset=abs(frame)-1
    return ''.join(CODON.get(seq[i:i+3],'X') for i in range(offset,len(seq)-2,3))


def orfs(seq:Seq,min_aa=60,include_partial=False):
    if sequence_type(seq.seq)!='ADN/ARN': raise ValueError('Détection d’ORF : fichier nucléotidique requis')
    found=[]; n=len(seq.seq.replace('-',''))
    for strand in (1,-1):
        dna=seq.seq.replace('-','').replace('U','T') if strand==1 else reverse_complement(seq.seq.replace('-',''))
        for frame in range(3):
            start=None
            for pos in range(frame,len(dna)-2,3):
                codon=dna[pos:pos+3]
                if start is None and codon=='ATG': start=pos
                if codon in ('TAA','TAG','TGA'):
                    if start is not None:
                        aa=translate_nt(dna[start:pos],1)
                        if len(aa)>=min_aa:
                            left,right=(start+1,pos+3) if strand==1 else (n-pos-2,n-start)
                            found.append(dict(seqid=seq.id,start=left,end=right,frame=(frame+1)*strand,
                                strand='+' if strand==1 else '-',length_aa=len(aa),complete=True,protein=aa))
                    start=None
            if include_partial and start is not None:
                end=len(dna)-(len(dna)-start)%3
                aa=translate_nt(dna[start:end],1)
                if len(aa)>=min_aa:
                    left,right=(start+1,end) if strand==1 else (n-end+1,n-start)
                    found.append(dict(seqid=seq.id,start=left,end=right,frame=(frame+1)*strand,
                        strand='+' if strand==1 else '-',length_aa=len(aa),complete=False,protein=aa))
    return sorted(found,key=lambda o:(-o['length_aa'],o['start']))


def find_motifs(seq:Seq,motif:str,both_strands=True,regex=False):
    if len(motif)>200: raise ValueError('Motif trop long (200 caractères max)')
    if regex: expr=motif
    else:
        try: expr=''.join(IUPAC[c] for c in motif.upper())
        except KeyError as e: raise ValueError(f'Code IUPAC inconnu : {e}')
    try: pattern=re.compile(f'(?=({expr}))',re.I)
    except re.error as e: raise ValueError(f'Expression régulière invalide : {e}')
    out=[]; n=len(seq.seq)
    for strand,s in [('+',seq.seq),('-',reverse_complement(seq.seq))] if both_strands and sequence_type(seq.seq)=='ADN/ARN' else [('+',seq.seq)]:
        for m in pattern.finditer(s):
            length=len(m.group(1))
            if length==0: continue
            start=m.start()+1
            a,b=(start,start+length-1) if strand=='+' else (n-(start+length-1)+1,n-start+1)
            out.append(dict(seqid=seq.id,start=a,end=b,strand=strand,matched=m.group(1)))
            if len(out)>=100000: return out
    return out


def sliding_scores(seq:str,window=101,step=25,kind='gc'):
    if window<1 or step<1: raise ValueError('Fenêtre et pas doivent être positifs')
    scores=[]
    if len(seq)<window: window=len(seq)
    for i in range(0,len(seq)-window+1,step):
        part=seq[i:i+window]; c=collections.Counter(part)
        if kind=='gc':
            valid=sum(c[x] for x in 'ACGTU')
            val=(c['G']+c['C'])/valid*100 if valid else None
        elif kind=='entropy':
            valid=sum(v for k,v in c.items() if k not in '-.?')
            val=-sum((v/valid)*math.log2(v/valid) for k,v in c.items() if valid and k not in '-.?') if valid else None
        elif kind=='gap': val=100*sum(c[x] for x in '-.')/len(part)
        else: raise ValueError(f'Métrique inconnue : {kind}')
        scores.append(dict(start=i+1,end=i+window,center=i+1+(window-1)/2,value=round(val,5) if val is not None else None))
    return scores


def msa_metrics(seqs:list[Seq]):
    if not seqs: raise ValueError('MSA vide')
    size=len(seqs[0].seq)
    if any(len(r.seq)!=size for r in seqs): raise ValueError('Longueurs différentes : aligner avec MAFFT')
    consensus=[]; positions=[];variable=0;info=0;gapped=0;invariant=0
    for j in range(size):
        counter=collections.Counter(r.seq[j] for r in seqs)
        residues={k:v for k,v in counter.items() if k not in '-.?'}
        n=sum(residues.values()); gap=sum(v for k,v in counter.items() if k in '-.?')
        best=max(residues,key=lambda k:(residues[k],k)) if residues else '-'
        cons=best if residues else '-'; consensus.append(cons)
        f=residues[best]/n if n else 0
        entropy=-sum((c/n)*math.log2(c/n) for c in residues.values()) if n else 0
        variability=len(residues)>1; variable+=variability
        parsimony=sum(c>=2 for c in residues.values())>=2; info+=parsimony
        gapped+=gap>0; invariant+=(len(residues)==1 and gap==0)
        positions.append(dict(position=j+1,consensus=cons,conservation=round(100*f,3),
            gap_pct=round(100*gap/len(seqs),3),entropy=round(entropy,5),variable=variability,
            parsimony_informative=parsimony))
    return dict(n=len(seqs),length=size,variable=variable,parsimony_informative=info,
                gapped=gapped,invariant=invariant,consensus=''.join(consensus),positions=positions)


def distance_matrix(seqs:list[Seq]):
    if len(seqs)>250: raise ValueError('Matrice limitée à 250 séquences (coût quadratique)')
    if len({len(s.seq) for s in seqs})>1: raise ValueError('Alignement multiple requis')
    out=[]
    for i,a in enumerate(seqs):
        row=[]
        for j,b in enumerate(seqs):
            if j>i: row.append(None);continue
            good=[(u,v) for u,v in zip(a.seq,b.seq) if u not in '-.?' and v not in '-.?']
            row.append(sum(u!=v for u,v in good)/len(good) if good else None)
        out.append(row)
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)): out[i][j]=out[j][i]
    return out


def dotplot(seq_a,seq_b,k=9,limit=60000, max_repeat=30):
    """Points exacts par k-mers, limités pour éviter les répétitions pathologiques."""
    if not 2<=k<=60: raise ValueError('k compris entre 2 et 60')
    index=collections.defaultdict(list)
    for j in range(len(seq_b)-k+1):
        w=seq_b[j:j+k]
        if '-' not in w and '.' not in w: index[w].append(j+1)
    index={key:pos for key,pos in index.items() if len(pos)<=max_repeat}
    points=[]; truncated=False
    for i in range(len(seq_a)-k+1):
        for j in index.get(seq_a[i:i+k],[]):
            points.append([i+1,j])
            if len(points)>=limit: truncated=True;break
        if truncated: break
    return dict(points=points,k=k,truncated=truncated,xlen=len(seq_a),ylen=len(seq_b))


def sha256(path:Path):
    digest=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()


def export_table(path:Path,rows:list[dict],fieldnames=None,sep='\t'):
    if fieldnames is None:
        fieldnames=list(rows[0]) if rows else []
    with open(path,'w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames,delimiter=sep,extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)
