"""AlignView biological algorithms and models, using only the standard library."""
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

    def __post_init__(self):
        """Reject malformed evidence before it can enter a result or project."""
        if not isinstance(self.query,str) or not self.query or not isinstance(self.subject,str) or not self.subject:
            raise ValueError('HSP query and reference identifiers must be nonempty strings.')
        if any(c.isspace() for identifier in (self.query,self.subject) for c in identifier):
            raise ValueError('HSP identifiers cannot contain whitespace.')
        for name in ('align_length','qstart','qend','sstart','send'):
            value=getattr(self,name)
            if type(value) is not int or value<1:
                raise ValueError('HSP '+name+' must be a positive integer.')
        for name in ('identity','evalue','bitscore'):
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
                raise ValueError('HSP '+name+' must be a finite number.')
        if not 0<=self.identity<=100:raise ValueError('HSP identity must be between 0 and 100.')
        if self.evalue<0:raise ValueError('HSP E-value must be nonnegative.')
        for name in ('qframe','sframe'):
            value=getattr(self,name)
            if type(value) is not int or value not in range(-3,4):
                raise ValueError('HSP '+name+' must be an integer from -3 to 3.')
        if not isinstance(self.qseq,str) or not isinstance(self.sseq,str):
            raise ValueError('HSP aligned rows must be strings.')
        if bool(self.qseq)!=bool(self.sseq):
            raise ValueError('HSP aligned query and reference rows must be provided together.')
        if self.qseq and (len(self.qseq)!=self.align_length or len(self.sseq)!=self.align_length):
            raise ValueError('HSP aligned rows must have the reported alignment length.')
        if self.qseq and any(not re.fullmatch(r'[A-Za-z*?.-]+',row) for row in (self.qseq,self.sseq)):
            raise ValueError('HSP aligned rows contain invalid symbols.')

    @property
    def qspan(self): return tuple(sorted((self.qstart, self.qend)))

    @property
    def sspan(self): return tuple(sorted((self.sstart, self.send)))

    @property
    def gap_columns(self): return sum(a in '-.' or b in '-.' for a,b in zip(self.qseq,self.sseq))

    @property
    def mismatches(self): return sum(a.upper()!=b.upper() and a not in '-.' and b not in '-.' for a,b in zip(self.qseq,self.sseq))

    @property
    def matches(self): return sum(a.upper()==b.upper() and a not in '-.' for a,b in zip(self.qseq,self.sseq))

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
                if not description: raise ValueError(f'Empty FASTA header (line {line_no})')
                name=description.split()[0]; parts=[]
            elif name is None:
                raise ValueError(f'Content before the first FASTA header (line {line_no})')
            else:
                if re.search(r'\s',line): line=''.join(line.split())
                parts.append(line)
    if name is not None: yield Seq(name,description,''.join(parts).upper())


def read_fasta(path: Path | str, max_sequences: int | None = None) -> list[Seq]:
    seqs=[]; seen=set()
    for record in fasta_iter(path):
        if record.id in seen: raise ValueError(f'Duplicate FASTA ID: {record.id}')
        if not record.seq: raise ValueError(f'Empty sequence: {record.id}')
        if re.search(r'[^A-Z*?\-.]', record.seq): raise ValueError(f'Invalid characters in {record.id}')
        seqs.append(record);seen.add(record.id)
        if max_sequences and len(seqs)>max_sequences: raise ValueError(f'More than {max_sequences} sequences; use a subset')
    if not seqs: raise ValueError(f'No FASTA sequences: {path}')
    return seqs


def write_fasta(seqs: Iterable[Seq], path: Path, width: int = 80):
    if type(width) is not int or width<1:raise ValueError('FASTA line width must be a positive integer.')
    with open(path,'w',encoding='utf-8') as f:
        for r in seqs:
            if not r.id or any(c.isspace() for c in r.id):
                raise ValueError('FASTA identifiers must be nonempty and cannot contain whitespace.')
            if '\n' in r.description or '\r' in r.description:
                raise ValueError('FASTA descriptions cannot contain line breaks.')
            description=r.description.strip()
            header=description if description and description.split()[0]==r.id else r.id+(' '+description if description else '')
            f.write('>'+header+'\n')
            for i in range(0,len(r.seq),width): f.write(r.seq[i:i+width]+'\n')


def sequence_type(seq: str) -> str:
    a=set(seq.upper())-set('-.?')
    if not a: return 'unknown'
    if a <= NUCLEOTIDES: return 'DNA/RNA'
    if a <= set('ABCDEFGHIKLMNPQRSTVWXYZJUO*'): return 'protein'
    return 'unknown'


def sequence_stats(r: Seq) -> dict:
    s=r.seq.upper(); valid=[c for c in s if c not in '-.?']; total=len(s); n=len(valid)
    gc=sum(x in 'GC' for x in valid)
    count=collections.Counter(valid)
    probs=[v/n for v in count.values()] if n else []
    return dict(id=r.id,length=total,ungapped=n,gaps=total-n,kind=sequence_type(s),
                gc_pct=round(100*gc/n,3) if n and sequence_type(s)=='DNA/RNA' else None,
                n_pct=round(100*count['N']/n,3) if n else 0,
                entropy=round(-sum(p*math.log2(p) for p in probs),4),
                ambiguous=sum(c not in 'ACGTU' for c in valid) if sequence_type(s)=='DNA/RNA' else None)


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
    """NCBI: 15 columns; DIAMOND: 13, without sframe. Each line represents one HSP."""
    def frame(value):
        number=float(value)
        if not math.isfinite(number) or not number.is_integer():
            raise ValueError('reading frames must be integers')
        return int(number)
    result=[]
    with open(path,encoding='utf-8') as f:
        for index,line in enumerate(f,1):
            if not line.strip() or line.startswith('#'): continue
            a=line.rstrip('\n').split('\t')
            expected=15 if engine=='blast' else 13
            if len(a)!=expected: raise ValueError(f'HSP row {index}: {len(a)} fields instead of {expected}')
            try:
                if engine=='blast':
                    result.append(HSP(a[0],a[1],float(a[2]),int(a[3]),int(a[4]),int(a[5]),int(a[6]),int(a[7]),float(a[8]),float(a[9]),frame(a[11]),frame(a[12]),a[13],a[14]))
                else:
                    result.append(HSP(a[0],a[1],float(a[2]),int(a[3]),int(a[4]),int(a[5]),int(a[6]),int(a[7]),float(a[8]),float(a[9]),frame(a[10]),0,a[11],a[12]))
            except (ValueError,IndexError) as e: raise ValueError(f'Parsing HSP row {index}: {e}') from e
    return sorted(result,key=lambda h:(-h.bitscore,h.evalue,h.subject))


def parse_mmseqs_tab(path) -> list[HSP]:
    """MMseqs easy-search, pident is a percentage; positions are 1-based."""
    rows=[]
    with open(path,encoding='utf-8') as stream:
        for index,line in enumerate(stream,1):
            if not line.strip() or line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=12:raise ValueError(f'MMseqs2 row {index}: expected 12 fields, received {len(fields)}')
            try:
                q,t,pid,alen,qs,qe,ts,te,e,bits,qaln,taln=fields
                rows.append(HSP(q,t,float(pid),int(alen),int(qs),int(qe),int(ts),int(te),float(e),float(bits),0,0,qaln,taln))
            except ValueError as exc:raise ValueError(f'MMseqs2 row {index}: {exc}') from exc
    return sorted(rows,key=lambda h:(-h.bitscore,h.evalue,h.subject))

def group_hsps(hsps: list[HSP]) -> dict[str,list[HSP]]:
    by=collections.defaultdict(list)
    for h in hsps: by[h.subject].append(h)
    return dict(by)


def rank_hsps(hsps):
    """Rank local evidence deterministically, independently of input ordering."""
    return sorted(hsps,key=lambda h:(-h.bitscore,h.evalue,-h.identity,h.subject,
        h.qspan,h.sspan,h.qstart,h.sstart,h.qframe,h.sframe,h.qseq,h.sseq))


def normalize_search_record(record: Seq, kind=None) -> Seq:
    """Use ungapped biological coordinates without changing immutable input data.

    Gap symbols are removed. Unknown positions remain positions (N or X), and
    RNA U is represented as T for nucleotide search engines.
    """
    kind=kind or sequence_type(record.seq)
    if kind not in ('DNA/RNA','protein'):
        raise ValueError('Search input must have a nucleotide or protein alphabet: '+record.id)
    sequence=record.seq.upper().replace('-','').replace('.','')
    sequence=sequence.replace('?', 'N' if kind=='DNA/RNA' else 'X')
    if kind=='DNA/RNA':sequence=sequence.replace('U','T')
    if not sequence:raise ValueError('No biological residues remain after removing gaps: '+record.id)
    return Seq(record.id,record.description,sequence)


def validate_search_results(hsps, query: Seq, references: list[Seq], mode: str):
    """Validate a fresh engine result before replacing a successful search.

    Coordinates are 1-based in the native sequences; translated HSP rows use
    amino acids. Partial-codon/frameshift alignment is not enabled by AlignView.
    Empty searches are valid, but incomplete or foreign rows are not.
    """
    from av_units import units_for
    units=units_for(mode,sequence_type(query.seq))
    refmap={record.id:record for record in references}
    if len(refmap)!=len(references):raise ValueError('Reference identifiers must be unique.')
    result=list(hsps)
    for hsp in result:
        if not isinstance(hsp,HSP):raise ValueError('Search returned a row that is not an HSP.')
        hsp.__post_init__()
        if hsp.query!=query.id:raise ValueError('Search returned an unexpected query: '+hsp.query)
        if hsp.subject not in refmap:raise ValueError('Search returned an unknown reference: '+hsp.subject)
        if hsp.qspan[1]>len(query.seq) or hsp.sspan[1]>len(refmap[hsp.subject].seq):
            raise ValueError('Search coordinates exceed the input sequence length: '+hsp.subject)
        if not hsp.qseq:raise ValueError('Search returned HSPs without aligned query/reference rows.')
        for name,row,span,translated,frame,start,end,length in (
            ('query',hsp.qseq,hsp.qspan,units.translated_query,hsp.qframe,hsp.qstart,hsp.qend,len(query.seq)),
            ('reference',hsp.sseq,hsp.sspan,units.translated_subject,hsp.sframe,hsp.sstart,hsp.send,len(refmap[hsp.subject].seq))):
            residues=sum(c not in '-.' for c in row)
            expected=residues*(3 if translated else 1)
            if span[1]-span[0]+1!=expected:
                raise ValueError('Search '+name+' coordinates disagree with its aligned residues: '+hsp.subject)
            if translated and frame==0:
                raise ValueError('Translated search is missing the '+name+' reading frame: '+hsp.subject)
            if translated:
                phase=(start-1)%3+1 if frame>0 else (length-start)%3+1
                if abs(frame)!=phase or (frame>0 and start>end) or (frame<0 and start<end):
                    raise ValueError('Search '+name+' reading frame disagrees with its native coordinates: '+hsp.subject)
    return rank_hsps(result)


def validate_alignment_results(aligned, expected, *, nucleotide=False):
    """Accept a complete reordered MSA, never missing or changed sequences.

    Gap placement may change; biological residues may not. RNA U and DNA T are
    equivalent only when the caller explicitly declares nucleotide input.
    Input order is restored so the query and selected references remain stable.
    """
    aligned=list(aligned);expected=list(expected)
    if not expected or not aligned:raise ValueError('The alignment contains no sequences.')
    expected_ids=[record.id for record in expected]
    observed_ids=[record.id for record in aligned]
    if len(set(expected_ids))!=len(expected_ids):raise ValueError('Alignment input identifiers must be unique.')
    if len(set(observed_ids))!=len(observed_ids):raise ValueError('The alignment contains duplicate identifiers.')
    if set(observed_ids)!=set(expected_ids):
        missing=sorted(set(expected_ids)-set(observed_ids))
        extra=sorted(set(observed_ids)-set(expected_ids))
        detail=[]
        if missing:detail.append('missing: '+', '.join(missing[:5]))
        if extra:detail.append('unexpected: '+', '.join(extra[:5]))
        raise ValueError('The alignment changed its sequence set ('+'; '.join(detail)+').')
    lengths={len(record.seq) for record in aligned}
    if len(lengths)!=1 or 0 in lengths:raise ValueError('The alignment must have equal, nonempty sequence lengths.')
    def residues(sequence):
        value=sequence.upper().replace('-','').replace('.','')
        return value.replace('U','T') if nucleotide else value
    outmap={record.id:record for record in aligned}
    for record in expected:
        if residues(record.seq)!=residues(outmap[record.id].seq):
            raise ValueError('The alignment changed biological residues for '+record.id+'.')
    return [outmap[identifier] for identifier in expected_ids]


def select_hits(hsps,query_len,identity=0,coverage=0,max_hits=DEFAULT_HITS,min_bits=0):
    filtered=rank_hsps(h for h in hsps if h.identity>=identity and h.bitscore>=min_bits)
    groups=group_hsps(filtered); selected={}
    for name,rows in groups.items():
        spans=merge_intervals([h.qspan for h in rows],query_len)
        cov=interval_length(spans)/query_len*100 if query_len else 0
        if cov>=coverage: selected[name]=rows
    keys=sorted(selected,key=lambda k: (-selected[k][0].bitscore,selected[k][0].evalue,k))[:max_hits]
    return {k:selected[k] for k in keys}


def coverage_data(hsps, query_len, ref_lengths:dict):
    out=[]
    for subj,grp in group_hsps(rank_hsps(hsps)).items():
        qs=merge_intervals([h.qspan for h in grp],query_len)
        slen=ref_lengths.get(subj,0)
        ss=merge_intervals([h.sspan for h in grp],slen) if slen else []
        out.append(dict(subject=subj,hsps=len(grp),query_segments=qs,
            query_uncovered=missing_intervals(qs,query_len),query_coverage=interval_length(qs)/max(1,query_len)*100,
            reference_length=slen,reference_segments=ss,reference_uncovered=missing_intervals(ss,slen),
            reference_coverage=interval_length(ss)/slen*100 if slen else None,
            best_identity=grp[0].identity,best_bitscore=grp[0].bitscore,best_evalue=grp[0].evalue,
            max_identity=max(h.identity for h in grp),
            query_aligned_length=interval_length(qs),reference_aligned_length=interval_length(ss) if slen else None,
            query_strands=sorted(set('-' if h.qframe<0 or (not h.qframe and h.qstart>h.qend) else '+' for h in grp)),
            reference_strands=sorted(set('-' if h.sframe<0 or (not h.sframe and h.sstart>h.send) else '+' for h in grp)),
            query_frames=sorted(set(h.qframe for h in grp if h.qframe)),
            reference_frames=sorted(set(h.sframe for h in grp if h.sframe)),
            frames=sorted(set(h.qframe for h in grp if h.qframe))))
    return out


def reverse_complement(seq): return seq.translate(COMPLEMENT)[::-1].replace('U','T')


def translate_nt(seq, frame=1, table=1):
    from av_genetics import code_table, translate_codon
    code_table(table)
    if frame not in (1,2,3,-1,-2,-3): raise ValueError('Invalid translation frame')
    seq=seq.upper().replace('U','T').replace('-','').replace('.','')
    if frame<0: seq=reverse_complement(seq)
    offset=abs(frame)-1
    return ''.join(translate_codon(seq[i:i+3],table) for i in range(offset,len(seq)-2,3))


def orfs(seq:Seq,min_aa=60,include_partial=False,table=1,starts='atg'):
    from av_genetics import stop_codons, start_codons, context_dependent
    stops=stop_codons(table)
    if context_dependent(table):raise ValueError(f'Table {table}: context-dependent stops; simple ORF detection requires a CDS context that is unavailable here.')
    if starts not in ('atg','table'):raise ValueError('Invalid ORF start policy.')
    initiators={'ATG'} if starts=='atg' else start_codons(table)
    if sequence_type(seq.seq)!='DNA/RNA': raise ValueError('ORF detection requires nucleotide input')
    found=[]; normalized=seq.seq.upper().replace('-','').replace('.','').replace('U','T');n=len(normalized)
    for strand in (1,-1):
        dna=normalized if strand==1 else reverse_complement(normalized)
        for frame in range(3):
            start=None
            for pos in range(frame,len(dna)-2,3):
                codon=dna[pos:pos+3]
                if start is None and codon in initiators: start=pos
                if codon in stops:
                    if start is not None:
                        aa=translate_nt(dna[start:pos],1,table)
                        if starts=='table' and aa:aa='M'+aa[1:]
                        if len(aa)>=min_aa:
                            left,right=(start+1,pos+3) if strand==1 else (n-pos-2,n-start)
                            found.append(dict(seqid=seq.id,start=left,end=right,frame=(frame+1)*strand,
                                strand='+' if strand==1 else '-',length_aa=len(aa),complete=True,protein=aa))
                    start=None
            if include_partial and start is not None:
                end=len(dna)-(len(dna)-start)%3
                aa=translate_nt(dna[start:end],1,table)
                if starts=='table' and aa:aa='M'+aa[1:]
                if len(aa)>=min_aa:
                    left,right=(start+1,end) if strand==1 else (n-end+1,n-start)
                    found.append(dict(seqid=seq.id,start=left,end=right,frame=(frame+1)*strand,
                        strand='+' if strand==1 else '-',length_aa=len(aa),complete=False,protein=aa))
    return sorted(found,key=lambda o:(-o['length_aa'],o['start']))


def find_motifs(seq:Seq,motif:str,both_strands=True,regex=False):
    if len(motif)>200: raise ValueError('Motif too long (maximum 200 characters)')
    if regex: expr=motif
    else:
        try: expr=''.join(IUPAC[c] for c in motif.upper())
        except KeyError as e: raise ValueError(f'Unknown IUPAC code: {e}')
    try: pattern=re.compile(f'(?=({expr}))',re.I)
    except re.error as e: raise ValueError(f'Invalid regular expression: {e}')
    out=[]; n=len(seq.seq)
    for strand,s in [('+',seq.seq),('-',reverse_complement(seq.seq))] if both_strands and sequence_type(seq.seq)=='DNA/RNA' else [('+',seq.seq)]:
        for m in pattern.finditer(s):
            length=len(m.group(1))
            if length==0: continue
            start=m.start()+1
            a,b=(start,start+length-1) if strand=='+' else (n-(start+length-1)+1,n-start+1)
            out.append(dict(seqid=seq.id,start=a,end=b,strand=strand,matched=m.group(1)))
            if len(out)>=100000: return out
    return out


def sliding_scores(seq:str,window=101,step=25,kind='gc'):
    if window<1 or step<1: raise ValueError('Window and step must be positive')
    if not seq:return []
    seq=seq.upper()
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
        else: raise ValueError(f'Unknown metric: {kind}')
        scores.append(dict(start=i+1,end=i+window,center=i+1+(window-1)/2,value=round(val,5) if val is not None else None))
    return scores


def msa_metrics(seqs:list[Seq]):
    if not seqs or not seqs[0].seq: raise ValueError('The alignment must contain nonempty sequences.')
    size=len(seqs[0].seq)
    if any(len(r.seq)!=size for r in seqs): raise ValueError('Different lengths: align with MAFFT')
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
    if len(seqs)>250: raise ValueError('Matrix limited to 250 sequences (quadratic cost)')
    if len({len(s.seq) for s in seqs})>1: raise ValueError('Multiple alignment required')
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
    """Exact k-mer points, capped to avoid pathological repetitions."""
    if not 2<=k<=60: raise ValueError('k must be between 2 and 60')
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
