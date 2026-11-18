#!/usr/bin/env python3
import os,sys,re
from pathlib import Path
name=Path(sys.argv[0]).name
args=sys.argv[1:]
def opt(k,default=None):
    return args[args.index(k)+1] if k in args and args.index(k)+1<len(args) else default
def fasta_records(path):
    records=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):records.append([line[1:].split()[0], ''])
        elif records:records[-1][1]+=line.strip().upper().replace('-','').replace('.','')
    return records

if any(flag in args for flag in ('--version','-version','-h','-help')):
    print(name+' 1.0-test');raise SystemExit(0)
if name in ('makeblastdb','diamond') and ('makedb' in args or name=='makeblastdb'):
    dest=opt('--db') or opt('-out')
    Path(dest+'.stub').write_text(opt('--in') or opt('-in'))
elif name in ('blastn','blastp','blastx','tblastn','tblastx','diamond'):
    dest=opt('-out') or opt('--out')
    queries=fasta_records(opt('-query') or opt('--query'))
    database=opt('-db') or opt('--db')
    references=fasta_records(Path(database+'.stub').read_text())
    qid,qsequence=queries[0]
    mode=args[0] if name=='diamond' else name
    qfactor=3 if mode in ('blastx','tblastx') else 1
    sfactor=3 if mode in ('tblastn','tblastx') else 1
    codons=dict(zip([a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
    def alignment_row(sequence,start,length,factor):
        region=sequence[start-1:start-1+length*factor]
        return region if factor==1 else ''.join(codons.get(region[i:i+3],'X') for i in range(0,len(region),3))
    rows=[]
    chosen=[references[0],references[0],references[min(1,len(references)-1)]]
    for index,(rid,rsequence) in enumerate(chosen):
        start=([3,51,120][index] if qfactor==3 and len(qsequence)>=180 else [1,3,5][index])
        start=min(start,max(1,len(qsequence)-qfactor+1))
        length=min(20,(len(qsequence)-start+1)//qfactor,len(rsequence)//sfactor)
        if not length:continue
        qrow=alignment_row(qsequence,start,length,qfactor)
        srow=alignment_row(rsequence,1,length,sfactor)
        identity=100*sum(a==b for a,b in zip(qrow,srow))/length
        qframe=(start-1)%3+1 if qfactor==3 else 0
        sframe=1 if sfactor==3 else 0
        values=[qid,rid,f'{identity:.3f}',str(length),str(start),str(start+length*qfactor-1),
                '1',str(length*sfactor),['1e-40','1e-30','1e-20'][index],str([145,130,100][index])]
        if name=='diamond':values += [str(qframe),qrow,srow]
        else:values += [str(100*length*qfactor/len(qsequence)),str(qframe),str(sframe),qrow,srow]
        rows.append(values)
    Path(dest).write_text(''.join('\t'.join(row)+'\n' for row in rows))
elif name=='mafft':
    def fasta(p):
        d=[];nm=None;s=''
        for line in open(p):
            if line.startswith('>'):
                if nm is not None:d.append((nm,s))
                nm=line.strip();s=''
            else:s+=line.strip()
        if nm is not None:d.append((nm,s))
        return d
    data=fasta(args[-1]);
    if '--addfragments' in args:data+=fasta(opt('--addfragments'))
    size=max(len(s) for n,s in data)
    print(''.join(n+'\n'+s.ljust(size,'-')+'\n' for n,s in data),end='')
elif name in ('needle','water'):
    dest=Path(opt('-outfile'))
    if '-aformat3' in args:
        dest.write_text('>QUERY\nMKAA-A\n>REF\nMKAAEA\n')
    else:dest.write_text('Mock pairwise Q MKAA-A\n              |||| |\nReference MKAAEA\n')
elif name=='alv':
    print('MOCK ALV VIEW')
else:
    print(f'Unknown mock tool: {name}',file=sys.stderr)
    sys.exit(1)
