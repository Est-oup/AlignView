#!/usr/bin/env python3
import sys
from pathlib import Path
name=Path(sys.argv[0]).name
args=sys.argv[1:]
def get(k):
    return args[args.index(k)+1] if k in args and args.index(k)+1<len(args) else None

def read(path):
    records=[]
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):records.append([line,''])
        elif records:records[-1][1]+=line.strip()
    return records

if name=='mmseqs':
    dest=Path(args[3])
    q=read(args[1]);r=read(args[2]);a=q[0][0][1:].split()[0];b=r[0][0][1:].split()[0]
    dest.write_text(f'{a}\t{b}\t95.0\t10\t1\t10\t1\t10\t1e-25\t92\tMPLDIGDAVK\tMPLDIGDAVR\n')
elif name in ('muscle','muscle5','clustalo'):
    infile=get('-align') if name!='clustalo' else get('-i')
    outfile=get('-output') if name!='clustalo' else get('-o')
    records=read(infile);length=max(len(x[1]) for x in records)
    Path(outfile).write_text(''.join(x[0]+'\n'+x[1].ljust(length,'-')+'\n' for x in records))
elif name in ('FastTree','fasttree'):
    sequences=read(args[-1]);names=[r[0][1:].split()[0] for r in sequences]
    print('('+','.join(n+':0.1' for n in names)+');')
else:raise SystemExit('Unknown mock '+name)
