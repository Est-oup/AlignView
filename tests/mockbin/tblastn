#!/usr/bin/env python3
import os,sys,re
from pathlib import Path
name=Path(sys.argv[0]).name
args=sys.argv[1:]
def opt(k,default=None):
    return args[args.index(k)+1] if k in args and args.index(k)+1<len(args) else default
if name in ('makeblastdb','diamond') and ('makedb' in args or name=='makeblastdb'):
    dest=opt('--db') or opt('-out')
    Path(dest+'.stub').write_text('mock database')
elif name in ('blastn','blastp','blastx','tblastn','tblastx','diamond'):
    dest=opt('-out') or opt('--out')
    if name=='diamond':
        fmt=[
          ['query_ASV_test','ref_HSP_1', '97.5','20','3','62','1','20','1e-40','145','3','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGKV'],
          ['query_ASV_test','ref_HSP_1', '95.0','20','50','109','1','20','1e-30','130','3','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGAV'],
          ['query_ASV_test','ref_HSP_2', '93.0','20','120','179','1','20','1e-20','100','3','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGAV']]
    else:
        fmt=[
          ['query_ASV_test','ref_HSP_1','97.5','20','3','62','1','20','1e-40','145','22.6','3','0','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGKV'],
          ['query_ASV_test','ref_HSP_1','95.0','20','50','109','1','20','1e-30','130','22.6','3','0','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGAV'],
          ['query_ASV_test','ref_HSP_2','93.0','20','120','179','1','20','1e-20','100','22.6','3','0','MPLDIGDAVKRHITATDGKV','MPLDIGDAVKRHITATDGAV']]
    Path(dest).write_text(''.join('\t'.join(r)+'\n' for r in fmt))
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
