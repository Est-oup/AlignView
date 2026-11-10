#!/usr/bin/env python3
"""Integration smoke tests with actual installed engines, never tests/mockbin."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_core import read_fasta
from av_setup import probe_version


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--require',nargs='*',default=[])
    parser.add_argument('--tuned',action='store_true',help='Exercise nondefault scientific settings in every engine')
    args=parser.parse_args()
    tools=('mafft','blastn','diamond','muscle','clustalo','needle','FastTree','mmseqs')
    available={name:shutil.which(name) for name in tools}
    for name,path in available.items():
        if path and '/tests/mock' in path:raise SystemExit('Remove mock tools from PATH before running this script.')
    missing=[name for name in args.require if not available.get(name)]
    if missing:raise SystemExit('Required real tools missing: '+', '.join(missing))
    results=[]
    dna=('ATGGCAGACGTTGATGCTAAAGGTGTTGTTGCTGCTGACGGTGACGTTGCTGCTGTTGATGACGCT'*3)
    protein='MPLDIGDAVKRHITATDGKVMPLDIGDAVKRHITATDGKVMPLDIGDAVKRHITATDGKV'
    with tempfile.TemporaryDirectory(prefix='alignview_real_') as tmp:
        tmp=Path(tmp)
        query=tmp/'query.fna';query.write_text('>query\n'+dna+'\n')
        refs=tmp/'refs.fna';refs.write_text('>reference\n'+dna+'\n>variant\n'+dna[:90]+'AAA'+dna[93:]+'\n')
        proteins=tmp/'proteins.faa';proteins.write_text('>p1\n'+protein+'\n>p2\n'+protein[:-3]+'AAA\n>p3\n'+protein[:-6]+'GGGAAA\n')
        for name,path in available.items():
            if not path:
                results.append(dict(tool=name,status='not installed'));continue
            out=tmp/('result_'+name)
            options=[str(proteins),'--msa','align','--aligner',name] if name in ('mafft','muscle','clustalo') else []
            if name=='blastn':options=[str(query),str(refs),'--mode','blastn','--task','blastn']
            elif name=='diamond':options=[str(proteins),str(proteins),'--query-id','p1','--mode','diamond-blastp']
            elif name=='needle':options=[str(query),str(refs),'--mode','global-dna','--pairwise','global']
            elif name=='FastTree':
                if not available['mafft']:results.append(dict(tool=name,status='MAFFT required for fixture'));continue
                options=[str(proteins),'--msa','align','--tree']
            elif name=='mmseqs':options=[str(proteins),str(proteins),'--query-id','p1','--mode','mmseqs-blastp']
            if args.tuned:
                options+=dict(
                    mafft=['--mafft','local','--mafft-op','2','--mafft-ep','.1','--mafft-maxiterate','3'],
                    blastn=['--blast-word-size','7','--blast-reward','2','--blast-penalty','-3',
                        '--blast-gapopen','5','--blast-gapextend','2','--blast-masking','no',
                        '--blast-soft-masking','true','--blast-strand','plus','--blast-max-hsps','2','--search-candidates','37'],
                    diamond=['--diamond-matrix','BLOSUM62','--diamond-gapopen','11','--diamond-gapextend','1',
                        '--diamond-masking','0','--diamond-comp-stats','0','--diamond-max-hsps','2','--search-candidates','37'],
                    muscle=['--muscle-method','super5','--muscle-perm','abc','--muscle-perturb','2','--muscle-refineiters','3'],
                    clustalo=['--clustalo-iterations','1','--clustalo-full','--clustalo-full-iter'],
                    needle=['--matrix','EDNAFULL','--needle-endopen','3','--needle-endextend','.2'],
                    FastTree=['--fasttree-model','LG','--fasttree-gamma','--fasttree-boot','50','--fasttree-seed','42'],
                    mmseqs=['--mmseqs-sensitivity','6.5','--mmseqs-mask','0','--mmseqs-comp-bias','0','--search-candidates','37'],
                )[name]
            command=[sys.executable,str(ROOT/'alignview.py'),*options,'--module-policy','never','--threads','1','--no-menu','--save','--out',str(out)]
            started=time.perf_counter()
            run=subprocess.run(command,capture_output=True,text=True,timeout=180)
            status='passed' if run.returncode==0 else 'failed'
            if status=='passed':
                if not (out/'session.json').is_file():status='failed: missing project'
                if name in ('mafft','muscle','clustalo','FastTree'):
                    records=read_fasta(out/'alignment_active.fasta')
                    if len(records)!=3 or len({len(r.seq) for r in records})!=1:status='failed: invalid MSA'
                elif name in ('blastn','diamond','mmseqs'):
                    if not json.loads((out/'session.json').read_text())['hsps']:status='failed: no hits for identical sequences'
            row=dict(tool=name,executable=path,version=probe_version(path,os.environ),status=status,
                tuned=args.tuned,elapsed_s=round(time.perf_counter()-started,4),command=command)
            if status=='passed':
                state=json.loads((out/'session.json').read_text())
                row['calculation_runs']=state.get('calculation_runs',{})
            if status!='passed':row['output']=(run.stdout+'\n'+run.stderr)[-6000:]
            results.append(row)
    encoded=json.dumps(results,indent=2)+'\n'
    if args.output:args.output.write_text(encoded,encoding='utf-8')
    print(encoded,end='')
    return int(any(row['status'].startswith('failed') for row in results))


if __name__=='__main__':raise SystemExit(main())
