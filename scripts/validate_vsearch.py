#!/usr/bin/env python3
"""Exercise four real VSEARCH workflows; no mock or biological accuracy claim."""
from pathlib import Path
import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from av_setup import probe_version


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--engine',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    engine=args.engine or shutil.which('vsearch') or ROOT/'tools/bin/vsearch'
    engine=Path(engine).resolve()
    if not engine.is_file() or 'mock' in str(engine).lower():
        raise SystemExit('A real VSEARCH executable is required; use --engine PATH.')
    version=probe_version(str(engine),os.environ.copy())
    results=[]
    with tempfile.TemporaryDirectory(prefix='alignview_real_vsearch_') as temporary:
        temporary=Path(temporary);rng=random.Random(42)
        first=''.join(rng.choice('ACGT') for _ in range(240))
        second=''.join(rng.choice('ACGT') for _ in range(240))
        records=[('q1',first),('q1_copy',first),('q1_copy2',first),('q2',second),('q2_copy',second),('mosaic',first[:120]+second[120:])]
        query=temporary/'query.fna';query.write_text(''.join('>'+name+'\n'+sequence+'\n' for name,sequence in records))
        refs=temporary/'refs.fna';refs.write_text('>reference_1\n'+first+'\n>reference_2\n'+second+'\n')
        for action in ('search','dereplicate','cluster','chimeras'):
            target=temporary/action
            command=[sys.executable,str(ROOT/'alignview.py'),str(query)]
            if action=='search':command += [str(refs),'--query-id','q1']
            command+=['--vsearch',action,'--executable','vsearch='+str(engine),'--vsearch-masking','none',
                      '--vsearch-minwordmatches','0','--no-menu','--module-policy','never',
                      '--report','--export-profile','publication','--save','--out',str(target)]
            run=subprocess.run(command,capture_output=True,text=True,timeout=60)
            if run.returncode:raise SystemExit(action+': '+run.stdout+run.stderr)
            data=json.loads((target/'vsearch'/action/'result.json').read_text())
            assert data['commands'] and all(row['exit_code']==0 for row in data['commands'])
            assert (target/'exports/vsearch'/action/'result.json').is_file()
            assert 'id="vsearch"' in (target/'report.html').read_text()
            if action=='search':
                rows=list(csv.DictReader((target/'vsearch/search/matches.tsv').open(),delimiter='\t'))
                assert any(row['query']=='q1' and row['target']=='reference_1' and float(row['id'])==100 for row in rows)
                assert not any('evalue' in row or 'bitscore' in row for row in rows)
            if action in ('dereplicate','cluster'):
                assert data['summary']['representatives']==3
                assert data['summary']['total_abundance']==6
            resume=[sys.executable,str(ROOT/'alignview.py'),'--resume',str(target),
                    '--html-out',str(temporary/(action+'_resumed.html')),'--no-menu']
            resumed=subprocess.run(resume,capture_output=True,text=True,timeout=30)
            assert resumed.returncode==0,resumed.stdout+resumed.stderr
            assert 'id="vsearch"' in (temporary/(action+'_resumed.html')).read_text()
            results.append({'action':action,'status':'passed','version':version,'engine':str(engine),
                            'summary':data['summary'],'commands':data['commands'],
                            'checks':['native outputs','publication exports','English standalone HTML','resume without rerunning engines']})
            print(action+': passed · '+version)
    if args.output:
        args.output.write_text(json.dumps(results,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    return results


if __name__=='__main__':main()
