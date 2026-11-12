"""Reproducible execution of bioinformatics tools and HPC modules."""
from __future__ import annotations
import os
import shlex
import shutil
import sys
import subprocess
import time
import signal
from pathlib import Path
from dataclasses import dataclass
from av_parameters import validate
from av_core import Seq, HSP, fasta_iter, write_fasta, parse_blast_tab, translate_nt

GROUPS = {
    'BLAST': (['blastn','blastp','blastx','tblastn','tblastx','makeblastdb'],['blast','ncbi-blast','blast+','ncbi-blast+']),
    'DIAMOND': (['diamond'],['diamond']),
    'MAFFT': (['mafft'],['mafft']),
    'EMBOSS': (['needle','water','transeq'],['emboss']),
    'ALV': (['alv'],['alv']),
    'MMSEQS': (['mmseqs'],['mmseqs2','mmseqs']),
    'MUSCLE': (['muscle','muscle5'],['muscle','muscle5']),
    'CLUSTALO': (['clustalo'],['clustalo','clustal-omega']),
    'FASTTREE': (['FastTree','fasttree'],['fasttree','FastTree']),
    'VSEARCH': (['vsearch'],['vsearch']),
}

class CommandError(RuntimeError): pass

@dataclass(frozen=True)
class Command:
    argv: list[str]
    label: str
    log: str
    stdout: Path | None = None
    allow_empty: bool = False


class Toolchain:
    def __init__(self,workspace:Path, ui=None, auto_modules=True, executable_paths=None, tool_dirs=None):
        self.workspace=workspace; self.ui=ui; self.env=os.environ.copy()
        self.auto_modules=auto_modules and os.environ.get('ALIGNVIEW_AUTO_MODULES','1')!='0'
        self.loaded=[];self.history=[];self.failed_modules=set();self.version_cache={}
        self.overrides={}
        for spec in executable_paths or []:
            if '=' not in spec:raise ValueError(f'--executable expects COMMAND=/path: {spec}')
            name,filename=spec.split('=',1)
            known={cmd for cmds,_ in GROUPS.values() for cmd in cmds}
            if name not in known:raise ValueError(f'Unknown executable: {name}')
            path=Path(filename).expanduser().resolve()
            if not path.is_file() or not os.access(path,os.X_OK):
                raise ValueError(f'Executable missing or not executable: {path}')
            self.overrides[name]=str(path)
        paths=[]
        for folder in tool_dirs or []:
            path=Path(folder).expanduser().resolve()
            if not path.is_dir():raise ValueError(f'Missing tools directory: {path}')
            paths.append(str(path))
        if paths:self.env['PATH']=os.pathsep.join([*paths,self.env.get('PATH','')])
        self.refresh_managed()
        (workspace/'logs').mkdir(parents=True,exist_ok=True)

    def refresh_managed(self):
        from av_setup import managed_prefix
        folder=str(managed_prefix()/'bin')
        paths=self.env.get('PATH','').split(os.pathsep)
        if folder not in paths:paths.append(folder)
        self.env['PATH']=os.pathsep.join(paths)
        self.failed_modules.clear()

    def has(self,command):
        return command in self.overrides or shutil.which(command,path=self.env.get('PATH')) is not None

    def executable(self,command):return self.overrides.get(command) or shutil.which(command,path=self.env.get('PATH')) or command

    def source(self,command):
        if command in self.overrides:return 'explicit path'
        path=shutil.which(command,path=self.env.get('PATH'))
        if path:
            from av_setup import managed_prefix
            if Path(path).parent==managed_prefix()/'bin':return 'AlignView environment'
            return 'module HPC / PATH' if self.loaded else 'PATH / Conda'
        return 'not found'


    def available(self,command):
        if self.has(command): return True
        if not self.auto_modules: return False
        group=next((g for g,(names,_) in GROUPS.items() if command in names),None)
        if group is None: return False
        custom=os.environ.get('ALIGNVIEW_MODULE_'+group)
        candidates=custom.split(';') if custom else GROUPS[group][1]
        if self.ui: self.ui.info(f'{command} not on PATH: searching for an HPC module…')
        for candidate in candidates:
            if (command,candidate) in self.failed_modules:continue
            env=self._module_env(candidate.strip())
            if env and shutil.which(command,path=env.get('PATH')):
                self.env=env
                self.loaded.append((group,candidate))
                if self.ui: self.ui.success(f'Module loaded: {candidate} ({command})')
                return True
            self.failed_modules.add((command,candidate))
        return False

    def require(self,command):
        if not self.available(command):
            group=next((g for g,(names,_) in GROUPS.items() if command in names),'')
            hint=f' ; module load {GROUPS[group][1][0]}' if group else ''
            from av_setup import package_for
            raise CommandError(f'Missing command: {command}{hint}.\nInstall: alignview --setup {package_for(command)} --install-tools\nOr provide --tools-dir /path/bin, --executable {command}=/path/{command}, or ALIGNVIEW_MODULE_{group}.')
        return command

    def _module_env(self,name):
        """Read a module environment in an isolated bash process without modifying the user's shell."""
        if not name or any(x in name for x in '\n\r\0'): return None
        script='''
if ! type module >/dev/null 2>&1; then
  for f in /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash /usr/share/modules/init/bash; do
    if [[ -r "$f" ]]; then source "$f" >/dev/null 2>&1; break; fi
  done
fi
if type module >/dev/null 2>&1 && module load "$1" >/dev/null 2>&1; then env -0; fi
'''
        try:
            r=subprocess.run(['bash','-c',script,'alignview',name],env=self.env,stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,timeout=12,check=False)
            if r.returncode or not r.stdout: return None
            out={}
            for item in r.stdout.split(b'\0'):
                if b'=' in item:
                    k,v=item.split(b'=',1);out[k.decode('utf-8','replace')]=v.decode('utf-8','surrogateescape')
            return out if 'PATH' in out else None
        except (OSError,subprocess.TimeoutExpired): return None

    def run(self,args:list[str],label:str,log_name:str,stdout_path:Path | None=None,stdin_path:Path | None=None,
            allow_empty=False):
        self.require(args[0]); args=[self.executable(args[0]),*args[1:]]; logfile=self.next_log(log_name)
        if args[0] not in self.version_cache:
            from av_setup import probe_version
            self.version_cache[args[0]]=probe_version(args[0],self.env)
        if self.ui:
            self.ui.stage(label)
            self.ui.info('Tool: '+Path(args[0]).name+' · '+self.version_cache[args[0]])
            self.ui.info('Executable: '+args[0])
        start=time.time()
        with logfile.open('w',encoding='utf-8') as err:
            with stdout_path.open('wb') if stdout_path else open(os.devnull,'wb') as out:
                with stdin_path.open('rb') if stdin_path else open(os.devnull,'rb') as inp:
                    p=None
                    try:
                        p=subprocess.Popen(args,cwd=self.workspace,env=self.env,stdin=inp,stdout=out,stderr=err,start_new_session=True)
                        frames='◴◷◶◵';counter=0
                        while p.poll() is None:
                            if self.ui and self.ui.tty:
                                sys.stdout.write(f'  {frames[counter%4]} {label} · {time.time()-start:.0f} s\r')
                                sys.stdout.flush();counter+=1
                            time.sleep(.15)
                        if self.ui and self.ui.tty:
                            sys.stdout.write('\033[2K\r');sys.stdout.flush()
                    except KeyboardInterrupt:
                        if p is not None and p.poll() is None:
                            os.killpg(p.pid,signal.SIGTERM)
                            try:p.wait(timeout=3)
                            except subprocess.TimeoutExpired:
                                os.killpg(p.pid,signal.SIGKILL);p.wait()
                        self.history.append(dict(command=shlex.join(args),duration_s=round(time.time()-start,2),exit_code=p.returncode if p is not None else -2,interrupted=True))
                        raise
                    except OSError as exc: raise CommandError(f'Cannot launch {args[0]}: {exc}') from exc
        elapsed=time.time()-start
        self.history.append(dict(command=shlex.join(args),argv=args,cwd=str(self.workspace),
                                 stdout=str(stdout_path) if stdout_path else None,
                                 version=self.version_cache[args[0]],duration_s=round(elapsed,2),exit_code=p.returncode,
                                 stderr_log=str(logfile.relative_to(self.workspace))))
        if p.returncode:
            tail='\n'.join(logfile.read_text(errors='replace').splitlines()[-12:])
            raise CommandError(f'Failed {args[0]} (code {p.returncode})\n{tail}\nLog: {logfile}')
        if stdout_path and not allow_empty and (not stdout_path.exists() or stdout_path.stat().st_size==0):
            raise CommandError(f'{label}: sortie vide ({stdout_path.name})')
        if self.ui: self.ui.success(f'{label} · {elapsed:.1f} s')
        return logfile

    def execute(self,commands):
        for command in commands:
            self.run(command.argv,command.label,command.log,stdout_path=command.stdout,allow_empty=command.allow_empty)

    def next_log(self,name):
        path=self.workspace/'logs'/name
        index=2
        while path.exists():
            path=self.workspace/'logs'/(Path(name).stem+f'_{index}'+Path(name).suffix);index+=1
        return path

    def describe(self,commands):
        lines=['Commands prepared from the same arguments used for execution.',
               'This preview runs no analysis engine.',
               'Omitted options keep the installed version\'s defaults.',
               'HPC modules may resolve the executable path at run time.',
               '', 'Working directory: '+str(self.workspace)]
        for command in commands:
            argv=[self.executable(command.argv[0]),*command.argv[1:]]
            lines+=['',command.label,'Tool: '+command.argv[0]+' · '+self.source(command.argv[0]),
                    'Executable: '+argv[0],shlex.join(argv)+((' > '+shlex.quote(str(command.stdout))) if command.stdout else ''),
                    'Standard error log: '+str(self.next_log(command.log).relative_to(self.workspace))]
        return '\n'.join(lines)

    def engine_info(self,command):
        """Display-only discovery; no module load and no scientific calculation."""
        if not self.has(command):return command+' · not found'
        path=self.executable(command)
        if path not in self.version_cache:
            from av_setup import probe_version
            self.version_cache[path]=probe_version(path,self.env)
        return command+' · '+self.version_cache[path]+' · '+self.source(command)+' · '+path

    @staticmethod
    def _option(args,options,name,flag):
        value=getattr(options,name,None)
        if value is not None and value not in ('native',''):args.extend([flag,str(value)])

    def search_commands(self,mode,query:Path,refs:Path,threads=4,evalue=1e-5,max_target=100,sensitivity='sensitive',task='',options=None):
        if options is not None:validate(options,context='search',mode=mode)
        raw=self.workspace/'raw_hsps.tsv'; db=self.workspace/'database';out=self.workspace
        if mode=='mmseqs-blastp':
            fmt='query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits,qaln,taln'
            args=['mmseqs','easy-search',str(query),str(refs),str(raw),str(out/'mmseqs_tmp'),
                  '--format-output',fmt,'-a','--max-seqs',str(max_target),'-e',str(evalue),'--threads',str(threads)]
            for name,flag in (('mmseqs_sensitivity','-s'),('mmseqs_mask','--mask'),('mmseqs_comp_bias','--comp-bias-corr')):
                self._option(args,options,name,flag)
            return [Command(args,'MMseqs2 easy-search (protein)','mmseqs_search.log',allow_empty=True)]
        if mode.startswith('diamond-'):
            index=Command(['diamond','makedb','--in',str(refs),'--db',str(db)],'Index DIAMOND','diamond_db.log')
            m=mode.removeprefix('diamond-')
            args=['diamond',m,'--query',str(query),'--db',str(db),'--threads',str(threads),
                  '--evalue',str(evalue),'--max-target-seqs',str(max_target),f'--{sensitivity}',
                  '--outfmt','6','qseqid','sseqid','pident','length','qstart','qend','sstart','send',
                  'evalue','bitscore','qframe',('qseq_translated' if m=='blastx' else 'qseq'),'sseq', '--out',str(raw)]
            if m=='blastx': args+=['--min-orf',str(getattr(options,'diamond_min_orf',1)),'--query-gencode',str(getattr(options,'genetic_code',1))]
            for name,flag in (('matrix','--matrix'),('gapopen','--gapopen'),('gapextend','--gapextend'),
                              ('masking','--masking'),('comp_stats','--comp-based-stats'),('max_hsps','--max-hsps')):
                self._option(args,options,'diamond_'+name,flag)
            return [index,Command(args,f'Recherche DIAMOND {m.upper()}','diamond_search.log',allow_empty=True)]
        dbtype='nucl' if mode in ('blastn','tblastn','tblastx') else 'prot'
        index=Command(['makeblastdb','-in',str(refs),'-dbtype',dbtype,'-parse_seqids','-out',str(db)],
                 'Index BLAST local','blast_db.log')
        args=[mode,'-query',str(query),'-db',str(db),'-num_threads',str(threads),'-evalue',str(evalue),
              '-max_target_seqs',str(max_target),'-outfmt',
              '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qcovhsp qframe sframe qseq sseq',
              '-out',str(raw)]
        if task: args+=['-task',task]
        if mode in ('blastx','tblastx'):args+=['-query_gencode',str(getattr(options,'genetic_code',1))]
        if mode in ('tblastn','tblastx'):args+=['-db_gencode',str(getattr(options,'reference_genetic_code',None) or getattr(options,'genetic_code',1))]
        for name,flag in (('word_size','-word_size'),('max_hsps','-max_hsps'),('soft_masking','-soft_masking')):
            self._option(args,options,'blast_'+name,flag)
        self._option(args,options,'blast_masking','-dust' if mode=='blastn' else '-seg')
        if mode in ('blastn','blastx','tblastx'):self._option(args,options,'blast_strand','-strand')
        if mode=='blastn':
            self._option(args,options,'blast_reward','-reward');self._option(args,options,'blast_penalty','-penalty')
        else:self._option(args,options,'blast_matrix','-matrix')
        if mode!='tblastx':
            self._option(args,options,'blast_gapopen','-gapopen');self._option(args,options,'blast_gapextend','-gapextend')
        if mode in ('blastp','blastx','tblastn'):self._option(args,options,'blast_comp_stats','-comp_based_stats')
        return [index,Command(args,f'Recherche {mode.upper()}','blast_search.log',allow_empty=True)]

    def search(self,mode,query:Path,refs:Path,threads=4,evalue=1e-5,max_target=100,sensitivity='sensitive',task='',options=None):
        self.execute(self.search_commands(mode,query,refs,threads,evalue,max_target,sensitivity,task,options))
        raw=self.workspace/'raw_hsps.tsv'
        if mode=='mmseqs-blastp':
            from av_core import parse_mmseqs_tab
            return parse_mmseqs_tab(raw)
        if mode.startswith('diamond-'):return parse_blast_tab(raw,'diamond')
        return parse_blast_tab(raw,'blast')

    def mafft_command(self,input:Path,output:Path,strategy='auto',threads=4,fragments:Path|None=None,options=None):
        iterations=getattr(options,'mafft_maxiterate',None)
        settings={'auto':['--auto'],'global':['--globalpair','--maxiterate','1000'],
                  'local':['--localpair','--maxiterate','1000'],'fast':['--retree','2','--maxiterate','0']}
        if strategy not in settings:raise ValueError('MAFFT strategy: auto, local, global, fast')
        if not fragments and strategy=='auto' and iterations is not None:
            raise ValueError('Explicit MAFFT refinement iterations require global, local or fast instead of auto.')
        flags=['--auto'] if fragments else settings[strategy][:]
        if not fragments and iterations is not None:flags[flags.index('--maxiterate')+1]=str(iterations)
        args=['mafft',*flags,'--thread',str(threads)]
        self._option(args,options,'mafft_op','--op');self._option(args,options,'mafft_ep','--ep')
        if fragments:args+=['--addfragments',str(fragments)]
        args.append(str(input))
        return Command(args,f'MAFFT {strategy if not fragments else "addfragments (auto)"}',f'mafft_{output.stem}.log',output)

    def mafft(self,input:Path,output:Path,strategy='auto',threads=4,fragments:Path|None=None,options=None):
        self.execute([self.mafft_command(input,output,strategy,threads,fragments,options)])
        return output

    def msa_command(self,input:Path,output:Path,engine='mafft',strategy='auto',threads=4,options=None):
        if options is not None:validate(options,context='msa')
        if engine=='mafft':return self.mafft_command(input,output,strategy,threads,options=options)
        if engine=='muscle':
            program='muscle5' if self.has('muscle5') and not self.has('muscle') else 'muscle'
            method=getattr(options,'muscle_method','align')
            args=[program,'-'+method,str(input),'-output',str(output),'-threads',str(threads)]
            for name,flag in (('muscle_perturb','-perturb'),('muscle_perm','-perm'),('muscle_refineiters','-refineiters')):
                self._option(args,options,name,flag)
            return Command(args,'MUSCLE 5: '+method,'muscle.log')
        elif engine=='clustalo':
            args=['clustalo','-i',str(input),'-o',str(output),'--force','--outfmt=fasta','--threads',str(threads)]
            self._option(args,options,'clustalo_iterations','--iterations')
            if getattr(options,'clustalo_full',False):args+=['--full']
            if getattr(options,'clustalo_full_iter',False):args+=['--full-iter']
            return Command(args,'Clustal Omega: MSA','clustalo.log')
        else:raise ValueError(f'Unknown MSA engine: {engine}')

    def msa(self,input:Path,output:Path,engine='mafft',strategy='auto',threads=4,options=None):
        if engine=='muscle' and not self.available('muscle'):self.require('muscle5')
        self.execute([self.msa_command(input,output,engine,strategy,threads,options)])
        if not output.is_file() or output.stat().st_size==0:raise CommandError(f'{engine}: missing MSA file')
        return output

    def fasttree_command(self,alignment:Path,output:Path,nucleotide=False,options=None):
        cmd='fasttree' if self.has('fasttree') and not self.has('FastTree') else 'FastTree'
        args=[cmd,'-nt'] if nucleotide else [cmd]
        model=getattr(options,'fasttree_model','auto')
        if model not in (('auto','JC','GTR') if nucleotide else ('auto','JTT','WAG','LG')):
            raise ValueError(f'FastTree model {model} incompatible with '+('DNA' if nucleotide else 'proteins'))
        if model in ('GTR','WAG','LG'):args+=['-'+model.lower()]
        if getattr(options,'fasttree_gamma',False):args+=['-gamma']
        support=getattr(options,'fasttree_boot',None)
        if support==0:args+=['-nosupport']
        elif support is not None:args+=['-boot',str(support)]
        self._option(args,options,'fasttree_seed','-seed')
        args.append(str(alignment))
        return Command(args,'FastTree: approximate tree','fasttree.log',output)

    def fasttree(self,alignment:Path,output:Path,nucleotide=False,options=None):
        if not self.available('FastTree'):self.require('fasttree')
        self.execute([self.fasttree_command(alignment,output,nucleotide,options)])
        if not output.read_text().strip().endswith(';'):
            raise CommandError('FastTree: incomplete Newick file')
        return output

    def emboss(self,algorithm,query:Seq,reference:Seq,outprefix:Path,mode='blastp',
              gapopen=10,gapextend=0.5,free_ends=False,matrix='',frame=1,sframe=1,options=None,preview=False):
        if algorithm not in ('needle','water'): raise ValueError('EMBOSS: needle or water')
        q=query;r=reference
        if mode in ('blastx','diamond-blastx'):
            q=Seq('QUERY','QUERY',translate_nt(query.seq,frame,getattr(options,'genetic_code',1)))
        if mode=='tblastn': r=Seq(reference.id,reference.id,translate_nt(reference.seq,sframe,getattr(options,'reference_genetic_code',None) or getattr(options,'genetic_code',1)))
        if mode=='tblastx':
            q=Seq('QUERY','QUERY',translate_nt(query.seq,frame,getattr(options,'genetic_code',1)))
            r=Seq(reference.id,reference.id,translate_nt(reference.seq,sframe,getattr(options,'reference_genetic_code',None) or getattr(options,'genetic_code',1)))
        if mode in ('blastx','diamond-blastx','tblastn','tblastx') and ('*' in q.seq or '*' in r.seq):
            if self.ui: self.ui.warn('The full translation contains stop codons: check the coding region and frame.')
        qfile=self.workspace/'emboss_query.fasta'; rfile=self.workspace/'emboss_ref.fasta'
        write_fasta([q],qfile);write_fasta([r],rfile)
        is_dna=mode in ('blastn','global-dna')
        base=[algorithm,'-asequence',str(qfile),'-bsequence',str(rfile),'-gapopen',str(gapopen),
              '-gapextend',str(gapextend),'-auto','-snucleotide1' if is_dna else '-sprotein1',
              '-snucleotide2' if is_dna else '-sprotein2']
        if algorithm=='needle' and not free_ends:
            endopen=getattr(options,'needle_endopen',None)
            endextend=getattr(options,'needle_endextend',None)
            base+=['-endweight','-endopen',str(gapopen if endopen is None else endopen),
                   '-endextend',str(gapextend if endextend is None else endextend)]
        elif algorithm=='needle':base+=['-endweight','N']
        if matrix: base+=['-datafile',matrix]
        text=outprefix.with_suffix('.txt');fasta=outprefix.with_suffix('.fasta')
        commands=[Command([*base,'-outfile',str(text)],f'{algorithm.upper()} pairwise','emboss_pairwise.log'),
                  Command([*base,'-aformat3','fasta','-outfile',str(fasta)],f'{algorithm.upper()} export FASTA','emboss_pairwise_fasta.log')]
        if preview:return self.describe(commands)
        self.execute(commands)
        return text,fasta

    def alv(self,aligned:Path,variable=False):
        self.require('alv')
        args=['alv','--only-variable','-k',str(aligned)] if variable else ['alv','-k',str(aligned)]
        env=self.env.copy()
        # ALV normally does not need LD_LIBRARY_PATH; avoid affecting less.
        p=subprocess.run([self.executable(args[0]),*args[1:]],env=env,cwd=self.workspace,capture_output=True,text=True,check=False)
        if p.returncode: raise CommandError('ALV: '+p.stderr.strip())
        return p.stdout
