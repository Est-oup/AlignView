"""Discover scientific engines and prepare isolated Conda/Bioconda installs."""
from __future__ import annotations

import json
import os
import re
from pathlib import Path
import shlex
import shutil
import subprocess

PACKAGES = {
    'blast': ['blast'], 'diamond': ['diamond'], 'mafft': ['mafft'],
    'emboss': ['emboss'], 'muscle': ['muscle>=5'], 'clustalo': ['clustalo'],
    'mmseqs': ['mmseqs2'], 'fasttree': ['fasttree'], 'alv': ['alv'],
    'minimal': ['blast', 'diamond', 'mafft'],
    'vsearch': ['vsearch'],
    'full': ['blast', 'diamond', 'mafft', 'emboss', 'muscle>=5', 'clustalo', 'mmseqs2', 'fasttree', 'vsearch', 'alv'],
}
ALIASES = {'muscle5':'muscle', 'FastTree':'fasttree', 'makeblastdb':'blast',
           'needle':'emboss', 'water':'emboss', 'transeq':'emboss'}
VERSION_FLAGS = {'mafft':['--version'], 'diamond':['version'], 'muscle':['-version'],
                 'vsearch':['--version'],
                 'muscle5':['-version'], 'clustalo':['--version'], 'mmseqs':['version'],
                 'FastTree':['-help'], 'fasttree':['-help'], 'alv':['--version'],
                 'needle':['-version'], 'water':['-version'], 'transeq':['-version']}


def managed_prefix():
    override = os.environ.get('ALIGNVIEW_TOOLS_PREFIX')
    base = Path(os.environ.get('XDG_DATA_HOME', str(Path.home()/'.local/share')))
    return Path(override).expanduser().resolve() if override else base/'alignview'/'tools'


def package_for(command):
    if command.startswith(('blast','tblast')):
        return 'blast'
    return ALIASES.get(command, command)


def install_plan(profile='minimal', prefix=None, manager=None):
    if profile not in PACKAGES:
        raise ValueError(f'Unknown profile: {profile}')
    prefix = Path(prefix or managed_prefix()).expanduser().resolve()
    program = manager or next((shutil.which(x) for x in ('micromamba','mamba','conda') if shutil.which(x)), None)
    if not program:
        raise ValueError('Install Micromamba, Mamba or Conda, then rerun --setup. Guide: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html')
    if prefix.exists() and not (prefix/'conda-meta').is_dir():
        raise ValueError(f'The prefix exists without a Conda environment: {prefix}. Choose --tools-prefix with an unused directory.')
    operation = 'install' if (prefix/'conda-meta').is_dir() else 'create'
    return [str(program), operation, '--prefix', str(prefix), '--override-channels',
            '--channel', 'conda-forge', '--channel', 'bioconda', '--strict-channel-priority',
            '--yes', *PACKAGES[profile]]


def setup_tools(ui, profile, prefix=None, install=False):
    from av_tools import CommandError
    command = install_plan(profile, prefix)
    ui.heading('INSTALL TOOLS / '+profile.upper())
    ui.info('Dedicated environment: '+str(Path(prefix or managed_prefix())))
    print('  '+shlex.join(command))
    if not install:
        ui.info('Plan only. Add --install-tools to run this command.')
        return command
    ui.info('Installation in progress; package-manager messages remain visible.')
    try:
        result = subprocess.run(command, check=False)
    except OSError as exc:
        raise CommandError(f'Installation impossible: {exc}') from exc
    if result.returncode:
        raise CommandError(f'Installation interrupted or failed (exit code {result.returncode}). Run --setup again with {profile} --install-tools ; check network access, disk space and channels.')
    ui.success('Tools installed. AlignView detects them without activating an environment.')
    return command


def diagnostic(tools, versions=True, load_modules=False):
    from av_tools import GROUPS
    rows = []
    for group,(commands,_) in GROUPS.items():
        for command in commands:
            found = tools.available(command) if load_modules else tools.has(command)
            path = tools.executable(command) if found else ''
            version = ''
            if versions and found:
                version = probe_version(path,tools.env)
            rows.append(dict(group=group, command=command, available=found, path=path,
                             source=tools.source(command), version=version,
                             setup_profile=package_for(command)))
    return rows


def probe_version(path,env):
    command=Path(path).name
    try:
        result=subprocess.run([path,*VERSION_FLAGS.get(command,['-version'])],
            env=env,capture_output=True,text=True,timeout=4,errors='replace')
        if result.returncode:return f'version unavailable (code {result.returncode})'
        lines=[line.strip() for line in (result.stdout+'\n'+result.stderr).splitlines() if line.strip()]
        if command=='vsearch':
            return next((line for line in lines if re.match(r'vsearch\s+v?\d',line,re.I)),'version unavailable')[:240]
        return next(iter(lines),'version unavailable')[:240]
    except (OSError,subprocess.TimeoutExpired):return 'version unavailable'


def software_hub(ui, tools):
    from av_tools import CommandError
    while True:
        rows = diagnostic(tools, versions=False)
        found = sum(row['available'] for row in rows)
        choice = ui.menu('LOGICIELS / INSTALLATION', [
            ('1','Availability and versions','Executables, source and setup command','DIAGNOSTIC'),
            ('2','Install the recommended core','BLAST + DIAMOND + MAFFT in a dedicated environment','INSTALLER'),
            ('3','Install one tool','Add only what is missing','INSTALLER'),
            ('4','Install all engines','Adds EMBOSS, MUSCLE, Clustal Omega, MMseqs2, FastTree, VSEARCH and ALV','INSTALLER'),
            ('m','Look for HPC modules','Load on demand according to the session policy','HPC'),
            ('q','Back','Resume the analysis','NAVIGATION')],
            f'{found}/{len(rows)} available commands · prefix {managed_prefix()}')
        if choice in (None,'q'):
            return
        if choice in ('1','m'):
            rows = diagnostic(tools, load_modules=choice=='m')
            text = '\n\n'.join(f'{"✓" if r["available"] else "○"} {r["command"]} · {r["source"]}\n  '+
                (f'{r["path"]}\n  {r["version"]}' if r['available'] else f'alignview --setup {r["setup_profile"]} --install-tools') for r in rows)
            ui.page(text, 'Scientific tools / diagnostics')
            continue
        profile = 'minimal' if choice=='2' else 'full'
        if choice=='3':
            options=[(key,key,', '.join(values),'OUTILS') for key,values in PACKAGES.items() if key not in ('minimal','full')]
            profile=ui.menu('CHOOSE A TOOL',options)
            if profile is None:
                continue
        try:
            setup_tools(ui,profile)
            if ui.yesno('Run this installation in the dedicated prefix?',False):
                setup_tools(ui,profile,install=True)
                tools.refresh_managed()
        except (ValueError,CommandError,OSError) as exc:
            ui.error(str(exc))
        ui.pause()
