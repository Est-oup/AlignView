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
        raise ValueError(f'Profil inconnu : {profile}')
    prefix = Path(prefix or managed_prefix()).expanduser().resolve()
    program = manager or next((shutil.which(x) for x in ('micromamba','mamba','conda') if shutil.which(x)), None)
    if not program:
        raise ValueError('Installer Micromamba, Mamba ou Conda, puis relancer --setup. Guide : https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html')
    if prefix.exists() and not (prefix/'conda-meta').is_dir():
        raise ValueError(f'Le préfixe existe sans environnement Conda : {prefix}. Choisir --tools-prefix avec un dossier libre.')
    operation = 'install' if (prefix/'conda-meta').is_dir() else 'create'
    return [str(program), operation, '--prefix', str(prefix), '--override-channels',
            '--channel', 'conda-forge', '--channel', 'bioconda', '--strict-channel-priority',
            '--yes', *PACKAGES[profile]]


def setup_tools(ui, profile, prefix=None, install=False):
    from av_tools import CommandError
    command = install_plan(profile, prefix)
    ui.heading('INSTALLER LES OUTILS / '+profile.upper())
    ui.info('Environnement dédié : '+str(Path(prefix or managed_prefix())))
    print('  '+shlex.join(command))
    if not install:
        ui.info('Plan uniquement. Ajouter --install-tools pour exécuter cette commande.')
        return command
    ui.info('Installation en cours ; les messages du gestionnaire restent visibles.')
    try:
        result = subprocess.run(command, check=False)
    except OSError as exc:
        raise CommandError(f'Installation impossible : {exc}') from exc
    if result.returncode:
        raise CommandError(f'Installation interrompue ou échouée (code {result.returncode}). Relancer --setup {profile} --install-tools ; vérifier réseau, espace disque et canaux.')
    ui.success('Outils installés. AlignView les détecte sans activer un environnement.')
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
        if result.returncode:return f'version indisponible (code {result.returncode})'
        lines=[line.strip() for line in (result.stdout+'\n'+result.stderr).splitlines() if line.strip()]
        if command=='vsearch':
            return next((line for line in lines if re.match(r'vsearch\s+v?\d',line,re.I)),'version indisponible')[:240]
        return next(iter(lines),'version indisponible')[:240]
    except (OSError,subprocess.TimeoutExpired):return 'version indisponible'


def software_hub(ui, tools):
    from av_tools import CommandError
    while True:
        rows = diagnostic(tools, versions=False)
        found = sum(row['available'] for row in rows)
        choice = ui.menu('LOGICIELS / INSTALLATION', [
            ('1','État et versions','Exécutables, origine et commande de réparation','DIAGNOSTIC'),
            ('2','Installer le socle conseillé','BLAST + DIAMOND + MAFFT dans un environnement dédié','INSTALLER'),
            ('3','Installer un outil précis','Ajouter uniquement ce qui manque','INSTALLER'),
            ('4','Installer tous les moteurs','Ajoute EMBOSS, MUSCLE, Clustal Omega, MMseqs2, FastTree, VSEARCH, ALV','INSTALLER'),
            ('m','Rechercher les modules HPC','Charger à la demande selon la politique de la session','HPC'),
            ('q','Retour','Reprendre l’analyse','NAVIGATION')],
            f'{found}/{len(rows)} commandes disponibles · préfixe {managed_prefix()}')
        if choice in (None,'q'):
            return
        if choice in ('1','m'):
            rows = diagnostic(tools, load_modules=choice=='m')
            text = '\n\n'.join(f'{"✓" if r["available"] else "○"} {r["command"]} · {r["source"]}\n  '+
                (f'{r["path"]}\n  {r["version"]}' if r['available'] else f'alignview --setup {r["setup_profile"]} --install-tools') for r in rows)
            ui.page(text, 'Outils scientifiques / diagnostic')
            continue
        profile = 'minimal' if choice=='2' else 'full'
        if choice=='3':
            options=[(key,key,', '.join(values),'OUTILS') for key,values in PACKAGES.items() if key not in ('minimal','full')]
            profile=ui.menu('CHOISIR UN OUTIL',options)
            if profile is None:
                continue
        try:
            setup_tools(ui,profile)
            if ui.yesno('Exécuter cette installation dans le préfixe dédié ?',False):
                setup_tools(ui,profile,install=True)
                tools.refresh_managed()
        except (ValueError,CommandError,OSError) as exc:
            ui.error(str(exc))
        ui.pause()
