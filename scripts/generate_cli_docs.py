#!/usr/bin/env python3
"""Generate the complete CLI reference from the application parser."""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import build_parser


def reference(locale='fr'):
    from av_i18n import set_language
    set_language(locale)
    parser=build_parser()
    lines=['# Référence de la ligne de commande','',
        'Cette page est générée depuis le parseur réel par `python3 scripts/generate_cli_docs.py`. Ne pas modifier ses tables à la main.','',
        'Dans les exemples, `alignview` désigne la commande installée. Sans installation, remplacer ce mot par `/chemin/vers/AlignView/alignview.sh`.','',
        '## Choisir une commande','',
        '```bash','alignview --demo','alignview --docs demarrage',
        'alignview aligne.fasta --msa view',
        'alignview sequences.fasta --msa align --report --save --out mon_projet --no-menu',
        'alignview requetes.fna references.faa --mode blastx --query-id mon_identifiant --save --out recherche --no-menu',
        'alignview aligne.fasta --msa view --tree --html-out arbre.html --no-menu',
        'alignview --resume mon_projet --html-out rapport_actualise.html --no-menu',
        'alignview --resume mon_projet --msa align --tree --report --save --out nouveau_calcul --no-menu','```','',
        '## Règles à retenir','',
        '- Les chemins relatifs partent du dossier depuis lequel la commande est lancée ; le lanceur ne change pas de dossier.',
        '- `--save --out DOSSIER` conserve un projet portable. `--out` seul ne déclenche pas la sauvegarde.',
        '- `--html-out FICHIER.html` conserve seulement le HTML autonome. `--export-dir DOSSIER` conserve les exports du profil choisi.',
        '- `--report` prépare le rapport en mémoire de session ; ajouter une destination ou une sauvegarde pour le conserver.',
        '- Un fichier ou dossier de destination existant est refusé en batch ; choisir un nouveau nom.',
        '- `--resume PROJET` ne relance aucun moteur par défaut. Seules les actions explicitement ajoutées (`--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif`) sont exécutées.',
        '- En reprise, les paramètres explicites deviennent les réglages futurs. Les paramètres historiques des résultats déjà calculés restent conservés. `--report` et `--html-out` reconstruisent le HTML sans relancer de moteur.',
        '- Un profil JSON fournit les valeurs de départ ; les options explicites sont prioritaires. Un profil de seuil (`--preset`) est distinct d’un profil de fichiers (`--export-profile`).','',
        'Les modes, les stratégies et leurs limites sont expliqués dans [les parcours](USAGE.md), [les paramètres](PARAMETERS.md) et [le guide des fichiers](FILES.md).','']
    if locale=='en':
        lines=['# Command-line reference','',
            'Generated from the actual parser by `python3 scripts/generate_cli_docs.py`. Do not edit the tables manually.','',
            'Use the installed `alignview` command, or replace it with `/path/to/AlignView/alignview.sh`.','',
            '## Common workflows','',
            '```bash','alignview --lang en --demo','alignview --lang en --docs quickstart',
            'alignview aligned.fasta --lang en --msa view --html-out report.html --no-menu',
            'alignview sequences.fasta --lang en --msa align --report --save --out project_01 --no-menu',
            'alignview queries.fna refs.fna --lang en --vsearch search --query-id chosen_id --report --save --out search_01 --no-menu',
            'alignview --resume project_01 --lang en --html-out updated.html --no-menu','```','',
            '## Storage and resuming','',
            '- Relative paths start from the launch directory; the launcher does not change it.',
            '- `--save --out DIRECTORY` retains a project. `--out` alone does not save.',
            '- `--html-out FILE.html` retains standalone HTML. `--export-dir DIRECTORY` retains exports of the selected profile.',
            '- `--report` alone prepares temporary HTML. Add a destination or saving to retain it.',
            '- Batch runs refuse an existing destination.',
            '- Resume alone reruns no scientific engine. Only explicit actions (`--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif`, `--vsearch`) request calculations.',
            '- Explicit parameters override profiles/saved future settings; historical settings remain associated with existing results.',
            '- `--genetic-code` selects the query/internal table; `--reference-genetic-code` selects or inherits the reference table.',
            '- `--lang en` selects English; `ALIGNVIEW_LANG=en` provides a default. Raw data and engine arguments are preserved.','',
            'See [Usage](USAGE.md), [Parameters](PARAMETERS.md), [VSEARCH](VSEARCH.md), [Genetic codes](GENETIC_CODES.md) and [Files](FILES.md).','']
    def cell(value):return str(value).replace('|','\\|').replace('\n',' ')
    for group in parser._action_groups:
        actions=[a for a in group._group_actions if a.help!=argparse.SUPPRESS]
        if not actions:continue
        title=('General help' if locale=='en' else 'Aide générale') if group.title=='options' else group.title
        lines+=['## '+title,'',('| Option | Default | Purpose and values |' if locale=='en' else '| Option | Défaut | Rôle et valeurs |'),'| --- | --- | --- |']
        for action in actions:
            flags=' / '.join(action.option_strings) or action.dest
            default='—' if action.default is None or action.default==argparse.SUPPRESS else repr(action.default)
            description=(action.help or '').replace('%%','%')
            if action.choices is not None:description+=(' Values: ' if locale=='en' else ' Valeurs : ')+', '.join(str(c) for c in action.choices)+'.'
            lines.append('| `'+cell(flags)+'` | '+cell(default)+' | '+cell(description)+' |')
        lines.append('')
    return '\n'.join(lines)


if __name__=='__main__':
    for locale,relative in (('fr','CLI.md'),('en','en/CLI.md')):
        target=ROOT/'docs'/relative;content=reference(locale)
        if '--check' in sys.argv:
            if not target.is_file() or target.read_text(encoding='utf-8')!=content:
                raise SystemExit('CLI references need regeneration: python3 scripts/generate_cli_docs.py')
        else:target.write_text(content,encoding='utf-8')
    print('French and English CLI references '+('are current' if '--check' in sys.argv else 'generated'))
