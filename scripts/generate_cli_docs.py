#!/usr/bin/env python3
"""Generate the complete CLI reference from the application parser."""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from alignview import build_parser


def reference():
    parser=build_parser()
    lines=['# Command-line reference','',
        'Generated from the actual parser by `python3 scripts/generate_cli_docs.py`. Do not edit the tables manually.','',
        'Use the installed `alignview` command, or replace it with `/path/to/AlignView/alignview.sh`.','',
        '## Common workflows','',
        '```bash','alignview --demo','alignview --docs quickstart',
        'alignview aligned.fasta --msa view --html-out report.html --no-menu',
        'alignview sequences.fasta --msa align --report --save --out project_01 --no-menu',
        'alignview queries.fna refs.fna --vsearch search --query-id chosen_id --report --save --out search_01 --no-menu',
        'alignview --resume project_01 --html-out updated.html --no-menu','```','',
        '## Storage and resuming','',
        '- Relative paths start from the launch directory; the launcher does not change it.',
        '- `--save --out DIRECTORY` retains a project. `--out` alone does not save.',
        '- `--html-out FILE.html` retains standalone HTML. `--export-dir DIRECTORY` retains exports of the selected profile.',
        '- `--report` alone prepares temporary HTML. Add a destination or saving to retain it.',
        '- Batch runs refuse an existing destination.',
        '- Resume alone reruns no scientific engine. Only explicit actions (`--msa`, `--tree`, `--pairwise`, `--dotplot`, `--motif`, `--vsearch`) request calculations.',
        '- Explicit parameters override profiles/saved future settings; historical settings remain associated with existing results.',
        '- `--genetic-code` selects the query/internal table; `--reference-genetic-code` selects or inherits the reference table.',
        'See [Usage](USAGE.md), [Parameters](PARAMETERS.md), [VSEARCH](VSEARCH.md), [Genetic codes](GENETIC_CODES.md) and [Files](FILES.md).','']
    def cell(value):return str(value).replace('|','\\|').replace('\n',' ')
    for group in parser._action_groups:
        actions=[a for a in group._group_actions if a.help!=argparse.SUPPRESS]
        if not actions:continue
        # argparse called this group "optional arguments" before Python 3.11.
        # Documentation must reflect our interface, not interpreter wording.
        title='General help' if group is parser._optionals else group.title
        lines+=['## '+title,'','| Option | Default | Purpose and values |','| --- | --- | --- |']
        for action in actions:
            flags=' / '.join(action.option_strings) or action.dest
            default='—' if action.default is None or action.default==argparse.SUPPRESS else repr(action.default)
            description=(action.help or '').replace('%%','%')
            if action.choices is not None:description+=' Values: '+', '.join(str(c) for c in action.choices)+'.'
            lines.append('| `'+cell(flags)+'` | '+cell(default)+' | '+cell(description)+' |')
        lines.append('')
    return '\n'.join(lines)


if __name__=='__main__':
    target=ROOT/'docs/CLI.md'
    content=reference()
    if '--check' in sys.argv:
        if not target.is_file() or target.read_text(encoding='utf-8')!=content:
            raise SystemExit('CLI reference needs regeneration: python3 scripts/generate_cli_docs.py')
    else:
        target.write_text(content,encoding='utf-8')
    print('CLI reference '+('is current' if '--check' in sys.argv else 'generated'))
