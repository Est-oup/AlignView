"""Offline manual shared by the terminal, wheel and zipapp."""
from __future__ import annotations

from importlib import resources
from pathlib import Path

TOPICS = {
    'index': ('README.md', 'Choose a guide', 'All guides, organized by task'),
    'quickstart': ('QUICKSTART.md', 'Getting started', 'A demo and your first results, step by step'),
    'usage': ('USAGE.md', 'Using AlignView', 'Menus, calculations and result exploration'),
    'files': ('FILES.md', 'Files and saving', 'Application files, outputs, temporary data and resuming'),
    'installation': ('INSTALLATION.md', 'Installing and moving AlignView', 'Python, scientific engines and Linux / HPC'),
    'parameters': ('PARAMETERS.md', 'Understanding parameters', 'Effects, defaults and profiles for each engine'),
    'html': ('HTML.md', 'Reading and sharing HTML', 'Figures, trees, filters and downloads'),
    'exports': ('EXPORTS.md', 'Reusing results', 'Formats, coordinates and export profiles'),
    'cli': ('CLI.md', 'Command-line reference', 'Generated from the actual argument parser'),
    'troubleshooting': ('TROUBLESHOOTING.md', 'Troubleshooting', 'Startup, engines, results and recovery'),
    'glossary': ('GLOSSARY.md', 'Understanding terminology', 'FASTA, HSP, coverage, MSA, supports and coordinates'),
    'vsearch': ('VSEARCH.md', 'VSEARCH analyses', 'Global search, dereplication, clustering and chimeras'),
    'translation': ('GENETIC_CODES.md', 'Translation tables', 'NCBI codes, frames, start codons and context-dependent stops'),
    'development': ('DEVELOPMENT.md', 'Development and testing', 'Architecture, tests and distribution builds'),
    'releasing': ('RELEASING.md', 'Preparing a release', 'Archives, validation and release limitations'),
    'audit': ('AUDIT.md', 'Validation scope', 'Measurements, real engines and known limitations'),
}
ALIASES = {'getting-started': 'quickstart', 'genetics': 'translation'}


def read_document(topic='index'):

    topic=ALIASES.get(topic,topic)
    if topic not in TOPICS:
        raise ValueError('Unknown guide: '+topic)
    filename = TOPICS[topic][0]
    local = Path(__file__).resolve().parent/'docs'/filename
    if local.is_file():
        return local.read_text(encoding='utf-8')
    try:
        return resources.files('alignview_docs').joinpath(filename).read_text(encoding='utf-8')
    except (ModuleNotFoundError, FileNotFoundError) as exc:
        raise ValueError('Missing documentation: '+filename+'. Use a complete AlignView distribution.') from exc


def documentation_menu(ui):
    with ui.section('Documentation'):
        while True:
            keys = list(TOPICS)
            options = [(str(i), title, description, '') for i, (_, title, description)
                       in enumerate(TOPICS.values(), 1)]
            choice = ui.menu('Complete manual', options,
                             hint='Offline manual · / find a guide · Esc back')
            if choice is None: return
            topic = keys[int(choice)-1]
            ui.page(read_document(topic), TOPICS[topic][1])
