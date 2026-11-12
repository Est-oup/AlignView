"""Offline manual shared by the terminal, wheel and zipapp."""
from __future__ import annotations

from importlib import resources
from pathlib import Path

TOPICS = {
    'index': ('README.md', 'Choisir un guide', 'Tous les guides, classés par besoin'),
    'demarrage': ('QUICKSTART.md', 'Premiers pas', 'Un exemple et les premiers résultats, étape par étape'),
    'parcours': ('USAGE.md', 'Utiliser toutes les fonctions', 'Menus, calculs et exploration des résultats'),
    'fichiers': ('FILES.md', 'Dossiers et sauvegardes', 'Application, sorties, temporaires et reprise'),
    'installation': ('INSTALLATION.md', 'Installer et déplacer AlignView', 'Python, moteurs scientifiques et Linux / HPC'),
    'parametres': ('PARAMETERS.md', 'Comprendre les paramètres', 'Effets, défauts et profils par moteur'),
    'html': ('HTML.md', 'Lire et partager le rapport HTML', 'Figures, arbre, filtres et téléchargements'),
    'exports': ('EXPORTS.md', 'Réutiliser les résultats', 'Formats, coordonnées et profils d’export'),
    'cli': ('CLI.md', 'Toutes les options de commande', 'Référence générée à partir du parseur réel'),
    'depannage': ('TROUBLESHOOTING.md', 'Résoudre un problème', 'Démarrage, moteurs, résultats et récupération'),
    'glossaire': ('GLOSSARY.md', 'Comprendre les termes', 'FASTA, HSP, couverture, MSA, supports et coordonnées'),
    'vsearch': ('VSEARCH.md','Analyses VSEARCH','Recherche globale, déduplication, groupes et chimères'),
    'translation': ('GENETIC_CODES.md','Tables de traduction','Codes NCBI, cadres, initiateurs et arrêts contextuels'),
    'developpement': ('DEVELOPMENT.md', 'Développer et tester', 'Architecture, tests et fabrication des distributions'),
    'publication': ('RELEASING.md', 'Préparer une publication', 'Archives, validations et limites de la release'),
    'audit': ('AUDIT.md', 'Portée des validations', 'Mesures, moteurs réels et limites connues'),
}
ALIASES={'getting-started':'demarrage','quickstart':'demarrage','usage':'parcours','files':'fichiers',
         'parameters':'parametres','troubleshooting':'depannage','glossary':'glossaire',
         'development':'developpement','releasing':'publication','genetics':'translation'}


def read_document(topic='index'):
    from av_i18n import language
    topic=ALIASES.get(topic,topic)
    if topic not in TOPICS:
        raise ValueError('Guide inconnu : '+topic)
    filename = TOPICS[topic][0]
    if language()=='en':
        english=Path(__file__).resolve().parent/'docs/en'/filename
        if english.is_file():return english.read_text(encoding='utf-8')
        try:return resources.files('alignview_docs').joinpath('en/'+filename).read_text(encoding='utf-8')
        except (ModuleNotFoundError,FileNotFoundError):pass
    local = Path(__file__).resolve().parent/'docs'/filename
    if local.is_file():
        return local.read_text(encoding='utf-8')
    try:
        return resources.files('alignview_docs').joinpath(filename).read_text(encoding='utf-8')
    except (ModuleNotFoundError, FileNotFoundError) as exc:
        raise ValueError('Documentation absente : '+filename+'. Utiliser une distribution AlignView complète.') from exc


def documentation_menu(ui):
    with ui.section('Documentation'):
        while True:
            keys = list(TOPICS)
            options = [(str(i), title, description, '') for i, (_, title, description)
                       in enumerate(TOPICS.values(), 1)]
            choice = ui.menu('Documentation complète', options,
                             hint='Documentation locale · / rechercher un guide · Échap retour')
            if choice is None: return
            topic = keys[int(choice)-1]
            ui.page(read_document(topic), TOPICS[topic][1])
