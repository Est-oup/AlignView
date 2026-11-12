# Contribuer

Lancer `python3 -m unittest discover -s tests -v`. Les tests ordinaires ne téléchargent aucun moteur. Voir [DEVELOPMENT.md](docs/DEVELOPMENT.md) pour les validations réelles et la construction.

Pour un problème, joindre une commande minimale, les versions d’AlignView et Python, le terminal, `alignview --diagnose-json` et un petit FASTA synthétique. Préciser : local, SSH ou Slurm.

Toute modification des coordonnées ou statistiques doit inclure un résultat connu attendu. Vérifier la navigation à 80×24, avec une longue liste et un long alignement. Conserver l’exécution sans dépendance Python obligatoire et le fonctionnement batch.

Les guides français et anglais sont consultables avec `--docs` et se trouvent dans `docs/` et `docs/en/`. Après une modification de CLI : `python3 scripts/generate_cli_docs.py`, puis `python3 scripts/check_docs.py`. Une nouvelle chaîne d’interface doit être traduite dans `av_locale_en.py`, sans modifier les identifiants scientifiques ni les arguments des moteurs. Une nouvelle option scientifique doit être reliée à sa commande réelle et à la provenance.
