# Contribuer

Lancer `python3 -m unittest discover -s tests -v`. Les tests ordinaires ne téléchargent aucun moteur. Voir [DEVELOPMENT.md](docs/DEVELOPMENT.md) pour les validations réelles et la construction.

Pour un problème, joindre une commande minimale, les versions d’AlignView et Python, le terminal, `alignview --diagnose-json` et un petit FASTA synthétique. Préciser : local, SSH ou Slurm.

Toute modification des coordonnées ou statistiques doit inclure un résultat connu attendu. Vérifier la navigation à 80×24, avec une longue liste et un long alignement. Conserver l’exécution sans dépendance Python obligatoire et le fonctionnement batch.
