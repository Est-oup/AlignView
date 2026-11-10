# Première publication 0.1.0

Le dossier racine à placer dans GitHub est celui contenant ce README et `pyproject.toml`, sans les données voisines ni les archives des prototypes.

## État préparé

Le nom public et la commande sont `AlignView` / `alignview`. La version publique démarre à 0.1.0. Les sources, le wheel, le zipapp, les environnements Conda, la recette Bioconda et la CI sont préparés. Aucune publication GitHub, PyPI ou Bioconda n’est effectuée ici.

## Avant le tag

1. Exécuter la CI du dépôt créé et les validations des moteurs réels.
2. Tester l’installation dans un environnement neuf, l’ouverture de projets transférés et les raccourcis dans les terminaux SSH utilisés.
3. Valider le profil minimal sur Linux et le cluster cible ; conserver versions et résultats. Vérifier aussi BLASTX et TBLASTN avec des séquences de référence connues.
4. Compléter URL du dépôt, auteur/mainteneur et contact. Ne pas inventer de DOI. Vérifier la disponibilité du nom sur les registres au moment de publier.
5. Actualiser `av_version.py`, le changelog et la recette ; construire et vérifier les artefacts.

```bash
python -m unittest discover -s tests -v
python scripts/validate_real_tools.py --require mafft blastn diamond
python -m build
python scripts/build_pyz.py
python -m twine check dist/*
python scripts/checksums.py
```

Le workflow `artifacts.yml`, déclenché par tag ou manuellement, construit les fichiers et les conserve comme artefacts Actions. Il ne crée pas de release publique et ne téléverse rien sur PyPI. Les fichiers sont `dist/alignview-0.1.0.tar.gz`, `dist/alignview-0.1.0-py3-none-any.whl`, `alignview.pyz` et `dist/SHA256SUMS`.

## PyPI

Une fois le compte, le nom et l’URL du dépôt finalisés, publier le wheel et le sdist validés. Un déploiement avec Trusted Publishing peut être configuré dans le vrai dépôt. Les exécutables scientifiques restent hors du wheel Python.

## Bioconda

La [recette modèle](../packaging/bioconda/meta.yaml.template) est un paquet `noarch: python`. Remplacer propriétaire, dépôt, URL de release stable, SHA256 réel et mainteneur GitHub. Tester puis soumettre la recette à `bioconda/bioconda-recipes` selon les [directives officielles](https://bioconda.github.io/contributor/guidelines.html).

Les moteurs restent facultatifs pour permettre une installation légère et les modules HPC. Les fichiers d’environnement et `--setup` fournissent les profils intégrés. Documenter cette distinction dans la soumission.

## Vers 1.0

La désignation 1.0 viendra après validation des parcours sur plusieurs postes et clusters, des moteurs annoncés, de l’interopérabilité des exports dans des outils tiers et de jeux biologiques de référence. Les limites connues sont maintenues dans l’audit. Une numérotation « 1.0 » seule ne remplace pas ces validations.
