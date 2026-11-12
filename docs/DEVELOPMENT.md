# Développer et vérifier AlignView

## Structure

- `alignview.py` : CLI, sessions et calculs ; `av_navigation.py` : parcours, prérequis et présentation.
- `av_defaults.py` : limite commune de 100 références et profils indépendants de cette limite.
- `av_parameters.py` : paramètres typés par moteur, validation, profils portables et instantanés de calcul.
- `av_vsearch.py` : commandes et résultats natifs VSEARCH ; `av_genetics.py` : tables NCBI hors ligne.
- `av_i18n.py`, `av_locale_en.py` : présentation français/anglais ; `av_docs.py` : guides hors ligne.
- `av_tree.py` : lecture Newick et figure SVG FastTree autonome, sans dépendance Python tierce.
- `av_core.py`, `av_units.py` : séquences, HSP, statistiques et unités.
- `av_ui.py`, `av_terminal.py`, `av_start.py`, `av_msa_tui.py` : menus, saisie UTF-8/CSI, sélecteurs, pager et visionneuse.
- `av_tools.py`, `av_setup.py` : moteurs, provenance, modules HPC, diagnostic et installation dédiée.
- `av_project.py` : instantanés portables, vérification des entrées, reprise et sauvegarde atomique.
- `av_visual.py`, `av_report_v4.py`, `av_report_v5.py`, `av_report_v6.py` : HTML/SVG ; les noms historiques sont conservés pour éviter de réécrire les calculs pendant la refonte des parcours.
- `av_exports.py` : formats et provenance ; `av_version.py` : version unique ; `av_demo.py` : exemples embarqués.

## Tests

```bash
python3 -m unittest discover -s tests -v
python3 scripts/generate_cli_docs.py --check
python3 scripts/check_docs.py
python3 scripts/validate_vsearch.py --output docs/vsearch-tools-results.json
python3 scripts/validate_real_tools.py --require mafft
python3 scripts/validate_real_tools.py --tuned --require mafft blastn diamond mmseqs muscle clustalo needle FastTree
python3 scripts/benchmark.py --rows 100 --columns 10000 --repeats 100
python3 scripts/check_report_browser.py  # optionnel : Chrome/Chromium installé
```

Les tests hérités couvrent les calculs et exports avec des exécutables simulés. `tests/test_release.py` ajoute la reprise, l’intégrité, les sauvegardes, les plans d’installation, les longues listes et de vrais terminaux PTY sans dépendance tierce. Le test historique utilisant `pexpect` est facultatif.

`tests/test_navigation.py` couvre les parcours complets, les retours F2, les prérequis, les écrans imbriqués, la séquence active et la limite de 100 dans la sélection et les exports. Les profils ne modifient plus le nombre de références. `scripts/preview_navigation.py` génère les aperçus de documentation à partir des vrais cadres de l’interface.

`tests/test_parameters.py` vérifie les réglages réellement transmis, la fidélité de l’aperçu, les profils invalides, les priorités CLI, les paramètres historiques après reprise, les scripts Slurm et l’accès terminal aux contrôles. Un nouveau paramètre doit être relié au constructeur de commande et à sa sélection contextuelle, pas seulement ajouté au parseur.

`tests/test_tree.py` vérifie la géométrie des distances, les supports, les identifiants cités, les arbres déséquilibrés, le repli explicite en cladogramme, les Newick invalides et l’intégration HTML/exports/reprise/invalidation. Le contrôle navigateur vérifie l’accès « Voir l’arbre », le filtrage Phylogénie et le contenu du téléchargement SVG.

`scripts/validate_real_tools.py` refuse le dossier des mocks. Un moteur absent est marqué `not installed`, jamais `passed`. La fixture vérifie des hits sur des séquences identiques et la structure des sorties ; elle ne mesure pas la justesse de tous les algorithmes sur des données difficiles.

Le benchmark peut comparer l’ancienne visionneuse : `--baseline /chemin/AlignView_Pro_v6`. Les temps dépendent de la machine ; ne pas extrapoler une mesure locale à une supériorité générale.

Le contrôle navigateur ouvre un rapport synthétique hors ligne et vérifie le maintien du défilement vertical, le pavé tactile, les bornes de colonnes, les pages de séquences, le filtre et une vue sans colonne variable. Il est distinct des tests sans navigateur.

## Construire

```bash
python -m pip install '.[dev]'
python -m build
python scripts/build_pyz.py
python alignview.pyz --demo --no-menu
python -m twine check dist/*.whl dist/*.tar.gz
python scripts/build_portable.py
python scripts/checksums.py
```

Sans le module `build`, avec les backends déjà présents :

```bash
python -m pip wheel --no-build-isolation --no-deps . -w dist
python -c "from setuptools.build_meta import build_sdist; build_sdist('dist')"
```

Tester le wheel dans un environnement vide et le `.pyz` depuis un dossier extérieur. La version des paquets provient de `av_version.VERSION` ; celle de la recette Bioconda doit être actualisée lors d’une release.

## Principes

Les paramètres scientifiques restent distincts des filtres appliqués aux résultats. Les commandes externes reçoivent une liste d’arguments, sans `shell=True`. Les moteurs interrompus sont terminés avec leur groupe de processus. Le shell HPC de l’utilisateur n’est pas modifié.

Une sauvegarde prépare le projet dans un dossier voisin, puis le renomme après une copie complète. Elle ne fusionne pas avec un dossier existant. La reprise valide les chemins du manifeste et les empreintes des entrées. Les anciens chemins exécutables ne sont pas réutilisés automatiquement sur la nouvelle machine.

Les projets chargent les FASTA et les alignements en mémoire. Un stockage indexé et une analyse incrémentale de jeux très volumineux restent un chantier distinct ; cette version n’annonce pas une capacité illimitée.

## Manuels et langues

Les guides source sont dans `docs/` et `docs/en/`, avec un index dans chaque langue. `scripts/generate_cli_docs.py` reconstruit les deux références CLI depuis le parseur. `scripts/check_docs.py` vérifie liens locaux, blocs de code et présence des ressources attendues. La construction du wheel installe ces guides dans `alignview_docs` ; le zipapp embarque les mêmes fichiers.

`tests/test_workflows.py` couvre notamment les chemins depuis le dossier de lancement, le HTML seul sans écrasement, la reprise avec paramètres explicites, les tables NCBI et leurs initiateurs/arrêts, les commandes VSEARCH natives et le sélecteur de requête dans un vrai terminal anglais. La vérification Chrome couvre aussi la traduction dynamique du HTML, la vue VSEARCH et la conservation des identifiants. Les logs bruts et arguments des moteurs ne sont pas traduits.

L’archive portable contient un seul dossier `AlignView/` avec sources, guides, exemples et outils de développement. Elle exclut archives historiques, préfixes scientifiques, Git, caches et sorties. Tester son extraction et son lancement depuis un autre dossier après toute modification du lanceur.
