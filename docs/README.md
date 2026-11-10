# Documentation d’AlignView

AlignView sert à explorer des séquences, rechercher des correspondances, créer des alignements et partager les résultats. Il fonctionne dans un terminal Linux, y compris par SSH. Le rapport HTML se consulte dans un navigateur.

## Partir de son besoin

| Je veux… | Guide | Accès depuis la commande |
| --- | --- | --- |
| Essayer sans installer de moteur | [Premiers pas](QUICKSTART.md) | `alignview --docs demarrage` |
| Trouver une fonction et comprendre ses résultats | [Tous les parcours](USAGE.md) | `alignview --docs parcours` |
| Savoir où sont l’application et mes sorties | [Dossiers et sauvegardes](FILES.md) | `alignview --docs fichiers` |
| Installer, déplacer ou utiliser sur HPC | [Installation](INSTALLATION.md) | `alignview --docs installation` |
| Choisir des paramètres et comprendre leurs effets | [Paramètres scientifiques](PARAMETERS.md) | `alignview --docs parametres` |
| Lire le HTML, l’arbre et les figures | [Rapport HTML](HTML.md) | `alignview --docs html` |
| Exporter vers R, Python ou un autre logiciel | [Formats et coordonnées](EXPORTS.md) | `alignview --docs exports` |
| Automatiser ou retrouver une option | [Référence CLI complète](CLI.md) | `alignview --docs cli` |
| Résoudre une erreur | [Dépannage](TROUBLESHOOTING.md) | `alignview --docs depannage` |
| Comprendre FASTA, HSP, MSA, E-value… | [Glossaire](GLOSSARY.md) | `alignview --docs glossaire` |

Dans ces exemples, `alignview` est la commande installée. Avec le dossier portable, la remplacer par `/chemin/vers/AlignView/alignview.sh`. Aucune installation n’est nécessaire pour consulter la documentation ou la démonstration.

## Depuis l’application

Au démarrage, choisir **Lire la documentation**. Dans un projet ouvert : **Réglages et aide → Documentation complète**. Les guides sont embarqués dans le dossier portable, le wheel et le fichier `.pyz` ; la consultation ne dépend pas d’Internet. Dans la page de lecture, les flèches déplacent le texte, `/` recherche et Échap revient au choix des guides.

## Quatre distinctions utiles

1. Un **FASTA** contient des séquences. Un **MSA** est un alignement de ces séquences, avec des colonnes comparables et éventuellement des gaps.
2. Une **recherche de similarités** trouve des segments locaux. Elle ne construit pas automatiquement un alignement multiple de séquences entières.
3. Un **projet** permet de reprendre les calculs ; un **HTML** présente les résultats ; des **exports** permettent de les réutiliser ailleurs.
4. Les **réglages actuels** préparent le prochain calcul. Les **paramètres historiques** décrivent le calcul ayant réellement produit un résultat existant.

## Pour contribuer ou publier

[Développement et tests](DEVELOPMENT.md) · [Publication](RELEASING.md) · [Audit et limites](AUDIT.md) · [Revue de la navigation](UX_REVIEW.md) · [Revue des calculs](CALCULATION_REVIEW.md) · [Contribuer](../CONTRIBUTING.md).

La version locale est 0.1.0. Les résultats des moteurs réels et les limites de validation sont décrits dans l’audit ; cette documentation ne présente pas la préparation locale comme une publication PyPI ou Bioconda.
