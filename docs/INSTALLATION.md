# Installer, lancer et déplacer AlignView

## Démarrer sans installation Python supplémentaire

La cible principale est Linux avec Python 3.10+. Copier **tout le dossier `AlignView/`**, puis lancer :

```bash
/chemin/vers/AlignView/alignview.sh --version
/chemin/vers/AlignView/alignview.sh --demo
/chemin/vers/AlignView/alignview.sh --docs demarrage
```

Depuis l’intérieur du dossier, `./alignview.sh` suffit. Le lanceur n’a besoin d’aucun script dans son dossier parent et ne change pas le dossier courant. Pour séparer application et résultats :

```bash
cd /chemin/vers/mes_donnees
/chemin/vers/AlignView/alignview.sh
```

Si le bit exécutable a été perdu, utiliser `bash /chemin/vers/AlignView/alignview.sh` ou le rétablir avec `chmod +x`. Il est aussi possible de lancer `python3 /chemin/vers/AlignView/alignview.py`.

L’application ne nécessite pas de bibliothèque Python tierce pour fonctionner. Les moteurs scientifiques sont facultatifs ; ils deviennent nécessaires seulement pour les calculs correspondants. Sous Windows, utiliser un environnement Linux/WSL ; cette release vise Linux. macOS et les terminaux cibles restent à valider séparément.

## Choisir un mode de distribution

| Mode | Usage | Ce qu’il faut conserver |
| --- | --- | --- |
| Dossier portable | Utiliser, déplacer ou développer sans pip | Le dossier `AlignView/` entier |
| Archive `AlignView-0.1.0-portable.zip` | Distribuer le dossier prêt à lancer | Extraire son dossier `AlignView/` |
| `alignview.pyz` | Un seul fichier Python autonome | Le `.pyz`, Python et les moteurs nécessaires |
| Wheel | Installer une commande `alignview` dans un environnement Python | L’environnement Python installé |

Le `.pyz` contient aussi les guides. Il ne contient pas d’exécutables BLAST/MAFFT/FastTree. Exemple : `python3 /chemin/alignview.pyz --docs fichiers`.

Pour une commande installée avec un wheel :

```bash
python3 -m venv /chemin/vers/venv_alignview
source /chemin/vers/venv_alignview/bin/activate
python -m pip install /chemin/vers/AlignView/dist/alignview-0.1.0-py3-none-any.whl
alignview --demo
```

Depuis les sources, `python -m pip install /chemin/vers/AlignView` fonctionne également avec un environnement de construction approprié. Cette installation est facultative pour utiliser le lanceur du dossier.

Le script `bash scripts/install_user.sh` copie volontairement le `.pyz` sous `~/.local/bin/alignview`. Il refuse un fichier existant sans `--force`. C’est une commodité optionnelle qui installe une copie hors du dossier ; ne pas l’utiliser si l’on veut conserver uniquement le dossier portable.

## De quels moteurs ai-je besoin ?

| Fonction | Moteur | Profil `--setup` |
| --- | --- | --- |
| Lire un FASTA/MSA, statistiques, motifs, traductions, HTML et SVG | Aucun | Aucun |
| Recherche ADN/ADN, protéine/protéine, recherches traduites NCBI | BLAST+ | `blast` |
| Recherche protéique / ADN traduit contre protéines | DIAMOND | `diamond` |
| Recherche protéine/protéine alternative | MMseqs2 | `mmseqs` |
| Alignement multiple et ajout de fragments | MAFFT | `mafft` |
| Alignement multiple alternatif | MUSCLE 5 ou Clustal Omega | `muscle`, `clustalo` |
| Comparaison globale/locale de deux séquences | EMBOSS Needle/Water | `emboss` |
| Calcul d’arbre | FastTree | `fasttree` |
| Visionneuse externe facultative | ALV | `alv` |

Le profil `minimal` regroupe BLAST, DIAMOND et MAFFT. Le profil `full` ajoute les autres moteurs. La visionneuse interne suffit sans ALV. MUSCLE 3 n’est pas compatible avec les commandes prévues pour MUSCLE 5.

## Installer les moteurs dans l’environnement dédié

Il faut déjà disposer de **Micromamba, Mamba ou Conda**. AlignView choisit le premier disponible et affiche la commande avant installation. Le plan sans `--install-tools` n’installe rien.

```bash
./alignview.sh --setup mafft
./alignview.sh --setup mafft --install-tools
./alignview.sh --setup minimal --install-tools
./alignview.sh --setup fasttree --install-tools
./alignview.sh --diagnose
```

Par défaut, le préfixe est `$XDG_DATA_HOME/alignview/tools`, ou `~/.local/share/alignview/tools`. Un préfixe existant n’est étendu que s’il contient `conda-meta`. Le shell de l’utilisateur n’est pas modifié et il n’est pas nécessaire d’activer cet environnement pour que l’application le détecte.

## Garder aussi les moteurs dans le dossier AlignView

Depuis le dossier de l’application :

```bash
./alignview.sh --setup full --tools-prefix "$PWD/tools"
./alignview.sh --setup full --tools-prefix "$PWD/tools" --install-tools
./alignview.sh --diagnose
```

Les sources, la documentation et cet environnement se trouvent alors tous sous `AlignView/`. Le lanceur détecte automatiquement `tools/conda-meta` lors des lancements suivants. Un `ALIGNVIEW_TOOLS_PREFIX` explicite ou `--tools-prefix` reste prioritaire. Le Python lancé directement et une commande wheel utilisent le préfixe explicite ou le préfixe utilisateur : fournir `--tools-prefix /chemin/AlignView/tools` pour choisir ce dossier.

Un environnement Conda contient des exécutables et parfois des chemins liés au poste d’installation. Le ranger dans l’application ne le transforme pas en environnement universellement transportable. Après un déplacement ou transfert, vérifier `--diagnose` et réinstaller le préfixe si nécessaire. Les distributions Python ne l’embarquent pas.

Pour un autre préfixe durable :

```bash
export ALIGNVIEW_TOOLS_PREFIX=/scratch/mon_compte/alignview-tools
/chemin/vers/AlignView/alignview.sh --setup minimal --install-tools
/chemin/vers/AlignView/alignview.sh --diagnose
```

## Utiliser les programmes déjà installés

Ordre de résolution : `--executable` explicite, dossiers `--tools-dir`, PATH existant, préfixe dédié, puis modules HPC à la demande lorsque cette politique est activée.

```bash
alignview --diagnose --tools-dir /opt/bio/bin
alignview seqs.fasta --msa align --executable mafft=/opt/mafft/bin/mafft --save --no-menu
alignview --diagnose-json --module-policy never
```

On peut répéter `--tools-dir` et `--executable`. Un chemin invalide est signalé avant le calcul. Les fichiers `environment-minimal.yml` et `environment-full.yml` décrivent aussi les profils Conda fournis.

## Linux distant et HPC

Les modules peuvent être chargés à la demande lorsqu’un outil manque avec `--module-policy auto`. Pour utiliser seulement les chemins, l’environnement dédié et les modules déjà chargés : `--module-policy never`. Le diagnostic normal ne charge pas de nouveaux modules ; l’écran **Logiciels et installation → Rechercher les modules HPC** permet une recherche explicite.

Les ressources de calcul se règlent avec `--threads`. FastTree standard reste monothread. Choisir `TMPDIR` et une destination de résultats adaptés au cluster ; voir [Dossiers et sauvegardes](FILES.md).

**Enregistrer et partager → Préparer un calcul sur HPC** produit un script Slurm dans le projet. Il n’est jamais soumis automatiquement. Vérifier partition, mémoire, durée, chemin de l’application et disponibilité des moteurs avant `sbatch`. Les entrées du script sont celles du projet sauvegardé.

Le HTML n’a pas besoin de serveur distant. Le transférer vers son poste avec `scp` ou l’explorateur distant et l’ouvrir dans un navigateur. Voir [Dépannage](TROUBLESHOOTING.md) pour les problèmes de moteurs, terminal ou chemins.
