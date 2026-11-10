# Dépannage

## Commencer par ces informations

```bash
pwd
/chemin/vers/AlignView/alignview.sh --version
/chemin/vers/AlignView/alignview.sh --diagnose
/chemin/vers/AlignView/alignview.sh --docs fichiers
```

`pwd` explique la base des chemins relatifs. Le diagnostic liste les exécutables détectés, leur provenance, leur version et leur profil d’installation. `--diagnose-json` donne un JSON exploitable automatiquement. Il n’installe rien et ne charge pas de module HPC pour ce diagnostic.

## L’application ne démarre pas

| Message ou symptôme | Cause à vérifier | Action |
| --- | --- | --- |
| `No such file or directory` dans un chemin de lancement | Dossier déplacé, mauvais nom ou casse | Utiliser le lanceur **dans `AlignView/`**, avec la majuscule réelle ; mettre à jour son alias éventuel |
| `Permission denied` sur le lanceur | Bit exécutable perdu lors d’une copie | `bash /chemin/vers/AlignView/alignview.sh --demo` ou `chmod +x /chemin/vers/AlignView/alignview.sh` |
| Python absent ou trop ancien | Il faut Python 3.10+ | Vérifier `python3 --version` ; sur HPC, utiliser le module Python disponible |
| `ModuleNotFoundError: av_…` | Un seul fichier Python a été copié | Copier tout `AlignView/`, ou utiliser le `.pyz` complet |
| Pas de terminal interactif | Sortie redirigée, tâche batch ou IDE sans TTY | Fournir les FASTA et `--no-menu`, ou ouvrir un terminal ; `--docs` fonctionne aussi sans TTY |
| Affichage de touches étrange | Terminal ou transport des touches | Essayer `--simple-ui`, agrandir le terminal et consulter `?` dans la visionneuse |

## Le moteur est absent ou échoue

1. Lire la commande manquante dans le message, puis le diagnostic.
2. Si les programmes sont déjà installés, fournir `--tools-dir /chemin/bin` ou `--executable mafft=/chemin/mafft`.
3. Sur HPC, charger le module correspondant, ou laisser `--module-policy auto` le rechercher à la demande. `never` évite les tentatives de chargement.
4. Pour une installation dédiée, examiner `--setup mafft`, puis lancer `--setup mafft --install-tools`. Il faut Micromamba, Mamba ou Conda et un accès aux canaux.

Un moteur lancé mais en échec n’est pas un moteur absent. Lire son journal stderr dans `logs/`, puis ses paramètres réellement exécutés dans la provenance. Une matrice ou des coûts de gaps incompatibles sont refusés par le moteur ; modifier les paramètres dans l’écran de préparation.

Pour MAFFT, un nombre explicite de raffinements demande `--mafft global`, `local` ou `fast`, et non `auto`. Pour MUSCLE, utiliser une version 5. Pour les modèles FastTree, GTR est réservé à l’ADN et WAG/LG aux protéines.

## Problèmes de données et de résultats

| Situation | Explication et suite |
| --- | --- |
| Plusieurs requêtes sans `--query-id` en batch | Choisir le premier mot de l’en-tête FASTA voulu ; une recherche porte sur une seule requête active |
| Un FASTA unique et un mode BLASTX/BLASTN… | Une recherche demande un second FASTA de références ; pour un seul fichier, utiliser l’exploration ou le MSA |
| Longueurs différentes avec `--msa view` | Le fichier n’est pas un MSA rectangulaire ; créer un alignement avec `--msa align` |
| ADN et protéines mélangés | Séparer les types ou préparer une traduction biologique appropriée ; ne pas aligner arbitrairement ces alphabets |
| Recherche terminée, aucune correspondance | Vérifier le mode, les séquences, l’E-value et les filtres identité/couverture/bitscore ; absence de hit n’est pas absence prouvée d’homologie |
| Des hits disparaissent après filtre | La sélection a changé ; le MSA, l’arbre et les vues dérivées peuvent être invalidés. Refaire les analyses nécessaires |
| Modifier E-value ne change pas les hits affichés | C’est un paramètre de calcul : relancer la recherche. Les filtres de résultats réutilisent les HSP déjà calculés |
| Moins de 100 références | 100 est un maximum de références retenues ; les données, le moteur, les candidats et les filtres peuvent en retenir moins |
| Limite conservée à 40 après reprise d’un ancien projet | La reprise conserve les réglages du projet. Choisir 100 dans les réglages ou ajouter explicitement `--hits 100` |
| Coordonnées BLASTX apparemment trois fois trop longues | QUERY est en nt et l’alignement HSP en aa ; consulter `coordinate_units.tsv` et `coordinate_audit.tsv` |

## Le HTML ou l’arbre n’apparaît pas

- **Aucun fichier après `--report --no-menu`** : le rapport était temporaire. Ajouter `--html-out rapport.html`, ou `--report --save --out projet`.
- **Un projet, mais pas de HTML** : générer le rapport explicitement ou choisir le profil d’export publication/complet.
- **Pas de figure d’arbre** : vérifier qu’un arbre a été calculé, qu’il existe au moins trois séquences alignées et que `tree_fasttree.nwk` est présent. Le HTML ne calcule pas lui-même l’arbre.
- **Projet ancien avec Newick seulement** : `alignview --resume projet --html-out actualise.html --no-menu` crée la figure sans relancer FastTree.
- **Figure d’arbre indisponible** : lire le message ; un Newick incomplet, plusieurs arbres concaténés ou des longueurs invalides empêchent son dessin.
- **Vieille version d’un résultat dans le navigateur** : reconstruire le HTML, ouvrir le nouveau fichier et recharger la page. Un HTML déjà exporté est un instantané.
- **Onglets ou boutons inactifs** : ouvrir le fichier dans un navigateur avec JavaScript activé, plutôt que dans un aperçu d’IDE qui désactive les scripts.
- **Pas d’ouverture automatique sur SSH** : récupérer le fichier avec `scp` ou l’explorateur distant, puis l’ouvrir localement. Aucun serveur temporaire n’est nécessaire.

## Sauvegarde et reprise

Un dossier de destination existant n’est pas écrasé : choisir `projet_02` ou une autre destination. `--out` seul ne conserve pas le projet ; ajouter `--save`. Un HTML isolé n’est pas un dossier de reprise.

Si la reprise indique qu’un fichier a changé, les empreintes d’entrée ne correspondent plus. Reprendre une copie intacte ou ouvrir les nouvelles données comme une nouvelle analyse. Ne pas modifier arbitrairement les SHA256 du manifeste pour masquer le changement.

Après une erreur, chercher le chemin de récupération affiché et reprendre `alignview --resume /chemin/alignview_recovery_…/project`. Après un arrêt brutal, la récupération peut ne pas exister. Un enregistrement régulier reste le moyen de conserver le travail.

## Informations utiles pour signaler un problème

Conserver la version AlignView/Python, la commande ou le parcours suivi, le message exact, le journal du moteur, les paramètres réellement exécutés et un petit FASTA reproduisant le problème lorsque ses données sont partageables. Voir [contribuer](../CONTRIBUTING.md). Les versions et chemins se trouvent dans le projet et dans le diagnostic ; ne pas attribuer une erreur de moteur à une comparaison scientifique non testée.
