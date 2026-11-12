# Analyses VSEARCH

VSEARCH est proposé pour les **séquences nucléotidiques** dans **Analyser → Analyses VSEARCH**. Le même moteur apparaît dans le choix des méthodes de recherche lorsque requête et références sont nucléotidiques. Les protéines sont refusées avant lancement.

Le moteur, sa version, son origine et son chemin sont visibles dans l’écran de préparation. **Voir les commandes et la préparation** ne lance aucun calcul scientifique. **Paramètres VSEARCH** permet de régler identité, brins, couvertures, masquage, heuristiques, longueurs et abondances. L’installation se prépare avec `alignview --setup vsearch` et s’exécute avec `--install-tools`.

## Quatre parcours

| Analyse | Séquences utilisées | Fichiers principaux |
| --- | --- | --- |
| Recherche globale | La requête active contre le second FASTA de références | `vsearch/search/matches.tsv` avec alignements et métriques VSEARCH |
| Déduplication | Toutes les séquences du premier FASTA | `vsearch/dereplicate/unique.fasta`, `dereplication.uc` |
| Regroupement | Toutes les séquences du premier FASTA | `vsearch/cluster/centroids.fasta`, `clusters.uc` |
| Chimères de novo | Toutes les séquences du premier FASTA, dédupliquées avant UCHIME | `vsearch/chimeras/unique.fasta`, `chimeras.fasta`, `nonchimeras.fasta`, `borderline.fasta`, `uchime.tsv` |

Ces analyses ne remplacent pas le FASTA source, ne changent pas le MSA courant et ne convertissent pas les résultats en HSP BLAST. Plusieurs analyses VSEARCH peuvent coexister dans le même projet. Refaire une action remplace son résultat après réussite ; les journaux des tentatives sont conservés.

## Recherche : sélectionner la bonne requête

Si le premier FASTA contient plusieurs séquences, l’interface demande maintenant laquelle utiliser, y compris dans le parcours guidé **Comparer à des références**. `/` filtre la liste par identifiant. **Données → Changer la requête active** permet d’en choisir une autre ensuite.

En batch, le choix est explicite :

```bash
alignview queries.fna refs.fna --vsearch search --query-id ma_sequence --vsearch-id 0.97 --report --save --out recherche_vs_01 --no-menu
```

La recherche est **globale entre séquences**, avec un préfiltre de mots. Elle ne cherche pas toutes les occurrences d’un motif à l’intérieur d’une longue référence. `--vsearch-id` est une fraction entre 0 et 1 ; l’identité rapportée est un pourcentage selon `--vsearch-iddef` (2 par défaut). Les colonnes qcov et tcov sont calculées par VSEARCH et ne sont pas une union de HSP locaux.

VSEARCH ne calcule pas d’E-value ou de bitscore pour ces alignements nucléotidiques. AlignView n’affiche donc pas de valeurs fictives. Les filtres BLAST d’E-value/bitscore ne s’appliquent pas à ces résultats. `--hits` fixe le maximum de références rapportées, indépendamment de l’arrêt sur candidats acceptés/rejetés.

Pour une petite recherche exhaustive explicite, désactiver les trois arrêts/préfiltres :

```bash
alignview query.fna refs.fna --vsearch search --vsearch-id 0.9 --vsearch-maxaccepts 0 --vsearch-maxrejects 0 --vsearch-minwordmatches 0 --vsearch-masking none --save --out vs_exhaustif_01 --no-menu
```

Cette combinaison peut coûter beaucoup plus cher sur une grande base. Abaisser seulement l’identité ne désactive pas le préfiltre de mots. Les valeurs par défaut et leur effet sont dans [CLI](CLI.md).

## Séquences uniques, groupes et abondances

```bash
alignview reads.fna --vsearch dereplicate --save --out uniques_01 --no-menu
alignview reads.fna --vsearch cluster --vsearch-id 0.97 --report --save --out groupes_01 --no-menu
alignview uniques_avec_sizes.fna --vsearch cluster --vsearch-sizein --save --out groupes_02 --no-menu
```

Les abondances sont exportées sous la forme `;size=N;`. Par défaut, chaque entrée FASTA compte pour une séquence. Activer **Lire les abondances des en-têtes** / `--vsearch-sizein` si les entrées possèdent déjà des abondances à additionner. Ne pas présenter un groupe à 97 % comme une espèce identifiée : c’est un regroupement selon l’algorithme et les réglages choisis.

## Chimères

```bash
alignview reads.fna --vsearch chimeras --report --save --out chimeres_01 --no-menu
```

AlignView déduplique d’abord les séquences, produit leurs abondances puis lance UCHIME de novo. Le rapport sépare séquences signalées chimériques, non chimériques et indéterminées. Les abondances, les parents présents, le masquage et les seuils influencent le résultat. Des références uniques d’abondance uniforme ne remplacent pas un jeu d’amplicons avec ses comptages ; aucune chimère détectée ne signifie pas une garantie d’absence de chimère. UCHIME de novo travaille sur le brin plus, indépendamment du choix de brin appliqué aux autres étapes.

## Préparation, limites et conservation

Les gaps d’un FASTA aligné sont retirés et U est normalisé en T pour préparer l’entrée VSEARCH ; cette préparation est annoncée et enregistrée. Les limites de longueur sont explicites (`--vsearch-min-length`, `--vsearch-max-length`) et les exclusions sont comptées. Les entrées internes sont normalisées en majuscules : seuls `dust` et `none` sont proposés, sans prétendre conserver un masquage préexistant en minuscules.

Les résumés, paramètres, commandes, versions et fichiers sont conservés dans le projet, le HTML et le dossier `vsearch/` des exports. Le HTML montre au plus 200 lignes par action ; les TSV/UC/FASTA du projet et des exports restent complets. Le choix de table génétique ne traduit pas les données VSEARCH : ce moteur reçoit des nucléotides.

Documentation du moteur : [recherche globale](https://torognes.github.io/vsearch/commands/vsearch-usearch_global.1.html), [champs de sortie](https://torognes.github.io/vsearch/misc/vsearch-userfields.7.html), [VSEARCH et ses commandes](https://torognes.github.io/vsearch/).
