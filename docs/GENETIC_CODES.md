# Choisir une table de traduction

AlignView embarque les **27 tables actives du fichier NCBI `gc.prt`** consulté le 10 novembre 2026. La table **1, Standard**, reste la valeur de départ. Le numéro doit correspondre à l’organisme, au compartiment et à l’annotation biologique étudiés ; il n’est pas déduit automatiquement d’une séquence.

## Choisir dans l’interface

Dans **Analyser → Autres analyses → Traduction et régions codantes → Choisir la table de traduction**, sélectionner le numéro **et son nom NCBI**. Le même réglage est proposé avant les calculs traduits dans **Paramètres du moteur / Paramètres de calcul**.

- **Table QUERY** : traductions internes de la requête, ORF, BLASTX/DIAMOND BLASTX et côté requête de TBLASTX.
- **Table des références** : héritée de QUERY par défaut ; peut être différente pour TBLASTN/TBLASTX, Needle/Water traduits et dotplot traduit côté référence.
- **Codons de départ ORF** : `atg` conserve ATG uniquement ; `table` utilise les initiateurs listés par NCBI pour la table choisie.

```bash
alignview --list-genetic-codes
alignview query.fna refs.faa --mode blastx --genetic-code 4 --save --out traduction_04 --no-menu
alignview query.faa refs.fna --mode tblastn --reference-genetic-code 11 --save --out recherche_11 --no-menu
alignview --write-profile code_bacterien.json --genetic-code 11 --orf-starts table
```

Les tableaux locaux sont disponibles hors ligne dans toutes les distributions. Les profils, projets, commandes et blocs de provenance conservent le choix. Les moteurs reçoivent leurs options natives `-query_gencode`, `-db_gencode` ou `--query-gencode` lorsqu’elles s’appliquent. Une version de moteur incompatible avec une table peut la refuser : son journal donne alors la cause ; le logiciel ne remplace pas silencieusement la table par le code standard.

## Effet concret

`TGA` produit un STOP avec la table 1 et un tryptophane `W` avec les tables 2 ou 4. `ATA` produit une isoleucine avec la table 1 et une méthionine avec la table 2. La traduction de six cadres, les protéines utilisées pour les MSA traduits, les comparaisons Needle/Water et les dotplots traduits utilisent tous le réglage correspondant.

Les codons ambigus IUPAC sont résolus lorsque toutes leurs interprétations donnent le même résidu ; sinon `X` est produit. Un codon incomplet en fin de séquence est ignoré. Les traductions six cadres conservent les STOP `*` ; les préparations de MSA et de dotplot traduits les remplacent par `X` en l’annonçant.

## Initiation et contexte

La traduction d’un fragment ne réécrit pas arbitrairement le premier résidu en méthionine : ce n’est pas une traduction CDS validée. Dans un ORF détecté avec `--orf-starts table`, un initiateur alternatif est rendu `M` en tête de la protéine, conformément au rôle d’initiation choisi. Le détecteur reste un outil simple entre initiateur et STOP ; il ne traite pas les introns, l’édition ARN ou les règles biologiques contextuelles.

Les tables 27, 28 et 31 comportent des codons dont le rôle d’arrêt dépend du contexte. Leur traduction interne sans contexte utilise l’acide aminé de `ncbieaa`. La détection simple d’ORF les refuse explicitement, car elle ne dispose pas de l’annotation CDS nécessaire pour décider de ces arrêts. Ce refus ne concerne pas l’affichage de ces tables ni la traduction simple de fragments.

Une modification de table prépare les calculs suivants. Elle ne réinterprète pas rétroactivement les résultats existants ; relancer le calcul concerné et reconstruire le HTML. Consulter les paramètres historiques pour savoir quelle table a réellement produit un résultat.

Sources : [codes génétiques NCBI](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), [données `gc.prt`](https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt). La source, la date et son SHA256 sont conservés dans `av_genetics.py`.

Les analyses internes de traduction et d’ORF retirent les gaps avant traitement : leurs coordonnées correspondent aux séquences préparées, pas aux colonnes d’un MSA. Conserver le FASTA source et la préparation pour remettre ces positions en contexte.
