# Glossaire et exemples

| Terme | Sens dans AlignView |
| --- | --- |
| FASTA | Fichier texte contenant des en-têtes `>identifiant description` puis des séquences. Il peut contenir une ou plusieurs séquences. |
| Identifiant | Premier mot après `>` ; c’est lui qu’utilise `--query-id`. Les identifiants d’un même fichier doivent distinguer les séquences. |
| Requête / QUERY | Séquence active comparée aux références. Dans un FASTA unique, elle est aussi la cible des analyses interactives de motifs, traduction et extraction. |
| Références | Séquences auxquelles la requête est comparée. Elles ne constituent pas forcément une base taxonomique validée. |
| Hit / correspondance | Dans l’interface et `--hits`, référence retenue après la recherche et ses filtres. Une référence peut avoir plusieurs HSP. |
| HSP | Segment d’alignement local rapporté par le moteur. Il possède ses coordonnées, son identité, son score et éventuellement son cadre. |
| Identité locale | Pourcentage de caractères identiques dans un HSP ; ce n’est pas un pourcentage d’identité sur toute la requête. |
| Couverture de la requête | Part des positions de la requête touchées par l’union des HSP retenus pour une référence. Les chevauchements ne sont comptés qu’une fois. |
| Couverture de la référence | Même calcul, avec la longueur et les coordonnées de cette référence comme dénominateur et repère. |
| E-value | Quantité attendue de correspondances de score comparable dues au hasard dans l’espace de recherche du moteur. Un seuil plus petit est plus strict ; il ne s’agit pas d’une probabilité d’appartenance taxonomique. |
| Bitscore | Score local normalisé rapporté par le moteur. Une valeur plus élevée indique un meilleur score dans ce cadre de calcul. |
| MSA | Alignement multiple : chaque colonne met en regard les caractères des séquences ; les gaps sont conservés. |
| Gap | Caractère d’alignement, généralement `-`, représentant une absence de résidu dans une colonne donnée. |
| Consensus | Résumé par colonne calculé à partir des séquences alignées. Il ne constitue pas nécessairement une séquence biologique observée. |
| Conservation | Similarité des caractères observés dans une colonne. Vérifier les gaps et le nombre de séquences comparables. |
| Entropie | Mesure de la diversité des caractères. Une valeur élevée peut décrire de la variation, de l’ambiguïté ou un alignement hétérogène. |
| Distance p | Fraction de différences dans les positions comparables d’une paire, en excluant gaps/caractères manquants. Ce n’est pas une distance corrigée par un modèle évolutif. |
| nt / aa | Nucléotides / acides aminés. Les unités doivent rester distinctes en recherche traduite. |
| Cadre de lecture | Décalage et orientation de traduction d’une séquence nucléotidique : +1, +2, +3, −1, −2 ou −3. |
| ORF | Ici, segment simple ATG → STOP complet dans l’un des six cadres, selon la table NCBI choisie (1 par défaut). Ce n’est pas une annotation de gène. |
| K-mer | Mot de longueur k prélevé dans une séquence ; les dotplots cherchent des mots exactement identiques. |
| Newick | Format textuel d’arbre avec parenthèses, noms, longueurs après `:` et éventuellement étiquettes internes. |
| Support SH-like | Support local calculé par FastTree. Ne pas le présenter comme un bootstrap d’arbres. |
| Phylogramme | Arbre dont les longueurs horizontales sont proportionnelles aux longueurs de branches. |
| Cladogramme | Dessin de la topologie sans échelle de distance ; utilisé si les longueurs ne permettent pas un phylogramme. |
| Provenance | Paramètres, entrées, versions, commandes et journaux qui décrivent comment les résultats ont été produits. |
| SHA256 | Empreinte d’un fichier permettant de vérifier qu’il n’a pas changé. Ce n’est pas une signature d’auteur. |
| Profil de seuil | Standard, Plus permissif ou Plus strict : règle notamment l’E-value et la sensibilité DIAMOND. |
| Profil de calcul JSON | Ensemble réutilisable de paramètres ; ne contient pas les séquences. |
| Profil d’export | Standard, Publication ou Complet : détermine les fichiers préparés pour la sauvegarde ou l’export. |

## Exemple : forte identité, faible couverture

Une requête de 1 000 nt possède un HSP de 100 nt avec 99 % d’identité. La couverture est proche de 10 %, malgré une très forte identité sur ce segment. Cela ne suffit pas à conclure que les deux séquences sont identiques à 99 % sur toute leur longueur.

## Exemple : deux segments qui se chevauchent

Les HSP couvrent les positions 1–100 et 81–180 d’une requête de 200 nt. Leur union comporte 180 positions, donc une couverture de 90 %. Additionner 100 + 100 aurait compté deux fois les 20 positions du chevauchement.

## Exemple : coordonnées et traduction

Les positions 3–62 incluses couvrent 60 nt, ce qui peut correspondre à 20 aa traduits. En BLASTX, les coordonnées QUERY sont en nt, celles de la référence en aa, et les colonnes du HSP en aa. Les coordonnées d’un MSA restent ses colonnes, gaps compris.

Un intervalle 1-based inclusif `3–62` devient l’intervalle BED 0-based half-open `[2, 62)`. Un résultat sur le brin inverse peut avoir un début numérique supérieur à la fin dans le HSP ; consulter le brin et les fichiers d’audit avant toute conversion.
