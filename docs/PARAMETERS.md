# Régler et comprendre un calcul

Pour trouver une option rapidement, consulter la [référence CLI complète](CLI.md) ou lancer `alignview --docs cli`. Cette page explique leur sens ; les parcours complets sont dans [Utiliser AlignView](USAGE.md).

Avant de lancer une analyse, ouvrir **Paramètres de calcul** (alignement, Needle/Water ou arbre) ou **Paramètres du moteur** (recherche). Seuls les réglages applicables sont présentés. **F1** donne leur effet ; `/` retrouve un paramètre. La commande peut être inspectée depuis l’écran de préparation, sans lancer de moteur.

## Calcul et filtrage

| Réglage | Effet | Recalcul nécessaire ? |
| --- | --- | --- |
| Graine, tâche, matrice, gaps, masquage, sensibilité | Change les alignements que le moteur peut détecter et leurs scores | Oui |
| E-value de recherche | Seuil transmis au moteur | Oui pour récupérer des résultats exclus |
| Identité, couverture et bitscore AlignView | Filtre les HSP déjà rapportés ; couverture par union des HSP retenus | Non |
| Maximum de références retenues | Classe et limite la sélection après recherche ; 100 par défaut | Non, dans la limite des résultats déjà disponibles |
| Limite de candidats du moteur | Limite en amont, propre au moteur ; peut changer les résultats | Oui |
| Stratégie MSA, coûts de gap, raffinements | Change les colonnes de l’alignement multiple | Oui |
| Modèle et support FastTree | Change l’inférence ou l’évaluation de l’arbre | Oui |

Les filtres ne récupèrent pas les correspondances que le moteur n’a pas produites. La limite automatique de candidats vaut `max(100, 10 × hits)` : 1000 pour 100 références retenues. `--search-candidates` la remplace. BLAST/DIAMOND reçoivent leur limite de cibles ; MMseqs2 reçoit `--max-seqs`, un plafond de candidats du préfiltre, pas une garantie sur le nombre de résultats. Ces contrôles n’ont pas une signification algorithmique identique.

Les profils de seuil `standard`, `sensible` et `strict` règlent E-value à `1e-5`, `1e-2` et `1e-10`. Ils choisissent aussi la sensibilité DIAMOND si ce moteur est utilisé. Ils ne règlent pas les graines de BLAST ou `-s` de MMseqs2. Le nombre de hits reste indépendant.

## Valeurs de départ

**Défaut du moteur installé** signifie que l’option est omise. Le résultat dépend alors des défauts de cette version. Les versions et les commandes sont enregistrées lors de l’exécution ; AlignView ne prétend pas développer toutes les décisions internes du moteur.

Les décisions d’AlignView sont explicites : MAFFT `--auto`, G-INS-i/L-INS-i avec 1000 raffinements si choisis, DIAMOND `--sensitive`, BLASTX DIAMOND `--min-orf 1`, Needle avec extrémités pénalisées, 4 CPU pour les outils acceptant les threads. FastTree standard reste monothread. Pour rétablir les réglages d’un calcul, utiliser **Rétablir les valeurs de départ** dans son écran de paramètres.

Pour MAFFT, `--maxiterate` explicite exige une stratégie `global`, `local` ou `fast` : `auto` choisit sa stratégie. Le mode `fast` utilise FFT-NS-2 (`--retree 2 --maxiterate 0`, sauf raffinement explicite). Les coûts `--op` et `--ep` sont propres à MAFFT ; ils ne correspondent pas aux coûts EMBOSS ou BLAST.

L’ajout de fragments impose **MAFFT en deux étapes** : alignement des références avec la stratégie réglée, puis `--addfragments --auto`. Les réglages de gaps s’appliquent aux deux étapes ; le nombre de raffinements ne s’applique qu’au premier alignement. Le MSA des HSP respecte désormais le moteur, la stratégie et le recouvrement minimal de l’ancre choisis.

Needle et Water ont leurs propres coûts et matrice. Les pénalités terminales sont réservées à Needle. En l’absence de coûts terminaux explicites, AlignView reprend `gapopen` et `gapextend`. **Extrémités libres** désactive leur pénalisation. Ce choix diffère du défaut EMBOSS natif et peut fortement changer l’alignement.

FastTree choisit JC pour l’ADN et JTT pour les protéines en mode `auto`. GTR s’applique à l’ADN ; JTT, WAG et LG aux protéines. Gamma20 réévalue et redimensionne les branches après l’optimisation CAT ; ce n’est pas une nouvelle recherche de topologie sous Gamma. Les rééchantillonnages du support local SH-like ne constituent pas un bootstrap d’arbres.

Les traductions utilisent uniquement le **code génétique standard (1)**. `--frame` choisit le cadre de QUERY dans les analyses complètes traduites, pas les cadres recherchés par BLASTX. `--sframe` choisit celui de la référence pour Needle/Water. Sans cadre forcé : meilleur HSP, sinon +1. Pour un autre code génétique, fournir des protéines traduites avec le code approprié. Le détecteur d’ORF reste un détecteur simple ATG→STOP.

## Profils portables

Dans **Profils de calcul**, enregistrer un nouveau JSON ou en charger un existant. Le profil contient les réglages des moteurs et les filtres ; les paramètres des moteurs inactifs restent disponibles pour une utilisation ultérieure. Il ne contient ni séquences ni chemins d’exécutables. Les versions des moteurs doivent être contrôlées sur chaque poste.

```bash
# Préparer un profil sans lancer un calcul
alignview --write-profile recherche.json --blast-word-size 7 --blast-masking yes

# Les options explicites gagnent sur le profil
alignview query.fna refs.fna --mode blastn --task blastn \
  --load-profile recherche.json --hits 100 --save --no-menu

# Reprendre un projet avec de nouveaux réglages, sans relancer ses calculs
alignview --resume projet --load-profile recherche.json
```

Les fichiers de profil inconnus, les options inconnues, les valeurs hors plage et les types invalides sont rejetés avant application. Un enregistrement n’écrase pas de fichier existant. Charger un profil ne recalcule pas les résultats ; si les filtres changent la sélection, les dérivés devenus invalides sont retirés. Les projets repris sans profil conservent leurs réglages enregistrés, y compris une ancienne limite de hits.

## Provenance

**Réglages → Calculs effectués et provenance** distingue les paramètres d’un calcul réussi des valeurs préparées pour le prochain lancement. Une recherche dont les paramètres changent est signalée comme à relancer. Un ancien alignement reste consultable avec son ancienne configuration identifiée.

`parameters.json`, `session.json` et les exports conservent les paramètres actuels et les blocs `search_run` / `calculation_runs`. Chaque commande exécutée conserve ses arguments, sa version, son statut, sa durée, le journal stderr et, si utilisé, le fichier stdout. Les projets et rapports gardent ces informations ; les scripts Slurm reprennent les réglages scientifiques supplémentaires. Les chemins temporaires des commandes originales sont des traces d’exécution ; les entrées sauvegardées et le script Slurm servent à préparer une nouvelle exécution portable.

L’aperçu utilise les mêmes constructeurs de commandes que l’exécution. Il peut préparer des FASTA intermédiaires dans le dossier temporaire, sans calcul externe et sans remplacer l’alignement actif. La résolution d’un exécutable via module HPC peut encore changer son chemin au lancement. Les matrices/coûts incompatibles sont refusés par le moteur dans son journal : l’aperçu ne remplace pas ses validations.

## Référence CLI

Les options d’un autre moteur ne sont pas transmises au calcul courant. `--help` décrit aussi chaque effet. Les options booléennes supplémentaires acceptent une forme `--no-…`.

| Option | Valeur de départ | Effet |
| --- | --- | --- |
| `--threads` | 4 | Ressources des moteurs qui acceptent des threads ; FastTree standard reste monothread. |
| `--evalue` | 1e-05 | Seuil transmis au moteur. Le modifier nécessite une nouvelle recherche. |
| `--search-candidates` | Défaut natif / automatique | Auto = max(100, 10 × références retenues). Ce plafond peut changer les résultats ; distinct des 100 hits affichés. |
| `--task` | défaut natif | Change graines et score. Sans option : BLASTN utilise megablast ; autres programmes : tâche native. |
| `--sensitivity` | sensitive | Recherche plus sensible = davantage de calcul ; ce réglage est propre à DIAMOND. |
| `--blast-word-size` | Défaut natif / automatique | Petites graines : détection plus sensible, calcul plus coûteux ; valeurs autorisées selon la tâche. |
| `--blast-reward` | Défaut natif / automatique | Gain pour une base identique ; dépend des coûts de mismatch et de gap. |
| `--blast-penalty` | Défaut natif / automatique | Score négatif pour une substitution ; change le score et les statistiques. |
| `--blast-matrix` | Défaut natif / automatique | Ex. BLOSUM62, BLOSUM45, PAM30. Les coûts de gap doivent être compatibles. |
| `--blast-gapopen` | Défaut natif / automatique | Coût propre à BLAST, indépendant de Needle/Water. TBLASTX est non gappé. |
| `--blast-gapextend` | Défaut natif / automatique | Coût propre à BLAST. Seules certaines combinaisons matrice/coûts sont acceptées. |
| `--blast-masking` | natif | Masque la faible complexité. Le désactiver peut révéler des correspondances peu spécifiques. |
| `--blast-soft-masking` | natif | true : masque les graines ; false : masque aussi l’extension. Défaut dépendant du programme. |
| `--blast-strand` | both | both : deux brins ; plus/minus : limite la recherche à un brin. |
| `--blast-max-hsps` | Défaut natif / automatique | Limite les segments rapportés, donc potentiellement la couverture calculée ensuite. |
| `--blast-comp-stats` | natif | 0 : désactivées ; 1/2/3 : corrections NCBI selon le moteur. Change scores et E-values. |
| `--diamond-matrix` | Défaut natif / automatique | Matrice protéique ; les coûts de gap doivent être compatibles. |
| `--diamond-gapopen` | Défaut natif / automatique | Coût propre à DIAMOND, indépendant des autres moteurs. |
| `--diamond-gapextend` | Défaut natif / automatique | Pénalise l’allongement des gaps. |
| `--diamond-masking` | natif | 0 : désactivé ; 1 : tantan ; seg : SEG (DIAMOND ≥ 2.0.12). |
| `--diamond-comp-stats` | natif | 0 : désactivée ; 1 : correction compositionnelle ; 3/4 : ajustements de matrice, selon version. |
| `--diamond-max-hsps` | Défaut natif / automatique | 0 : tous ; plusieurs HSP peuvent améliorer la couverture et augmenter le coût. |
| `--diamond-min-orf` | 1 | AlignView utilisait déjà 1, contrairement au choix automatique de DIAMOND. 0 : automatique natif. |
| `--mmseqs-sensitivity` | Défaut natif / automatique | De 1 (rapide) à 7.5 (sensible). Indépendante du réglage DIAMOND. |
| `--mmseqs-mask` | natif | 1 : faible complexité masquée ; 0 : désactivé. |
| `--mmseqs-comp-bias` | natif | 0/1 : désactive/active la correction du biais de composition. |
| `--aligner` | mafft | Les algorithmes et leurs paramètres ne sont pas interchangeables. |
| `--mafft` | auto | auto : choix natif ; global : G-INS-i ; local : L-INS-i ; fast : FFT-NS-2. |
| `--mafft-op` | Défaut natif / automatique | Pénalité de gap en alignement groupe à groupe ; indépendante d’EMBOSS. |
| `--mafft-ep` | Défaut natif / automatique | Paramètre d’offset de MAFFT ; ne se confond pas avec gapextend d’EMBOSS. |
| `--mafft-maxiterate` | Défaut natif / automatique | Global/local : 1000 par défaut AlignView. Pour régler ce nombre, choisir global, local ou fast plutôt que auto. |
| `--muscle-method` | align | align : PPP ; super5 : prévu pour des ensembles plus grands. MUSCLE 3 non pris en charge. |
| `--muscle-perturb` | 0 | 0 : HMM standard ; une autre graine explore la variabilité de l’alignement. |
| `--muscle-perm` | none | Change le parcours de l’arbre. Une seule permutation produit un seul MSA. |
| `--muscle-refineiters` | Défaut natif / automatique | Nombre d’itérations de raffinement ; coût et résultat peuvent changer. |
| `--clustalo-iterations` | 0 | Cycles combinés arbre guide / HMM. 0 : aucun raffinement itératif. |
| `--clustalo-full` | False | Matrice de distances complète plutôt que mBed ; coût quadratique. |
| `--clustalo-full-iter` | False | Matrice complète pendant les raffinements ; requiert des itérations. |
| `--hsp-overlap` | 0.5 | Fraction de la région requête ancre recouverte pour retenir un segment dans le MSA exploratoire. |
| `--frame` | Défaut natif / automatique | Vide : cadre du meilleur HSP, sinon +1. Ne réduit pas les six cadres recherchés par BLASTX. |
| `--sframe` | Défaut natif / automatique | Vide : cadre du meilleur HSP, sinon +1. Applicable à TBLASTN/TBLASTX. |
| `--gapopen` | 10 | Pénalité Needle/Water uniquement ; dépend de la matrice utilisée. |
| `--gapextend` | 0.5 | Pénalité Needle/Water uniquement. |
| `--matrix` | défaut natif | Vide : EDNAFULL en ADN, EBLOSUM62 en protéines (défauts EMBOSS). |
| `--free-ends` | False | AlignView pénalise les extrémités par défaut ; EMBOSS natif les laisse libres. Sans effet sur Water. |
| `--needle-endopen` | Défaut natif / automatique | Vide : même coût que gapopen AlignView ; ignoré si extrémités libres. |
| `--needle-endextend` | Défaut natif / automatique | Vide : même coût que gapextend AlignView ; ignoré si extrémités libres. |
| `--fasttree-model` | auto | auto : JC pour ADN, JTT pour protéines. GTR : ADN ; JTT/WAG/LG : protéines. |
| `--fasttree-gamma` | False | Après optimisation CAT, réévalue et redimensionne les branches sous Gamma20 ; ne refait pas la topologie. |
| `--fasttree-boot` | Défaut natif / automatique | Support SH-like, pas un bootstrap d’arbres. 0 : désactive les supports locaux. |
| `--fasttree-seed` | Défaut natif / automatique | Fixe la graine aléatoire pour les calculs stochastiques. |

## Sources des moteurs

## Modifier les réglages avec méthode

1. Vérifier les séquences, les alphabets et le résultat attendu avant de modifier un score.
2. Conserver une première analyse et sa provenance.
3. Modifier un paramètre ou un petit groupe cohérent, consulter la commande préparée, puis relancer le calcul concerné.
4. Comparer les résultats sur les mêmes entrées : nombre de références, couverture, régions alignées, gaps, qualité du MSA. Davantage de hits ne signifie pas automatiquement une meilleure analyse biologique.
5. Enregistrer le nouveau projet sous un autre nom et garder les versions des moteurs.

| Besoin | Où intervenir | Ce que cela ne fait pas seul |
| --- | --- | --- |
| Garder davantage de références | `--hits` ; vérifier aussi `--search-candidates` | Ne recrée pas des candidats que le moteur n’avait pas rapportés |
| Assouplir le seuil d’une recherche | E-value, graines/masquage, paramètres de sensibilité du moteur | Ne modifie pas une recherche déjà calculée sans relance |
| Filtrer une recherche existante | Identité, couverture, bitscore, maximum de références | Ne change pas le score des HSP d’origine |
| Modifier un MSA | Stratégie et paramètres du moteur MSA utilisé | Ne règle pas les gaps de BLAST ou d’EMBOSS |
| Changer la comparaison globale | Matrice/coûts EMBOSS, extrémités Needle | Ne produit pas un nouveau MSA |
| Modifier l’arbre | Modèle, Gamma20, supports et graine FastTree | Ne corrige pas un mauvais MSA et ne crée pas d’enracinement biologique |

En reprise, une option explicite telle que `--mafft-op 2` modifie les réglages futurs même sans profil JSON. Ajouter explicitement `--msa align` pour recalculer un MSA ; les paramètres historiques du résultat précédent restent conservés jusqu’à son remplacement. `--export-profile` concerne les fichiers, tandis que `--preset` concerne le seuil de recherche.

## Sources des moteurs

[BLAST+ NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279684/) · [DIAMOND](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options) · [MMseqs2](https://github.com/soedinglab/MMseqs2/wiki) · [MAFFT](https://mafft.ddbj.nig.ac.jp/alignment/software/manual/manual.html) · [MUSCLE align](https://www.drive5.com/muscle5/manual/cmd_align.html) · [MUSCLE Super5](https://www.drive5.com/muscle5/manual/cmd_super5.html) · [Clustal Omega](https://github.com/GSLBiotech/clustal-omega/blob/master/README) · [EMBOSS Needle](https://www.bioinformatics.nl/cgi-bin/emboss/help/needle) · [FastTree](https://morgannprice.github.io/fasttree/).

Le périmètre et les limites de validation sont décrits dans la [revue des calculs](CALCULATION_REVIEW.md).
