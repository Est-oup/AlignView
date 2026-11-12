# Revue critique des calculs et de leurs réglages

Révision du 10 novembre 2026. Relecture des appels de moteurs, comparaison des documentations officielles, corrections du code et essais locaux. Il s’agit d’une revue fonctionnelle et d’intégration ; aucun classement de précision biologique entre logiciels n’est déduit de ces essais.

## Ce qui donnait un effet de boîte noire

| Constat dans le code | Conséquence | Correction |
| --- | --- | --- |
| E-value et quelques options générales accessibles ; nombreux paramètres natifs absents | Impossible de tester l’effet des graines, du masquage ou des coûts de gaps | Registre typé commun aux menus, CLI, profils et provenance ; réglages contextualisés par moteur |
| Sensibilité générique affichée alors que seul DIAMOND recevait ce choix | Réglage sans effet sur MMseqs2, confusion avec BLAST | Sensibilité numérique MMseqs2 dédiée, champ DIAMOND limité à DIAMOND, profils de seuil renommés |
| Limite de candidats calculée comme `max(100, hits × 10)` et ignorée par MMseqs2 | Plafond scientifique caché, comportement variable selon le moteur | Plafond visible et réglable ; transmission à `--max-seqs` MMseqs2, avec sa signification propre |
| Extraction des alignements MMseqs2 sans demande explicite du traceback | Calcul des séquences alignées insuffisamment explicite | `-a` demandé pour produire les alignements avec gaps exploités par AlignView |
| Paramètres `gapopen`, `gapextend`, `matrix` uniquement EMBOSS, sans distinction assez visible | Risque de croire avoir réglé BLAST ou le MSA | Noms séparés `--blast-*`, `--diamond-*`, `--mafft-*` et libellés explicites |
| G-INS-i/L-INS-i codés avec 1000 raffinements sans réglage dédié | Coût et stratégie peu contrôlables | Stratégie accessible, raffinements et gaps modifiables ; conflit `auto`/raffinements explicites refusé |
| MSA des HSP forçant MAFFT auto et seuil d’ancre fixé à 50 % | Stratégie choisie ignorée, sélection des régions implicite | Stratégie respectée, seuil d’ancre réglable, préparation inspectable |
| Fragments forçant MAFFT auto pour les deux étapes, même avec un autre moteur affiché | Moteur et méthode réels ambigus | MAFFT annoncé et imposé pour cette tâche ; références alignées avec la stratégie choisie, ajout automatique distinct |
| DIAMOND BLASTX avec `--min-orf 1` fixé | Écart au défaut natif peu visible | Choix explicité et modifiable ; 0 restaure l’automatisme natif |
| Needle pénalisant les extrémités avec les mêmes coûts que les gaps internes | Choix différent d’EMBOSS natif, important pour des séquences partielles | Extrémités libres et coûts terminaux séparés ; aucun réglage d’extrémités proposé à Water |
| FastTree sans choix du modèle ou de support | Difficile de savoir quel arbre a été inféré | JC/GTR ou JTT/WAG/LG selon alphabet, Gamma20, support SH-like et graine accessibles |
| Paramètres actuels pouvant remplacer le contexte d’un ancien résultat | Attribution scientifique ambiguë après modification | Instantané par calcul réussi, commandes/versions conservées, réglages modifiés signalés |
| Journaux stderr réécrits lors d’une relance | Une ancienne commande pouvait pointer vers un nouveau journal | Noms de journaux distincts à chaque exécution |
| Nouvelle sortie MSA/arbre écrite directement sur la précédente | Un échec de relance pouvait endommager un résultat encore actif | Sorties préparées séparément puis remplacées après succès et validation |

## Enseignements d’autres interfaces

| Documentation consultée | Fonction documentée | Application dans AlignView |
| --- | --- | --- |
| [Jalview, dialogue de paramètres](https://www.jalview.org/help/html/webServices/webServicesParams.html) | Modifier avant lancement, aide par paramètre, jeux de paramètres réutilisables | Accès aux paramètres dans la préparation du calcul, F1, profils JSON |
| [AliView, documentation officielle](https://www.ormbunkar.se/aliview/) | Préréglages de programmes d’alignement et configuration de commandes externes | Moteurs distincts, profils réutilisables et inspection des arguments générés |
| [UGENE, alignement MUSCLE](https://ugene.net/docs/extensions/umuscle/muscle-aligning/) | Dialogue avec modes adaptés aux tâches, dont grands alignements et raffinement | Séparer question biologique, portée des séquences et réglages du moteur ; MUSCLE 5 expose `align` et `super5` selon sa propre documentation |

Ces fonctions sont établies par leurs documentations. Elles n’ont pas fait l’objet d’un test utilisateur comparatif local. Les réglages d’un MUSCLE intégré à une autre application ne sont pas supposés compatibles avec la CLI de MUSCLE 5.

Le choix d’AlignView est une interface terminal avec champs validés et aperçu issu des mêmes constructeurs que l’exécution. Le périmètre est plus restreint que des commandes externes entièrement libres. Il couvre les contrôles exposés dans le [guide des paramètres](PARAMETERS.md), sans prétendre exposer toutes les options natives.

## Validation locale

Deux passages ont été réalisés sur les huit vrais moteurs : valeurs de départ et paramètres modifiés. Les résultats conservent versions, commandes, statuts et provenance dans [real-tools-results.json](real-tools-results.json) et [tuned-tools-results.json](tuned-tools-results.json).

| Moteur | Version locale du passage complet | Paramètres modifiés testés |
| --- | --- | --- |
| BLASTN | 2.17.0+ | Tâche, graine, reward/penalty, gaps, DUST, soft masking, brin, HSP, candidats |
| DIAMOND | 2.2.8 | Matrice, gaps, masquage, composition, HSP et candidats |
| MMseqs2 | 18.8cc5c | Sensibilité 6.5, masquage et correction désactivés, 37 candidats et traceback |
| MAFFT | 7.526 | L-INS-i, ouverture de gap 2, offset 0.1, 3 raffinements |
| MUSCLE | 5.3 | Super5, permutation abc, perturbation 2, 3 raffinements |
| Clustal Omega | 1.2.4 | Une itération et matrices de distances complètes |
| EMBOSS Needle | 6.6.0.0 | EDNAFULL et pénalités terminales distinctes |
| FastTree | 2.2.0 | LG, Gamma20, 50 rééchantillonnages locaux et graine 42 |

Les fixtures vérifient des hits sur séquences identiques, des alignements valides et des projets reprenables. Les tests de non-régression couvrent aussi la fidélité aperçu/exécution, l’isolation des réglages de chaque moteur, l’effet du seuil HSP sur les séquences retenues, les profils invalides, les priorités CLI, la provenance après reprise, les scripts Slurm et l’accès aux paramètres dans un vrai pseudo-terminal. Le contrôle Chrome vérifie notamment que le résumé de recherche montre le seuil réellement exécuté plutôt qu’une nouvelle valeur préparée.

Pour reproduire les passages avec les moteurs installés :

```bash
python3 scripts/validate_real_tools.py \
  --require mafft blastn diamond mmseqs muscle clustalo needle FastTree
python3 scripts/validate_real_tools.py --tuned \
  --require mafft blastn diamond mmseqs muscle clustalo needle FastTree
python3 -m unittest discover -s tests -v
```

Les moteurs du passage complet ont été installés dans un préfixe temporaire, sans remplacer les logiciels du poste. Ces tests ne couvrent pas toutes les combinaisons de matrices/coûts, toutes les combinaisons de codes génétiques alternatifs, toutes les versions de moteurs, les données difficiles ou le cluster cible. La révision initiale utilisait le code standard ; la révision suivante ajoute les tables NCBI choisies explicitement. Les journaux d’erreur des moteurs restent l’autorité sur les incompatibilités natives.

## Révision complémentaire : utilisation et moteurs

La révision complémentaire ajoute les guides français/anglais embarqués, le sélecteur de requête avant comparaison guidée, les actions VSEARCH natives, les 27 tables NCBI hors ligne et un export HTML explicite. Les rapports historiques de moteurs ci-dessus décrivent leurs passages d’origine ; ils ne prétendent pas tester toutes les tables ni tous les nouveaux réglages. Voir [VSEARCH](VSEARCH.md), [codes de traduction](GENETIC_CODES.md) et [règles de fichiers](FILES.md).
