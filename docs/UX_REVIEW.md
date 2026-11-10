# Revue critique de l’interface

Révision du 10 novembre 2026, après retour d’utilisation sur la première préparation 0.1.0.

## Diagnostic

La première refonte avait ajouté des fonctions et raccourcis sans réduire suffisamment les décisions. Les tests validaient calculs et touches ; ils ne démontraient pas la clarté du parcours complet.

| Problème observé | Conséquence | Modification |
| --- | --- | --- |
| Jusqu’à 13 entrées à l’accueil, destinations dupliquées | Difficile de savoir par où commencer | Cinq rubriques stables, un contexte et une étape conseillée |
| Nombreuses méthodes de menus aux vocabulaires variables | Profondeur et repères incohérents | Contrôleur unique `av_navigation.py` |
| Moteurs et profils proposés trop tôt | Il faut connaître les logiciels avant de poser une question | Démarrage par intention, proposition compatible et aperçu avant calcul |
| Actions proposées sans données nécessaires | Commandes sans effet, erreurs tardives | Choix indisponibles avec leur prérequis expliqué |
| Libellé et description dans une ligne étroite | Texte tronqué difficile à comparer | Une colonne, description sous le libellé |
| Retours locaux uniquement | Perte de repères dans les sous-menus | Fil de navigation, Échap/← et F2 accueil |
| Touches vidées pendant les transitions | Une touche rapide semble sans effet | Mode terminal changé sans vider les touches en attente |
| Aide dans un écran alternatif imbriqué | Écran ou souris mal restaurés | Gestion des écrans imbriqués vérifiée dans un PTY |
| Rafraîchissement permanent au repos | Sortie SSH inutile, scintillement possible | Rendu après modification ou redimensionnement |
| Profils changeant aussi le maximum de hits | Nombre de résultats modifié de façon peu visible | Limite 100 indépendante du profil ; `--hits` prioritaire |
| Séquence active implicite dans un FASTA unique | Analyses ciblant la première séquence sans choix clair | Sélection explicite, FASTA complet conservé pour le MSA |
| Historique temporaire assimilé à une sauvegarde | Risque de quitter en pensant le travail conservé | Reprise d’historique toujours marquée à enregistrer |
| Réglages futurs confondus avec le calcul existant | Interprétation ambiguë après modification | Paramètres calculés conservés séparément, état de relance visible |
| Rapport avec cartes vides et toutes les figures ouvertes | Même surcharge dans le navigateur | Guide de lecture, vues selon les données, figures dépliables |

## Parcours retenus

L’accueil suit les tâches : **Explorer, Analyser, Données, Enregistrer et partager, Réglages et aide**. Le moteur reste accessible sans imposer sa connaissance au premier écran.

Un alignement existant est chargé sans moteur. Un FASTA brut reste chargé avant calcul ; l’aperçu précise les séquences et le moteur. Une comparaison commence par une recherche, dont les correspondances et couvertures sont regroupées. Créer ensuite un MSA reste une action explicite.

Projet, rapport et exports sont décrits par leur usage : reprendre, présenter, réutiliser. Un enregistrement est visible et ne provoque pas de deuxième copie à la fermeture.

## Vérification

Les tests parcourent de vrais pseudo-terminaux : démo → visionneuse → F2 ; réglages imbriqués → F2 ; prérequis manquant ; alignement guidé ; sélecteur ; changement de requête et historique ; sauvegarde et fermeture.

Avec 120 références, les tests vérifient que 100 sont retenues et exportées. Profils, limites explicites et reprise de paramètres existants sont testés séparément. Les cinq rubriques et leurs libellés tiennent dans un terminal 80×24.

Le rapport est contrôlé dans Chrome, en complément des tests de formats. Ces vérifications portent sur des comportements observables ; elles ne remplacent pas les retours d’usage sur les données et terminaux de travail.

La [revue des calculs](CALCULATION_REVIEW.md) complète cette révision : les contrôles scientifiques sont accessibles avant lancement, propres au moteur sélectionné, accompagnés de leur effet et réutilisables dans des profils. L’aperçu des commandes et la provenance distinguent configuration future et calcul effectué.

## Limites

Les grands FASTA restent en mémoire, les rapports gardent des modules historiques et toutes les versions de moteurs/plateformes ne sont pas validées. Les huit moteurs ont désormais un essai réel local, avec des paramètres de départ et modifiés. L’audit scientifique conserve les limites de portée. Aucune satisfaction utilisateur ou supériorité générale n’est déduite du nombre de tests.
