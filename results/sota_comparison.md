# État de l'art externe — pour situer nos scores

> Nos chiffres viennent de **classifieurs légers (LogReg/RF/MLP) sur embeddings GELÉS**,
> sans fine-tuning. Les SOTA ci-dessous utilisent des modèles dédiés, souvent fine-tunés
> end-to-end et entraînés sur des jeux bien plus grands. **L'écart est attendu** : notre
> objectif est de *sonder l'information contenue dans les embeddings*, pas de battre les
> outils spécialisés. À expliciter dans la discussion.

## Transmembrane (TM vs non-TM / topologie)

| Outil | Tâche | Performance rapportée | Référence |
|---|---|---|---|
| **DeepTMHMM** | topologie TM (per-residue) | ~0.98 acc détection TM/globulaire | Hallgren et al. 2022 (bioRxiv) |
| **TMbed** | détection + segments TM | ~0.98 recall protéines TM | Bernhofer & Rost 2022 (BMC Bioinformatics) |
| Notre meilleur (ProtT5) | TM binaire (per-protéine) | **MCC 0.37** | ce travail |

→ Écart majeur : tâche binaire simplifiée + classes très déséquilibrées (8 % positifs) + probe linéaire sur embedding gelé.

## Localisation subcellulaire

| Outil | Tâche | Performance | Référence |
|---|---|---|---|
| **DeepLoc 2.0** | 10 compartiments | ~0.75 acc / Q10 ; MCC variable par classe | Thumuluri et al. 2022 (NAR) |
| **DeepLoc 1.0** | 10 compartiments | ~0.78 acc (jeu d'origine) | Almagro Armenteros et al. 2017 |
| Notre meilleur (Ankh-L) | 5 classes regroupées | **macro-F1 0.57** | ce travail |

→ Nos 5 classes regroupées ≠ les 10 de DeepLoc ; comparaison qualitative seulement.

## Classification de fold / structure

| Outil | Tâche | Performance | Référence |
|---|---|---|---|
| SOTA fold recognition (ex. DeepFRI, profilage) | SCOP fold (centaines de classes) | très variable selon granularité | Gligorijević et al. 2021 (DeepFRI, Nat. Commun.) |
| Probe ESM2/ProtT5 (littérature) | classe/superfamille | F1 élevé sur 4 classes top-niveau | divers benchmarks PLM |
| Notre meilleur (ProtT5) | SCOP 4 classes (a/b/c/d) | **macro-F1 0.76** | ce travail |

→ Nous restons sur les **4 classes top-niveau** (pas les folds fins), tâche plus facile : nos 0.76 sont cohérents avec un probe linéaire sur embedding.

## Espèce (4 organismes ATLAS)

Pas de SOTA standard (tâche non conventionnelle, restreinte aux 4 espèces ≥30 prot d'ATLAS).
Nos scores élevés (macro-F1 ~0.82, ProtT5/Ankh) suggèrent que l'organisme d'origine est
fortement encodé — biais de composition en acides aminés probablement. À discuter.

## Flexibilité / désordre / accessibilité (per-protéine)

Pas d'outil SOTA per-protéine directement comparable (ces propriétés sont usuellement
prédites per-résidu). Nos résultats per-résidu (`dt_results_*.csv`) se rapprochent des
ordres de grandeur de Soufir et al.

---

*Placeholders à confirmer avant rédaction finale : valeurs exactes DeepLoc 2.0 (acc/MCC par
compartiment), TMbed/DeepTMHMM (métrique précise selon le benchmark cité).*
