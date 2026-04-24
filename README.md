# Decoding Protein Embeddings — ESM-C & Ankh2

**M2 Biologie-Informatique · Université Paris Cité · INSERM UMR_S1134**

> Extension du travail de Soufir et al. sur les embeddings de protéines aux modèles de deuxième génération ESM-C et Ankh2, avec comparaison aux baselines 1G (ESM2, ProtT5, Ankh1G).

---

## 🎯 Objectif

Ce projet applique la méthodologie d'analyse des embeddings protéiques développée par Emma Soufir aux modèles de **deuxième génération** :
- **ESM-C** (300M et 600M) — EvolutionaryScale, déc. 2024
- **Ankh2-Large** — Elnaggar / Rost Lab

Les embeddings sont générés sur les 1390 protéines de la base **ATLAS** (simulations de dynamique moléculaire), puis utilisés pour prédire des propriétés structurales et dynamiques par résidu (RMSF, Neq, Bfactor, structure secondaire, accessibilité).

---

## 📊 Résultats — Decision Tree (test set, F1 / MCC)

| Tâche | ESM-C 300M | ESM-C 600M | Ankh2-Large |
|---|---|---|---|
| RMSF (flexibilité) | 0.63 / 0.22 | 0.64 / 0.24 | **0.65 / 0.29** |
| Neq (désordre) | 0.75 / 0.35 | 0.76 / 0.39 | **0.79 / 0.46** |
| Bfactor | 0.69 / 0.17 | 0.68 / 0.17 | **0.70 / 0.17** |
| Accessibilité solvant | 0.82 / 0.53 | 0.82 / 0.53 | **0.84 / 0.58** |
| Structure secondaire 3 classes | 0.60 / 0.42 | 0.64 / 0.49 | **0.70 / 0.57** |
| Structure secondaire 8 classes | 0.20 / 0.30 | 0.22 / 0.36 | **0.25 / 0.42** |

**Ankh2-Large domine sur toutes les tâches. ESM-C 600M > 300M surtout sur SS3/SS8.**

---

## 🗂️ Structure du repo

```
├── analysis/
│   ├── full_prot_emb_2g.py     # Adaptation pipeline Soufir pour modèles 2G
│   └── notebooks/              # Notebooks Jupyter (analyses interactives)
├── scripts/
│   ├── fetch_atlas_fasta.py    # Récupération FASTA depuis RCSB PDB
│   ├── download_atlas_data.py  # Labels ATLAS (RMSF/Bfactor/Neq)
│   ├── prepare_atlas_labels.py # Conversion labels
│   ├── compute_dssp.py         # Calcul DSSP v4 sur PDB ATLAS
│   ├── run_embeddings_m4.sh    # Génération embeddings sur Apple M4 (MPS)
│   ├── build_all_datasets.sh   # CSV dynamics (rmsf, neq, bfact)
│   ├── build_dssp_datasets.sh  # CSV structure (acc, sec3, sec8)
│   └── run_analyses.sh         # DT sur toutes les variables
├── results/
│   ├── figures/                # Confusion matrices, arbres de décision
│   └── dt_results_2g.csv       # Tableau métriques complet
├── pixi.toml                   # Environnement reproductible
└── README.md
```

---

## 🚀 Reproduire les résultats

### Prérequis

- Python 3.11 + [pixi](https://pixi.sh)
- [PLM-API](https://gitlab.dsimb.inserm.fr/cretin/plm-api) pour la génération d'embeddings
- Accès à la base [ATLAS](https://www.dsimb.inserm.fr/ATLAS)

### 1. Installer l'environnement

```bash
git clone https://github.com/assadiab/projet_long
cd projet_long
pixi install
```

### 2. Télécharger les données ATLAS

```bash
pixi run python scripts/download_atlas_data.py --workers 4
pixi run python scripts/prepare_atlas_labels.py
```

### 3. Générer les embeddings

```bash
# Sur Apple Silicon M4 (MPS)
source plm-api/.venv/bin/activate
bash scripts/run_embeddings_m4.sh

# Les 3 modèles tournent en ~8 min total sur M4
```

### 4. Calculer les DSSP

```bash
# dssp est déjà dans pixi.toml
pixi run python scripts/compute_dssp.py
# Génère Datasets/ATLAS/data/{protein_id}_dssp.tsv pour les 1390 protéines
```

### 5. Construire les datasets et lancer les analyses

```bash
pixi shell
bash scripts/build_all_datasets.sh      # dynamics : rmsf, neq, bfact
bash scripts/build_dssp_datasets.sh     # structure : acc, sec3, sec8
bash scripts/run_analyses.sh            # DT sur toutes les variables
```

---

## 🧬 Modèles

| Modèle | Génération | Dimension | Statut |
|---|---|---|---|
| ESM-C 300M | 2G — cible | 960 | ✅ généré (1390 protéines) |
| ESM-C 600M | 2G — cible | 1152 | ✅ généré (1390 protéines) |
| Ankh2-Large | 2G — cible | 1536 | ✅ généré (1390 protéines) |
| ESM2 650M | 1G — baseline | 1280 | réf. Soufir et al. |
| Ankh-Large (1G) | 1G — baseline | 1536 | réf. Soufir et al. |
| ProtT5-XL | 1G — baseline | 1024 | réf. Soufir et al. |

---

## 📚 Références

- Vander Meersche et al. (2023). **ATLAS: protein flexibility description from atomistic MD simulations.** *Nucleic Acids Research.* [doi:10.1093/nar/gkad1084](https://doi.org/10.1093/nar/gkad1084)
- Hayes et al. (2024). **Simulating 500 million years of evolution with a language model.** *Science.* (ESM-C / ESM3)
- Elnaggar et al. (2023). **Ankh: Optimized Protein Language Model Unlocks General-Purpose Modelling.** *arXiv.*
- Lin et al. (2023). **Evolutionary-scale prediction of atomic-level protein structure with a language model.** *Science.* (ESM2)

---

## 👩‍💻 Auteure

**Assa Diabira** — M2 Biologie-Informatique, Université Paris Cité  
Encadrant : Pr. Jean-Christophe Gelly · Co-encadrant : Gabriel Cretin  
INSERM UMR_S1134 — Hôpital Necker, Paris

---

*Dernière mise à jour : 23 avril 2026 — fin S3*