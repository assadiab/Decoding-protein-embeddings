# Decoding Protein Embeddings — ESM-C & Ankh2

**M2 Biologie-Informatique · Université Paris Cité · INSERM UMR_S1134**

> Extending the work of Soufir et al. on protein language model embeddings to second-generation models ESM-C and Ankh2, with comparison against 1G baselines (ESM2, ProtT5, Ankh1G).

---

## 🎯 Objectif

Ce projet applique la méthodologie d'analyse des embeddings protéiques développée par Emma Soufir aux modèles de **deuxième génération** :
- **ESM-C** (300M et 600M) — EvolutionaryScale, déc. 2024
- **Ankh2-Large** — Elnaggar / Rost Lab

Les embeddings sont générés sur les protéines de la base **ATLAS** (1390 protéines avec simulations de dynamique moléculaire), puis utilisés pour prédire des propriétés structurales et dynamiques par résidu (RMSF, Neq, structure secondaire, accessibilité).

---

## 📊 Résultats (mis à jour au fil de l'avancement)

### État d'avancement

| Tâche | Statut | Semaine |
|---|---|---|
| Setup environnement + pipeline | ✅ Complet | S1 |
| API ATLAS — scripts de téléchargement | ✅ Complet | S1 |
| Validation 3 modèles 2G (subset 20 protéines) | ✅ Complet | S1 |
| Génération embeddings 1390 protéines (cluster SFBI) | 🔄 En cours | S2 |
| Premières analyses DT/NN | ⏳ À faire | S2 |
| Benchmark comparatif 2G vs 1G | ⏳ À faire | S2–S3 |
| SHAP + interprétabilité | ⏳ À faire | S3 |
| Rapport final | ⏳ À faire | S3–S4 |

### Performances — à compléter (S2/S3)

| Tâche | ESM-C 300M | ESM-C 600M | Ankh2 | ESM2 650M | ProtT5 | Ankh1G |
|---|---|---|---|---|---|---|
| RMSF (Spearman) | — | — | — | ref. Soufir | ref. Soufir | ref. Soufir |
| Neq (Spearman) | — | — | — | ref. Soufir | ref. Soufir | ref. Soufir |
| Sec. structure (F1) | — | — | — | ref. Soufir | ref. Soufir | ref. Soufir |
| Solvent acc. (F1) | — | — | — | ref. Soufir | ref. Soufir | ref. Soufir |

---

## 🗂️ Structure du repo

```
├── scripts/                    # Pipeline de données
│   ├── fetch_atlas_fasta.py    # Récupération séquences FASTA depuis RCSB PDB
│   ├── download_atlas_data.py  # Téléchargement labels ATLAS (RMSF/Bfactor/Neq)
│   ├── prepare_atlas_labels.py # Conversion labels vers format pipeline
│   ├── cluster_setup.sh        # Setup PLM-API sur cluster SFBI
│   └── slurm/                  # Jobs SLURM pour génération embeddings
│       ├── run_esmc_300m.sh
│       ├── run_esmc_600m.sh
│       ├── run_ankh2.sh
│       └── run_all.sh
├── analysis/                   # Scripts et notebooks d'analyse
│   ├── full_prot_emb_2g.py     # Adaptation du script Soufir pour modèles 2G
│   └── notebooks/              # Jupyter notebooks (S2/S3)
├── results/                    # Figures et tableaux (S3)
│   └── figures/
├── pixi.toml                   # Configuration environnement reproductible
└── README.md
```

---

## 🚀 Reproduire les résultats

### Prérequis

- Python 3.11
- [pixi](https://pixi.sh) pour l'environnement d'analyse
- [PLM-API](https://gitlab.dsimb.insern.fr/cretin/plm-api) pour la génération d'embeddings
- Accès à la base [ATLAS](https://www.dsimb.insern.fr/ATLAS)

### 1. Installer l'environnement

```bash
git clone https://github.com/assadiabo/plm-decoding-esmc-ankh2
cd plm-decoding-esmc-ankh2

# Environnement d'analyse
pixi install

# PLM-API (cloner séparément)
git clone https://gitlab.dsimb.insern.fr/cretin/plm-api
cd plm-api && pip install . && cd ..
```

### 2. Télécharger les données ATLAS

```bash
# Séquences FASTA depuis RCSB PDB
pixi run python scripts/fetch_atlas_fasta.py \
    --id-files path/to/id_train.txt path/to/id_test.txt \
    --out Datasets/ATLAS/atlas_sequences.fasta

# Labels dynamiques depuis l'API ATLAS
pixi run python scripts/download_atlas_data.py --workers 4
pixi run python scripts/prepare_atlas_labels.py
```

### 3. Générer les embeddings

```bash
# Local (test subset)
plm-api/.venv/bin/plm-embeddings \
    -i Datasets/ATLAS/atlas_sequences.fasta \
    -o Datasets/embeddings/ -m esmc_300M

# Cluster SFBI (production complète)
bash scripts/cluster_setup.sh
sbatch scripts/slurm/run_all.sh
```

### 4. Construire les datasets et lancer les analyses

```bash
pixi shell
python analysis/full_prot_emb_2g.py id_train.txt \
    Datasets/embeddings/esmc_300M.safetensors \
    -emb esmc_300m --y rmsf \
    -dssp_dir Datasets/labels/dssp/ \
    -md_dir Datasets/labels/md/
```

---

## 🧬 Modèles

| Modèle | Génération | Dimension | Statut |
|---|---|---|---|
| ESM-C 300M | 2G — cible | 960 | ✅ validé |
| ESM-C 600M | 2G — cible | 1152 | ✅ validé |
| Ankh2-Large | 2G — cible | 1536 | ✅ validé |
| ESM2 650M | 1G — baseline | 1280 | ref. Soufir |
| ESM2 3B | 1G — baseline | 2560 | ref. Soufir |
| Ankh-Large | 1G — baseline | 1536 | ref. Soufir |
| ProtT5-XL | 1G — baseline | 1024 | ref. Soufir |
| Ankh3-Large | 3G — bonus | TBD | ⏳ S3 |

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

*Dernière mise à jour : 16 avril 2026 — fin S1*