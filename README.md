# Decoding Protein Embeddings — ESM-C & Ankh2

**M2 Bioinformatics · Université Paris Cité · INSERM UMR_S1134**

Do second-generation protein language models (ESM-C, Ankh2) better encode global protein properties than established baselines? We benchmark 6 PLMs on 5 global prediction tasks across 1,390 proteins from the ATLAS database — and find a surprising answer.

---

## Key Finding

> **Second-generation models (ESM-C, Ankh2) do not systematically outperform first-generation baselines on global protein properties.** ProtT5 and Ankh-Large (1G) remain competitive across fold, transmembrane, localization, disorder and accessibility prediction. However, Ankh2-Large shows a clear advantage on local per-residue tasks (secondary structure, flexibility).

This suggests that architectural improvements in 2G models primarily benefit residue-level tasks rather than protein-level aggregated representations.

---

## Results — Global Properties (best classifier per task, test set)

Metric: **macro-F1** for multi-class tasks (fold, localization), **MCC** for binary tasks.

| Task | ESM-C 300M | ESM-C 600M | Ankh2-Large | ESM2-650M | Ankh-L 1G | ProtT5 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| Fold classification (F1 ↑) | 0.677 | 0.674 | 0.747 | 0.756 | 0.749 | **0.763** |
| Subcellular localization (F1 ↑) | 0.484 | 0.521 | 0.569 | 0.481 | 0.572 | **0.558** |
| Transmembrane (MCC ↑) | 0.252 | 0.241 | 0.340 | 0.323 | 0.257 | **0.370** |
| Global disorder (MCC ↑) | 0.265 | 0.261 | **0.371** | 0.298 | 0.364 | 0.349 |
| Mean accessibility (MCC ↑) | 0.742 | 0.664 | 0.767 | 0.741 | 0.736 | **0.764** |

## Results — Local Properties / Per-residue (Decision Tree, test set, F1 / MCC)

| Task | ESM-C 300M | ESM-C 600M | Ankh2-Large | ESM2-650M | Ankh-L 1G | ProtT5 |
|---|:-:|:-:|:-:|:-:|:-:|:-:|
| RMSF / flexibility | 0.63 / 0.22 | 0.64 / 0.24 | **0.65 / 0.29** | 0.645 | 0.675 | 0.660 |
| Neq / disorder | 0.75 / 0.35 | 0.76 / 0.39 | **0.79 / 0.46** | 0.735 | 0.785 | 0.755 |
| B-factor | 0.69 / 0.17 | 0.68 / 0.17 | **0.70 / 0.17** | 0.700 | 0.710 | 0.705 |
| Solvent accessibility | 0.82 / 0.53 | 0.82 / 0.53 | **0.84 / 0.58** | 0.825 | 0.845 | 0.825 |
| Secondary structure (3) | 0.60 / 0.42 | 0.64 / 0.49 | **0.70 / 0.57** | 0.560 | 0.705 | 0.615 |
| Secondary structure (8) | 0.20 / 0.30 | 0.22 / 0.36 | **0.25 / 0.42** | 0.190 | 0.245 | 0.210 |

On per-residue tasks, **Ankh2-Large dominates** across all 6 variables.

---

## Methods

**Dataset** — 1,390 non-redundant proteins from [ATLAS](https://www.dsimb.inserm.fr/ATLAS) (atomistic MD simulations). Train/test split: 973 / 417 proteins.

**Embeddings** — Generated via [PLM-API](https://gitlab.dsimb.inserm.fr/cretin/plm-api) on Apple M4 (MPS acceleration). All 6 models run on the full 1,390-protein set.

**Global pipeline** — Per-residue embeddings `[L × D]` aggregated by **mean pooling** → fixed-size protein representation `[D]`. Labels (fold/TM/localization) retrieved from ATLAS API + UniProt REST API. Classifiers: Logistic Regression, Random Forest, MLP with 15-fold stratified cross-validation.

**Local pipeline** — Per-residue embeddings matched to DSSP-derived labels (SS3, SS8, accessibility) and MD simulation labels (RMSF, Neq, B-factor). Classifier: Decision Tree following Soufir et al. (2024) methodology.

**Interpretability** — SHAP values (top-20 dimensions per task) + drop curves (progressive dimension masking 0→100%).

---

## Models

| Model | Generation | Dim | Architecture |
|---|---|---|---|
| ESM-C 300M | 2G | 960 | Transformer encoder, UniRef + metagenomics |
| ESM-C 600M | 2G | 1152 | Transformer encoder, UniRef + metagenomics |
| Ankh2-Large | 2G | 1536 | T5 encoder + SiLU, ankh2-ext2 variant |
| ESM2-650M | 1G baseline | 1280 | BERT encoder, UniRef50 |
| Ankh-Large | 1G baseline | 1536 | T5 encoder + ReLU |
| ProtT5-XL | 1G baseline | 1024 | T5 encoder, UniRef50 |

---

## Repository Structure

```
├── analysis/
│   ├── global/
│   │   ├── build_global_dataset.py      # Mean pooling → per-protein CSV
│   │   ├── train_global_classifiers.py  # LogReg + RF + MLP, 15-fold CV
│   │   └── shap_drop_global.py          # SHAP + drop curves
│   ├── full_prot_emb_2g.py              # Per-residue dataset builder
│   └── notebooks/
│       └── results_comparison.ipynb     # Full comparison notebook
├── scripts/
│   ├── fetch_global_labels.py           # ATLAS API + UniProt → global labels
│   ├── compute_dssp.py                  # DSSP v4 on ATLAS PDB files
│   ├── run_embeddings_m4.sh             # Embedding generation (Apple M4 / MPS)
│   ├── run_embeddings_1g.sh             # 1G baselines
│   └── run_analyses.sh                  # Per-residue DT analyses
├── results/
│   ├── global_results.csv               # Global tasks — all models × classifiers
│   ├── dt_results_2g.csv                # Per-residue — 2G models
│   ├── dt_results_1g.csv                # Per-residue — 1G baselines
│   └── figures/                         # SHAP, drop curves, heatmaps, PCA
├── pixi.toml                            # Reproducible environment (Python 3.11)
└── README.md
```

---

## Reproduce

### Requirements

- Python 3.11 + [pixi](https://pixi.sh)
- [PLM-API](https://gitlab.dsimb.inserm.fr/cretin/plm-api) for embedding generation
- Access to [ATLAS](https://www.dsimb.inserm.fr/ATLAS) database

### Setup

```bash
git clone https://github.com/assadiab/protein-embeddings-atlas
cd protein-embeddings-atlas
pixi install
```

### 1. Fetch data

```bash
# ATLAS protein sequences
pixi run python scripts/download_atlas_data.py --workers 4

# Global labels (fold, TM, localization via ATLAS API + UniProt)
pixi run python scripts/fetch_global_labels.py
# → Datasets/ATLAS/global_labels/atlas_global_labels.tsv

# DSSP (secondary structure + accessibility per residue)
pixi run python scripts/compute_dssp.py \
    --pdb_dir Datasets/ATLAS/data/ \
    --output_dir Datasets/ATLAS/data/
```

### 2. Generate embeddings

```bash
# Requires PLM-API (see plm-api/ setup)
source plm-api/.venv/bin/activate
bash scripts/run_embeddings_m4.sh   # 2G models (~8 min on M4)
bash scripts/run_embeddings_1g.sh   # 1G baselines
```

### 3. Run global analyses (WP2)

```bash
pixi shell
python analysis/global/build_global_dataset.py      # mean pooling → CSV
python analysis/global/train_global_classifiers.py  # LogReg + RF + MLP
python analysis/global/shap_drop_global.py          # SHAP + drop curves
# Results → results/global_results.csv + results/figures/global_*
```

### 4. Run per-residue analyses (extension)

```bash
bash scripts/build_all_datasets.sh
bash scripts/build_dssp_datasets.sh
bash scripts/run_analyses.sh
```

---

## References

- Vander Meersche Y. et al. (2023). **ATLAS: protein flexibility description from atomistic MD simulations.** *Nucleic Acids Research.* [doi:10.1093/nar/gkad1084](https://doi.org/10.1093/nar/gkad1084)
- Hayes T. et al. (2025). **Simulating 500 million years of evolution with a language model.** *Science*, 387(6736). (ESM-C)
- Elnaggar A. et al. (2023). **Ankh: Optimized Protein Language Model Unlocks General-Purpose Modelling.** *arXiv:2301.06568.*
- Lin Z. et al. (2023). **Evolutionary-scale prediction of atomic-level protein structure with a language model.** *Science*, 379, 1123–1130. (ESM2)
- Elnaggar A. et al. (2022). **ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning.** *IEEE TPAMI*, 44(10).

---

**Assa Diabira** · M2 Bioinformatics, Université Paris Cité  
Supervisor: Pr. Jean-Christophe Gelly · Co-supervisor: Gabriel Cretin · INSERM UMR_S1134
