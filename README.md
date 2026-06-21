# Protein Embedding Benchmark — 1G vs 2G PLMs on ATLAS

A reproducible benchmark that measures **how much structural and dynamic information is encoded in protein language model (PLM) embeddings**, and whether second-generation models (ESM-C, Ankh2) actually beat the first generation (ESM2, ProtT5, Ankh) at it.

Six PLMs are embedded over the **1390 ATLAS proteins** and probed on two fronts:

- **Global properties** (one label per protein, mean-pooled embeddings): fold class, subcellular localization, transmembrane status, global disorder, mean solvent accessibility.
- **Per-residue properties** (one label per residue): RMSF, Neq, B-factor, 3/8-state secondary structure, solvent accessibility.

The probing classifiers are deliberately simple (Logistic Regression / Random Forest / MLP, Decision Tree for per-residue) so that **scores reflect the embeddings, not the head**.

---

## TL;DR — what the benchmark says

- **Bigger generation ≠ better embeddings.** On *global* tasks the 1G models (ProtT5, Ankh-Large) match or beat every 2G model. The "500M-years-of-evolution" jump does not translate into more linearly-decodable global structure.
- **2G wins only on fine-grained local structure.** Ankh2-Large leads clearly on per-residue secondary structure (SS3/SS8), where 600M > 300M too — extra capacity helps resolution, not abstraction.
- **A linear probe is usually enough.** Logistic Regression is the best head on most global tasks: the information is already linearly separable in embedding space.
- **ESM-C is the weak link here.** Both 300M and 600M trail the others on global tasks despite being the newest models.

---

## Results

### Global properties (per-protein, best classifier, test set)

Metric: macro-F1 for fold/localization, MCC for tm/disorder/accessibility. **Bold = best per row.**

| Task | Metric | ESM-C 300M | ESM-C 600M | Ankh2-Large | ESM2-650M | Ankh-Large | ProtT5 |
|---|---|---|---|---|---|---|---|
| Fold (SCOP a/b/c/d) | F1 | 0.677 | 0.674 | 0.747 | 0.756 | 0.749 | **0.763** |
| Localization (5-class) | F1 | 0.484 | 0.521 | 0.569 | 0.481 | **0.572** | 0.558 |
| Transmembrane | MCC | 0.252 | 0.241 | 0.340 | 0.323 | 0.257 | **0.370** |
| Global disorder | MCC | 0.265 | 0.261 | **0.371** | 0.298 | 0.364 | 0.349 |
| Mean accessibility | MCC | 0.742 | 0.664 | **0.767** | 0.741 | 0.736 | 0.764 |

→ 1G models top 3 of 5 global tasks. 2G only leads on disorder and accessibility, and only via Ankh2.

### Per-residue properties (Decision Tree, test set, F1)

| Task | ESM-C 300M | ESM-C 600M | Ankh2-Large | ESM2-650M | Ankh-Large | ProtT5 |
|---|---|---|---|---|---|---|
| RMSF (flexibility) | 0.63 | 0.64 | 0.65 | 0.63 | **0.66** | 0.65 |
| Neq (disorder) | 0.75 | 0.76 | **0.79** | 0.73 | 0.78 | 0.75 |
| B-factor | 0.69 | 0.68 | **0.70** | 0.69 | **0.70** | 0.69 |
| Accessibility | 0.82 | 0.82 | **0.84** | 0.82 | **0.84** | 0.82 |
| SS3 | 0.60 | 0.64 | **0.70** | 0.55 | **0.70** | 0.61 |
| SS8 | 0.20 | 0.22 | **0.25** | 0.19 | 0.24 | 0.21 |

→ Ankh-family models (2G *and* 1G) lead the local tasks; the 2G advantage is real but concentrated on secondary structure.

Full per-classifier numbers: [`results/global_results.csv`](results/global_results.csv), [`results/dt_results_2g.csv`](results/dt_results_2g.csv), [`results/dt_results_1g.csv`](results/dt_results_1g.csv).

---

## Models benchmarked

| Model | Generation | Dim | Source |
|---|---|---|---|
| ESM-C 300M | 2G | 960 | EvolutionaryScale |
| ESM-C 600M | 2G | 1152 | EvolutionaryScale |
| Ankh2-Large | 2G | 1536 | ElnaggarLab |
| ESM2-650M | 1G | 1280 | Meta AI |
| Ankh-Large | 1G | 1536 | ElnaggarLab |
| ProtT5-XL | 1G | 1024 | Rost Lab |

All embeddings are generated with [PLM-API](https://gitlab.dsimb.inserm.fr/cretin/plm-api) over the **1390 ATLAS** proteins (MD-derived flexibility labels), with an Apple Silicon (MPS) patch for local runs.

---

## How it works

```
ATLAS sequences (1390)
        │  PLM-API  ->  Datasets/embeddings/full/*.safetensors  [seq_len, dim]
        ▼
  ┌─────────────────────────────┬──────────────────────────────┐
  │ GLOBAL pipeline             │ PER-RESIDUE pipeline          │
  │ mean pooling -> [dim]       │ keep all positions            │
  │ + ATLAS/UniProt labels      │ + DSSP / MD labels            │
  │ LogReg / RF / MLP, 15-fold  │ Decision Tree                 │
  │ SHAP + drop curves          │ confusion matrices, PCA       │
  └─────────────────────────────┴──────────────────────────────┘
        ▼
  results/*.csv + results/figures/*
```

**Labels.** Global labels are assembled by `scripts/fetch_global_labels.py`: fold from the ATLAS metadata API (SCOP class), transmembrane + localization from UniProt (batched), global disorder = mean RMSF, mean accessibility from DSSP. Continuous targets are binarized at the **train-set median only** (no leakage). The train/test split follows ATLAS (973 / 417 proteins).

**Interpretability.** For the best classifier per (model, task), `analysis/global/shap_drop_global.py` produces SHAP top-20 dimension importances and dimension-dropout curves (mask 0→100% of dims, 3 repeats) — see `results/figures/global_shap_*` and `global_drop_*`.

---

## Repository layout

```
analysis/
  global/
    build_global_dataset.py        # mean pooling + label join -> per-protein CSVs
    train_global_classifiers.py    # LogReg/RF/MLP, 15-fold CV, test eval
    shap_drop_global.py            # SHAP + dimension-dropout curves
  full_prot_emb_2g.py              # per-residue dataset builder (2G adaptation)
  notebooks/results_comparison.ipynb
scripts/
  fetch_global_labels.py           # ATLAS metadata + UniProt -> global labels
  fetch_atlas_fasta.py / download_atlas_data.py / prepare_atlas_labels.py
  compute_dssp.py                  # DSSP v4 over ATLAS PDBs
  run_embeddings_m4.sh / run_embeddings_1g.sh
  build_*_datasets.sh / run_analyses*.sh / run_pca.sh
results/
  global_results.csv               # global benchmark (90 rows)
  dt_results_2g.csv / dt_results_1g.csv
  figures/                         # global_* (SHAP/drop/heatmaps) + per-residue
pixi.toml                          # reproducible environment
```

---

## Reproduce

```bash
git clone https://github.com/assadiab/projet_long && cd projet_long
pixi install

# 1. data
pixi run python scripts/download_atlas_data.py --workers 4
pixi run python scripts/prepare_atlas_labels.py
pixi run python scripts/compute_dssp.py

# 2. embeddings (separate PLM-API venv; MPS-ready on Apple Silicon)
source plm-api/.venv/bin/activate
bash scripts/run_embeddings_m4.sh      # 2G
bash scripts/run_embeddings_1g.sh      # 1G baselines

# 3. global benchmark
pixi run python scripts/fetch_global_labels.py
pixi run python analysis/global/build_global_dataset.py
pixi run python analysis/global/train_global_classifiers.py
pixi run python analysis/global/shap_drop_global.py

# 4. per-residue benchmark
pixi shell
bash scripts/build_all_datasets.sh && bash scripts/build_dssp_datasets.sh
bash scripts/run_analyses.sh
```

---

## Notes & limitations

- **Aggregation propensity** is not benchmarked (would require a local Aggrescan3D install).
- Fold and localization labels cover ~51% of ATLAS (no SCOP / no UniProt annotation for the rest), so those two tasks run on ~700 proteins; transmembrane, disorder and accessibility cover the full set.
- RMSF is read in Å (ATLAS update of 2025-11-21).

---

## References

- Vander Meersche et al. (2024). *ATLAS: protein flexibility description from atomistic MD simulations.* Nucleic Acids Research. [doi:10.1093/nar/gkad1084](https://doi.org/10.1093/nar/gkad1084)
- Hayes et al. (2024). *Simulating 500 million years of evolution with a language model.* (ESM-C / ESM3)
- Elnaggar et al. (2023). *Ankh: Optimized Protein Language Model.* arXiv.
- Lin et al. (2023). *Evolutionary-scale prediction of atomic-level protein structure with a language model.* Science. (ESM2)

---

*Built on the Soufir et al. embedding-decoding methodology. Author: Assa Diabira — M2 Bio-Informatics, Université Paris Cité — INSERM UMR_S1134.*
