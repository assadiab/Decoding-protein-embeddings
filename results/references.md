# References

## Dataset

- **ATLAS** — Vander Meersche, Y., Cretin, G., Gheeraert, A., Gelly, J.-C., Galochkina, T. (2024).
  *ATLAS: protein flexibility description from atomistic molecular dynamics simulations.*
  **Nucleic Acids Research**, 52(D1), D384–D392. https://doi.org/10.1093/nar/gkad1084

## Reference methodology

- **Soufir, E. et al.** (2024). *Decoding protein language model embeddings* —
  per-residue probing of structural properties from PLM embeddings.

## Protein language models

### 2G (targets)
- **ESM-C / ESM Cambrian** — EvolutionaryScale Team (2024). *ESM Cambrian.*
  https://www.evolutionaryscale.ai/blog/esm-cambrian (300M / 600M models).
  ESM3 lineage: Hayes, T. et al. (2024). *Simulating 500 million years of evolution with a
  language model.* **Science**, eads0018. https://doi.org/10.1126/science.ads0018
- **Ankh2 / Ankh-Large** — Elnaggar, A. et al. (2023). *Ankh: Optimized Protein Language Model
  Unlocks General-Purpose Modelling.* **arXiv:2301.06568**. https://arxiv.org/abs/2301.06568

### 1G (baselines)
- **ESM2** — Lin, Z. et al. (2023). *Evolutionary-scale prediction of atomic-level protein
  structure with a language model.* **Science**, 379(6637), 1123–1130.
  https://doi.org/10.1126/science.ade2574
- **ProtT5 (ProtTrans)** — Elnaggar, A. et al. (2022). *ProtTrans: Toward Understanding the
  Language of Life Through Self-Supervised Learning.* **IEEE TPAMI**, 44(10), 7112–7127.
  https://doi.org/10.1109/TPAMI.2021.3095381

## Annotations & tools

- **DSSP** — Kabsch, W. & Sander, C. (1983). *Dictionary of protein secondary structure.*
  **Biopolymers**, 22(12), 2577–2637. (mkdssp v4: Joosten et al. 2011, NAR.)
- **SCOP** — Murzin, A.G. et al. (1995). *SCOP: a structural classification of proteins database.*
  **J. Mol. Biol.**, 247(4), 536–540.
- **UniProt** — The UniProt Consortium (2023). *UniProt: the Universal Protein Knowledgebase in
  2023.* **Nucleic Acids Research**, 51(D1), D523–D531. https://doi.org/10.1093/nar/gkac1052

## SOTA tools (comparison)

- **DeepLoc 2.0** — Thumuluri, V. et al. (2022). *DeepLoc 2.0: multi-label subcellular
  localization prediction.* **Nucleic Acids Research**, 50(W1), W228–W234.
  https://doi.org/10.1093/nar/gkac278
- **TMbed** — Bernhofer, M. & Rost, B. (2022). *TMbed: transmembrane proteins predicted through
  language model embeddings.* **BMC Bioinformatics**, 23, 326.
  https://doi.org/10.1186/s12859-022-04873-x
- **DeepTMHMM** — Hallgren, J. et al. (2022). *DeepTMHMM predicts alpha and beta transmembrane
  proteins using deep neural networks.* **bioRxiv**. https://doi.org/10.1101/2022.04.08.487609

## ML / interpretability

- **SHAP** — Lundberg, S.M. & Lee, S.-I. (2017). *A Unified Approach to Interpreting Model
  Predictions.* **NeurIPS 30**.
- **scikit-learn** — Pedregosa, F. et al. (2011). *Scikit-learn: Machine Learning in Python.*
  **JMLR**, 12, 2825–2830.
