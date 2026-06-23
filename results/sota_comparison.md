# External state of the art

> Our scores come from **lightweight classifiers (LogReg/RF/MLP) on FROZEN embeddings**,
> without fine-tuning. The SOTA tools below use dedicated, often end-to-end fine-tuned
> models trained on much larger datasets. The gap is expected: the goal here is to
> *probe the information contained in the embeddings*, not to beat specialized tools.

## Transmembrane (TM vs non-TM / topology)

| Tool | Task | Reported performance | Reference |
|---|---|---|---|
| **DeepTMHMM** | TM topology (per-residue) | ~0.98 acc TM/globular detection | Hallgren et al. 2022 (bioRxiv) |
| **TMbed** | TM detection + segments | ~0.98 recall on TM proteins | Bernhofer & Rost 2022 (BMC Bioinformatics) |
| Best here (ProtT5) | binary TM (per-protein) | **MCC 0.37** | this work |

Large gap: simplified binary task + strong class imbalance (8% positives) + linear probe on frozen embeddings.

## Subcellular localization

| Tool | Task | Performance | Reference |
|---|---|---|---|
| **DeepLoc 2.0** | 10 compartments | ~0.75 acc / Q10 | Thumuluri et al. 2022 (NAR) |
| **DeepLoc 1.0** | 10 compartments | ~0.78 acc | Almagro Armenteros et al. 2017 |
| Best here (Ankh-L) | 5 grouped classes | **macro-F1 0.57** | this work |

Our 5 grouped classes differ from DeepLoc's 10; qualitative comparison only.

## Fold / structure classification

| Tool | Task | Performance | Reference |
|---|---|---|---|
| DeepFRI | SCOP fold (hundreds of classes) | varies with granularity | Gligorijević et al. 2021 (Nat. Commun.) |
| ESM2/ProtT5 probes | class / superfamily | high F1 on 4 top-level classes | various PLM benchmarks |
| Best here (ProtT5) | SCOP 4 classes (a/b/c/d) | **macro-F1 0.76** | this work |

We use the 4 top-level classes (not fine folds), an easier task; 0.76 is consistent with a linear probe.

## Species (4 ATLAS organisms)

No standard SOTA (non-conventional task, restricted to the 4 ATLAS species with >=30 proteins).
High scores (macro-F1 ~0.82, ProtT5/Ankh) suggest the source organism is strongly encoded,
likely an amino-acid composition bias.

## Aggregation propensity

Aggrescan3D (native, Python 2.7) is not installable on Apple Silicon. We use a reproducible
Python-3 structural proxy (Kyte-Doolittle hydrophobicity + Chou-Fasman beta propensity,
modulated by solvent accessibility). Not directly comparable to native A3D scores.

## Flexibility / disorder / accessibility (per-protein)

No directly comparable per-protein SOTA (these are usually predicted per-residue). Our
per-residue results (`dt_results_*.csv`) are in line with Soufir et al.
