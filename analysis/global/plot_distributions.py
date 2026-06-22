#!/usr/bin/env python3
"""FIX 4 — Distributions des classes cibles (train vs test) pour les tâches globales.

Justifie le choix de MCC sur les tâches déséquilibrées.
Lit atlas_global_labels.tsv + le split deciphering, applique les mêmes
préparations qu'à l'entraînement (filtre fold a/b/c/d, fusion mitochondrion,
binarisation disorder/acc à la médiane TRAIN, filtre species 4 classes).

Sortie :
  - results/figures/global_class_distributions.png  (multi-panneaux)
  - texte : % par classe (train/test) pour le rapport
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
LABELS_TSV = ROOT / "Datasets" / "ATLAS" / "global_labels" / "atlas_global_labels.tsv"
ID_TRAIN = ROOT / "deciphering" / "id_train.txt"
ID_TEST = ROOT / "deciphering" / "id_test.txt"
FIG = ROOT / "results" / "figures" / "global_class_distributions.png"

TASKS = ["fold_label", "localization_class", "species_label",
         "tm_label", "disorder_global", "acc_mean", "aggregation_score"]
FOLD_KEEP = {"a", "b", "c", "d"}


def load_ids(p):
    return [l.strip() for l in p.read_text().splitlines() if l.strip()]


def prep(df_tr, df_te, task):
    tr, te = df_tr.copy(), df_te.copy()
    if task == "fold_label":
        tr = tr[tr.fold_label.isin(FOLD_KEEP)]; te = te[te.fold_label.isin(FOLD_KEEP)]
        return tr.fold_label, te.fold_label
    if task == "localization_class":
        f = lambda v: "other" if v == "mitochondrion" else v
        tr = tr[tr.localization_class != "NA"]; te = te[te.localization_class != "NA"]
        return tr.localization_class.map(f), te.localization_class.map(f)
    if task == "species_label":
        tr = tr[tr.species_label != "NA"]; te = te[te.species_label != "NA"]
        return tr.species_label, te.species_label
    if task == "tm_label":
        tr = tr[tr.tm_label != "NA"]; te = te[te.tm_label != "NA"]
        return tr.tm_label.astype(float).astype(int), te.tm_label.astype(float).astype(int)
    # disorder / acc : binarisation à la médiane TRAIN
    tr = tr[tr[task] != "NA"]; te = te[te[task] != "NA"]
    med = tr[task].astype(float).median()
    return ((tr[task].astype(float) > med).astype(int),
            (te[task].astype(float) > med).astype(int))


def main():
    df = pd.read_csv(LABELS_TSV, sep="\t", keep_default_na=False)
    df = df.set_index("pdb_chain")
    tr_ids = [i for i in load_ids(ID_TRAIN) if i in df.index]
    te_ids = [i for i in load_ids(ID_TEST) if i in df.index]
    df_tr, df_te = df.loc[tr_ids].reset_index(), df.loc[te_ids].reset_index()

    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    print("=== Distributions des classes (train / test) ===")
    for ax, task in zip(axes.ravel(), TASKS):
        ytr, yte = prep(df_tr, df_te, task)
        ctr = ytr.value_counts().sort_index()
        cte = yte.value_counts().sort_index()
        classes = sorted(set(ctr.index) | set(cte.index), key=str)
        x = np.arange(len(classes))
        w = 0.4
        ptr = [ctr.get(c, 0) for c in classes]
        pte = [cte.get(c, 0) for c in classes]
        ax.bar(x - w/2, ptr, w, label=f"train (n={len(ytr)})", color="#4C72B0")
        ax.bar(x + w/2, pte, w, label=f"test (n={len(yte)})", color="#C44E52")
        ax.set_title(task)
        ax.set_xticks(x); ax.set_xticklabels([str(c) for c in classes],
                                             rotation=20, ha="right", fontsize=8)
        ax.legend(fontsize=8)

        print(f"\n[{task}] train={len(ytr)} test={len(yte)}")
        for c in classes:
            ptr_pct = 100*ctr.get(c, 0)/len(ytr) if len(ytr) else 0
            pte_pct = 100*cte.get(c, 0)/len(yte) if len(yte) else 0
            print(f"  {str(c):14s} train {ctr.get(c,0):4d} ({ptr_pct:4.1f}%) | "
                  f"test {cte.get(c,0):4d} ({pte_pct:4.1f}%)")

    # masquer les panneaux inutilisés (7 tâches sur une grille 2x4)
    for ax in axes.ravel()[len(TASKS):]:
        ax.axis("off")

    plt.suptitle("Distribution des classes cibles — tâches globales", fontsize=14)
    plt.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(FIG, dpi=160, bbox_inches="tight")
    print(f"\n[OK] -> {FIG}")


if __name__ == "__main__":
    main()
