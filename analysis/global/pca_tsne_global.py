#!/usr/bin/env python3
"""FIX 5 — PCA & t-SNE 2D des embeddings mean-pooled (espace global).

Montre que les embeddings par protéine s'organisent en clusters biologiques.
Pour 2-3 modèles clés, projette les embeddings train en 2D (PCA puis t-SNE),
coloriés par fold_label puis par tm_label.

Sortie :
  results/figures/global_pca_{model}_{label}.png
  results/figures/global_tsne_{model}_{label}.png

Usage : python pca_tsne_global.py [model ...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = ROOT / "datasets_global"
FIG_DIR = ROOT / "results" / "figures"

DEFAULT_MODELS = ["ankh2_large", "esmc_600M", "prot_t5_xl_uniref50"]
LABELS = ["fold_label", "tm_label"]
FOLD_KEEP = {"a", "b", "c", "d"}


def scatter(emb2d, labels, title, out, discrete=True):
    plt.figure(figsize=(7, 6))
    cats = sorted(set(labels))
    cmap = plt.get_cmap("tab10")
    for i, c in enumerate(cats):
        m = labels == c
        plt.scatter(emb2d[m, 0], emb2d[m, 1], s=14, alpha=0.7,
                    color=cmap(i % 10), label=str(c))
    plt.title(title)
    plt.xlabel("comp 1"); plt.ylabel("comp 2")
    plt.legend(title="classe", fontsize=8, markerscale=1.3)
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    -> {out.name}")


def run_model(model):
    f = DATA_DIR / f"{model}_train.csv"
    if not f.exists():
        print(f"[SKIP] {model} : CSV absent")
        return
    df = pd.read_csv(f, keep_default_na=False)
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    X = df[dim_cols].to_numpy(dtype=np.float32)
    Xs = StandardScaler().fit_transform(X)
    print(f"\n### {model} ({len(df)} prot, {len(dim_cols)} dims) ###")

    pca = PCA(n_components=2, random_state=0).fit(Xs)
    X_pca = pca.transform(Xs)
    var = pca.explained_variance_ratio_ * 100

    print("  t-SNE…")
    X_tsne = TSNE(n_components=2, perplexity=30, init="pca",
                  random_state=0).fit_transform(Xs)

    for label in LABELS:
        if label == "fold_label":
            mask = df["fold_label"].isin(FOLD_KEEP).to_numpy()
        else:
            mask = (df["tm_label"] != "NA").to_numpy()
        lab = df.loc[mask, label].to_numpy().astype(str)

        scatter(X_pca[mask], lab,
                f"PCA — {model} · {label}  (PC1 {var[0]:.0f}%, PC2 {var[1]:.0f}%)",
                FIG_DIR / f"global_pca_{model}_{label}.png")
        scatter(X_tsne[mask], lab,
                f"t-SNE — {model} · {label}",
                FIG_DIR / f"global_tsne_{model}_{label}.png")


def main():
    models = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_MODELS
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for m in models:
        run_model(m)
    print("\n[OK] FIX 5 terminé.")


if __name__ == "__main__":
    main()
