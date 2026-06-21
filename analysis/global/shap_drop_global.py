#!/usr/bin/env python3
"""ÉTAPE 4 — SHAP + drop curves sur le meilleur classifieur par (modèle, tâche).

Lit results/global_results.csv pour identifier le meilleur clf (score_test) par
(model, task), ré-entraîne, puis :
  - SHAP : top-20 dimensions les plus importantes -> results/figures/global_shap_{model}_{task}.png
  - Drop curves : masque 0/25/50/75/90/95/100% des dimensions (3 répétitions),
    trace le score vs % supprimé -> results/figures/global_drop_{model}_{task}.png

Pour limiter le coût, SHAP utilise un sous-échantillon du test set.

Usage :
    python shap_drop_global.py [model ...] [--tasks t1,t2]
    # défaut : tous les modèles, toutes les tâches
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from train_global_classifiers import (
    ALL_MODELS, TASKS, load_split, prepare_task, make_clfs, score_fn)

ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS = ROOT / "results"
FIG_DIR = RESULTS / "figures"
GLOBAL_CSV = RESULTS / "global_results.csv"

DROP_FRACTIONS = [0.0, 0.25, 0.50, 0.75, 0.90, 0.95, 1.0]
N_REPEAT = 3
SHAP_SAMPLE = 150          # taille sous-échantillon pour SHAP


def best_clf_per_task(df):
    """Retourne {(model, task): clf_name} d'après score_test max."""
    best = {}
    for (model, task), grp in df.groupby(["model", "task"]):
        row = grp.loc[grp["score_test"].idxmax()]
        best[(model, task)] = row["clf"]
    return best


def shap_plot(clf, clf_name, Xtr, Xte, model, task):
    import shap
    rng = np.random.default_rng(0)
    idx = rng.choice(len(Xte), size=min(SHAP_SAMPLE, len(Xte)), replace=False)
    Xs = Xte[idx]
    try:
        if clf_name == "RF":
            explainer = shap.TreeExplainer(clf)
            sv = explainer.shap_values(Xs)
            # multiclasse -> liste ; moyenne sur classes
            if isinstance(sv, list):
                arr = np.mean([np.abs(s) for s in sv], axis=0)
            else:
                arr = np.abs(sv)
                if arr.ndim == 3:
                    arr = arr.mean(axis=2)
        else:
            bg = shap.sample(Xtr, 50, random_state=0)
            explainer = shap.KernelExplainer(clf.predict_proba, bg)
            sv = explainer.shap_values(Xs, nsamples=100, silent=True)
            arr = np.mean([np.abs(s) for s in sv], axis=0) if isinstance(sv, list) \
                else np.abs(sv)
            if arr.ndim == 3:
                arr = arr.mean(axis=2)
    except Exception as e:
        print(f"    [SHAP ERR] {model}/{task}: {e}")
        return

    mean_abs = arr.mean(axis=0)               # importance moyenne par dimension
    top = np.argsort(mean_abs)[::-1][:20]
    plt.figure(figsize=(7, 6))
    plt.barh([f"dim_{i}" for i in top][::-1], mean_abs[top][::-1], color="#4C72B0")
    plt.xlabel("Importance SHAP moyenne (|valeur|)")
    plt.title(f"SHAP top-20 — {model} · {task} ({clf_name})")
    plt.tight_layout()
    out = FIG_DIR / f"global_shap_{model}_{task}.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"    SHAP -> {out.name}")


def drop_curve(clf_factory, Xtr, ytr, Xte, yte, scorer, model, task):
    rng = np.random.default_rng(0)
    dim = Xtr.shape[1]
    means, stds = [], []
    for frac in DROP_FRACTIONS:
        reps = []
        n_drop = int(round(frac * dim))
        for r in range(N_REPEAT):
            if n_drop == 0:
                Xtr2, Xte2 = Xtr, Xte
            elif n_drop >= dim:
                # tout masqué -> features à 0 (baseline aléatoire)
                Xtr2 = np.zeros_like(Xtr)
                Xte2 = np.zeros_like(Xte)
            else:
                drop = rng.choice(dim, size=n_drop, replace=False)
                Xtr2, Xte2 = Xtr.copy(), Xte.copy()
                Xtr2[:, drop] = 0.0
                Xte2[:, drop] = 0.0
            clf = clf_factory()
            try:
                clf.fit(Xtr2, ytr)
                reps.append(scorer(yte, clf.predict(Xte2)))
            except Exception:
                reps.append(np.nan)
        means.append(np.nanmean(reps))
        stds.append(np.nanstd(reps))

    plt.figure(figsize=(7, 5))
    xs = [f * 100 for f in DROP_FRACTIONS]
    plt.errorbar(xs, means, yerr=stds, marker="o", capsize=3, color="#C44E52")
    plt.xlabel("% dimensions masquées")
    plt.ylabel(TASKS[task])
    plt.title(f"Drop curve — {model} · {task}")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = FIG_DIR / f"global_drop_{model}_{task}.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"    drop -> {out.name}")


def main():
    argv = sys.argv[1:]
    tasks = list(TASKS.keys())
    if "--tasks" in argv:
        i = argv.index("--tasks")
        tasks = argv[i + 1].split(",")
        del argv[i:i + 2]
    models = argv if argv else ALL_MODELS

    if not GLOBAL_CSV.exists():
        print("global_results.csv absent — lancer l'Étape 3 d'abord.")
        sys.exit(1)
    df = pd.read_csv(GLOBAL_CSV)
    best = best_clf_per_task(df)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    for model in models:
        df_train, dim_cols = load_split(model, "train")
        df_test, _ = load_split(model, "test")
        print(f"\n########## {model} ##########")
        for task in tasks:
            clf_name = best.get((model, task))
            if clf_name is None:
                print(f"  [SKIP] {task} : pas de résultat")
                continue
            Xtr, ytr, Xte, yte = prepare_task(task, df_train, df_test, dim_cols)
            if len(ytr) < 30 or len(np.unique(ytr)) < 2:
                print(f"  [SKIP] {task} : trop peu d'exemples")
                continue
            print(f"  --- {task} (meilleur={clf_name}) ---")
            clf = make_clfs()[clf_name]
            clf.fit(Xtr, ytr)
            shap_plot(clf, clf_name, Xtr, Xte, model, task)
            scorer = score_fn(TASKS[task])
            drop_curve(lambda: make_clfs()[clf_name], Xtr, ytr, Xte, yte,
                       scorer, model, task)

    print("\n[OK] Étape 4 terminée.")


if __name__ == "__main__":
    main()
