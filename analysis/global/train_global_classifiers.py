#!/usr/bin/env python3
"""ÉTAPE 3 — Entraîne LogReg + RandomForest + MLP sur les datasets globaux.

Pour chaque (modèle × tâche × classifieur) :
  - 15-fold stratified CV sur le train set
  - évaluation finale sur le test set
  - métrique : macro-F1 (fold, localization) ou MCC (tm, disorder, acc)

Préparation par tâche (sans fuite — seuils calculés sur TRAIN uniquement) :
  - fold_label        : garder a/b/c/d, dropper e/f/g et NA  (macro-F1)
  - tm_label          : binaire 0/1, dropper NA              (MCC)
  - localization_class: fusionner mitochondrion -> other, dropper NA (macro-F1)
  - disorder_global   : binariser à la médiane du TRAIN      (MCC)
  - acc_mean          : binariser à la médiane du TRAIN      (MCC)

Usage :
    python train_global_classifiers.py [model ...] [--tasks t1,t2] [--quick]
    --quick : 3-fold au lieu de 15 (test rapide)

Sortie : results/global_results.csv
    model | task | clf | metric | score_cv_mean | score_cv_std | score_test | n_train | n_test
"""
import sys
import csv
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import f1_score, matthews_corrcoef

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = ROOT / "datasets_global"
OUT_CSV = ROOT / "results" / "global_results.csv"

ALL_MODELS = [
    "esmc_300M", "esmc_600M", "ankh2_large",
    "esm2_t33_650M_UR50D", "ankh_large", "prot_t5_xl_uniref50",
]

# task -> (metric_name, scoring_for_cv)
TASKS = {
    "fold_label": "f1_macro",
    "localization_class": "f1_macro",
    "species_label": "f1_macro",
    "tm_label": "mcc",
    "disorder_global": "mcc",
    "acc_mean": "mcc",
    "aggregation_score": "mcc",
}

LABEL_COLS = ["fold_label", "tm_label", "localization_class",
              "disorder_global", "acc_mean", "species_label",
              "aggregation_score"]
FOLD_KEEP = {"a", "b", "c", "d"}


def make_clfs():
    return {
        "LogReg": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, class_weight="balanced")),
        "RF": RandomForestClassifier(
            n_estimators=300, class_weight="balanced", n_jobs=-1, random_state=0),
        "MLP": make_pipeline(
            StandardScaler(),
            MLPClassifier(hidden_layer_sizes=(256, 64), max_iter=400,
                          early_stopping=True, random_state=0)),
    }


def load_split(model, split):
    # keep_default_na=False : garder "NA" comme chaîne (sinon pandas -> NaN float)
    df = pd.read_csv(DATA_DIR / f"{model}_{split}.csv", keep_default_na=False)
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    return df, dim_cols


def prepare_task(task, df_train, df_test, dim_cols):
    """Retourne X_train, y_train, X_test, y_test préparés pour la tâche.
    Les seuils continus sont calculés sur le TRAIN uniquement."""
    tr, te = df_train.copy(), df_test.copy()

    if task == "fold_label":
        tr = tr[tr["fold_label"].isin(FOLD_KEEP)]
        te = te[te["fold_label"].isin(FOLD_KEEP)]
        ytr, yte = tr["fold_label"].astype(str), te["fold_label"].astype(str)

    elif task == "localization_class":
        def fuse(v):
            return "other" if v == "mitochondrion" else v
        tr = tr[tr["localization_class"] != "NA"]
        te = te[te["localization_class"] != "NA"]
        ytr = tr["localization_class"].map(fuse).astype(str)
        yte = te["localization_class"].map(fuse).astype(str)

    elif task == "species_label":
        # 4 espèces ATLAS >= 30 prot ; "NA" (autres organismes) droppé
        tr = tr[tr["species_label"] != "NA"]
        te = te[te["species_label"] != "NA"]
        ytr, yte = tr["species_label"].astype(str), te["species_label"].astype(str)

    elif task == "tm_label":
        tr = tr[tr["tm_label"] != "NA"]
        te = te[te["tm_label"] != "NA"]
        ytr = tr["tm_label"].astype(float).astype(int)
        yte = te["tm_label"].astype(float).astype(int)

    elif task in ("disorder_global", "acc_mean", "aggregation_score"):
        # Seuil de binarisation = MÉDIANE DU TRAIN (pas de fuite).
        # FIX 6 — seuils mesurés sur ce dataset :
        #   disorder_global  : médiane train ≈ 1.27 Å (RMSF moyen)
        #   acc_mean         : médiane train ≈ 57.5 Å² (ASA ABSOLUE DSSP)
        #   aggregation_score: médiane train ≈ -0.20 (proxy A3D, sans unité)
        # NB : acc_mean est l'ASA absolue, PAS la RSA relative — le seuil
        # littérature 0.35 RSA (Benhamouche) n'est donc pas directement
        # transposable. La médiane donne des classes ~50/50 (équilibre neutre).
        # NB2 : aggregation_score est un proxy structural (Python 3), pas A3D
        # natif — le seuil -0.7 d'A3D (Benhamouche) n'est pas transposable.
        tr = tr[tr[task] != "NA"]
        te = te[te[task] != "NA"]
        med = tr[task].astype(float).median()       # seuil = médiane TRAIN
        ytr = (tr[task].astype(float) > med).astype(int)
        yte = (te[task].astype(float) > med).astype(int)
    else:
        raise ValueError(task)

    Xtr = tr[dim_cols].to_numpy(dtype=np.float32)
    Xte = te[dim_cols].to_numpy(dtype=np.float32)
    # Encodage entier des labels (fit sur train) : évite le bug MLP early_stopping
    # avec labels string et homogénéise le typage.
    ytr_s = np.asarray(ytr).astype(str)
    yte_s = np.asarray(yte).astype(str)
    le = LabelEncoder().fit(ytr_s)
    # garder en test uniquement les classes vues en train
    seen = np.isin(yte_s, le.classes_)
    Xte, yte_s = Xte[seen], yte_s[seen]
    return Xtr, le.transform(ytr_s), Xte, le.transform(yte_s)


def score_fn(metric):
    if metric == "f1_macro":
        return lambda yt, yp: f1_score(yt, yp, average="macro")
    return lambda yt, yp: matthews_corrcoef(yt, yp)


def cv_scoring(metric):
    return "f1_macro" if metric == "f1_macro" else "matthews_corrcoef"


def run(models, tasks, n_splits):
    rows = []
    fold_rows = []        # 1 ligne par (model, task, clf, fold)
    metric_score = {t: score_fn(m) for t, m in TASKS.items()}

    for model in models:
        df_train, dim_cols = load_split(model, "train")
        df_test, _ = load_split(model, "test")
        print(f"\n########## {model} ({len(dim_cols)} dims) ##########")

        for task in tasks:
            metric = TASKS[task]
            Xtr, ytr, Xte, yte = prepare_task(task, df_train, df_test, dim_cols)
            nclass = len(np.unique(ytr))
            print(f"\n--- {task} [{metric}] | train={len(ytr)} test={len(yte)} "
                  f"| {nclass} classes ---")
            if len(ytr) < 30 or nclass < 2:
                print("  [SKIP] trop peu d'exemples / classes")
                continue

            # nb de folds limité par la plus petite classe
            min_class = np.min(np.bincount(
                pd.factorize(ytr)[0]))
            k = max(2, min(n_splits, min_class))
            skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=0)

            for clf_name, clf in make_clfs().items():
                try:
                    cv = cross_val_score(clf, Xtr, ytr, cv=skf,
                                         scoring=cv_scoring(metric), n_jobs=-1)
                    clf.fit(Xtr, ytr)
                    yp = clf.predict(Xte)
                    test_score = metric_score[task](yte, yp)
                    print(f"  {clf_name:7s} CV={cv.mean():.3f}±{cv.std():.3f} "
                          f"test={test_score:.3f}")
                    rows.append({
                        "model": model, "task": task, "clf": clf_name,
                        "metric": metric,
                        "score_cv_mean": round(float(cv.mean()), 4),
                        "score_cv_std": round(float(cv.std()), 4),
                        "score_test": round(float(test_score), 4),
                        "n_train": len(ytr), "n_test": len(yte),
                    })
                    for fi, sc in enumerate(cv):
                        fold_rows.append({
                            "model": model, "task": task, "clf": clf_name,
                            "metric": metric, "fold": fi,
                            "score_cv": round(float(sc), 4), "k": k,
                        })
                except Exception as e:
                    print(f"  [ERR] {clf_name}: {e}")

    return rows, fold_rows


def main():
    argv = [a for a in sys.argv[1:]]
    quick = "--quick" in argv
    argv = [a for a in argv if a != "--quick"]
    tasks = list(TASKS.keys())
    if "--tasks" in argv:
        i = argv.index("--tasks")
        tasks = argv[i + 1].split(",")
        del argv[i:i + 2]
    models = argv if argv else ALL_MODELS
    n_splits = 3 if quick else 15

    print(f"Modèles: {models}\nTâches: {tasks}\nCV: {n_splits}-fold")
    rows, fold_rows = run(models, tasks, n_splits)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["model", "task", "clf", "metric", "score_cv_mean",
            "score_cv_std", "score_test", "n_train", "n_test"]
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[OK] {len(rows)} lignes -> {OUT_CSV}")

    fold_csv = OUT_CSV.parent / "global_results_per_fold.csv"
    fcols = ["model", "task", "clf", "metric", "fold", "score_cv", "k"]
    with fold_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fcols)
        w.writeheader()
        w.writerows(fold_rows)
    print(f"[OK] {len(fold_rows)} lignes -> {fold_csv}")


if __name__ == "__main__":
    main()
