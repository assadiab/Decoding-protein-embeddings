#!/usr/bin/env python3
"""Per-residue pipeline with 3 classifiers (LogReg / RF / MLP).

Le pipeline local d'origine n'utilisait que le Decision Tree (dt_results_*.csv,
kept for comparability with Soufir). Here we add LogReg + RF + MLP on the
same per-residue datasets (datasets_emb_*/), to align the methodology with the
pipeline global.

The per-residue datasets are large (~230k residues/task). To keep the
runtime reasonable we **subsample the TRAIN to SUBSAMPLE residues** for
CV and training. The final evaluation uses the **full test set**.

Metric: F1 (binary for rmsf/neq/bfact/acc, macro for sec3/sec8) - same
convention as dt_results_*.csv, for direct comparability.

Output: results/local_results_full.csv
    model | task | clf | metric | score_cv_mean | score_cv_std | score_test | n_train_sub | n_test

Usage : python local_train_full.py [--quick] [--subsample N] [--folds K]
"""
import sys
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import f1_score
import warnings
warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parent.parent
OUT_CSV = ROOT / "results" / "local_results_full.csv"

# datasets_emb_* directory -> clean model name (aligned with global)
MODELS = {
    "esmc_300m": "esmc_300M",
    "esmc_600m": "esmc_600M",
    "ankh2_large": "ankh2_large",
    "esm33": "esm2_t33_650M_UR50D",
    "ankhl": "ankh_large",
    "t5": "prot_t5_xl_uniref50",
}
TASKS = ["rmsf", "neq", "bfact", "acc", "sec3", "sec8"]
MULTICLASS = {"sec3", "sec8"}

SUBSAMPLE = 30000
FOLDS = 3
RNG = 0


def make_clfs():
    return {
        "LogReg": make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=1000, class_weight="balanced", n_jobs=-1)),
        "RF": RandomForestClassifier(
            n_estimators=100, class_weight="balanced", n_jobs=-1, random_state=RNG),
        "MLP": make_pipeline(
            StandardScaler(),
            MLPClassifier(hidden_layer_sizes=(128,), max_iter=150,
                          early_stopping=True, random_state=RNG)),
    }


def load_xy(path, task, subsample=None):
    df = pd.read_csv(path)
    emb_cols = [c for c in df.columns if c.startswith("emb")]
    y = df[f"{task}_categ"].to_numpy()
    X = df[emb_cols].to_numpy(dtype=np.float32)
    if subsample and len(y) > subsample:
        rng = np.random.default_rng(RNG)
        idx = rng.choice(len(y), size=subsample, replace=False)
        X, y = X[idx], y[idx]
    return X, y


def f1_metric(task):
    avg = "macro" if task in MULTICLASS else "binary"
    return lambda yt, yp: f1_score(yt, yp, average=avg)


def cv_scoring(task):
    return "f1_macro" if task in MULTICLASS else "f1"


def main():
    argv = sys.argv[1:]
    quick = "--quick" in argv
    sub = SUBSAMPLE
    folds = FOLDS
    if "--subsample" in argv:
        sub = int(argv[argv.index("--subsample") + 1])
    if "--folds" in argv:
        folds = int(argv[argv.index("--folds") + 1])
    if quick:
        sub, folds = 5000, 2
    sel_models = [a for a in argv if a in MODELS]
    model_dirs = sel_models if sel_models else list(MODELS.keys())

    print(f"subsample={sub} folds={folds} models={model_dirs}")
    rows = []
    for mdir in model_dirs:
        mname = MODELS[mdir]
        base = ROOT / f"datasets_emb_{mdir}"
        print(f"\n########## {mname} ##########")
        for task in TASKS:
            ftr = base / f"emb_all_positions_{task}_train.csv"
            fte = base / f"emb_all_positions_{task}_test.csv"
            if not ftr.exists() or not fte.exists():
                print(f"  [SKIP] {task}: missing CSV")
                continue
            Xtr, ytr = load_xy(ftr, task, subsample=sub)
            Xte, yte = load_xy(fte, task)              # full test set

            # encode for type homogeneity (multiclass string)
            le = LabelEncoder().fit(np.concatenate([ytr.astype(str), yte.astype(str)]))
            ytr_e, yte_e = le.transform(ytr.astype(str)), le.transform(yte.astype(str))

            metric = "f1_macro" if task in MULTICLASS else "f1"
            scorer = f1_metric(task)
            print(f"\n--- {task} [{metric}] | train_sub={len(ytr_e)} test={len(yte_e)} ---")
            skf = StratifiedKFold(n_splits=folds, shuffle=True, random_state=RNG)
            for clf_name, clf in make_clfs().items():
                try:
                    cv = cross_val_score(clf, Xtr, ytr_e, cv=skf,
                                         scoring=cv_scoring(task), n_jobs=-1)
                    clf.fit(Xtr, ytr_e)
                    ts = scorer(yte_e, clf.predict(Xte))
                    print(f"  {clf_name:7s} CV={cv.mean():.3f}±{cv.std():.3f} test={ts:.3f}")
                    rows.append({
                        "model": mname, "task": task, "clf": clf_name,
                        "metric": metric,
                        "score_cv_mean": round(float(cv.mean()), 4),
                        "score_cv_std": round(float(cv.std()), 4),
                        "score_test": round(float(ts), 4),
                        "n_train_sub": len(ytr_e), "n_test": len(yte_e),
                    })
                except Exception as e:
                    print(f"  [ERR] {clf_name}: {e}")

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["model", "task", "clf", "metric", "score_cv_mean", "score_cv_std",
            "score_test", "n_train_sub", "n_test"]
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[OK] {len(rows)} lignes -> {OUT_CSV}")


if __name__ == "__main__":
    main()
