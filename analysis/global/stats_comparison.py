#!/usr/bin/env python3
"""Statistical test 2G vs 1G per global task.

For each task, gather the 15-fold CV scores of the 3 2G models
on one side and the 3 1G baselines on the other (best classifier per model),
and test whether the difference is significant (Mann-Whitney U, two-sided).

Source: results/global_results_per_fold.csv (+ global_results.csv to
pick the best clf per model x task).

Output: results/global_stats_comparison.csv
    task | mean_2g | mean_1g | delta | p_value | significant
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent.parent
RESULTS = ROOT / "results"
AGG = RESULTS / "global_results.csv"
PERFOLD = RESULTS / "global_results_per_fold.csv"
OUT = RESULTS / "global_stats_comparison.csv"

MODELS_2G = ["esmc_300M", "esmc_600M", "ankh2_large"]
MODELS_1G = ["esm2_t33_650M_UR50D", "ankh_large", "prot_t5_xl_uniref50"]
TASK_ORDER = ["fold_label", "localization_class", "species_label",
              "tm_label", "disorder_global", "acc_mean", "aggregation_score"]


def best_clf(agg, model, task):
    sub = agg[(agg.model == model) & (agg.task == task)]
    if sub.empty:
        return None
    return sub.loc[sub.score_test.idxmax(), "clf"]


def fold_scores(perfold, agg, models, task):
    """Concatenate the per-fold scores of the best clf of each model."""
    vals = []
    for m in models:
        clf = best_clf(agg, m, task)
        if clf is None:
            continue
        s = perfold[(perfold.model == m) & (perfold.task == task) &
                    (perfold.clf == clf)]["score_cv"].to_numpy()
        vals.append(s)
    return np.concatenate(vals) if vals else np.array([])


def main():
    agg = pd.read_csv(AGG)
    perfold = pd.read_csv(PERFOLD)

    rows = []
    print("=== Test 2G vs 1G (Mann-Whitney U, scores par fold du meilleur clf) ===\n")
    for task in TASK_ORDER:
        s2 = fold_scores(perfold, agg, MODELS_2G, task)
        s1 = fold_scores(perfold, agg, MODELS_1G, task)
        if len(s2) == 0 or len(s1) == 0:
            print(f"[SKIP] {task}")
            continue
        try:
            u, p = mannwhitneyu(s2, s1, alternative="two-sided")
        except ValueError:
            p = np.nan
        delta = s2.mean() - s1.mean()
        sig = "yes" if (p == p and p < 0.05) else "no"
        arrow = "2G>1G" if delta > 0 else "1G>2G"
        print(f"{task:20s} 2G={s2.mean():.3f} 1G={s1.mean():.3f} "
              f"delta={delta:+.3f} ({arrow}) p={p:.4f} sig={sig}")
        rows.append({
            "task": task,
            "mean_2g": round(float(s2.mean()), 4),
            "mean_1g": round(float(s1.mean()), 4),
            "delta": round(float(delta), 4),
            "p_value": round(float(p), 5) if p == p else "NA",
            "significant": sig,
            "n_folds_2g": len(s2), "n_folds_1g": len(s1),
        })

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print(f"\n[OK] -> {OUT}")


if __name__ == "__main__":
    main()
