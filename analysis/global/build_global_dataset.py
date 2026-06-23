#!/usr/bin/env python3
"""Build the global (per-protein) datasets by mean pooling.

For each model (safetensors), reduce each embedding [seq_len, dim] -> [dim] by
mean pooling, join with atlas_global_labels.tsv on pdb_chain, and write a train +
test CSV following the deciphering/id_{train,test}.txt split.

Usage:
    python build_global_dataset.py [model_name ...]
    # no argument: all models
    # e.g. python build_global_dataset.py esmc_300M   (subset test)

Output: datasets_global/{model}_train.csv and {model}_test.csv
Columns: pdb_chain, dim_0..dim_{N-1}, then the label columns.
"""
import sys
import csv
from pathlib import Path

import numpy as np
from safetensors import safe_open

ROOT = Path(__file__).resolve().parent.parent.parent
EMB_DIR = ROOT / "Datasets" / "embeddings" / "full"
LABELS_TSV = ROOT / "Datasets" / "ATLAS" / "global_labels" / "atlas_global_labels.tsv"
ID_TRAIN = ROOT / "deciphering" / "id_train.txt"
ID_TEST = ROOT / "deciphering" / "id_test.txt"
OUT_DIR = ROOT / "datasets_global"

LABEL_COLS = ["fold_label", "tm_label", "localization_class",
              "disorder_global", "acc_mean", "species_label",
              "aggregation_score"]

# All available models (file name without extension)
ALL_MODELS = [
    "esmc_300M", "esmc_600M", "ankh2_large",       # 2G
    "esm2_t33_650M_UR50D", "ankh_large", "prot_t5_xl_uniref50",  # 1G
]


def load_ids(path):
    return [l.strip() for l in path.read_text().splitlines() if l.strip()]


def load_labels():
    labels = {}
    with LABELS_TSV.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            labels[row["pdb_chain"]] = {c: row[c] for c in LABEL_COLS}
    return labels


def build_model(model, labels, train_ids, test_ids):
    emb_path = EMB_DIR / f"{model}.safetensors"
    if not emb_path.exists():
        print(f"[SKIP] {model}: missing safetensors")
        return

    print(f"\n=== {model} ===")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with safe_open(str(emb_path), "numpy") as f:
        emb_keys = set(f.keys())
        dim = f.get_tensor(next(iter(emb_keys))).shape[1]
        dim_cols = [f"dim_{i}" for i in range(dim)]
        header = ["pdb_chain"] + dim_cols + LABEL_COLS

        for split_name, ids in (("train", train_ids), ("test", test_ids)):
            out_csv = OUT_DIR / f"{model}_{split_name}.csv"
            written, skipped = 0, []
            with out_csv.open("w", newline="") as out:
                w = csv.writer(out)
                w.writerow(header)
                for pid in ids:
                    if pid not in emb_keys:
                        skipped.append((pid, "no_emb"))
                        continue
                    if pid not in labels:
                        skipped.append((pid, "no_label"))
                        continue
                    vec = f.get_tensor(pid).mean(axis=0)  # [dim]
                    row = [pid] + [f"{x:.6f}" for x in vec] + \
                          [labels[pid][c] for c in LABEL_COLS]
                    w.writerow(row)
                    written += 1
            print(f"  {split_name}: {written} written, {len(skipped)} skipped "
                  f"-> {out_csv.name}")
            for pid, reason in skipped:
                print(f"    [SKIP] {pid} ({reason})")


def main():
    models = sys.argv[1:] if len(sys.argv) > 1 else ALL_MODELS
    labels = load_labels()
    train_ids = load_ids(ID_TRAIN)
    test_ids = load_ids(ID_TEST)
    print(f"Labels: {len(labels)} | train: {len(train_ids)} | test: {len(test_ids)}")
    for model in models:
        build_model(model, labels, train_ids, test_ids)
    print("\n[OK] done.")


if __name__ == "__main__":
    main()
