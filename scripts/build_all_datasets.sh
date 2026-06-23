#!/bin/bash
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/build_datasets_$(date +%Y%m%d_%H%M).log"

cd "$ROOT"

MODELS=("esmc_300m:Datasets/embeddings/full/esmc_300M.safetensors"
        "esmc_600m:Datasets/embeddings/full/esmc_600M.safetensors"
        "ankh2_large:Datasets/embeddings/full/ankh2_large.safetensors")

# Dynamics only (no DSSP files available for acc/sec3/sec8)
DYNAMICS_VARS="rmsf neq bfact"

MD_DIR="Datasets/ATLAS/labels_md"
TRAIN_IDS="deciphering/id_train.txt"
TEST_IDS="deciphering/id_test.txt"

echo "[$(date)] Starting dataset build - dynamics (rmsf, neq, bfact)" | tee -a "$LOG"

for entry in "${MODELS[@]}"; do
    EMB="${entry%%:*}"
    EMB_FILE="${entry##*:}"

    for Y in $DYNAMICS_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"

        python analysis/full_prot_emb_2g.py \
            "$TRAIN_IDS" "$EMB_FILE" \
            -emb "$EMB" --y "$Y" \
            -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"

        python analysis/full_prot_emb_2g.py \
            "$TEST_IDS" "$EMB_FILE" \
            -emb "$EMB" --y "$Y" \
            -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"

        echo "[$(date)] $EMB x $Y done" | tee -a "$LOG"
    done
done

echo "[$(date)] All datasets generated" | tee -a "$LOG"
