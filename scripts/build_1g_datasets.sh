#!/bin/bash
# Build per-residue CSV datasets for 1G baseline models.
# Requires embeddings in Datasets/embeddings/full/{model}.safetensors
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/build_1g_datasets_$(date +%Y%m%d_%H%M).log"

cd "$ROOT"

DYNAMICS_VARS="rmsf neq bfact"
DSSP_VARS="acc sec3 sec8"
MD_DIR="Datasets/ATLAS/labels_md"
DSSP_DIR="Datasets/ATLAS/data"
TRAIN_IDS="deciphering/id_train.txt"
TEST_IDS="deciphering/id_test.txt"

echo "[$(date)] Starting 1G dataset build" | tee -a "$LOG"

# esm2_t33_650M_UR50D → emb key esm33
PLM_NAME="esm2_t33_650M_UR50D"; EMB="esm33"
EMB_FILE="Datasets/embeddings/full/${PLM_NAME}.safetensors"
if [ -f "$EMB_FILE" ]; then
    for Y in $DYNAMICS_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
    for Y in $DSSP_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
else
    echo "[$(date)] [SKIP] $PLM_NAME - missing safetensors" | tee -a "$LOG"
fi

# ankh_large → emb key ankhl
PLM_NAME="ankh_large"; EMB="ankhl"
EMB_FILE="Datasets/embeddings/full/${PLM_NAME}.safetensors"
if [ -f "$EMB_FILE" ]; then
    for Y in $DYNAMICS_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
    for Y in $DSSP_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
else
    echo "[$(date)] [SKIP] $PLM_NAME - missing safetensors" | tee -a "$LOG"
fi

# prot_t5_xl_uniref50 → emb key t5
PLM_NAME="prot_t5_xl_uniref50"; EMB="t5"
EMB_FILE="Datasets/embeddings/full/${PLM_NAME}.safetensors"
if [ -f "$EMB_FILE" ]; then
    for Y in $DYNAMICS_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
    for Y in $DSSP_VARS; do
        echo "[$(date)] === $EMB x $Y ===" | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TRAIN_IDS" "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
        python analysis/full_prot_emb_2g.py "$TEST_IDS"  "$EMB_FILE" -emb "$EMB" --y "$Y" -dssp_dir "$DSSP_DIR" -md_dir "$MD_DIR" 2>&1 | tee -a "$LOG"
    done
else
    echo "[$(date)] [SKIP] $PLM_NAME - missing safetensors" | tee -a "$LOG"
fi

echo "[$(date)] All 1G datasets generated" | tee -a "$LOG"
