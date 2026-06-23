#!/bin/bash
# SHAP feature importance analysis for all 2G models x all tasks.
# Launch from project root: bash scripts/run_shap.sh

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/shap_$(date +%Y%m%d_%H%M).log"
FIG_DIR="$ROOT/results/figures"
mkdir -p "$FIG_DIR"

cd "$ROOT/deciphering/scripts"

MODELS="esmc_300m esmc_600m ankh2_large"
VARS="rmsf neq bfact acc sec3 sec8"

echo "[$(date)] Starting SHAP" | tee -a "$LOG"

for EMB in $MODELS; do
    EMB_UPPER=$(echo "$EMB" | tr '[:lower:]' '[:upper:]' | tr -d '_')
    for Y in $VARS; do
        TRAIN="../../datasets_emb_${EMB}/emb_all_positions_${Y}_train.csv"
        TEST="../../datasets_emb_${EMB}/emb_all_positions_${Y}_test.csv"
        [ -f "$TRAIN" ] || continue

        echo "[$(date)] === SHAP $EMB x $Y ===" | tee -a "$LOG"
        python pca_full_prot.py "$TRAIN" "$TEST" --analysis shap 2>&1 | tee -a "$LOG" \
            || echo "[WARN] SHAP $EMB x $Y failed" | tee -a "$LOG"

        OUT_DIR="output_${EMB_UPPER}"
        for f in "$OUT_DIR"/shap_*.png; do
            [ -f "$f" ] || continue
            mv "$f" "$FIG_DIR/${EMB}_$(basename "$f")" 2>/dev/null || true
        done
        echo "[$(date)] $EMB x $Y done" | tee -a "$LOG"
    done
done

echo "[$(date)] SHAP done" | tee -a "$LOG"
