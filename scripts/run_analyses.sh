#!/bin/bash
# Run pca_full_prot.py analyses on all 2G model × variable combinations.
# Launch from project root: bash scripts/run_analyses.sh
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/analyses_$(date +%Y%m%d_%H%M).log"
SCRIPT="$ROOT/deciphering/scripts/pca_full_prot.py"
FIG_DIR="$ROOT/results/figures"
mkdir -p "$FIG_DIR"

cd "$ROOT/deciphering/scripts"

MODELS="esmc_300m esmc_600m ankh2_large"
VARS="rmsf neq bfact acc sec3 sec8"

echo "[$(date)] Démarrage analyses DT" | tee -a "$LOG"

for EMB in $MODELS; do
    for Y in $VARS; do
        TRAIN="../../datasets_emb_${EMB}/emb_all_positions_${Y}_train.csv"
        TEST="../../datasets_emb_${EMB}/emb_all_positions_${Y}_test.csv"

        if [ ! -f "$TRAIN" ] || [ ! -f "$TEST" ]; then
            echo "[$(date)] [SKIP] $EMB × $Y — CSV manquant" | tee -a "$LOG"
            continue
        fi

        echo "[$(date)] === $EMB × $Y (dt) ===" | tee -a "$LOG"
        python pca_full_prot.py "$TRAIN" "$TEST" --analysis dt 2>&1 | tee -a "$LOG"
        echo "[$(date)] $EMB × $Y terminé" | tee -a "$LOG"

        # Move figures to results/figures/ with model prefix
        EMB_TYPE=$(echo "$EMB" | tr '[:lower:]' '[:upper:]' | tr '_' '')
        OUT_DIR="output_${EMB_TYPE}"
        if [ -d "$OUT_DIR" ]; then
            for f in "$OUT_DIR"/*.png "$OUT_DIR"/*.pdf; do
                [ -f "$f" ] || continue
                base=$(basename "$f")
                mv "$f" "$FIG_DIR/${EMB}_${base}" 2>/dev/null || true
            done
        fi
    done
done

echo "[$(date)] Toutes les analyses terminées" | tee -a "$LOG"
