#!/bin/bash
# Run DT analyses on all 1G baseline model x variable combinations.
# Launch from project root: bash scripts/run_analyses_1g.sh
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/analyses_1g_$(date +%Y%m%d_%H%M).log"
FIG_DIR="$ROOT/results/figures"
mkdir -p "$FIG_DIR"

cd "$ROOT/deciphering/scripts"

VARS="rmsf neq bfact acc sec3 sec8"

echo "[$(date)] Starting 1G DT analyses" | tee -a "$LOG"

run_model() {
    local EMB="$1"
    local PREFIX="$2"
    local EMB_UPPER="$3"

    for Y in $VARS; do
        TRAIN="../../datasets_emb_${EMB}/emb_all_positions_${Y}_train.csv"
        TEST="../../datasets_emb_${EMB}/emb_all_positions_${Y}_test.csv"

        if [ ! -f "$TRAIN" ] || [ ! -f "$TEST" ]; then
            echo "[$(date)] [SKIP] $EMB x $Y - missing CSV" | tee -a "$LOG"
            continue
        fi

        echo "[$(date)] === $EMB x $Y (dt) ===" | tee -a "$LOG"
        python pca_full_prot.py "$TRAIN" "$TEST" --analysis dt 2>&1 | tee -a "$LOG"
        echo "[$(date)] $EMB x $Y done" | tee -a "$LOG"

        OUT_DIR="output_${EMB_UPPER}"
        if [ -d "$OUT_DIR" ]; then
            for f in "$OUT_DIR"/*.png "$OUT_DIR"/*.pdf; do
                [ -f "$f" ] || continue
                mv "$f" "$FIG_DIR/${PREFIX}_$(basename "$f")" 2>/dev/null || true
            done
        fi
    done
}

run_model "esm33"  "esm2_650m"  "ESM33"
run_model "ankhl"  "ankh_large" "ANKHL"
run_model "t5"     "prot_t5"    "T5"

echo "[$(date)] All 1G analyses done" | tee -a "$LOG"
