#!/bin/bash
# Run PCA analysis for one representative task (rmsf) across all 3 models.
# Generates PCA scatter plots + explained variance in results/figures/
ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/results/pca_$(date +%Y%m%d_%H%M).log"
SCRIPT="$ROOT/deciphering/scripts/pca_full_prot.py"
FIG_DIR="$ROOT/results/figures"

cd "$ROOT/deciphering/scripts"

MODELS="esmc_300m esmc_600m ankh2_large"
# Run PCA on rmsf (dynamics) and sec3 (structure) as representative tasks
VARS="rmsf sec3"

echo "[$(date)] Démarrage PCA" | tee -a "$LOG"

for EMB in $MODELS; do
    for Y in $VARS; do
        TRAIN="../../datasets_emb_${EMB}/emb_all_positions_${Y}_train.csv"
        TEST="../../datasets_emb_${EMB}/emb_all_positions_${Y}_test.csv"

        [ -f "$TRAIN" ] || continue

        echo "[$(date)] === PCA $EMB × $Y ===" | tee -a "$LOG"
        python pca_full_prot.py "$TRAIN" "$TEST" --analysis pca 2>&1 | tee -a "$LOG" || echo "[$(date)] [WARN] PCA $EMB × $Y échec" | tee -a "$LOG"

        # Move + prefix figures
        EMB_UPPER=$(echo "$EMB" | tr '[:lower:]' '[:upper:]' | tr -d '_')
        OUT_DIR="output_${EMB_UPPER}"
        if [ -d "$OUT_DIR" ]; then
            for f in "$OUT_DIR"/*.png; do
                [ -f "$f" ] || continue
                mv "$f" "$FIG_DIR/${EMB}_$(basename $f)" 2>/dev/null || true
            done
        fi

        echo "[$(date)] $EMB × $Y PCA terminé" | tee -a "$LOG"
    done
done

echo "[$(date)] PCA terminée" | tee -a "$LOG"
