#!/bin/bash
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/Datasets/embeddings/full/run_1g_$(date +%Y%m%d_%H%M).log"

source "$ROOT/plm-api/.venv/bin/activate"

echo "[$(date)] Démarrage génération embeddings 1G" | tee -a "$LOG"

for MODEL in esm2_t33_650M_UR50D ankh_large prot_t5_xl_uniref50; do
    echo "[$(date)] === $MODEL ===" | tee -a "$LOG"
    plm-embeddings \
        -i "$ROOT/Datasets/ATLAS/atlas_sequences.fasta" \
        -o "$ROOT/Datasets/embeddings/full/" \
        -m "$MODEL" \
        --toks-per-batch 512 2>&1 | tee -a "$LOG"
    echo "[$(date)] $MODEL terminé" | tee -a "$LOG"
done

echo "[$(date)] Tous les embeddings 1G générés" | tee -a "$LOG"
