#!/bin/bash
set -e

ROOT="/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code"
LOG="$ROOT/Datasets/embeddings/full/run_$(date +%Y%m%d_%H%M).log"

source "$ROOT/plm-api/.venv/bin/activate"

echo "[$(date)] Starting 2G embedding generation - 1390 proteins" | tee -a "$LOG"

for MODEL in esmc_300M esmc_600M ankh2_large; do
    echo "[$(date)] === $MODEL ===" | tee -a "$LOG"
    plm-embeddings \
        -i "$ROOT/Datasets/ATLAS/atlas_sequences.fasta" \
        -o "$ROOT/Datasets/embeddings/full/" \
        -m "$MODEL" \
        --toks-per-batch 512 2>&1 | tee -a "$LOG"
    echo "[$(date)] $MODEL done" | tee -a "$LOG"
done

echo "[$(date)] All embeddings generated" | tee -a "$LOG"
