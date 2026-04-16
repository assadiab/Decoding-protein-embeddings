#!/bin/bash
#SBATCH --job-name=esmc_600M
#SBATCH --output=logs/esmc_600M_%j.out
#SBATCH --error=logs/esmc_600M_%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4

# ── ESM-C 600M — ATLAS complet (1390 protéines) ──────────────────────────
# Modèle plus lourd : 1152 dims, plus de mémoire GPU requise.
# Réduire --toks-per-batch si OOM (essayer 1024 puis 512).
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_DIR="$HOME/plm-project"
FASTA="$PROJECT_DIR/Datasets/ATLAS/atlas_sequences.fasta"
OUT_DIR="$PROJECT_DIR/Datasets/embeddings"
VENV="$PROJECT_DIR/plm-api/.venv/bin/plm-embeddings"

mkdir -p "$OUT_DIR" logs/

echo "[$(date)] Starting ESM-C 600M on $(grep -c '^>' $FASTA) sequences"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"

$VENV \
  -i "$FASTA" \
  -o "$OUT_DIR" \
  -m esmc_600M \
  --toks-per-batch 2048 \
  --device cuda

echo "[$(date)] Done → $OUT_DIR/esmc_600M.safetensors"
