#!/bin/bash
#SBATCH --job-name=esmc_300M
#SBATCH --output=logs/esmc_300M_%j.out
#SBATCH --error=logs/esmc_300M_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

# ESM-C 300M - full ATLAS (1390 proteins)
# Run from the project root:
#   sbatch scripts/slurm/run_esmc_300M.sh
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_DIR="$HOME/plm-project"   # adjust to the path on the cluster
FASTA="$PROJECT_DIR/Datasets/ATLAS/atlas_sequences.fasta"
OUT_DIR="$PROJECT_DIR/Datasets/embeddings"
VENV="$PROJECT_DIR/plm-api/.venv/bin/plm-embeddings"

mkdir -p "$OUT_DIR" logs/

echo "[$(date)] Starting ESM-C 300M on $(grep -c '^>' $FASTA) sequences"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"

$VENV \
  -i "$FASTA" \
  -o "$OUT_DIR" \
  -m esmc_300M \
  --toks-per-batch 2048 \
  --device cuda

echo "[$(date)] Done → $OUT_DIR/esmc_300M.safetensors"
