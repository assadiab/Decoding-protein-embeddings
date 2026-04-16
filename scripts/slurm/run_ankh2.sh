#!/bin/bash
#SBATCH --job-name=ankh2_large
#SBATCH --output=logs/ankh2_large_%j.out
#SBATCH --error=logs/ankh2_large_%j.err
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4

# ── Ankh2-Large (ext2) — ATLAS complet (1390 protéines) ──────────────────
# Modèle : ElnaggarLab/ankh2-ext2 — 1536 dims
# ⚠️ C'est la variante ext2, pas le ankh2-large standard — à documenter
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_DIR="$HOME/plm-project"
FASTA="$PROJECT_DIR/Datasets/ATLAS/atlas_sequences.fasta"
OUT_DIR="$PROJECT_DIR/Datasets/embeddings"
VENV="$PROJECT_DIR/plm-api/.venv/bin/plm-embeddings"

mkdir -p "$OUT_DIR" logs/

echo "[$(date)] Starting Ankh2-Large on $(grep -c '^>' $FASTA) sequences"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"

$VENV \
  -i "$FASTA" \
  -o "$OUT_DIR" \
  -m ankh2_large \
  --toks-per-batch 2048 \
  --device cuda

echo "[$(date)] Done → $OUT_DIR/ankh2_large.safetensors"
