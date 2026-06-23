#!/bin/bash
# Setup PLM-API on the SFBI cluster
# Run ONCE from your $HOME on the cluster.
#
# Cluster prerequisites:
#   - Python 3.11 available (module load python/3.11 or conda)
#   - Git available
#   - Internet access (for pip + HuggingFace)
#
# Usage :
#   ssh user@sfbi-cluster
#   bash cluster_setup.sh
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_DIR="$HOME/plm-project"
echo "=== Setup PLM project dans $PROJECT_DIR ==="

# 1. Create directory structure
mkdir -p "$PROJECT_DIR/Datasets/ATLAS"
mkdir -p "$PROJECT_DIR/Datasets/embeddings"
mkdir -p "$PROJECT_DIR/logs"

# 2. Copy the repo (or clone from GitHub once public)
#    Option A - copy from your local machine:
#    rsync -av --exclude='.git' --exclude='.venv' --exclude='Datasets/embeddings' \
#        /local/path/to/Code/ user@sfbi-cluster:$PROJECT_DIR/
#
#    Option B - clone directly if the repo is on GitLab/GitHub:
#    git clone https://gitlab.com/... "$PROJECT_DIR"

echo ""
echo ">>> Manual step: copy files via rsync or git clone"
echo ">>> rsync command from your Mac:"
echo "    rsync -av --exclude='.venv' --exclude='Datasets/embeddings' \\"
echo "        '/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code/' \\"
echo "        user@sfbi-cluster:$PROJECT_DIR/"
echo ""

# 3. Check Python
echo "=== Python version ==="
python3 --version || { echo "ERROR: python3 not available - load the module"; exit 1; }

# 4. Create venv for PLM-API
echo "=== Creating PLM-API venv ==="
python3 -m venv "$PROJECT_DIR/plm-api/.venv"
source "$PROJECT_DIR/plm-api/.venv/bin/activate"

pip install --upgrade pip
pip install -e "$PROJECT_DIR/plm-api/"

echo "=== PLM-API installed ==="
plm-embeddings --list-models

# 5. Check that CUDA is available
python3 -c "import torch; print('CUDA:', torch.cuda.is_available(), '| GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

echo ""
echo "=== Setup done ==="
echo "Pour lancer les jobs :"
echo "  cd $PROJECT_DIR && bash scripts/slurm/run_all.sh"
echo ""
echo "Check PROJECT_DIR in the SLURM scripts before submitting:"
echo "  grep PROJECT_DIR scripts/slurm/run_esmc_300M.sh"
