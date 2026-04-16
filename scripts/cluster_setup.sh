#!/bin/bash
# ── Setup PLM-API sur le cluster SFBI ────────────────────────────────────
# À lancer UNE SEULE FOIS depuis votre $HOME sur le cluster.
#
# Pré-requis cluster :
#   - Python 3.11 disponible (module load python/3.11 ou conda)
#   - Git disponible
#   - Accès internet (pour pip + HuggingFace)
#
# Usage :
#   ssh user@sfbi-cluster
#   bash cluster_setup.sh
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_DIR="$HOME/plm-project"
echo "=== Setup PLM project dans $PROJECT_DIR ==="

# 1. Créer structure de dossiers
mkdir -p "$PROJECT_DIR/Datasets/ATLAS"
mkdir -p "$PROJECT_DIR/Datasets/embeddings"
mkdir -p "$PROJECT_DIR/logs"

# 2. Copier le repo (ou cloner depuis GitHub une fois public)
#    Option A — copier depuis votre machine locale :
#    rsync -av --exclude='.git' --exclude='.venv' --exclude='Datasets/embeddings' \
#        /local/path/to/Code/ user@sfbi-cluster:$PROJECT_DIR/
#
#    Option B — cloner directement si le repo est sur GitLab/GitHub :
#    git clone https://gitlab.com/... "$PROJECT_DIR"

echo ""
echo ">>> Étape manuelle : copier les fichiers via rsync ou git clone"
echo ">>> Commande rsync depuis votre Mac :"
echo "    rsync -av --exclude='.venv' --exclude='Datasets/embeddings' \\"
echo "        '/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code/' \\"
echo "        user@sfbi-cluster:$PROJECT_DIR/"
echo ""

# 3. Vérifier Python
echo "=== Python version ==="
python3 --version || { echo "ERREUR: python3 non disponible — charger le module"; exit 1; }

# 4. Créer venv pour PLM-API
echo "=== Création venv PLM-API ==="
python3 -m venv "$PROJECT_DIR/plm-api/.venv"
source "$PROJECT_DIR/plm-api/.venv/bin/activate"

pip install --upgrade pip
pip install -e "$PROJECT_DIR/plm-api/"

echo "=== PLM-API installé ==="
plm-embeddings --list-models

# 5. Vérifier que CUDA est disponible
python3 -c "import torch; print('CUDA:', torch.cuda.is_available(), '| GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'N/A')"

echo ""
echo "=== Setup terminé ==="
echo "Pour lancer les jobs :"
echo "  cd $PROJECT_DIR && bash scripts/slurm/run_all.sh"
echo ""
echo "⚠️  Vérifier PROJECT_DIR dans les scripts SLURM avant de soumettre :"
echo "  grep PROJECT_DIR scripts/slurm/run_esmc_300M.sh"
