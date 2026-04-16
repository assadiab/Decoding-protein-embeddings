#!/bin/bash
# ── Lance tous les jobs embeddings en séquence (dépendances SLURM) ───────
# Utiliser si le cluster a une seule GPU disponible.
# Sinon supprimer --dependency pour lancer en parallèle.
#
# Usage (depuis la racine du projet sur le cluster) :
#   bash scripts/slurm/run_all.sh
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "=== Soumission des jobs SLURM ==="

JOB1=$(sbatch --parsable scripts/slurm/run_esmc_300M.sh)
echo "ESM-C 300M → job $JOB1"

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 scripts/slurm/run_esmc_600M.sh)
echo "ESM-C 600M → job $JOB2 (après $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 scripts/slurm/run_ankh2.sh)
echo "Ankh2-Large → job $JOB3 (après $JOB2)"

echo ""
echo "Suivi : squeue -u \$USER"
echo "Logs  : tail -f logs/esmc_300M_\${JOB1}.out"
