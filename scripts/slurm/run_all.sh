#!/bin/bash
# Launch all embedding jobs in sequence (SLURM dependencies)
# Use if the cluster has a single available GPU.
# Otherwise remove --dependency to run in parallel.
#
# Usage (from the project root on the cluster):
#   bash scripts/slurm/run_all.sh
# ─────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "=== Soumission des jobs SLURM ==="

JOB1=$(sbatch --parsable scripts/slurm/run_esmc_300M.sh)
echo "ESM-C 300M → job $JOB1"

JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 scripts/slurm/run_esmc_600M.sh)
echo "ESM-C 600M -> job $JOB2 (after $JOB1)"

JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 scripts/slurm/run_ankh2.sh)
echo "Ankh2-Large -> job $JOB3 (after $JOB2)"

echo ""
echo "Suivi : squeue -u \$USER"
echo "Logs  : tail -f logs/esmc_300M_\${JOB1}.out"
