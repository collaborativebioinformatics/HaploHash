#!/bin/bash
#SBATCH --job-name=haplo-gtex-analysis
#SBATCH --account=standard
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=annotate_atlas/gtex/logs/analysis-%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_CACHE_DIR=/tmp/haplohash-uv-cache
export MPLCONFIGDIR=/tmp/haplohash-mpl
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
uv run --project annotate_atlas --frozen --offline python annotate_atlas/gtex_eqtl.py "${1:-tissue}" --index "${SLURM_ARRAY_TASK_ID:-0}"
