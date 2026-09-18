#!/bin/bash
#SBATCH --job-name=haplo-gtex-download
#SBATCH --account=standard
#SBATCH --partition=medium
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=annotate_atlas/gtex/logs/download-%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
export UV_CACHE_DIR=/tmp/haplohash-uv-cache
uv run --project annotate_atlas --frozen --offline python annotate_atlas/gtex_download.py "${SLURM_ARRAY_TASK_ID:-common}" --connections "${GTEX_CONNECTIONS:-1}"
