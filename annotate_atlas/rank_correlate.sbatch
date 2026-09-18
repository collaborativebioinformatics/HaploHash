#!/bin/bash
#SBATCH --job-name=rank_correlate
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/rank_correlate_%j.log

set -euo pipefail
REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
cd "$REPO/annotate_atlas"
PY="$REPO/.venv_avi/bin/python3"

"$PY" rank_blocks.py "$REPO/block_scores.tsv" block_avi_rank.tsv
"$PY" correlate_avi_eqtl.py block_avi_rank.tsv "$REPO/eqtl_by_block.tsv" avi_vs_eqtl
echo "[$(date)] rank + correlate done"
