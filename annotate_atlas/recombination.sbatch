#!/bin/bash
#SBATCH --job-name=recombination
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/recombination_%j.log

set -euo pipefail
REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
cd "$REPO/annotate_atlas"
PY="$REPO/.venv_avi/bin/python3"

MAP=/dcai/projects/iu_0142/team11/data/genetic_map/genetic_map_hg38_withX.txt.gz
mkdir -p /dcai/projects/iu_0142/team11/data/genetic_map
if [[ ! -s "$MAP" ]]; then
  echo "[$(date)] downloading hg38 genetic map"
  wget -q -O "$MAP.part" "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
  mv "$MAP.part" "$MAP"
fi

echo "[$(date)] computing per-block recombination rate"
"$PY" recomb_rate.py blocks.bed "$MAP" recomb_rate_by_block.tsv

echo "[$(date)] correlating against AVI, eQTL, diversity"
"$PY" correlate_recombination.py recomb_rate_by_block.tsv block_avi_rank.tsv "$REPO/eqtl_by_block.tsv" \
  "$REPO/intersect-constraint-AVI/block_constraint_avi_joined.tsv" recombination_vs_avi_eqtl_diversity

echo "[$(date)] done"
