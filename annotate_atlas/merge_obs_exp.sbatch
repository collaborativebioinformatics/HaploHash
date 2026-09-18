#!/bin/bash
#SBATCH --job-name=merge_obs_exp
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/merge_obs_exp_%j.log

set -euo pipefail
REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
cd "$REPO/annotate_atlas"
PY="$REPO/.venv_avi/bin/python3"
OUT=obs_exp_by_block.tsv

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X)
first=1
> "$OUT.part"
for c in "${CHROMS[@]}"; do
  f="obs_exp/chr${c}.tsv"
  if [[ ! -s "$f" ]]; then
    echo "MISSING: $f" >&2
    exit 1
  fi
  if [[ $first -eq 1 ]]; then
    cat "$f" >> "$OUT.part"
    first=0
  else
    tail -n +2 "$f" >> "$OUT.part"
  fi
done
mv "$OUT.part" "$OUT"
echo "[$(date)] merged: $(wc -l < "$OUT") rows"

"$PY" correlate_obs_exp_eqtl.py "$OUT" "$REPO/block_scores.tsv" "$REPO/eqtl_by_block.tsv" obs_exp_vs_eqtl
echo "[$(date)] obs/exp vs eQTL correlation done"
