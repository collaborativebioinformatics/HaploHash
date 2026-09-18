#!/bin/bash
#SBATCH --job-name=obs_exp
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:30:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/obs_exp_%A_%a.log

set -euo pipefail
source /etc/profile.d/z00_lmod.sh
module load GCC/13.3.0 BCFtools/1.21

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X)
CHROM=chr${CHROMS[$SLURM_ARRAY_TASK_ID]}

REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
BIN=$REPO/annotate_atlas/bin
VCF_DIR=/dcai/projects/iu_0142/team11/data/1000G
AVI=/dcai/projects/iu_0142/team11/data/avi/alphagenome_variant_impact_score_snvs.tsv.gz
OUT=$REPO/annotate_atlas/obs_exp/${CHROM}.tsv

mkdir -p "$REPO/annotate_atlas/obs_exp"

echo "[$(date)] ${CHROM}: computing obs/exp on $(hostname)"
"$REPO/.venv_avi/bin/python3" "$REPO/annotate_atlas/genome_obs_exp.py" \
  "$REPO/annotate_atlas/blocks.bed" "$VCF_DIR" "$AVI" "$OUT" \
  --tabix "$BIN/tabix" --bcftools bcftools --chrom "$CHROM"
echo "[$(date)] ${CHROM}: done"
