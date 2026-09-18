#!/bin/bash
#SBATCH --job-name=dl_1000g
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/dl_1000g_%A_%a.log

set -euo pipefail

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 X)
CHROM=${CHROMS[$SLURM_ARRAY_TASK_ID]}

BASE="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
NAME="ALL.chr${CHROM}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
OUT_DIR=/dcai/projects/iu_0142/team11/data/1000G

mkdir -p "$OUT_DIR"

if [[ -s "$OUT_DIR/$NAME" && -s "$OUT_DIR/$NAME.tbi" ]]; then
  echo "[$(date)] chr${CHROM}: already present, skipping"
  exit 0
fi

echo "[$(date)] chr${CHROM}: downloading $NAME"
wget -q -O "$OUT_DIR/$NAME.part" "$BASE/$NAME"
mv "$OUT_DIR/$NAME.part" "$OUT_DIR/$NAME"
wget -q -O "$OUT_DIR/$NAME.tbi.part" "$BASE/$NAME.tbi"
mv "$OUT_DIR/$NAME.tbi.part" "$OUT_DIR/$NAME.tbi"
echo "[$(date)] chr${CHROM}: done, $(du -h "$OUT_DIR/$NAME" | cut -f1)"
