#!/bin/bash
#SBATCH --job-name=gene_density
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/gene_density_%j.log

set -euo pipefail
REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
cd "$REPO/annotate_atlas"
PY="$REPO/.venv_avi/bin/python3"

GTF=/dcai/projects/iu_0142/team11/data/gencode/gencode.v46.annotation.gtf.gz
mkdir -p /dcai/projects/iu_0142/team11/data/gencode
if [[ ! -s "$GTF" ]]; then
  echo "[$(date)] downloading GENCODE v46 GTF"
  wget -q -O "$GTF.part" "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"
  mv "$GTF.part" "$GTF"
fi

echo "[$(date)] counting genes per block"
"$PY" gene_density.py blocks.bed "$GTF" gene_density_by_block.tsv

echo "[$(date)] correlating against AVI + eQTL"
"$PY" correlate_gene_density.py gene_density_by_block.tsv block_avi_rank.tsv "$REPO/eqtl_by_block.tsv" gene_density_vs_avi_eqtl

echo "[$(date)] done"
