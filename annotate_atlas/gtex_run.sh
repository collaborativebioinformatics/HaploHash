#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
gtex_storage=$(realpath -m "${1:?Usage: bash annotate_atlas/gtex_run.sh SHARED_STORAGE_DIRECTORY}")
mkdir -p "$gtex_storage/data" "$gtex_storage/work"
for gtex_subdir in data work; do
    gtex_link="annotate_atlas/gtex/$gtex_subdir"
    if [[ -e "$gtex_link" || -L "$gtex_link" ]]; then
        [[ $(realpath "$gtex_link") == "$gtex_storage/$gtex_subdir" ]] || { echo "Storage mismatch: $gtex_link" >&2; exit 1; }
    else
        ln -s "$gtex_storage/$gtex_subdir" "$gtex_link"
    fi
done
mkdir -p annotate_atlas/gtex/logs
gtex_connections=${GTEX_CONNECTIONS:-1}
[[ $gtex_connections =~ ^([1-9]|1[0-6])$ ]] || { echo "GTEX_CONNECTIONS must be 1–16" >&2; exit 1; }
gtex_parallel=$((20 / gtex_connections))
uv sync --project annotate_atlas
gtex_last_index=$(($(wc -l < annotate_atlas/gtex/tissues.tsv) - 2))
gtex_common_job=$(sbatch --parsable annotate_atlas/gtex_download.sbatch)
gtex_download_job=$(sbatch --parsable --array="0-${gtex_last_index}%${gtex_parallel}" annotate_atlas/gtex_download.sbatch)
gtex_prepare_job=$(sbatch --parsable --dependency="afterok:${gtex_common_job}" annotate_atlas/gtex_analyze.sbatch prepare)
gtex_analysis_job=$(sbatch --parsable --array="0-${gtex_last_index}%8" --dependency="aftercorr:${gtex_download_job},afterok:${gtex_prepare_job}" annotate_atlas/gtex_analyze.sbatch)
gtex_combine_job=$(sbatch --parsable --dependency="afterok:${gtex_analysis_job}" annotate_atlas/gtex_analyze.sbatch combine)
echo "GTEx jobs: downloads ${gtex_common_job},${gtex_download_job}; preparation ${gtex_prepare_job}; analysis ${gtex_analysis_job}; results ${gtex_combine_job}"
