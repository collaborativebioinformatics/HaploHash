#!/bin/bash
#SBATCH --job-name=setup_venv
#SBATCH --account=iu_0142
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs/logs/setup_venv_%j.log

set -euo pipefail
REPO=/dcai/projects/iu_0142/team11/Pathway_analysis_into_graphs
"$REPO/.venv_avi/bin/pip" install --no-input matplotlib scipy
echo "[$(date)] venv ready"
