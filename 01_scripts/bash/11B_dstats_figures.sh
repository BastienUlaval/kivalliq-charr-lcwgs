#!/bin/bash
# =============================================================================
# 11B_dstats_figures.sh — D-stats publication figure + table (Fig. 5, Table 3)
# Run AFTER all 11_dstats.sh array tasks complete.
#
#   plot_dstats.R — forest/bar plots (observed vs transversion-only) +
#   Table_Dstats_JAY.tsv
#
# Submit: sbatch 01_scripts/bash/11B_dstats_figures.sh
# =============================================================================
#SBATCH -J "11B_dstats_figs"
#SBATCH -o logs/11B_dstats_figs_%j.out
#SBATCH -e logs/11B_dstats_figs_%j.err
#SBATCH -c 2
#SBATCH -p small
#SBATCH --time=0-06:00
#SBATCH --mem=20G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

DSTAT_DIR="${BASE_DIR:-$(pwd)}/14_dstats"
FIG_DIR="${FIG_DIR:-99_figures}/main"

log_msg "=== D-STATS FIGURES (H1=JAY) ==="

"$RSCRIPT" 01_scripts/R/plot_dstats.R \
    "$DSTAT_DIR" "$FIG_DIR" "JAY"

log_msg "=== D-STATS FIGURES COMPLETE ==="
