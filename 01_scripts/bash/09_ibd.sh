#!/bin/bash
# =============================================================================
# 09_ibd.sh — Isolation by Distance (Mantel test, FST/(1-FST) vs hydro distance)
#
# PREREQUISITE: run 01_scripts/R/compute_hydro_distances.R first (separate
# conda env with 'marmap' — see that script's header) to produce
# ${IBD_DIR}/hydro_distances.csv, which ibd_analysis.R reads.
#
# Submit: sbatch 01_scripts/bash/09_ibd.sh
# =============================================================================
#SBATCH -J "09_ibd"
#SBATCH -o logs/09_ibd_%j.out
#SBATCH -e logs/09_ibd_%j.err
#SBATCH -c 4
#SBATCH -p small
#SBATCH --time=0-12:00
#SBATCH --mem=20G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

log_msg "=== ISOLATION BY DISTANCE ==="

mkdir -p "${IBD_DIR}"

"$RSCRIPT" 01_scripts/R/ibd_analysis.R \
    "${FST_DIR}/fst_results.tsv" \
    "${INFO_DIR}" \
    "${IBD_DIR}" \
    "${FIG_DIR}/main" \
    "${TABLE_DIR}/main"

log_msg "=== IBD COMPLETE ==="
