#!/bin/bash
# =============================================================================
# 02_genotype_likelihoods.sh — Compute GL and MAF on canonical SNPs
# Runs per chromosome via SLURM array, then concatenates
# Submit: sbatch 01_scripts/bash/02_genotype_likelihoods.sh
# =============================================================================
#SBATCH -J "02_gl"
#SBATCH -o logs/02_gl_%A_%a.out
#SBATCH -e logs/02_gl_%A_%a.err
#SBATCH -c 4
#SBATCH -p small
#SBATCH --time=1-00:00
#SBATCH --mem=60G
#SBATCH --array=1-41%20

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load angsd
ulimit -S -n 2048

NB_CPU=4
REGION_NUM="${SLURM_ARRAY_TASK_ID}"
REGION=$(sed -n "${REGION_NUM}p" "${INFO_DIR}/regions.txt")
BAMLIST="${INFO_DIR}/bam.filelist"
SITES="${INFO_DIR}/sites_by_chr/sites_${SUFFIX}_chr${REGION_NUM}_canonical"

read -r N_IND MIN_IND MAX_DEPTH <<< "$(compute_filters "$BAMLIST")"

require_file "$SITES" "Canonical sites for chr${REGION_NUM}"

log_msg "=== GL — CHROMOSOME ${REGION_NUM} ==="

angsd -P $NB_CPU -nQueueSize 50 \
  -doMaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 \
  -anc "${INFO_DIR}/genome.fasta" \
  -sites "$SITES" \
  -remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
  -uniqueOnly 1 -only_proper_pairs 1 \
  -r "$REGION" \
  -b "$BAMLIST" \
  -out "${GL_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}"

log_msg "=== GL — CHROMOSOME ${REGION_NUM} COMPLETE ==="
