#!/bin/bash
# =============================================================================
# 01_ngsparalog.sh — Identify paralogous (deviant) sites with ngsParalog
# Runs per chromosome via SLURM array
# Submit: sbatch 01_scripts/bash/01_ngsparalog.sh
# =============================================================================
#SBATCH -J "01_ngsparalog"
#SBATCH -o logs/01_ngsparalog_%A_%a.out
#SBATCH -e logs/01_ngsparalog_%A_%a.err
#SBATCH -c 4
#SBATCH -p small
#SBATCH --time=1-00:00
#SBATCH --mem=60G
#SBATCH --array=1-41%10

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load angsd/0.931
module load samtools
ulimit -S -n 2048

# ─── Parameters ──────────────────────────────────────────────────────────────
NB_CPU=4
REGION_NUM="${SLURM_ARRAY_TASK_ID}"
REGION=$(sed -n "${REGION_NUM}p" "${INFO_DIR}/regions.txt")
BAMLIST="${INFO_DIR}/bam.filelist"

read -r N_IND MIN_IND MAX_DEPTH <<< "$(compute_filters "$BAMLIST")"

log_msg "=== CHROMOSOME ${REGION_NUM} : ${REGION} ==="
log_msg "N=${N_IND}, MIN_IND=${MIN_IND}, MAX_DEPTH=${MAX_DEPTH}"

# ─── Step 1: Call SNPs with ANGSD (MAF + quality filters) ───────────────────
log_msg "Step 1: ANGSD SNP calling..."

angsd -P $NB_CPU -nQueueSize 20 \
  -doMaf 1 -GL 2 -doMajorMinor 1 -doCounts 1 \
  -remove_bads 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
  -uniqueOnly 1 -only_proper_pairs 1 \
  -minInd "$MIN_IND" -minMaf "$MIN_MAF" \
  -setMaxDepth "$MAX_DEPTH" -setMinDepthInd "$MIN_DEPTH" \
  -b "$BAMLIST" -r "$REGION" \
  -out "${NGSPARALOG_DIR}/all_${SUFFIX}_chr${REGION_NUM}"

# ─── Step 2: Extract site positions ─────────────────────────────────────────
log_msg "Step 2: Extracting site positions..."

gunzip -f "${NGSPARALOG_DIR}/all_${SUFFIX}_chr${REGION_NUM}.mafs.gz"

MAFS="${NGSPARALOG_DIR}/all_${SUFFIX}_chr${REGION_NUM}.mafs"
SITES="${INFO_DIR}/sites_by_chr/sites_${SUFFIX}_chr${REGION_NUM}"
BED="${INFO_DIR}/sites_by_chr/sites_${SUFFIX}_chr${REGION_NUM}.bed"

"$RSCRIPT" 01_scripts/R/make_sites_list.R "$MAFS" "$SITES"
angsd sites index "$SITES"
awk '{print $1"\t"$2-1"\t"$2}' "$SITES" > "$BED"

# ─── Step 3: Run ngsParalog ─────────────────────────────────────────────────
log_msg "Step 3: Running ngsParalog..."

samtools mpileup -b "$BAMLIST" -l "$BED" -r "$REGION" -q 0 -Q 0 --ff UNMAP,DUP | \
  "${NGSPARALOG_PATH}/ngsParalog" calcLR \
    -infile - \
    -outfile "${NGSPARALOG_DIR}/all_${SUFFIX}_chr${REGION_NUM}.ngsparalog" \
    -minQ 20 -minind "$MIN_IND" -mincov "$MIN_DEPTH" -allow_overwrite 1

# ─── Step 4: Separate canonical / deviant sites ─────────────────────────────
log_msg "Step 4: Classifying canonical vs deviant..."

"$RSCRIPT" 01_scripts/R/classify_paralog.R \
  "${NGSPARALOG_DIR}/all_${SUFFIX}_chr${REGION_NUM}.ngsparalog" \
  "$SITES" "$PVAL_THRESHOLD"

angsd sites index "${SITES}_canonical"
angsd sites index "${SITES}_deviant"

N_CANON=$(wc -l < "${SITES}_canonical")
N_DEV=$(wc -l < "${SITES}_deviant")
log_msg "  Canonical: ${N_CANON}, Deviant: ${N_DEV}"

log_msg "=== CHROMOSOME ${REGION_NUM} COMPLETE ==="
