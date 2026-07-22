#!/bin/bash
# =============================================================================
# 08_mask_deviants.sh -- Mask deviant (paralogous) regions in the genome
# Creates a masked reference genome for demographic inference
# Submit: sbatch 01_scripts/bash/08_mask_deviants.sh
# =============================================================================
#SBATCH -J "08_mask"
#SBATCH -o logs/08_mask_%j.out
#SBATCH -e logs/08_mask_%j.err
#SBATCH -c 1
#SBATCH -p small
#SBATCH --time=1-00:00
#SBATCH --mem=10G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load bedtools
module load samtools

MASK_LENGTH=150

log_msg "=== DEVIANT MASKING (${MASK_LENGTH}bp windows) ==="

require_file "${INFO_DIR}/genome.fasta"
require_file "${INFO_DIR}/genome.fasta.fai"

mkdir -p "${INFO_DIR}/mask_by_chr"

# --- Step 1: R script to detect deviants and create BED masks ----------------
"$RSCRIPT" 01_scripts/R/deviant_masking.R "$MASK_LENGTH" "${NGSPARALOG_DIR}" "${INFO_DIR}"

# --- Step 2: Combine all BED files ------------------------------------------
COMBINED_BED="${INFO_DIR}/mask_deviants_combined.bed"
cat "${INFO_DIR}/mask_by_chr/"mask_deviant_chr*.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge > "$COMBINED_BED"

TOTAL_MASKED=$(awk '{sum += $3 - $2} END {print sum}' "$COMBINED_BED")
GENOME_SIZE=$(awk '{sum += $2} END {print sum}' "${INFO_DIR}/genome.fasta.fai")
PCT=$(echo "scale=2; $TOTAL_MASKED * 100 / $GENOME_SIZE" | bc)

log_msg "  Masked: ${TOTAL_MASKED} bp (${PCT}% of genome)"

# --- Step 3: Mask the genome FASTA ------------------------------------------
log_msg "Masking genome FASTA..."

MASKED_GENOME="${INFO_DIR}/genome_masked_deviants.fasta"

bedtools maskfasta \
    -fi "${INFO_DIR}/genome.fasta" \
    -bed "$COMBINED_BED" \
    -fo "$MASKED_GENOME"

samtools faidx "$MASKED_GENOME"

log_msg "  Masked genome: ${MASKED_GENOME}"

log_msg "=== MASKING COMPLETE ==="
