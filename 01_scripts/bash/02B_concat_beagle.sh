#!/bin/bash
# =============================================================================
# 02B_concat_beagle.sh — Concatenate per-chromosome Beagle files
# Also builds the combined canonical sites file
# Run AFTER all 02_genotype_likelihoods.sh array jobs complete
# Submit: sbatch 01_scripts/bash/02B_concat_beagle.sh
# =============================================================================
#SBATCH -J "02B_concat"
#SBATCH -o logs/02B_concat_%j.out
#SBATCH -e logs/02B_concat_%j.err
#SBATCH -c 2
#SBATCH -p small
#SBATCH --time=0-12:00
#SBATCH --mem=20G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

log_msg "=== CONCATENATING BEAGLE FILES ==="

# ─── 1. Combine canonical site lists ────────────────────────────────────────
log_msg "Step 1: Combining canonical site lists..."

COMBINED_SITES="${INFO_DIR}/sites_all_${SUFFIX}_canonical"
> "$COMBINED_SITES"

for NUM in $(cat "${INFO_DIR}/regions_number.txt"); do
    SITES_CHR="${INFO_DIR}/sites_by_chr/sites_${SUFFIX}_chr${NUM}_canonical"
    if [[ -f "$SITES_CHR" ]]; then
        cat "$SITES_CHR" >> "$COMBINED_SITES"
    else
        log_msg "  WARNING: Missing sites for chr${NUM}"
    fi
done

module load angsd
angsd sites index "$COMBINED_SITES"

N_SITES=$(wc -l < "$COMBINED_SITES")
log_msg "  Total canonical sites: ${N_SITES}"

# ─── 2. Create nosex version (exclude sex chromosomes) ──────────────────────
log_msg "Step 2: Creating autosomal-only sites..."

NOSEX_SITES="${INFO_DIR}/sites_all_${SUFFIX}_canonical_nosex"
grep -v -E "$(echo $SEX_CHRS | tr ' ' '|')" "$COMBINED_SITES" > "$NOSEX_SITES"
angsd sites index "$NOSEX_SITES"

N_NOSEX=$(wc -l < "$NOSEX_SITES")
log_msg "  Autosomal canonical sites: ${N_NOSEX}"

# ─── 3. Concatenate Beagle files (all chromosomes) ──────────────────────────
log_msg "Step 3: Concatenating Beagle files..."

OUT_BEAGLE="${GL_DIR}/all_${SUFFIX}.beagle.gz"
TEMP_BEAGLE="${GL_DIR}/all_${SUFFIX}.beagle"

# Expected column count from first chromosome
FIRST_CHR=$(head -1 "${INFO_DIR}/regions_number.txt")
FIRST_FILE="${GL_DIR}/by_chr/all_${SUFFIX}_chr${FIRST_CHR}.beagle.gz"
require_file "$FIRST_FILE" "First chromosome Beagle"

EXPECTED_COLS=$(zcat "$FIRST_FILE" | head -1 | awk '{print NF}')
log_msg "  Expected columns: ${EXPECTED_COLS}"

# Write header from first file
zcat "$FIRST_FILE" | head -1 > "$TEMP_BEAGLE"

# Append data from all chromosomes with column validation
TOTAL_SNPS=0
SKIPPED=0

for NUM in $(cat "${INFO_DIR}/regions_number.txt"); do
    CHR_FILE="${GL_DIR}/by_chr/all_${SUFFIX}_chr${NUM}.beagle.gz"
    if [[ -f "$CHR_FILE" ]]; then
        CHR_SNPS=$(zcat "$CHR_FILE" | tail -n +2 | \
            awk -v cols="$EXPECTED_COLS" '{
                if (NF == cols) { print; c++ }
            } END { print c+0 > "/dev/stderr" }' >> "$TEMP_BEAGLE" 2>&1 | tail -1)
        TOTAL_SNPS=$((TOTAL_SNPS + CHR_SNPS))
    else
        log_msg "  WARNING: Missing Beagle for chr${NUM}"
        SKIPPED=$((SKIPPED + 1))
    fi
done

gzip -f "$TEMP_BEAGLE"

log_msg "  Total SNPs in concatenated Beagle: ${TOTAL_SNPS}"
[[ $SKIPPED -gt 0 ]] && log_msg "  WARNING: ${SKIPPED} chromosomes missing"

# ─── 4. Validate ────────────────────────────────────────────────────────────
log_msg "Step 4: Validation..."

N_COLS=$(zcat "$OUT_BEAGLE" | head -1 | awk '{print NF}')
N_IND=$(( (N_COLS - 3) / 3 ))
N_SNPS=$(zcat "$OUT_BEAGLE" | tail -n +2 | wc -l)

log_msg "  Beagle: ${N_IND} individuals × ${N_SNPS} SNPs"

log_msg "=== CONCATENATION COMPLETE ==="
