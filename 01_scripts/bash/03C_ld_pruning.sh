#!/bin/bash
# =============================================================================
# 03C_ld_pruning.sh — LD prune per chromosome with ngsLD
# Submit: cat regions_number.txt | parallel -j20 \
#   srun -p medium -c 1 --mem=20G --time=7-00:00 \
#   01_scripts/bash/03C_ld_pruning.sh {}
# =============================================================================

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

REGION_NUM="$1"
CHR_NUM=$(printf "%02d" $((10#$REGION_NUM)))

INPUT_LD="${LD_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}.ld"
INPUT_LD_GZ="${INPUT_LD}.gz"
BEAGLE_IN="${GL_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}.beagle.gz"
OUTPUT_PRUNED="${LD_DIR}/by_chr/chr${CHR_NUM}_${SUFFIX}.pruned"
BEAGLE_OUT="${LD_DIR}/beagle_pruned/chr${CHR_NUM}_${SUFFIX}.pruned.beagle.gz"

# Decompress LD if needed
if [[ -f "$INPUT_LD_GZ" ]] && [[ ! -f "$INPUT_LD" ]]; then
    gunzip -k "$INPUT_LD_GZ"
fi

require_file "$INPUT_LD" "LD file chr${REGION_NUM}"
require_file "$BEAGLE_IN" "Beagle chr${REGION_NUM}"

log_msg "=== LD PRUNING — chr${CHR_NUM} ==="

# Check data lines
DATA_LINES=$(tail -n +2 "$INPUT_LD" | wc -l)
if [[ "$DATA_LINES" -eq 0 ]]; then
    log_msg "  No LD data for chr${CHR_NUM}. Creating empty outputs."
    touch "$OUTPUT_PRUNED"
    zcat "$BEAGLE_IN" | head -1 | gzip > "$BEAGLE_OUT"
    exit 0
fi

# Clean LD file (remove header + NaN)
TEMP_LD="${LD_DIR}/by_chr/temp_chr${CHR_NUM}.ld"
tail -n +2 "$INPUT_LD" | grep -vi nan > "$TEMP_LD"

CLEAN_LINES=$(wc -l < "$TEMP_LD")
log_msg "  Clean LD pairs: ${CLEAN_LINES}"

# Python LD pruning
python "$PRUNE_SCRIPT" \
    --input "$TEMP_LD" \
    --max_dist "$LD_PRUNE_DIST" \
    --min_weight "$LD_PRUNE_WEIGHT" \
    --field_dist 3 \
    --field_weight 7 \
    --weight_precision 4 \
    --output "$OUTPUT_PRUNED"

rm -f "$TEMP_LD"

# Reformat output
sed -i 's/:/\t/g' "$OUTPUT_PRUNED"

PRUNED_N=$(wc -l < "$OUTPUT_PRUNED")
log_msg "  Pruned SNPs retained: ${PRUNED_N}"

# Extract pruned SNPs from Beagle
python 01_scripts/utils/beagle_extract_snps.py \
    "$BEAGLE_IN" "$OUTPUT_PRUNED" "$BEAGLE_OUT"

# Clean up decompressed LD
rm -f "$INPUT_LD"

log_msg "=== LD PRUNING chr${CHR_NUM} COMPLETE ==="
