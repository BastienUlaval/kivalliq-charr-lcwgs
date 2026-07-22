#!/bin/bash
# =============================================================================
# 03D_concat_pruned.sh -- Concatenate pruned Beagles; create nosex + regional subsets
# Run AFTER all 03C array jobs complete
# Submit: sbatch 01_scripts/bash/03D_concat_pruned.sh
# =============================================================================
#SBATCH -J "03D_concat_pruned"
#SBATCH -o logs/03D_concat_%j.out
#SBATCH -e logs/03D_concat_%j.err
#SBATCH -c 2
#SBATCH -p small
#SBATCH --time=0-12:00
#SBATCH --mem=30G

set -x
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

log_msg "=== CONCATENATE PRUNED BEAGLE FILES ==="

BAMLIST="${INFO_DIR}/bam.filelist"
TEMP="${LD_DIR}/all_${SUFFIX}.pruned.beagle"
OUT_ALL="${LD_DIR}/all_${SUFFIX}.pruned.beagle.gz"

# --- 1. Concatenate all chromosomes -----------------------------------------
FIRST=$(ls "${LD_DIR}/beagle_pruned/"*.pruned.beagle.gz 2>/dev/null | sort | head -1)
require_file "$FIRST" "First pruned beagle"

EXPECTED_COLS=$(zcat "$FIRST" | head -1 | awk '{print NF}')
zcat "$FIRST" | head -1 > "$TEMP"

TOTAL=0
for f in $(ls "${LD_DIR}/beagle_pruned/"*.pruned.beagle.gz | sort); do
    N=$(zcat "$f" | tail -n +2 | awk -v c="$EXPECTED_COLS" 'NF==c' | tee -a "$TEMP" | wc -l)
    TOTAL=$((TOTAL + N))
done

gzip -f "$TEMP"
log_msg "  All chromosomes: ${TOTAL} SNPs"

# --- 2. Create nosex version (exclude sex-linked chromosomes) ----------------
log_msg "Filtering sex chromosomes..."

OUT_NOSEX="${LD_DIR}/all_${SUFFIX}.prunednosex.beagle.gz"

# Build grep pattern from sex chromosome accessions
SEX_PATTERN=$(echo $SEX_CHRS | tr ' ' '|')

zcat "$OUT_ALL" | head -1 | gzip > "$OUT_NOSEX"
zcat "$OUT_ALL" | tail -n +2 | grep -vE "$SEX_PATTERN" | gzip >> "$OUT_NOSEX"

N_NOSEX=$(zcat "$OUT_NOSEX" | tail -n +2 | wc -l)
log_msg "  Autosomal pruned SNPs: ${N_NOSEX}"

# --- 3. Create regional subsets ---------------------------------------------
log_msg "Creating regional Beagle subsets..."

# Rankin Inlet
python 01_scripts/utils/filter_beagle_samples.py \
    "$BAMLIST" \
    "${INFO_DIR}/bamlists/rankin.bamlist" \
    "$OUT_NOSEX"
mv "${OUT_NOSEX%.beagle.gz}.subset.beagle.gz" \
   "${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz"

N_RANKIN=$(zcat "${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
log_msg "  Rankin: ${N_RANKIN} individuals"

# Naujaat
python 01_scripts/utils/filter_beagle_samples.py \
    "$BAMLIST" \
    "${INFO_DIR}/bamlists/naujaat.bamlist" \
    "$OUT_NOSEX"
mv "${OUT_NOSEX%.beagle.gz}.subset.beagle.gz" \
   "${LD_DIR}/naujaat_${SUFFIX}.prunednosex.beagle.gz"

N_NAUJAAT=$(zcat "${LD_DIR}/naujaat_${SUFFIX}.prunednosex.beagle.gz" | head -1 | awk '{print (NF-3)/3}')
log_msg "  Naujaat: ${N_NAUJAAT} individuals"

# --- 4. Create LD-pruned sites list (for downstream ANGSD calls) ------------
log_msg "Extracting pruned site positions..."

PRUNED_SITES="${INFO_DIR}/sites_ldpruned_${SUFFIX}_nosex"
zcat "$OUT_NOSEX" | tail -n +2 | awk -F'\t' '{split($1, a, "_"); print a[1]"\t"a[2]}' \
    > "$PRUNED_SITES"

module load angsd
angsd sites index "$PRUNED_SITES"

N_PRUNED_SITES=$(wc -l < "$PRUNED_SITES")
log_msg "  LD-pruned autosomal sites: ${N_PRUNED_SITES}"

log_msg "=== CONCATENATION COMPLETE ==="
echo ""
echo "Key outputs:"
echo "  All pruned:        ${OUT_ALL}"
echo "  Autosomal pruned:  ${OUT_NOSEX}"
echo "  Rankin subset:     ${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz"
echo "  Naujaat subset:    ${LD_DIR}/naujaat_${SUFFIX}.prunednosex.beagle.gz"
echo "  Sites for ANGSD:   ${PRUNED_SITES}"
