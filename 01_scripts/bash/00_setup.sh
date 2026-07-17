#!/bin/bash
# =============================================================================
# 00_setup.sh — Initialize angsd_pipeline2 directory structure
# Run ONCE before launching any analysis
# Usage: bash 01_scripts/bash/00_setup.sh
# =============================================================================

set -euo pipefail
source config/00_config.sh

log_msg "=== SETTING UP ANGSD PIPELINE 2 ==="

# ─────────────────────────────────────────────────────────────────────────────
# 1. Create directory tree
# ─────────────────────────────────────────────────────────────────────────────
log_msg "Creating directory structure..."

dirs=(
    "$INFO_DIR" "${INFO_DIR}/sites_by_chr" "${INFO_DIR}/bamlists"
    "$NGSPARALOG_DIR"
    "$GL_DIR" "${GL_DIR}/by_chr"
    "$LD_DIR" "${LD_DIR}/by_chr" "${LD_DIR}/ld_decay" "${LD_DIR}/beagle_pruned"
    "$PCA_DIR"
    "$ADMIX_DIR" "${ADMIX_DIR}/global" "${ADMIX_DIR}/naujaat" "${ADMIX_DIR}/rankin"
    "$FST_DIR"
    "$THETA_DIR"
    "$IBD_DIR"
    "$IBE_DIR" "${IBE_DIR}/env_data" "${IBE_DIR}/results"
    "$FIG_DIR" "${FIG_DIR}/main" "${FIG_DIR}/supplementary"
    "$TABLE_DIR" "${TABLE_DIR}/main" "${TABLE_DIR}/supplementary"
    "$LOG_DIR"
)

for d in "${dirs[@]}"; do
    mkdir -p "$d"
done

log_msg "  Created $(echo ${#dirs[@]}) directories."

# ─────────────────────────────────────────────────────────────────────────────
# 2. Symlink the genome reference
# ─────────────────────────────────────────────────────────────────────────────
log_msg "Symlinking genome reference..."

for ext in "" ".fai" ".dict"; do
    src="${GENOME}${ext}"
    dst="${INFO_DIR}/genome.fasta${ext}"
    if [[ -f "$src" ]] && [[ ! -e "$dst" ]]; then
        ln -s "$src" "$dst"
    fi
done

# Masked genome (deviant regions)
if [[ -f "$GENOME_MASKED" ]] && [[ ! -e "${INFO_DIR}/genome_masked.fasta" ]]; then
    ln -s "$GENOME_MASKED" "${INFO_DIR}/genome_masked.fasta"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 3. Copy essential metadata from your existing BAM/metadata source
# ─────────────────────────────────────────────────────────────────────────────
log_msg "Copying metadata from SOURCE_INFO_DIR..."

# BAM filelist
if [[ -f "${SOURCE_INFO_DIR}/bam.filelist" ]]; then
    cp "${SOURCE_INFO_DIR}/bam.filelist" "${INFO_DIR}/bam.filelist"
fi

# Regions file
if [[ -f "${SOURCE_INFO_DIR}/regions.txt" ]]; then
    cp "${SOURCE_INFO_DIR}/regions.txt" "${INFO_DIR}/regions.txt"
fi

# Regions number file
if [[ -f "${SOURCE_INFO_DIR}/regions_number.txt" ]]; then
    cp "${SOURCE_INFO_DIR}/regions_number.txt" "${INFO_DIR}/regions_number.txt"
fi

# Population file
if [[ -f "${SOURCE_INFO_DIR}/pop.txt" ]]; then
    cp "${SOURCE_INFO_DIR}/pop.txt" "${INFO_DIR}/pop.txt"
fi

# Info file (sample metadata)
if [[ -f "${SOURCE_INFO_DIR}/info.txt" ]]; then
    cp "${SOURCE_INFO_DIR}/info.txt" "${INFO_DIR}/info.txt"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 4. Generate BAM filelists by population
# ─────────────────────────────────────────────────────────────────────────────
log_msg "Generating per-population BAM filelists..."

if [[ -f "${INFO_DIR}/bam.filelist" ]] && [[ -f "${INFO_DIR}/pop.txt" ]]; then
    while read -r pop; do
        grep "${pop}s_" "${INFO_DIR}/bam.filelist" > "${INFO_DIR}/bamlists/${pop}.bamlist" 2>/dev/null || true
        n=$(wc -l < "${INFO_DIR}/bamlists/${pop}.bamlist")
        echo "  ${pop}: ${n} individuals"
    done < "${INFO_DIR}/pop.txt"
fi

# Regional BAM lists
log_msg "Generating regional BAM filelists..."
> "${INFO_DIR}/bamlists/rankin.bamlist"
for pop in $RANKIN_POPS; do
    cat "${INFO_DIR}/bamlists/${pop}.bamlist" >> "${INFO_DIR}/bamlists/rankin.bamlist" 2>/dev/null || true
done

> "${INFO_DIR}/bamlists/naujaat.bamlist"
for pop in $NAUJAAT_POPS; do
    cat "${INFO_DIR}/bamlists/${pop}.bamlist" >> "${INFO_DIR}/bamlists/naujaat.bamlist" 2>/dev/null || true
done

# ─────────────────────────────────────────────────────────────────────────────
# 5. Create regions file without sex chromosomes
# ─────────────────────────────────────────────────────────────────────────────
log_msg "Creating autosomal regions file..."

if [[ -f "${INFO_DIR}/regions.txt" ]]; then
    grep -v -E "$(echo $SEX_CHRS | tr ' ' '|')" "${INFO_DIR}/regions.txt" \
        > "${INFO_DIR}/regions_nosex.txt" || true
    n_auto=$(wc -l < "${INFO_DIR}/regions_nosex.txt")
    n_all=$(wc -l < "${INFO_DIR}/regions.txt")
    log_msg "  Autosomal regions: ${n_auto} / ${n_all} total"
fi

# ─────────────────────────────────────────────────────────────────────────────
# 6. Summary
# ─────────────────────────────────────────────────────────────────────────────
log_msg "=== SETUP COMPLETE ==="
echo ""
echo "Pipeline root:  ${BASE_DIR}"
echo "BAM files:      $(wc -l < "${INFO_DIR}/bam.filelist" 2>/dev/null || echo 'N/A') individuals"
echo "Regions:        $(wc -l < "${INFO_DIR}/regions.txt" 2>/dev/null || echo 'N/A') chromosomes/scaffolds"
echo "Populations:    $(wc -l < "${INFO_DIR}/pop.txt" 2>/dev/null || echo 'N/A') populations"
echo ""
echo "NOTE: bamlists for the reference populations (JAY, LLS, DV) are not"
echo "generated here — place them manually at \${INFO_DIR}/bamlists/{JAY,LLS,DV}.bamlist"
echo "before running the triangle-plot or D-stats scripts (12*, 13*)."
echo ""
echo "Next step: submit 01_ngsparalog.sh"
