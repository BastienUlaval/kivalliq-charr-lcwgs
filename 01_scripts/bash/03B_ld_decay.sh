#!/bin/bash
# =============================================================================
# 03B_ld_decay.sh -- Fit LD decay curve
# Submit: sbatch 01_scripts/bash/03B_ld_decay.sh
# =============================================================================
#SBATCH -J "03B_ld_decay"
#SBATCH -o logs/03B_ld_decay_%j.out
#SBATCH -e logs/03B_ld_decay_%j.err
#SBATCH -c 1
#SBATCH -p medium
#SBATCH --time=7-00:00
#SBATCH --mem=200G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

BAMLIST="${INFO_DIR}/bam.filelist"
read -r N_IND _ _ <<< "$(compute_filters "$BAMLIST")"

log_msg "=== LD DECAY ==="

# Decompress LD files needed for decay (first 10 chromosomes)
for NUM in $(head -10 "${INFO_DIR}/regions_number.txt"); do
    GZ="${LD_DIR}/by_chr/all_${SUFFIX}_chr${NUM}.ld.gz"
    PLAIN="${LD_DIR}/by_chr/all_${SUFFIX}_chr${NUM}.ld"
    if [[ -f "$GZ" ]] && [[ ! -f "$PLAIN" ]]; then
        gunzip -k "$GZ"
    fi
done

# Create list of LD files
ls -1 "${LD_DIR}/by_chr/"*.ld 2>/dev/null | head -10 > "${LD_DIR}/ld_decay/list.ldfiles"

N_FILES=$(wc -l < "${LD_DIR}/ld_decay/list.ldfiles")
log_msg "  Using ${N_FILES} chromosomes for LD decay"

require_file "$FIT_LDDECAY" "fit_LDdecay.R script"

"$RSCRIPT" "$FIT_LDDECAY" \
    --ld_files "${LD_DIR}/ld_decay/list.ldfiles" \
    --ld r2 \
    --n_ind "$N_IND" \
    --max_kb_dist 500 \
    --fit_level 100 \
    --fit_boot 100 \
    --plot_data \
    --plot_no_legend \
    -o "${LD_DIR}/ld_decay/ld_decay.pdf"

# Clean up decompressed files
for NUM in $(head -10 "${INFO_DIR}/regions_number.txt"); do
    rm -f "${LD_DIR}/by_chr/all_${SUFFIX}_chr${NUM}.ld"
done

log_msg "=== LD DECAY COMPLETE ==="
