#!/bin/bash
# =============================================================================
# 03A_ngsLD.sh -- Estimate LD per chromosome with ngsLD
# Submit: cat regions_number.txt | parallel -j10 \
#   sbatch -c 4 --mem 20G -p small --time 1-00:00 \
#   01_scripts/bash/03A_ngsLD.sh {}
# =============================================================================

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

REGION_NUM="$1"
REGION=$(sed -n "${REGION_NUM}p" "${INFO_DIR}/regions.txt")
NB_CPU=4
BAMLIST="${INFO_DIR}/bam.filelist"

read -r N_IND MIN_IND MAX_DEPTH <<< "$(compute_filters "$BAMLIST")"

BEAGLE="${GL_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}.beagle.gz"
SITES="${INFO_DIR}/sites_by_chr/sites_${SUFFIX}_chr${REGION_NUM}_canonical"

require_file "$BEAGLE" "Beagle chr${REGION_NUM}"
require_file "$SITES" "Sites chr${REGION_NUM}"

# Prepare tab-separated positions
SITES_TAB="${SITES}.tab"
sed 's/ /\t/g' "$SITES" > "$SITES_TAB"

# Count sites for this chromosome
NSITES=$(zcat "$BEAGLE" | tail -n +2 | wc -l)

log_msg "=== ngsLD -- chr${REGION_NUM}: ${NSITES} SNPs ==="

"$NGSLD_PATH" \
    --geno "$BEAGLE" \
    --probs \
    --n_ind "$N_IND" \
    --n_sites "$NSITES" \
    --min_maf "$MIN_MAF" \
    --pos "$SITES_TAB" \
    --max_kb_dist "$LD_MAX_KB" \
    --rnd_sample "$LD_RND_SAMPLE" \
    --n_threads "$NB_CPU" \
    --out "${LD_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}.ld"

gzip -f "${LD_DIR}/by_chr/all_${SUFFIX}_chr${REGION_NUM}.ld"

log_msg "=== ngsLD chr${REGION_NUM} COMPLETE ==="
