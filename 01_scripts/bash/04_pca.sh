#!/bin/bash
# =============================================================================
# 04_pca.sh -- PCA via PCAngsd (global, Rankin, Naujaat)
# Submit: sbatch 01_scripts/bash/04_pca.sh
# =============================================================================
#SBATCH -J "04_pca"
#SBATCH -o logs/04_pca_%j.out
#SBATCH -e logs/04_pca_%j.err
#SBATCH -c 8
#SBATCH -p medium
#SBATCH --time=7-00:00
#SBATCH --mem=200G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

NB_CPU=8

# Setup PCAngsd python environment
PYTHON_ENV="${PCANGSD_ENV}/bin/python"
require_file "$PYTHON_ENV" "PCAngsd Python"
export PATH="${PCANGSD_ENV}/bin:$PATH"

# --- Helper function --------------------------------------------------------
run_pcangsd() {
    local beagle="$1"
    local out_prefix="$2"
    local bamlist="$3"
    local label="$4"

    log_msg "--- PCAngsd: ${label} ---"
    require_file "$beagle" "${label} Beagle"

    local n_snps
    n_snps=$(zcat "$beagle" | tail -n +2 | wc -l)
    log_msg "  SNPs: ${n_snps}"

    "$PYTHON_ENV" "${PCANGSD_PATH}/pcangsd.py" \
        -beagle "$beagle" \
        -o "$out_prefix" \
        -threads $NB_CPU \
        -iter 200 \
        -tole 1e-5 \
        -sites_save

    if [[ -f "${out_prefix}.cov" ]]; then
        log_msg "  Covariance matrix created."
        "$RSCRIPT" 01_scripts/R/pca_eigen.R "${out_prefix}.cov" "$bamlist"
        log_msg "  PCA eigendecomposition done."
    else
        log_msg "  ERROR: Covariance matrix not created!"
        return 1
    fi
}

# --- Run three analyses -----------------------------------------------------
GLOBAL_BEAGLE="${LD_DIR}/all_${SUFFIX}.prunednosex.beagle.gz"
RANKIN_BEAGLE="${LD_DIR}/rankin_${SUFFIX}.prunednosex.beagle.gz"
NAUJAAT_BEAGLE="${LD_DIR}/naujaat_${SUFFIX}.prunednosex.beagle.gz"

run_pcangsd "$GLOBAL_BEAGLE" \
    "${PCA_DIR}/global_${SUFFIX}_prunednosex" \
    "${INFO_DIR}/bam.filelist" \
    "GLOBAL"

run_pcangsd "$RANKIN_BEAGLE" \
    "${PCA_DIR}/rankin_${SUFFIX}_prunednosex" \
    "${INFO_DIR}/bamlists/rankin.bamlist" \
    "RANKIN"

run_pcangsd "$NAUJAAT_BEAGLE" \
    "${PCA_DIR}/naujaat_${SUFFIX}_prunednosex" \
    "${INFO_DIR}/bamlists/naujaat.bamlist" \
    "NAUJAAT"

# --- Produce publication figures ---------------------------------------------
log_msg "Generating publication-quality PCA figures..."

"$RSCRIPT" 01_scripts/R/plot_pca.R \
    "${PCA_DIR}" \
    "${INFO_DIR}/info.txt" \
    "${FIG_DIR}/main"

log_msg "=== PCA COMPLETE ==="
