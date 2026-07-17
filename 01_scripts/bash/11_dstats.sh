#!/bin/bash
# =============================================================================
# 11_dstats.sh — D-statistics (ABBA-BABA) via ANGSD doAbbababa2
# Tests Atlantic introgression (LLS) -> Kivalliq Canadian populations,
# using DV (Dolly Varden, Babbage River) as the distant outgroup (H4).
#
# D topology tested: D(H1, H2 ; H3, H4)
#   H1 = JAY (pure Beringian Arctic reference, Kitikmeot)
#   H2 = pop_test (each Kivalliq population in turn)
#   H3 = LLS (Sweden, Atlantic — potential introgression source)
#   H4 = DV  (Dolly Varden, S. malma — distant outgroup)
#
# Positive D + |Z| > 3 : pop_test shares more alleles with LLS than JAY does
#                        => Atlantic introgression in pop_test
# D ~ 0                : no asymmetric gene flow detected
#
# Output table: Table_Dstats_JAY.tsv (via plot_dstats.R downstream)
#
# Submit: sbatch 01_scripts/bash/11_dstats.sh
# =============================================================================
#SBATCH -J "11_dstats"
#SBATCH -o logs/11_dstats_%A_%a.out
#SBATCH -e logs/11_dstats_%A_%a.err
#SBATCH -c 8
#SBATCH -p medium
#SBATCH --time=2-00:00
#SBATCH --mem=80G
#SBATCH --array=0-12%8

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(pwd)}"
source config/00_config.sh

module load angsd
ulimit -S -n 4096

NB_CPU="${SLURM_CPUS_PER_TASK:-8}"

# --- Helper: resolve the bamlist path for a population ----------------------
resolve_bamlist() {
    local pop="$1"
    shift
    local candidates=(
        "${INFO_DIR}/bamlists/${pop}.bamlist"
        "${INFO_DIR}/${pop}bam.filelist"
        "$@"
    )
    for p in "${candidates[@]}"; do
        [[ -n "$p" && -f "$p" && -s "$p" ]] && { echo "$p"; return 0; }
    done
    return 1
}

# --- D-stats configuration --------------------------------------------------
CAN_POPS_ARRAY=(AKL AUL CRB DIA MEL ITI KGJ NOP PAM SUP TIN WHI HOR)

H1_POP="JAY"           # pure Arctic reference
H3_POP="LLS"           # Atlantic source (Sweden)
H4_POP="DV"            # distant outgroup (Dolly Varden)

H2_POP="${CAN_POPS_ARRAY[$SLURM_ARRAY_TASK_ID]}"

# Skip the degenerate case D(JAY, JAY, ...)
if [[ "$H2_POP" == "$H1_POP" ]]; then
    log_msg "Skip: H1 == H2 == $H1_POP (degenerate test)"
    exit 0
fi

# --- Output folder and bamlists ---------------------------------------------
DSTAT_DIR="${BASE_DIR:-$(pwd)}/14_dstats"
mkdir -p "$DSTAT_DIR"

BAMLIST_H1=$(resolve_bamlist "$H1_POP") || { echo "ERROR: no bamlist for $H1_POP"; exit 1; }
BAMLIST_H2=$(resolve_bamlist "$H2_POP") || { echo "ERROR: no bamlist for $H2_POP"; exit 1; }
BAMLIST_H3=$(resolve_bamlist "$H3_POP") || { echo "ERROR: no bamlist for $H3_POP"; exit 1; }
BAMLIST_H4=$(resolve_bamlist "$H4_POP") || { echo "ERROR: no bamlist for $H4_POP (expected ${INFO_DIR}/bamlists/${H4_POP}.bamlist)"; exit 1; }

# --- Prepare ANGSD input files -----------------------------------------------
TRIO_PREFIX="${DSTAT_DIR}/${H2_POP}_vs_${H1_POP}_H3-${H3_POP}_H4-${H4_POP}"
COMBINED_BAMLIST="${TRIO_PREFIX}.bamlist"
SIZE_FILE="${TRIO_PREFIX}.sizeFile"

# Skip if already done
if [[ -f "${TRIO_PREFIX}.abbababa2" ]]; then
    log_msg "D-stat already computed for H2=$H2_POP, skip"
    exit 0
fi

# 1. Concatenate robustly (handles files missing a trailing newline)
> "$COMBINED_BAMLIST"
for f in "$BAMLIST_H1" "$BAMLIST_H2" "$BAMLIST_H3" "$BAMLIST_H4"; do
    cat "$f" >> "$COMBINED_BAMLIST"
    if [[ -n $(tail -c1 "$f" 2>/dev/null) ]]; then
        echo "" >> "$COMBINED_BAMLIST"
    fi
done
sed -i '/^$/d' "$COMBINED_BAMLIST"

# 2. Count individuals (grep -c "^" is more robust here than wc -l)
N_H1=$(grep -c "^" "$BAMLIST_H1")
N_H2=$(grep -c "^" "$BAMLIST_H2")
N_H3=$(grep -c "^" "$BAMLIST_H3")
N_H4=$(grep -c "^" "$BAMLIST_H4")

# 3. Generate the sizeFile
printf "%d\n%d\n%d\n%d\n" "$N_H1" "$N_H2" "$N_H3" "$N_H4" > "$SIZE_FILE"

log_msg "=== D-STAT: D($H1_POP, $H2_POP ; $H3_POP, $H4_POP) ==="
log_msg "  N individuals: H1=$N_H1  H2=$N_H2  H3=$N_H3  H4=$N_H4"

# --- Reference and site parameters ------------------------------------------
REF="${INFO_DIR}/genome.fasta"
require_file "$REF"

SITES_OPT=""
if [[ -f "${INFO_DIR}/sites_ldpruned_${SUFFIX}_nosex" ]]; then
    SITES="${INFO_DIR}/sites_ldpruned_${SUFFIX}_nosex"
    SITES_OPT="-sites $SITES -rf ${INFO_DIR}/regions_nosex.txt"
    log_msg "  Using LD-pruned sites: $SITES"
else
    SITES_OPT="-rf ${INFO_DIR}/regions_nosex.txt"
    log_msg "  No pruned sites available, using all autosomal regions"
fi

# --- Run doAbbababa2 --------------------------------------------------------
angsd -P "$NB_CPU" -nQueueSize 50 \
    -bam "$COMBINED_BAMLIST" \
    -ref "$REF" \
    $SITES_OPT \
    -doAbbababa2 1 \
    -sizeFile "$SIZE_FILE" \
    -doCounts 1 \
    -useLast 1 \
    -blockSize 5000000 \
    -minMapQ 30 -minQ 20 \
    -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 \
    -out "$TRIO_PREFIX"

log_msg "=== doAbbababa2 COMPLETE for H2=$H2_POP ==="

# --- Compute D and jackknife Z-score ----------------------------------------
# ANGSD ships an R helper (R/estAvgError.R) for this; point ANGSD_R_SCRIPT at
# it via 00_config.sh (ANGSD_DIR) if your install layout differs.
ANGSD_R_SCRIPT="${ANGSD_DIR:-$(dirname "$(command -v angsd)")}/R/estAvgError.R"

if [[ -f "$ANGSD_R_SCRIPT" ]]; then
    log_msg "Computing D, Z-score via jackknife..."
    "$RSCRIPT" "$ANGSD_R_SCRIPT" \
        angsdFile="${TRIO_PREFIX}" \
        out="${TRIO_PREFIX}_Dstat" \
        sizeFile="$SIZE_FILE" \
        nameFile=<(printf "%s\n%s\n%s\n%s\n" "$H1_POP" "$H2_POP" "$H3_POP" "$H4_POP")
    log_msg "  D-stat written: ${TRIO_PREFIX}_Dstat.txt"
else
    log_msg "  WARNING: ANGSD R script not found at $ANGSD_R_SCRIPT — set ANGSD_DIR in 00_config.sh"
fi

log_msg "=== D-STAT $H2_POP COMPLETE ==="
